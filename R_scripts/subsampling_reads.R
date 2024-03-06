#!/usr/bin/env Rscript
### load input variables ###
args = commandArgs(trailingOnly=TRUE)

# to run the script from the terminal: Rscript subsampling_reads.R [list of arguments]
# Arguments: path to the fastq file of each fraction (path_fastq_chr_ass, path_fastq_nucleo, path_fastq_cyto), path to the reference genome (optional), 
# path to the gtf file (optional), number of threads for the mapping (optional), minimum number of reads mapping on each gene (optional), 
# number of gene-level samplings to do (optional), path to minimap2, samtools and seqtk tools (optional) 
# and the different conditions with the same order of the corresponding fastq files (cond1, cond2, cond3)

for(v in args)
{
  vTmp <- strsplit(v,"=")[[1]]
  assign(vTmp[[1]],vTmp[[2]])
}

# default parameters
if (!exists('path_reference_genome')) {
  path_reference_genome <- '/path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
}

name_ref_genome <- strsplit(strsplit(path_reference_genome,'/')[[1]][length(strsplit(path_reference_genome,'/')[[1]])], '\\.fa')[[1]][1]

if (!exists('path_gtf_file')) {
  path_gtf_file <- '/path/to/Homo_sapiens.GRCh38.104.gtf'
}

if (!exists('threads')) {
  threads <- '17'
}

if (!exists('num_reads')) {
  num_reads <- '20'
}

if (!exists('num_subsampling')) {
  num_subsampling <- '5'
}

# minimap2 v2.26
if (!exists('minimap2')) {
  minimap2 <- '/path/to/minimap2'
}

# samtools v1.17
if (!exists('samtools')) {
  samtools <- '/path/to/samtools'
}

# seqtk v1.4
if (!exists('seqtk')) {
  seqtk <- '/path/to/seqtk'
}

suppressMessages(library("GenomicAlignments"))
suppressMessages(library("GenomicFeatures"))
suppressMessages(library("GenomicRanges"))

# list with the path of each fastq file
file_fastq <- list(path_fastq_chr_ass, path_fastq_cyto, path_fastq_nucleo)

# list with the conditions
order_conditions <- c(cond1,cond2,cond3)

subsampling <- function(file_fastq) {
  # for each fastq file perform the genome mapping and filtering if they have not been performed yet
  bam_filtered <- lapply(file_fastq, function(x) {
    path_bam <- paste0(strsplit(x,'\\.fastq')[[1]][1], '_spliced_mapped_to_', name_ref_genome, '.bam')
    path_bam_filtered <- paste0(strsplit(path_bam,'\\.bam')[[1]][1], '_filtered.bam')
    
    if (!file.exists(path_bam) & !file.exists(path_bam_filtered)) {
      system(command = paste(minimap2, "-ax splice -k 14 --seed 1 -t", threads, path_reference_genome, x, "|", samtools, "view -hSb |", samtools, "sort -o", path_bam))
      system(command = paste(samtools, "view -bh -q 0 -F2308", path_bam, "|", samtools, "sort -o", path_bam_filtered))
    } else if (file.exists(path_bam) & !file.exists(path_bam_filtered)) {
      system(command = paste(samtools, "view -bh -q 0 -F2308", path_bam, "|", samtools, "sort -o", path_bam_filtered))
    } else if (!file.exists(path_bam) & file.exists(path_bam_filtered)) {
      system(command = paste(minimap2, "-ax splice -k 14 --seed 1 -t", threads, path_reference_genome, x, "|", samtools, "view -hSb |", samtools, "sort -o", path_bam))
    }
    return(path_bam_filtered)
  })
  
  # read each filtered.bam file as table, each line reports one read 
  table_filtered_bam <- lapply(bam_filtered, function(x){
    table_bam <- readGAlignments(file = x, use.names = TRUE)
    return(table_bam)
  })
  
  txdb <- makeTxDbFromGFF(path_gtf_file)
  
  # genes
  genes_txdb <- GenomicFeatures::genes(txdb)
  
  # remove the reads mapping on multiple genes and return a table with the uniquely mapping reads
  table_filtered_bam_without_multiple_mapping_reads <- lapply(table_filtered_bam, function(x){
    genes_over <- findOverlaps(x,genes_txdb)
    table_bam <- x[queryHits(genes_over)][isUnique(queryHits(genes_over))]
    return(c(length(table_bam),table_bam))
  })
  
  # count the minimum number of reads mapping only on one gene across the three fractions, 
  # it will be used as threshold for the subsampling at library level
  threshold_library <- min(table_filtered_bam_without_multiple_mapping_reads[[1]][[1]], table_filtered_bam_without_multiple_mapping_reads[[2]][[1]], table_filtered_bam_without_multiple_mapping_reads[[3]][[1]])
  print(paste('Minimum number of mapping reads across the three fractions: ', as.character(threshold_library)))
  
  # for each fraction extract randomly a number of reads mapping only on one gene equal to the library level threshold
  table_subsampled_reads <- lapply(table_filtered_bam_without_multiple_mapping_reads, function(x){
    set.seed(1)
    selected_reads <- sample(names(x[[2]]), threshold_library)
    table_sub_reads <- x[[2]][selected_reads,]
    return(table_sub_reads)
  }) 
  
  # identify on which gene each read (after the library subsampling) maps
  overlapping_reads_gene <- lapply(table_subsampled_reads, function(x){
    genes_over2 <- findOverlaps(x,genes_txdb)
    return(genes_over2)
  })
  
  # count the number of reads on each gene
  num_reads_per_gene <- lapply(overlapping_reads_gene, function(x){
    gene_names <- names(genes_txdb[subjectHits(x)])
    genes_over_split <- split(queryHits(x), gene_names)
    number_reads_per_gene <- unlist(lapply(genes_over_split, length))
    return(number_reads_per_gene)
  })  
  
  # identify all the genes
  all_genes <- union(names(num_reads_per_gene[1][[1]]), names(num_reads_per_gene[2][[1]]))
  all_genes <- union(all_genes, names(num_reads_per_gene[3][[1]]))
  print(paste('Overall number of genes with at least one read in one of the fractions: ', as.character(length(all_genes))))
  
  # create a matrix with a row per gene and three columns, one for each fraction (initialize to 0)
  genes_counts <- matrix(data = 0, nrow = length(all_genes), ncol=3) 
  rownames(genes_counts) <- all_genes
  colnames(genes_counts) <- order_conditions
  
  # update the number of reads mapping on each gene in each fraction (letting 0 when there are no reads mapping on the gene)
  genes_counts[names(num_reads_per_gene[1][[1]]), order_conditions[1]] <- unname(num_reads_per_gene[1][[1]])
  genes_counts[names(num_reads_per_gene[2][[1]]), order_conditions[2]] <- unname(num_reads_per_gene[2][[1]])
  genes_counts[names(num_reads_per_gene[3][[1]]), order_conditions[3]] <- unname(num_reads_per_gene[3][[1]])
  
  # filter the genes with at least 20 mapping reads in all the fractions
  genes_counts_over_threshold <- genes_counts[which(genes_counts[,order_conditions[1]]>= as.numeric(num_reads) & genes_counts[,order_conditions[2]]>= as.numeric(num_reads) & genes_counts[,order_conditions[3]]>= as.numeric(num_reads)),]
  
  # save the minimum number of reads - across the three fractions - mapping on each gene with at least 20X in all the fractions
  genes_counts_over_threshold_min <- apply(genes_counts_over_threshold,1,min)
  path <- system(command = paste('dirname', path_fastq_cyto), intern=TRUE)
  save(genes_counts_over_threshold_min, file= paste0(path,'/genes_counts_over_threshold_min.Rda'))
  
  # save the names of the genes with at least 20X in all the fractions
  selected_genes <- rownames(genes_counts_over_threshold)
  print(paste('Number of genes with at least', num_reads, 'reads in all the fractions: ', as.character(length(selected_genes))))
  write(selected_genes, paste0(path,"/selected_genes_total_reads.txt"))
  
  # save a bed file with the coordinates of the selected genes
  bed <- cbind(as.data.frame(seqnames(genes_txdb)), start(genes_txdb), end(genes_txdb), genes_txdb$gene_id,'.', as.data.frame(strand(genes_txdb)))
  write.table(bed, paste0(strsplit(path_gtf_file, '\\.gtf')[[1]][1], '.genes.bed'), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  system(command=paste('cat', paste0(strsplit(path_gtf_file, '\\.gtf')[[1]][1], '.genes.bed'), '| grep -f', paste0(path,"/selected_genes_total_reads.txt"), '>', paste0(strsplit(path_gtf_file, '\\.gtf')[[1]][1], '.selected_genes_total_reads.bed')))
  
  path_actual_directory <- system(command='realpath .', intern = TRUE)
  
  for (i in 1:length(order_conditions)) {
    # for each fraction a dataframe is generated with the read IDs (reads after library level subsampling and only those mapping on one gene) 
    # and the gene on which they map
    read_gene <- data.frame(read_id = names(table_subsampled_reads[[i]][queryHits(overlapping_reads_gene[[i]])]),gene_name =names(genes_txdb[subjectHits(overlapping_reads_gene[[i]])]))
    
    for (j in 1:as.numeric(num_subsampling)) {      
      selected_reads_j <- c()
      
      # for each of the selected genes find which reads overlap with it in the analysed fraction and randomly
      # extract a number of reads equal to the minimum number of reads mapping on that gene between the fractions
      for (gene in selected_genes) {
        set.seed(j)
        reads_on_gene <- as.vector(read_gene[read_gene$gene_name==gene,]$read_id)
        selected_reads_j <- c(selected_reads_j, sample(reads_on_gene,genes_counts_over_threshold_min[gene]))
      }
      
      path_fastq_j <- paste0(strsplit(file_fastq[[i]],'\\.fastq')[[1]][1],"_selected_reads_",as.character(j),'.fastq')
      write(selected_reads_j, paste0(path_actual_directory, "/selected_reads_", order_conditions[i], "_",as.character(j),".txt"))      
      reads_j_txt <- paste0(path_actual_directory, "/selected_reads_", order_conditions[i], "_",as.character(j),".txt")
      system(command=paste(seqtk, "subseq", file_fastq[[i]], reads_j_txt, ">", path_fastq_j))
      print(path_fastq_j)
      system(command=paste("rm", reads_j_txt))
    }}   
  
}

subsampling(file_fastq)


