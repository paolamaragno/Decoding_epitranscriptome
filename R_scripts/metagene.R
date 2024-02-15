library('GenomicFeatures')
library('GenomicRanges')
library('pheatmap')
library('RColorBrewer')
library('cartography')
library('xlsx')
library('rtracklayer')
library('grid')

gtf_file <- "/path/to/Homo_sapiens.GRCh38.104.gtf"
gtf <- readGFF(gtf_file)
# extract the annotation of the genes
genes <- gtf[gtf$type=='gene',]

# GRangesList containing for each gene a GRanges with the coordinates of its 5'UTRs, 3'UTRs, stop codon, introns and exons
load('/path/to/R_data/protein_coding_genes_5UTR_3UTR_introns_exons_stop.RDa')

# identify the gene names of protein coding genes and how many of the annotated genes are protein coding
protein_coding <- genes[genes$gene_biotype=='protein_coding',]
protein_coding_names <- protein_coding$gene_id
length(names(protein_coding_genes_5UTR_3UTR_introns_exons_stop)[names(protein_coding_genes_5UTR_3UTR_introns_exons_stop) %in% protein_coding_names])

# extract the gene biotype of the genes annotated in protein_coding_genes_5UTR_3UTR_introns_exons_stop
type <- unlist(lapply(as.list(names(protein_coding_genes_5UTR_3UTR_introns_exons_stop)), function(x) {
  genes[genes$gene_id == x, 13]
}))

# create a data frame reporting for each annotated gene the corresponding biotype 
gene_biotype <- data.frame(gene_biotype=type)
rownames(gene_biotype) <- names(protein_coding_genes_5UTR_3UTR_introns_exons_stop)
save(gene_biotype, file='/path/to/R_data/gene_biotype.Rda')

#########################################

txdb <- makeTxDbFromGFF(gtf_file)
genes_txdb <- GenomicFeatures::genes(txdb)

# directory contains /hits_ELIGOS/ folder in which there are the RData objects with all ELIGOS hits or ELIGOS DRACH+ hits, 
# separately for each fraction, and /without_DRACH/ folder in which there are the RData objects with all ELIGOS hits in DRACH negative context, 
# separately for each fraction.
# create a folder in "directory" called /gene_counts/ where you save genes_counts_over_threshold_min.Rda object generated
# by subsampling_reads.R script during the library-level and gene-level subsamplings
count_mods <- function(directory) { 
  
  f2 <- function(R_objects) {
    
    # create three matrix that will be used to report the Spearman correlation between the number of hits
    # identified in each gene part and different gene features
    correlation_read_counts <- matrix(nrow=3, ncol=5)
    colnames(correlation_read_counts) <- c('5UTR', 'coding exon', 'intron', 'stop codon', '3UTR')
    rownames(correlation_read_counts) <- c('Chromatin', 'Nucleoplasm', 'Cytoplasm')
    correlation_gene_length <- matrix(nrow=3, ncol=5)
    colnames(correlation_gene_length) <- c('5UTR', 'coding exon', 'intron', 'stop codon', '3UTR')
    rownames(correlation_gene_length) <- c('Chromatin', 'Nucleoplasm', 'Cytoplasm')
    correlation_strand <- matrix(nrow=3, ncol=5)
    colnames(correlation_strand) <- c('5UTR', 'coding exon', 'intron', 'stop codon', '3UTR')
    rownames(correlation_strand) <- c('Chromatin', 'Nucleoplasm', 'Cytoplasm')
    
    # iterate on each RData contained in /hits_ELIGOS/ folder
    for (obj in R_objects) {
      n <- load(obj)
      hits <- get(n)
      genes <- unique(hits$gene_id)
      strand_genes <- as.vector(strand(genes_txdb[genes]))
      strand_genes[strand_genes=='+'] <- '1'
      strand_genes[strand_genes=='-'] <- '0'
      strand_genes_numeric <- as.numeric(strand_genes)
      names(strand_genes_numeric) <- genes
      
      # create a matrix in which, for each gene, the number of hits mapping on each gene part will be annotated
      m <- matrix(0, ncol = 5, nrow = length(genes)) 
      colnames(m) <- c('5UTR', 'coding exon', 'intron', 'stop codon', '3UTR')
      rownames(m) <- genes 
      # create a matrix to report for each gene whether it lacks the annotation of any gene part  
      m2 <- matrix(1, ncol = 5, nrow = length(genes)) 
      colnames(m2) <- c('5UTR', 'coding exon', 'intron', 'stop codon', '3UTR')
      rownames(m2) <- genes
      
      Strand <- rep(c('+','-'), 5) 
      Gene_part <- rep(c('5UTR', 'coding exon', 'intron', 'stop codon', '3UTR'), each = 2)
      Hits <- rep(0, 10)
      stats <- data.frame(Strand,Gene_part,Hits)
      
      # load the R object containing, for each gene, the minimum number of reads mapping on that gene across the three fractions
      load(paste0(directory,'/gene_counts/genes_counts_over_threshold_min.Rda'))
      read_counts <- genes_counts_over_threshold_min[rownames(m)] # read counts of the genes with at least one modification
      length_genes <- width(genes_txdb[rownames(m)]) # length of the genes with at least one modification
      strand_genes <- as.vector(strand(genes_txdb[genes])) # strand of the genes with at least one modification
      names(read_counts) <- rownames(m)
      names(length_genes) <- rownames(m)
      names(strand_genes) <- rownames(m)
      
      # remove the genes on which ELIGOS has identified hits but with no reads mapping on them (because minimap2 executed inside
      # ELIGOS pipeline may have assigned these reads to a gene different from the one chosen by minimap2 during the subsampling)
      if (length(names(read_counts[which(is.na(read_counts))])) != 0) {
        for (name in names(read_counts[which(is.na(read_counts))])) {
          read_counts <- read_counts[-which(names(read_counts) == name)]
          length_genes <- length_genes[-which(names(length_genes) == name)]
          strand_genes_numeric <- strand_genes_numeric[-which(names(strand_genes_numeric) == name)]
          strand_genes <- strand_genes[-which(names(strand_genes) == name)]
          m <- m[-which(rownames(m) == name),]
          m2 <- m2[-which(rownames(m2) == name),]
        }
      }
      
      # create a matrix to annotate the biotype of the genes that lack of one or more gene parts
      na_annotation <- matrix(0, nrow=3, ncol=6)
      rownames(na_annotation) <- c('protein_coding','Mt_rRNA', 'pseudogene')
      colnames(na_annotation) <- c('5UTR', 'coding exon', 'intron', 'stop codon', '3UTR', 'Tot')
      
      na_annotation[1,6] <- length(gene_biotype[rownames(m),][gene_biotype[rownames(m),]=='protein_coding'])
      na_annotation[2,6] <- length(gene_biotype[rownames(m),][gene_biotype[rownames(m),]=='Mt_rRNA'])
      na_annotation[3,6] <- length(rownames(m)) - sum(na_annotation[1,6],na_annotation[2,6])
      
      # iterate over the hits: considering the gene on which each hit maps, report 0 in m2 in correspondence of the gene
      # region(s) (if any) not present in that gene. Update m adding +1 in correspondence of the gene (row) and the gene part (column)
      # on which the hit maps
      for (hit in 1:length(hits)) {
        gene <- hits[hit]$gene_id
        if (gene %in% rownames(m)) {
          parts_gene <- protein_coding_genes_5UTR_3UTR_introns_exons_stop[gene][[1]]
          for (colname in colnames(m2)) {
            if (!colname %in% unique(parts_gene$feature)) {
              m2[gene,colname] <- 0
            }
          }
          over <- findOverlaps(hits[hit], parts_gene)
          if (length(subjectHits(over)) != 0) {
            which_features <- parts_gene[subjectHits(over)]$feature   
            for (colname in colnames(m)) { 
              if (colname %in% which_features) {
                strand_hit <- unique(as.vector(strand(hits[hit])))
                stats[(stats$Strand == strand_hit) & (stats$Gene_part == colname),3] <- stats[(stats$Strand == strand_hit) & (stats$Gene_part == colname),3] + 1
                m[gene, colname] <- m[gene, colname] + 1
              }
            }
          }
        }
      }
      
      # use the information in m2 to indicate in m with NA which are the genes that lack of one or more gene parts and which
      # they are 
      for (r in 1:nrow(m2)) {
        for (c in 1:ncol(m2)) {
          if (m2[r,c] == 0) {
            m[r,c] <- NA
            gene_type <- gene_biotype[rownames(m2)[r],]
            if ('pseudogene' %in% unlist(strsplit(gene_type, split='_'))) {
              na_annotation[3,c] <- na_annotation[3,c]+1
            } else {
              na_annotation[gene_type,c] <- na_annotation[gene_type,c]+1
            }
          }
        }
      }
      
      print(na_annotation)
      
      if ('chr' %in% unlist(strsplit(n, split='_'))) {
        fraction = c('Chromatin')
      } else if ('nucleo' %in% unlist(strsplit(n, split='_'))) {
        fraction = c('Nucleoplasm')
      } else {
        fraction = c('Cytoplasm')
      }
      
      # compute the Spearman correlation between the number of hits identified in each gene part and 
      # the number of reads mapping on the corresponding gene
      cor_read_counts <- unlist(lapply(1:ncol(m), function(x,y,i){
        round(cor(x[!is.na(y[,i])], y[,i][!is.na(y[,i])], method = 'spearman'),2)
      }, x = read_counts, y = m))
      cor_read_counts[is.na(cor_read_counts)] <- 0
      correlation_read_counts[fraction,] <- cor_read_counts
      
      # compute the Spearman correlation between the number of hits identified in each gene part and 
      # the length of the corresponding gene
      cor_gene_length <- unlist(lapply(1:ncol(m), function(x,y,i){
        round(cor(x[!is.na(y[,i])], y[,i][!is.na(y[,i])], method = 'spearman'),2)
      }, x = length_genes, y = m))
      cor_gene_length[is.na(cor_gene_length)] <- 0
      correlation_gene_length[fraction,] <- cor_gene_length
      
      # compute the Spearman correlation between the number of hits identified in each gene part and 
      # the strand of the corresponding gene
      cor_gene_strand <- unlist(lapply(1:ncol(m), function(x,y,i){
        round(cor(x[!is.na(y[,i])], y[,i][!is.na(y[,i])], method = 'spearman'),2)
      }, x = strand_genes_numeric, y = m))
      cor_gene_strand[is.na(cor_gene_strand)] <- 0
      correlation_strand[fraction,] <- cor_gene_strand
      
      # saturation of the number of hits mapping on each gene part
      if ('DRACH' %in% unlist(strsplit(n, split='_'))) {
        # saturation at 5
        for (col in 1:ncol(m)) {
          for (row in 1:nrow(m)) {
            if (! is.na(m[row,col])) {
              if (m[row,col] > 5) {
                m[row,col] = 5
              }
            }
          }
        }
      } else {
        # saturation at 10
        for (col in 1:ncol(m)) {
          for (row in 1:nrow(m)) {
            if (! is.na(m[row,col])) {
              if (m[row,col] > 10) {
                m[row,col] = 10
              }
            }
          }
        }
      }
      
      # evaluate if there is any kind of correlation between the length of the coding exons and the length of the 3'UTR
      # of the genes with at least 5 hits falling in the coding exons (when considering all the hits) or
      # at least 3 hits falling in the coding exons (when considering only DRACH+ hits)
      if ('DRACH' %in% unlist(strsplit(n, split='_'))) {
        genes_hits_coding <- rownames(m[m[,2]>3,])
      } else {
        genes_hits_coding <- rownames(m[m[,2]>5,])
      }
      
      length_coding_exons <- unlist(lapply(protein_coding_genes_5UTR_3UTR_introns_exons_stop[genes_hits_coding], function(x) {
        sum(width(x[x$feature == 'coding exon']))
      }))
      length_3UTR <- unlist(lapply(protein_coding_genes_5UTR_3UTR_introns_exons_stop[genes_hits_coding], function(x) {
        sum(width(x[x$feature == '3UTR']))
      }))
      
      if ('without' %in% unlist(strsplit(n, split='_'))) { 
        pdf(paste0(directory, '/without_DRACH/', n, '_hist.pdf'), width = 13, height = 7)
        par(mfrow=c(1,2))
        hist(length_coding_exons, main = as.character(cor(length_coding_exons,length_3UTR, method = 'spearman')))
        hist(length_3UTR)
        dev.off()
      } else {
        pdf(paste0(directory, '/hits_ELIGOS/', n, '_hist.pdf'), width = 13, height = 7)
        par(mfrow=c(1,2))
        hist(length_coding_exons, main = as.character(cor(length_coding_exons,length_3UTR, method = 'spearman')))
        hist(length_3UTR)
        dev.off()
      }
      
      ann <- data.frame(log10_gene_length = log10(length_genes)[rownames(m)], gene_strand = strand_genes[rownames(m)], log10_read_counts = log10(read_counts)[rownames(m)])
      rownames(ann) <- rownames(m)
      
      # plot the matrix m as heatmap
      if (!'DRACH' %in% unlist(strsplit(n, split='_'))) { 
        pdf(paste0(directory, '/hits_ELIGOS/', n, '.pdf'), height = 7, width = 6.5)
        pheatmap(m, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE, na_col = 'white',
                 annotation_row = ann, treeheight_row = 0, main = fraction, fontsize = 10, color = carto.pal(pal1 = 'red.pal',n1 = 11), 
                 legend_breaks = 0:10, legend_labels = c(as.character(0:9), '≥10'), fontsize_row = 10, fontsize_col = 10)
        grid.text('Transcriptional units', x=0.68, y=0.52, rot=270)
        dev.off()
      } else if (('DRACH' %in% unlist(strsplit(n, split='_'))) & (!'without' %in% unlist(strsplit(n, split='_')))) {
        pdf(paste0(directory, '/hits_ELIGOS/', n, '.pdf'), height = 7, width = 6.5)
        pheatmap(m, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE, na_col = 'white',
                 annotation_row = ann, treeheight_row = 0, main = fraction, fontsize = 10, color = brewer.pal(6, 'Reds'), 
                 legend_breaks = 0:5, legend_labels = c(as.character(0:4), '≥5'), fontsize_row = 10, fontsize_col = 10)
        grid.text('Transcriptional units', x=0.7, y=0.52, rot=270)
        dev.off()
      } else {
        pdf(paste0(directory, '/without_DRACH/', n, '.pdf'), height = 7, width = 6.5)
        pheatmap(m, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE, na_col = 'white',
                 annotation_row = ann, treeheight_row = 0, main = fraction, fontsize = 10, color = brewer.pal(6, 'Reds'), 
                 legend_breaks = 0:5, legend_labels = c(as.character(0:4), '≥5'), fontsize_row = 10, fontsize_col = 10)
        grid.text('Transcriptional units', x=0.7, y=0.52, rot=270)
        dev.off()
      }
    }

    if (!'DRACH' %in% unlist(strsplit(n, split='_'))) {
      write.xlsx(x = data.frame(correlation_read_counts),file = paste0(directory, '/hits_ELIGOS/correlation_read_counts.xlsx'),col.names = TRUE, row.names = TRUE)
      write.xlsx(x = data.frame(correlation_gene_length), file=paste0(directory, '/hits_ELIGOS/correlation_gene_length.xlsx'),col.names = TRUE, row.names = TRUE)
      write.xlsx(x = data.frame(correlation_strand),file=paste0(directory, '/hits_ELIGOS/correlation_strand.xlsx'),col.names = TRUE, row.names = TRUE)
    }  else if (('DRACH' %in% unlist(strsplit(n, split='_'))) & (!'without' %in% unlist(strsplit(n, split='_')))) {
      write.xlsx(x = data.frame(correlation_read_counts),file = paste0(directory, '/hits_ELIGOS/correlation_read_counts_DRACH.xlsx'),col.names = TRUE, row.names = TRUE)
      write.xlsx(x = data.frame(correlation_gene_length),file=paste0(directory, '/hits_ELIGOS/correlation_gene_length_DRACH.xlsx'),col.names = TRUE, row.names = TRUE)
      write.xlsx(x = data.frame(correlation_strand),file=paste0(directory, '/hits_ELIGOS/correlation_strand_DRACH.xlsx'),col.names = TRUE, row.names = TRUE)
    } else {
      write.xlsx(x = data.frame(correlation_read_counts),file = paste0(directory, '/without_DRACH/correlation_read_counts_without_DRACH.xlsx'),col.names = TRUE, row.names = TRUE)
      write.xlsx(x = data.frame(correlation_gene_length),file=paste0(directory, '/without_DRACH/correlation_gene_length_without_DRACH.xlsx'),col.names = TRUE, row.names = TRUE)
      write.xlsx(x = data.frame(correlation_strand),file=paste0(directory, '/without_DRACH/correlation_strand_without_DRACH.xlsx'),col.names = TRUE, row.names = TRUE)
    }
    }
  
  R_objects <- list.files(path = paste0(directory, '/hits_ELIGOS/'), pattern = 'Rda', full.names = TRUE)
  R_objects_all <- R_objects[unlist(lapply(strsplit(R_objects,'_'), function(x) {!'DRACH.Rda' %in% x}))]
  f2(R_objects_all)
  R_objects_DRACH <- R_objects[unlist(lapply(strsplit(R_objects,'_'), function(x) {'DRACH.Rda' %in% x}))]
  f2(R_objects_DRACH)
}

count_mods(directory= '/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/')
count_mods(directory='/path/to/fractions_eligos_4sU_library_gene_subsampling_rep1_min05_min05_mag1/')
count_mods(directory='/path/to/fractions_eligos_4sU_library_gene_subsampling_rep2_min05_min05_mag1/')
count_mods(directory='/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/')
count_mods(directory='/path/to/fractions_eligos_4sU_library_gene_subsampling_total_510645_min05_min05_mag1/')


###########

# function identical to the previous one but the path of input and output are different to adapt it 
# to the location in which ELIGOS DRACH+ and DRACH- hits only in nascent chromatin associated RNAs are stored
count_mods_nascent_only_chr <- function(directory) { 
  
  f2 <- function(R_objects) {
    
    # create three matrix that will be used to report the Spearman correlation between the number of hits
    # identified in each gene part and different gene features
    correlation_read_counts <- matrix(nrow=3, ncol=5)
    colnames(correlation_read_counts) <- c('5UTR', 'coding exon', 'intron', 'stop codon', '3UTR')
    rownames(correlation_read_counts) <- c('Chromatin', 'Nucleoplasm', 'Cytoplasm')
    correlation_gene_length <- matrix(nrow=3, ncol=5)
    colnames(correlation_gene_length) <- c('5UTR', 'coding exon', 'intron', 'stop codon', '3UTR')
    rownames(correlation_gene_length) <- c('Chromatin', 'Nucleoplasm', 'Cytoplasm')
    correlation_strand <- matrix(nrow=3, ncol=5)
    colnames(correlation_strand) <- c('5UTR', 'coding exon', 'intron', 'stop codon', '3UTR')
    rownames(correlation_strand) <- c('Chromatin', 'Nucleoplasm', 'Cytoplasm')
    
    # iterate on each RData contained in directory
    for (obj in R_objects) {
      n <- load(obj)
      hits <- get(n)
      genes <- unique(hits$gene_id)
      strand_genes <- as.vector(strand(genes_txdb[genes]))
      strand_genes[strand_genes=='+'] <- '1'
      strand_genes[strand_genes=='-'] <- '0'
      strand_genes_numeric <- as.numeric(strand_genes)
      names(strand_genes_numeric) <- genes
      
      # create a matrix in which, for each gene, the number of hits mapping on each gene part will be annotated
      m <- matrix(0, ncol = 5, nrow = length(genes)) 
      colnames(m) <- c('5UTR', 'coding exon', 'intron', 'stop codon', '3UTR')
      rownames(m) <- genes 
      # create a matrix to report for each gene whether it lacks the annotation of any gene part  
      m2 <- matrix(1, ncol = 5, nrow = length(genes)) 
      colnames(m2) <- c('5UTR', 'coding exon', 'intron', 'stop codon', '3UTR')
      rownames(m2) <- genes
      
      Strand <- rep(c('+','-'), 5) 
      Gene_part <- rep(c('5UTR', 'coding exon', 'intron', 'stop codon', '3UTR'), each = 2)
      Hits <- rep(0, 10)
      stats <- data.frame(Strand,Gene_part,Hits)
      
      # load the R object containing, for each gene, the minimum number of reads mapping on that gene across the three fractions
      load(paste0(directory,'/gene_counts/genes_counts_over_threshold_min.Rda'))
      read_counts <- genes_counts_over_threshold_min[rownames(m)] # read counts of the genes with at least one modification
      length_genes <- width(genes_txdb[rownames(m)]) # length of the genes with at least one modification
      strand_genes <- as.vector(strand(genes_txdb[genes])) # strand of the genes with at least one modification
      names(read_counts) <- rownames(m)
      names(length_genes) <- rownames(m)
      names(strand_genes) <- rownames(m)
      
      # remove the genes on which ELIGOS has identified hits but with no reads mapping on them (because minimap2 executed inside
      # ELIGOS pipeline may have assigned these reads to a gene different from the one chosen by minimap2 during the subsampling)
      if (length(names(read_counts[which(is.na(read_counts))])) != 0) {
        for (name in names(read_counts[which(is.na(read_counts))])) {
          read_counts <- read_counts[-which(names(read_counts) == name)]
          length_genes <- length_genes[-which(names(length_genes) == name)]
          strand_genes_numeric <- strand_genes_numeric[-which(names(strand_genes_numeric) == name)]
          strand_genes <- strand_genes[-which(names(strand_genes) == name)]
          m <- m[-which(rownames(m) == name),]
          m2 <- m2[-which(rownames(m2) == name),]
        }
      }
      
      # create a matrix to annotate the biotype of the genes that lack of one or more gene parts
      na_annotation <- matrix(0, nrow=3, ncol=6)
      rownames(na_annotation) <- c('protein_coding','Mt_rRNA', 'pseudogene')
      colnames(na_annotation) <- c('5UTR', 'coding exon', 'intron', 'stop codon', '3UTR', 'Tot')
      
      na_annotation[1,6] <- length(gene_biotype[rownames(m),][gene_biotype[rownames(m),]=='protein_coding'])
      na_annotation[2,6] <- length(gene_biotype[rownames(m),][gene_biotype[rownames(m),]=='Mt_rRNA'])
      na_annotation[3,6] <- length(rownames(m)) - sum(na_annotation[1,6],na_annotation[2,6])
      
      # iterate over the hits: considering the gene on which each hit maps, report 0 in m2 in correspondence of the gene
      # region(s) (if any) not present in that gene. Update m adding +1 in correspondence of the gene (row) and the gene part (column)
      # on which the hit maps
      for (hit in 1:length(hits)) {
        gene <- hits[hit]$gene_id
        if (gene %in% rownames(m)) {
          parts_gene <- protein_coding_genes_5UTR_3UTR_introns_exons_stop[gene][[1]]
          for (colname in colnames(m2)) {
            if (!colname %in% unique(parts_gene$feature)) {
              m2[gene,colname] <- 0
            }
          }
          over <- findOverlaps(hits[hit], parts_gene)
          if (length(subjectHits(over)) != 0) {
            which_features <- parts_gene[subjectHits(over)]$feature   
            for (colname in colnames(m)) { 
              if (colname %in% which_features) {
                strand_hit <- unique(as.vector(strand(hits[hit])))
                stats[(stats$Strand == strand_hit) & (stats$Gene_part == colname),3] <- stats[(stats$Strand == strand_hit) & (stats$Gene_part == colname),3] + 1
                m[gene, colname] <- m[gene, colname] + 1
              }
            }
          }
        }
      }
      
      # use the information in m2 to indicate in m with NA which are the genes that lack of one or more gene parts and which
      # they are 
      for (r in 1:nrow(m2)) {
        for (c in 1:ncol(m2)) {
          if (m2[r,c] == 0) {
            m[r,c] <- NA
            gene_type <- gene_biotype[rownames(m2)[r],]
            if ('pseudogene' %in% unlist(strsplit(gene_type, split='_'))) {
              na_annotation[3,c] <- na_annotation[3,c]+1
            } else {
              na_annotation[gene_type,c] <- na_annotation[gene_type,c]+1
            }
          }
        }
      }
      
      print(na_annotation)
      
      if ('chr' %in% unlist(strsplit(n, split='_'))) {
        fraction = c('Chromatin')
      } else if ('nucleo' %in% unlist(strsplit(n, split='_'))) {
        fraction = c('Nucleoplasm')
      } else {
        fraction = c('Cytoplasm')
      }
      
      # compute the Spearman correlation between the number of hits identified in each gene part and 
      # the number of reads mapping on the corresponding gene
      cor_read_counts <- unlist(lapply(1:ncol(m), function(x,y,i){
        round(cor(x[!is.na(y[,i])], y[,i][!is.na(y[,i])], method = 'spearman'),2)
      }, x = read_counts, y = m))
      cor_read_counts[is.na(cor_read_counts)] <- 0
      correlation_read_counts[fraction,] <- cor_read_counts
      
      # compute the Spearman correlation between the number of hits identified in each gene part and 
      # the length of the corresponding gene
      cor_gene_length <- unlist(lapply(1:ncol(m), function(x,y,i){
        round(cor(x[!is.na(y[,i])], y[,i][!is.na(y[,i])], method = 'spearman'),2)
      }, x = length_genes, y = m))
      cor_gene_length[is.na(cor_gene_length)] <- 0
      correlation_gene_length[fraction,] <- cor_gene_length
      
      # compute the Spearman correlation between the number of hits identified in each gene part and 
      # the strand of the corresponding gene
      cor_gene_strand <- unlist(lapply(1:ncol(m), function(x,y,i){
        round(cor(x[!is.na(y[,i])], y[,i][!is.na(y[,i])], method = 'spearman'),2)
      }, x = strand_genes_numeric, y = m))
      cor_gene_strand[is.na(cor_gene_strand)] <- 0
      correlation_strand[fraction,] <- cor_gene_strand
      
      # saturation of the number of hits mapping on each gene part
      for (col in 1:ncol(m)) {
          for (row in 1:nrow(m)) {
            if (! is.na(m[row,col])) {
              if (m[row,col] > 5) {
                m[row,col] = 5
              }
            }
          }
        }
      
      # evaluate if there is any kind of correlation between the length of the coding exons and the length of the 3'UTR
      # of the genes with 3 hits falling in the coding exons 
      genes_hits_coding <- rownames(m[m[,2]>3,])
      
      length_coding_exons <- unlist(lapply(protein_coding_genes_5UTR_3UTR_introns_exons_stop[genes_hits_coding], function(x) {
        sum(width(x[x$feature == 'coding exon']))
      }))
      length_3UTR <- unlist(lapply(protein_coding_genes_5UTR_3UTR_introns_exons_stop[genes_hits_coding], function(x) {
        sum(width(x[x$feature == '3UTR']))
      }))
      
      pdf(paste0(directory, '/', n, '_hist.pdf'), width = 13, height = 7)
      par(mfrow=c(1,2))
      hist(length_coding_exons, main = as.character(cor(length_coding_exons,length_3UTR, method = 'spearman')))
      hist(length_3UTR)
      dev.off()
      
      ann <- data.frame(log10_gene_length = log10(length_genes)[rownames(m)], gene_strand = strand_genes[rownames(m)], log10_read_counts = log10(read_counts)[rownames(m)])
      rownames(ann) <- rownames(m)
      
      # plot the matrix m as heatmap
      pdf(paste0(directory, '/', n, '.pdf'), height = 7, width = 6.5)
      pheatmap(m, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE, na_col = 'white',
                 annotation_row = ann, treeheight_row = 0, main = fraction, fontsize = 10, color = brewer.pal(6, 'Reds'), 
                 legend_breaks = 0:5, legend_labels = c(as.character(0:4), '≥5'), fontsize_row = 10, fontsize_col = 10)
        grid.text('Transcriptional units', x=0.7, y=0.52, rot=270)
        dev.off()
    }
    
    write.xlsx(x = data.frame(correlation_read_counts),file = paste0(directory, '/correlation_read_counts.xlsx'),col.names = TRUE, row.names = TRUE)
    write.xlsx(x = data.frame(correlation_gene_length),file=paste0(directory, '/correlation_gene_length.xlsx'),col.names = TRUE, row.names = TRUE)
    write.xlsx(x = data.frame(correlation_strand),file =paste0(directory, '/correlation_strand.xlsx'),col.names = TRUE, row.names = TRUE)
    
  }
  
  R_objects <- list.files(directory, pattern = 'Rda', full.names = TRUE)
  f2(R_objects)
}

# copy /gene_counts/ folder inside /only_chr/with_DRACH/ folder and /only_chr/without_DRACH/ 
# containing ELIGOS DRACH+ and DRACH- hits only in nascent chromatin associated RNAs
count_mods_nascent_only_chr(directory='/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/only_chr/with_DRACH/')
count_mods_nascent_only_chr(directory='/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/only_chr/without_DRACH/')
