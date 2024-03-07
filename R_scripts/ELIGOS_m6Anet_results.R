library('GenomicFeatures')
library('compEpiTools')
library('ggplot2')
library('GenomicAlignments')
library('pheatmap')
library('RColorBrewer')
library('seqinr')
library('xlsx')
library('bedr')

# path to GTF annotation file
gtf_file <- "/path/to/Homo_sapiens.GRCh38.104.gtf"
txdb <- makeTxDbFromGFF(gtf_file)
genes_txdb <- GenomicFeatures::genes(txdb)

# function to remove ELIGOS hits overlapping with the coordinates of the SNPs of SUM159 and k562 cell lines
rm_SNPs <- function(path_SNPs_SUM, path_SNPs_k562, hits) {
    # SNPs SUM159 on hg38
    bed_SUM <-read.table(path_SNPs_SUM)
    
    grange_bed_SUM <- GRanges(seqnames = bed_SUM$V1,
                            ranges = IRanges(start = bed_SUM$V2, end=bed_SUM$V2)) 
    grange_bed_SUM <- resize(grange_bed_SUM, 5, 'center')
  
    over_hits_eligos_bed_SUM <- suppressWarnings(findOverlaps(hits,grange_bed_SUM, type = 'any', ignore.strand=TRUE))
    hits_without_SNPs_SUM <- hits[-unique(queryHits(over_hits_eligos_bed_SUM))]
  
    # SNPs k562 on hg38
    vcf_k562<-read.vcf(path_SNPs_k562)
    bed_k562<-vcf2bed(vcf_k562, filename = NULL, header = FALSE, other = NULL, verbose = TRUE)

    grange_bed_k562 <- GRanges(seqnames = gsub('chr', '', bed_k562$V1),
                               ranges = IRanges(start = bed_k562$V2, end=bed_k562$V2))
    grange_bed_k562 <- resize(grange_bed_k562, 5, 'center')
    
    over_hits_eligos_bed_k562 <- suppressWarnings(findOverlaps(hits_without_SNPs_SUM,grange_bed_k562, type = 'any', ignore.strand=TRUE))
    hits_without_SNPs_SUM_k562 <- hits_without_SNPs_SUM[-unique(queryHits(over_hits_eligos_bed_k562))]
    
    return(hits_without_SNPs_SUM_k562)
}

# path_directory is the path to the directory containing three folders: /chr/, /nucleo/, /cyto/ each with the 5 txt files 
# produced by ELIGOS for each of the 5 samplings
# create a folder /hits_ELIGOS/ inside path_directory in which you will save the hits of each fraction 
# confirmed by 5 samplings
ELIGOS_results <- function(path_directory, 
                           name_pdf_overlap_5samplings_chr_ass, name_pdf_overlap_5samplings_nucleo, name_pdf_overlap_5samplings_cyto,
                           name_pdf_histogram_chr_ass, name_pdf_histogram_nucleo, name_pdf_histogram_cyto,
                           name_pdf_overlap_3fractions, p, ap, OR) {
  
  # load the txt files with ELIGOS output for all the 5 samplings of each fraction
  eligos_chr_ass <- list.files(path = paste0(path_directory,'chr'), pattern = "txt", full.names = TRUE)
  eligos_nucleo <- list.files(path = paste0(path_directory,'nucleo'), pattern = "txt", full.names = TRUE)
  eligos_cyto <- list.files(path = paste0(path_directory,'cyto'), pattern = "txt", full.names = TRUE)
  
  eligos_output <- function(dir, tol) {
    list_gr <- list()
    
    for (i in 1:length(dir)) {
      t <- read.table(dir[i], sep= '\t', header=TRUE)
      
      gr <- GRanges(seqnames = Rle(t$chrom),
                    ranges = IRanges(start = t$start_loc - tol, end=t$end_loc + tol),
                    strand = Rle(t$strand),
                    pval = t$pval,
                    pvalAdj = t$adjPval,
                    oddR = t$oddR,
                    rep = if (('chr' %in% strsplit(dir,split='/')[[1]]) & (!'nascent' %in% strsplit(dir,split='_')[[1]])) {
                      as.numeric(strsplit(strsplit(strsplit(dir[i], split = '\\.')[[1]][1], split='/')[[1]][8],split='_')[[1]][3])
                    } else if ((('nucleo' %in% strsplit(dir,split='/')[[1]]) | ('cyto' %in% strsplit(dir,split='/')[[1]])) & (!'nascent' %in% strsplit(dir,split='_')[[1]])) {
                      as.numeric(strsplit(strsplit(strsplit(dir[i], split = '\\.')[[1]][1], split='/')[[1]][8],split='_')[[1]][2])
                    } else if (('chr' %in% strsplit(dir,split='/')[[1]]) & ('nascent' %in% strsplit(dir,split='_')[[1]])) {
                      as.numeric(strsplit(strsplit(strsplit(dir[i], split = '\\.')[[1]][1], split='/')[[1]][8],split='_')[[1]][4])
                    } else if ((('nucleo' %in% strsplit(dir,split='/')[[1]]) | ('cyto' %in% strsplit(dir,split='/')[[1]])) & ('nascent' %in% strsplit(dir,split='_')[[1]])) {
                      as.numeric(strsplit(strsplit(strsplit(dir[i], split = '\\.')[[1]][1], split='/')[[1]][8],split='_')[[1]][3])
                    }
      )
      
      # filter the analysed sites on the base of the input parameters
      gr <- gr[(gr$pval <= p) & (gr$pvalAdj <= ap) & (gr$oddR >=OR)]
      
      list_gr <- c(list_gr, gr)
    }
    return(list_gr)
  }
  
  gr_eligos_chr_ass <- eligos_output(eligos_chr_ass, 2)
  # assign each hit to the gene on which it maps (this step is done since there are errors in ELIGOS results)
  for (i in 1:length(gr_eligos_chr_ass)) {
    gene_ids <- c()
    for (j in 1:length(gr_eligos_chr_ass[[i]])){
      over <- findOverlaps(gr_eligos_chr_ass[[i]][j],genes_txdb)
      gene_ids <- c(gene_ids, genes_txdb[subjectHits(over)[1]]$gene_id)
    }
    mcols(gr_eligos_chr_ass[[i]]) <- cbind(mcols(gr_eligos_chr_ass[[i]]), gene_id = gene_ids)
  }
  
  gr_eligos_nucleo <- eligos_output(eligos_nucleo, 2)
  for (i in 1:length(gr_eligos_nucleo)) {
    gene_ids <- c()
    for (j in 1:length(gr_eligos_nucleo[[i]])){
      over <- findOverlaps(gr_eligos_nucleo[[i]][j],genes_txdb)
      gene_ids <- c(gene_ids, genes_txdb[subjectHits(over)[1]]$gene_id)
    }
    mcols(gr_eligos_nucleo[[i]]) <- cbind(mcols(gr_eligos_nucleo[[i]]), gene_id = gene_ids)
  }
  
  gr_eligos_cyto <- eligos_output(eligos_cyto, 2)
  for (i in 1:length(gr_eligos_cyto)) {
    gene_ids <- c()
    for (j in 1:length(gr_eligos_cyto[[i]])){
      over <- findOverlaps(gr_eligos_cyto[[i]][j],genes_txdb)
      gene_ids <- c(gene_ids, genes_txdb[subjectHits(over)[1]]$gene_id)
    }
    mcols(gr_eligos_cyto[[i]]) <- cbind(mcols(gr_eligos_cyto[[i]]), gene_id = gene_ids)
  }
  
  # remove ELIGOS hits overlapping with the SNPs of SUM159 or K562 cells and create a heatmap that 
  # reports the overlap between the hits identified in one sampling vs the hits identified 
  # in each other sampling
  gr_eligos_chr_ass <- unlist(gr_eligos_chr_ass)
  gr_eligos_chr_ass_without_SNPs <- lapply(gr_eligos_chr_ass, function(x) rm_SNPs("/path/to/SNPs_SUM_hg38.bed","/path/to/IVT_k562_hg38.vcf",x))
  pdf(file = paste0(path_directory, 'chr/', name_pdf_overlap_5samplings_chr_ass), width = 5, height = 7)
  overlapOfGRanges(gr_eligos_chr_ass_without_SNPs,plot = TRUE)
  dev.off()
  
  gr_eligos_nucleo <- unlist(gr_eligos_nucleo)
  gr_eligos_nucleo_without_SNPs <- lapply(gr_eligos_nucleo, function(x) rm_SNPs("/path/to/SNPs_SUM_hg38.bed","/path/to/IVT_k562_hg38.vcf",x))
  pdf(file = paste0(path_directory, 'nucleo/', name_pdf_overlap_5samplings_nucleo), width = 5, height = 7)
  overlapOfGRanges(gr_eligos_nucleo_without_SNPs,plot = TRUE)
  dev.off()
  
  gr_eligos_cyto <- unlist(gr_eligos_cyto)
  gr_eligos_cyto_without_SNPs <- lapply(gr_eligos_cyto, function(x) rm_SNPs("/path/to/SNPs_SUM_hg38.bed","/path/to/IVT_k562_hg38.vcf",x))
  pdf(file = paste0(path_directory, 'cyto/', name_pdf_overlap_5samplings_cyto), width = 5, height = 7)
  overlapOfGRanges(gr_eligos_cyto_without_SNPs,plot = TRUE)
  dev.off()
  
  eligos_chr_ass_all_samplings_without_SNPs <- c(gr_eligos_chr_ass_without_SNPs[[1]],gr_eligos_chr_ass_without_SNPs[[2]],gr_eligos_chr_ass_without_SNPs[[3]],gr_eligos_chr_ass_without_SNPs[[4]],gr_eligos_chr_ass_without_SNPs[[5]])
  
  confirmed_by_chr <- data.frame(num_samplings=seq(1,5), num_hits=rep(0,5))
  confirmed_by_5_chr <- c()
  
  # identify the hits present in 5 samplings (any type of overlap)
  for (r in 1:length(eligos_chr_ass_all_samplings_without_SNPs)) {
    gr_r_rep <- queryHits(findOverlaps(eligos_chr_ass_all_samplings_without_SNPs, eligos_chr_ass_all_samplings_without_SNPs[r], type='any')) 
    confirmed_by_chr[length(unique(eligos_chr_ass_all_samplings_without_SNPs[gr_r_rep]$rep)),2] <- confirmed_by_chr[length(unique(eligos_chr_ass_all_samplings_without_SNPs[gr_r_rep]$rep)),2] +1
    if (length(unique(eligos_chr_ass_all_samplings_without_SNPs[gr_r_rep]$rep)) == 5) {
      confirmed_by_5_chr <- c(confirmed_by_5_chr, r)
    }
  }
  
  # the hits confirmed in 5 samplings are resized to the original coordinates returned by ELIGOS 
  hits_eligos_chr_ass_confirmed_5 <- resize(eligos_chr_ass_all_samplings_without_SNPs[confirmed_by_5_chr], 2, fix = 'center')
  hits_eligos_chr_ass_confirmed_5 <- reduce(hits_eligos_chr_ass_confirmed_5, ignore.strand = FALSE)
  print(table(width(hits_eligos_chr_ass_confirmed_5)))
  # resize the merged hits to a width at least of 10 nucleotides
  for (i in 1:length(hits_eligos_chr_ass_confirmed_5)) {
    if (width(hits_eligos_chr_ass_confirmed_5[i])<10) {
      hits_eligos_chr_ass_confirmed_5[i] <- resize(hits_eligos_chr_ass_confirmed_5[i], 10, fix='center')
    }
  }
  print(table(width(hits_eligos_chr_ass_confirmed_5)))
  
  # assign each hit to the gene on which it maps
  gene_ids <- c()
  for (i in 1:length(hits_eligos_chr_ass_confirmed_5)) {
    over <- findOverlaps(hits_eligos_chr_ass_confirmed_5[i],genes_txdb)
    gene_ids <- c(gene_ids, genes_txdb[subjectHits(over)[1]]$gene_id)
  }
  mcols(hits_eligos_chr_ass_confirmed_5) <- cbind(mcols(hits_eligos_chr_ass_confirmed_5), gene_id = gene_ids)
  
  eligos_nucleo_all_samplings_without_SNPs <- c(gr_eligos_nucleo_without_SNPs[[1]],gr_eligos_nucleo_without_SNPs[[2]],gr_eligos_nucleo_without_SNPs[[3]],gr_eligos_nucleo_without_SNPs[[4]],gr_eligos_nucleo_without_SNPs[[5]])
  
  confirmed_by_nucleo <- data.frame(num_samplings=seq(1,5), num_hits=rep(0,5))
  confirmed_by_5_nucleo <- c()
  
  for (r in 1:length(eligos_nucleo_all_samplings_without_SNPs)) {
    gr_r_rep <- queryHits(findOverlaps(eligos_nucleo_all_samplings_without_SNPs, eligos_nucleo_all_samplings_without_SNPs[r], type='any')) 
    confirmed_by_nucleo[length(unique(eligos_nucleo_all_samplings_without_SNPs[gr_r_rep]$rep)),2] <- confirmed_by_nucleo[length(unique(eligos_nucleo_all_samplings_without_SNPs[gr_r_rep]$rep)),2] +1
    if (length(unique(eligos_nucleo_all_samplings_without_SNPs[gr_r_rep]$rep)) == 5) {
      confirmed_by_5_nucleo <- c(confirmed_by_5_nucleo, r)
    }
  }
  
  hits_eligos_nucleo_confirmed_5 <- resize(eligos_nucleo_all_samplings_without_SNPs[confirmed_by_5_nucleo], 2, fix = 'center')
  hits_eligos_nucleo_confirmed_5 <- reduce(hits_eligos_nucleo_confirmed_5, ignore.strand = FALSE)
  print(table(width(hits_eligos_nucleo_confirmed_5)))
  for (i in 1:length(hits_eligos_nucleo_confirmed_5)) {
    if (width(hits_eligos_nucleo_confirmed_5[i])<10) {
      hits_eligos_nucleo_confirmed_5[i] <- resize(hits_eligos_nucleo_confirmed_5[i], 10, fix='center')
    }
  }
  print(table(width(hits_eligos_nucleo_confirmed_5)))
  
  gene_ids <- c()
  for (i in 1:length(hits_eligos_nucleo_confirmed_5)) {
    over <- findOverlaps(hits_eligos_nucleo_confirmed_5[i],genes_txdb)
    gene_ids <- c(gene_ids, genes_txdb[subjectHits(over)[1]]$gene_id)
  }
  mcols(hits_eligos_nucleo_confirmed_5) <- cbind(mcols(hits_eligos_nucleo_confirmed_5), gene_id = gene_ids)
  
  eligos_cyto_all_samplings_without_SNPs <- c(gr_eligos_cyto_without_SNPs[[1]],gr_eligos_cyto_without_SNPs[[2]],gr_eligos_cyto_without_SNPs[[3]],gr_eligos_cyto_without_SNPs[[4]],gr_eligos_cyto_without_SNPs[[5]])
  
  confirmed_by_cyto <- data.frame(num_samplings=seq(1,5), num_hits=rep(0,5))
  confirmed_by_5_cyto <- c()
  
  for (r in 1:length(eligos_cyto_all_samplings_without_SNPs)) {
    gr_r_rep <- queryHits(findOverlaps(eligos_cyto_all_samplings_without_SNPs, eligos_cyto_all_samplings_without_SNPs[r], type='any')) 
    confirmed_by_cyto[length(unique(eligos_cyto_all_samplings_without_SNPs[gr_r_rep]$rep)),2] <- confirmed_by_cyto[length(unique(eligos_cyto_all_samplings_without_SNPs[gr_r_rep]$rep)),2] +1
    if (length(unique(eligos_cyto_all_samplings_without_SNPs[gr_r_rep]$rep)) == 5) {
      confirmed_by_5_cyto <- c(confirmed_by_5_cyto, r)
    }
  }
  
  hits_eligos_cyto_confirmed_5 <- resize(eligos_cyto_all_samplings_without_SNPs[confirmed_by_5_cyto], 2, fix = 'center')
  hits_eligos_cyto_confirmed_5 <- reduce(hits_eligos_cyto_confirmed_5, ignore.strand = FALSE)
  print(table(width(hits_eligos_cyto_confirmed_5)))
  for (i in 1:length(hits_eligos_cyto_confirmed_5)) {
    if (width(hits_eligos_cyto_confirmed_5[i])<10) {
      hits_eligos_cyto_confirmed_5[i] <- resize(hits_eligos_cyto_confirmed_5[i], 10, fix='center')
    }
  }
  print(table(width(hits_eligos_cyto_confirmed_5)))
  
  gene_ids <- c()
  for (i in 1:length(hits_eligos_cyto_confirmed_5)) {
    over <- findOverlaps(hits_eligos_cyto_confirmed_5[i],genes_txdb)
    gene_ids <- c(gene_ids, genes_txdb[subjectHits(over)[1]]$gene_id)
  }
  mcols(hits_eligos_cyto_confirmed_5) <- cbind(mcols(hits_eligos_cyto_confirmed_5), gene_id = gene_ids)
  
  # barplot reporting the number of hits confirmed by n samplings
  p <- ggplot(confirmed_by_chr, aes(x=num_samplings, y = num_hits)) +
    geom_bar(stat="identity") +
    theme_classic(base_size = 16) +
    scale_x_continuous(breaks=seq(1,5)) +
    labs(title='Chromatin') +
    xlab('Num samplings with hit')+
    ylab('Number of hits')
  
  ggsave(paste0(path_directory, 'chr/', name_pdf_histogram_chr_ass), plot = p, height = 7, width = 7)
  
  p <- ggplot(confirmed_by_nucleo, aes(x=num_samplings, y = num_hits)) +
    geom_bar(stat="identity") +
    theme_classic(base_size = 16) +
    scale_x_continuous(breaks=seq(1,5)) +
    labs(title='Nucleoplasm') +
    xlab('Num samplings with hit')+
    ylab('Number of hits')
  
  ggsave(paste0(path_directory, 'nucleo/', name_pdf_histogram_nucleo), plot = p, height = 7, width = 7)
  
  p <- ggplot(confirmed_by_cyto, aes(x=num_samplings, y = num_hits)) +
    geom_bar(stat="identity") +
    theme_classic(base_size = 16) +
    scale_x_continuous(breaks=seq(1,5)) +
    labs(title='Cytoplasm') +
    xlab('Num samplings with hit')+
    ylab('Number of hits')
  
  ggsave(paste0(path_directory, 'cyto/', name_pdf_histogram_cyto), plot = p, height = 7, width = 7)
  
  mods <- matrix(0,ncol=3,nrow=3)
  colnames(mods) <- c('chr_ass','nucleo','cyto')
  rownames(mods) <- c('# hits >= 1 samplings','# hits in 5 samplings', '# hits in 5 samplings after reduce')
  
  mods[1,1] <- length(eligos_chr_ass_all_samplings)
  mods[2,1] <- confirmed_by_chr[5,2]
  mods[3,1] <- length(hits_eligos_chr_ass_confirmed_5)
  mods[1,2] <- length(eligos_nucleo_all_samplings)
  mods[2,2] <- confirmed_by_nucleo[5,2]
  mods[3,2] <- length(hits_eligos_nucleo_confirmed_5)
  mods[1,3] <- length(eligos_cyto_all_samplings)
  mods[2,3] <- confirmed_by_cyto[5,2]
  mods[3,3] <- length(hits_eligos_cyto_confirmed_5)
  
  write.xlsx(x = data.frame(mods),file = paste0(path_directory, 'number_hit_5_samplings.xlsx'),col.names = TRUE, row.names = TRUE)
  
  save(hits_eligos_chr_ass_confirmed_5, file=paste0(path_directory,'/hits_ELIGOS/hits_eligos_chr_ass_confirmed_5.Rda'))
  save(hits_eligos_nucleo_confirmed_5, file=paste0(path_directory,'/hits_ELIGOS/hits_eligos_nucleo_confirmed_5.Rda'))
  save(hits_eligos_cyto_confirmed_5, file=paste0(path_directory,'/hits_ELIGOS/hits_eligos_cyto_confirmed_5.Rda'))
  
  # heatmap with the pairwise overlap between ELIGOS hits in the different fractions
  overlap_hits_eligos_confirmed_5 <- list(hits_eligos_chr_ass_confirmed_5,hits_eligos_nucleo_confirmed_5,hits_eligos_cyto_confirmed_5)
  names(overlap_hits_eligos_confirmed_5) <- c('Chromatin', 'Nucleoplasm','Cytoplasm')
  
  pdf(file = paste0(path_directory, '/', name_pdf_overlap_3fractions), width = 7, height = 5)
  overlapOfGRanges(overlap_hits_eligos_confirmed_5,plot = TRUE)
  dev.off()
  
  return(overlap_hits_eligos_confirmed_5)
}

# path_directory is the path to the directory containing three folders: /chr/, /nucleo/, /cyto/ each with the 5 txt files 
# produced by ELIGOS for each of the 5 samplings
# create a folder /without_DRACH/ inside path_directory in which you will save the hits of each fraction 
# confirmed by 5 samplings that are DRACH-
DRACH_overlap_ELIGOS <- function(path_directory, hits_ELIGOS_chr, hits_ELIGOS_nucleo, hits_ELIGOS_cyto) {
  
  summary_table <- matrix(ncol=3, nrow=4)
  colnames(summary_table) <- c('Total number of hits', 'DRACH+ hits', 'DRACH- hits')
  rownames(summary_table) <- c('Chromatin', 'Nucleoplasm', 'Cytoplasm', 'Overlap')
  
  summary_table[1,1] <- as.character(length(hits_ELIGOS_chr))
  summary_table[2,1] <- as.character(length(hits_ELIGOS_nucleo))
  summary_table[3,1] <- as.character(length(hits_ELIGOS_cyto))
  
  confirmed_hits <- list(hits_ELIGOS_chr, hits_ELIGOS_nucleo, hits_ELIGOS_cyto)
  order <- order(c(length(hits_ELIGOS_chr), length(hits_ELIGOS_nucleo), length(hits_ELIGOS_cyto)))
  
  # compute the number of ELIGOS hits present in all the fractions
  overlap_eligos <- findOverlaps(confirmed_hits[[order[1]]],confirmed_hits[[order[2]]], type = 'any')
  overlap_eligos_3_sets <- confirmed_hits[[order[1]]][unique(queryHits(overlap_eligos))]
  overlap_eligos_chr_cyto_nucleo <- findOverlaps(overlap_eligos_3_sets,confirmed_hits[[order[3]]], type = 'any')
  hits_eligos_confirmed_by_all <- overlap_eligos_3_sets[unique(queryHits(overlap_eligos_chr_cyto_nucleo))]
  summary_table[4,1] <- as.character(length(hits_eligos_confirmed_by_all))
  
  load('/path/to/R_data/DRACH_forward_strand.Rda')
  # identify which and how many hits are DRACH+/DRACH-
  overlap_DRACH_chr <- findOverlaps(hits_ELIGOS_chr, DRACH, type='any', ignore.strand=FALSE)
  hits_eligos_chr_ass_confirmed_5_with_DRACH <- hits_ELIGOS_chr[unique(queryHits(overlap_DRACH_chr))]
  save(hits_eligos_chr_ass_confirmed_5_with_DRACH, file=paste0(path_directory,'/hits_ELIGOS/hits_eligos_chr_ass_confirmed_5_with_DRACH.Rda'))
  summary_table[1,2] <- paste0(as.character(length(hits_eligos_chr_ass_confirmed_5_with_DRACH)), ' - ', as.character(round(length(hits_eligos_chr_ass_confirmed_5_with_DRACH)*100/length(hits_ELIGOS_chr), 2)),'%')
  
  hits_eligos_chr_ass_confirmed_5_without_DRACH <- hits_ELIGOS_chr[-unique(queryHits(overlap_DRACH_chr))]
  save(hits_eligos_chr_ass_confirmed_5_without_DRACH, file=paste0(path_directory,'/without_DRACH/hits_eligos_chr_ass_confirmed_5_without_DRACH.Rda'))
  summary_table[1,3] <- as.character(length(hits_eligos_chr_ass_confirmed_5_without_DRACH))
  
  overlap_DRACH_nucleo <- findOverlaps(hits_ELIGOS_nucleo, DRACH, type='any', ignore.strand=FALSE)
  hits_eligos_nucleo_confirmed_5_with_DRACH <- hits_ELIGOS_nucleo[unique(queryHits(overlap_DRACH_nucleo))]
  save(hits_eligos_nucleo_confirmed_5_with_DRACH, file=paste0(path_directory,'/hits_ELIGOS/hits_eligos_nucleo_confirmed_5_with_DRACH.Rda'))
  summary_table[2,2] <- paste0(as.character(length(hits_eligos_nucleo_confirmed_5_with_DRACH)), ' - ', as.character(round(length(hits_eligos_nucleo_confirmed_5_with_DRACH)*100/length(hits_ELIGOS_nucleo), 2)),'%')
  
  hits_eligos_nucleo_confirmed_5_without_DRACH <- hits_ELIGOS_nucleo[-unique(queryHits(overlap_DRACH_nucleo))]
  save(hits_eligos_nucleo_confirmed_5_without_DRACH, file=paste0(path_directory,'/without_DRACH/hits_eligos_nucleo_confirmed_5_without_DRACH.Rda'))
  summary_table[2,3] <- as.character(length(hits_eligos_nucleo_confirmed_5_without_DRACH))
  
  overlap_DRACH_cyto <- findOverlaps(hits_ELIGOS_cyto, DRACH, type='any', ignore.strand=FALSE)
  hits_eligos_cyto_confirmed_5_with_DRACH <- hits_ELIGOS_cyto[unique(queryHits(overlap_DRACH_cyto))]
  save(hits_eligos_cyto_confirmed_5_with_DRACH, file=paste0(path_directory,'/hits_ELIGOS/hits_eligos_cyto_confirmed_5_with_DRACH.Rda'))
  summary_table[3,2] <- paste0(as.character(length(hits_eligos_cyto_confirmed_5_with_DRACH)), ' - ', as.character(round(length(hits_eligos_cyto_confirmed_5_with_DRACH)*100/length(hits_ELIGOS_cyto), 2)),'%')
  
  hits_eligos_cyto_confirmed_5_without_DRACH <- hits_ELIGOS_cyto[-unique(queryHits(overlap_DRACH_cyto))]
  save(hits_eligos_cyto_confirmed_5_without_DRACH, file=paste0(path_directory,'/without_DRACH/hits_eligos_cyto_confirmed_5_without_DRACH.Rda'))
  summary_table[3,3] <- as.character(length(hits_eligos_cyto_confirmed_5_without_DRACH))
  
  # merge ELIGOS DRACH- hits from the three fractions and resize the corresponding regions to 20 nucleotides
  union_hits_without_DRACH <- c(hits_eligos_chr_ass_confirmed_5_without_DRACH,hits_eligos_nucleo_confirmed_5_without_DRACH,hits_eligos_cyto_confirmed_5_without_DRACH)
  union_hits_without_DRACH <- reduce(union_hits_without_DRACH, ignore.strand = FALSE)
  print(table(width(union_hits_without_DRACH)))
  union_hits_without_DRACH <- resize(union_hits_without_DRACH,20,fix = 'center')
  save(union_hits_without_DRACH, file=paste0(path_directory,'/without_DRACH/union_hits_without_DRACH.Rda'))
  
  # save the fasta sequences corresponding to each of the 20-nucleotides DRACH- regions  
  chrs <- readDNAStringSet('/path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa', format="fasta")
  names(chrs) <- gsub(x = names(chrs), pattern = " .*", replacement = "")
  for (j in 1:length(union_hits_without_DRACH)) {
    write.fasta(file.out = paste0(path_directory, '/ELIGOS_without_DRACH.fa'), sequences = as.character(getSeq(chrs, union_hits_without_DRACH[j])), names = as.character(j),open = 'a')
  }
  
  # compute the number of ELIGOS DRACH+ hits present in all the fractions
  confirmed_hits_DRACH <- list(hits_eligos_chr_ass_confirmed_5_with_DRACH, hits_eligos_nucleo_confirmed_5_with_DRACH, hits_eligos_cyto_confirmed_5_with_DRACH)
  order <- order(c(length(hits_eligos_chr_ass_confirmed_5_with_DRACH), length(hits_eligos_nucleo_confirmed_5_with_DRACH), length(hits_eligos_cyto_confirmed_5_with_DRACH)))
  
  overlap_eligos_DRACH <- findOverlaps(confirmed_hits_DRACH[[order[1]]],confirmed_hits_DRACH[[order[2]]], type = 'any')
  overlap_eligos_3_sets_DRACH <- confirmed_hits_DRACH[[order[1]]][unique(queryHits(overlap_eligos_DRACH))]
  overlap_eligos_chr_cyto_nucleo_DRACH <- findOverlaps(overlap_eligos_3_sets_DRACH,confirmed_hits_DRACH[[order[3]]], type = 'any')
  hits_eligos_confirmed_by_all_DRACH <- overlap_eligos_3_sets_DRACH[unique(queryHits(overlap_eligos_chr_cyto_nucleo_DRACH))]
  summary_table[4,2] <- as.character(length(hits_eligos_confirmed_by_all_DRACH))
  
  # compute the number of ELIGOS DRACH- hits present in all the fractions
  confirmed_hits_without_DRACH <- list(hits_eligos_chr_ass_confirmed_5_without_DRACH, hits_eligos_nucleo_confirmed_5_without_DRACH, hits_eligos_cyto_confirmed_5_without_DRACH)
  order <- order(c(length(hits_eligos_chr_ass_confirmed_5_without_DRACH), length(hits_eligos_nucleo_confirmed_5_without_DRACH), length(hits_eligos_cyto_confirmed_5_without_DRACH)))
  
  overlap_eligos_without_DRACH <- findOverlaps(confirmed_hits_without_DRACH[[order[1]]],confirmed_hits_without_DRACH[[order[2]]], type = 'any')
  overlap_eligos_3_sets_without_DRACH <- confirmed_hits_without_DRACH[[order[1]]][unique(queryHits(overlap_eligos_without_DRACH))]
  overlap_eligos_chr_cyto_nucleo_without_DRACH <- findOverlaps(overlap_eligos_3_sets_without_DRACH,confirmed_hits_without_DRACH[[order[3]]], type = 'any')
  hits_eligos_confirmed_by_all_without_DRACH <- overlap_eligos_3_sets_without_DRACH[unique(queryHits(overlap_eligos_chr_cyto_nucleo_without_DRACH))]
  summary_table[4,3] <- as.character(length(hits_eligos_confirmed_by_all_without_DRACH))
  
  write.xlsx(x = data.frame(summary_table),file = paste0(path_directory, 'summary_table.xlsx'),col.names = TRUE, row.names = TRUE)
  
  # heatmap representing the spatial overlap between DRACH+ hits of the three fractions
  union_hits_eligos <- c(hits_eligos_chr_ass_confirmed_5_with_DRACH, hits_eligos_nucleo_confirmed_5_with_DRACH, hits_eligos_cyto_confirmed_5_with_DRACH)
  union_hits_eligos <- reduce(union_hits_eligos, ignore.strand = FALSE)
  
  heatmap_matrix_eligos <- matrix(NA,ncol = 3)
  for (i in 1:length(union_hits_eligos)) { 
    if ((length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_chr_ass_confirmed_5_with_DRACH, type = 'any')))!=0) &
        (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_nucleo_confirmed_5_with_DRACH, type = 'any')))!=0) &
        (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_cyto_confirmed_5_with_DRACH, type = 'any')))!=0))  {
      heatmap_matrix_eligos <- rbind(heatmap_matrix_eligos, c(1,1,1))
    } else if ((length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_chr_ass_confirmed_5_with_DRACH, type = 'any')))!=0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_nucleo_confirmed_5_with_DRACH, type = 'any')))!=0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_cyto_confirmed_5_with_DRACH, type = 'any')))==0)) {
      heatmap_matrix_eligos <- rbind(heatmap_matrix_eligos, c(1,1,0))
    } else if ((length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_chr_ass_confirmed_5_with_DRACH, type = 'any')))!=0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_nucleo_confirmed_5_with_DRACH, type = 'any')))==0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_cyto_confirmed_5_with_DRACH, type = 'any')))!=0)) {
      heatmap_matrix_eligos <- rbind(heatmap_matrix_eligos, c(1,0,1))
    } else if ((length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_chr_ass_confirmed_5_with_DRACH, type = 'any')))!=0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_nucleo_confirmed_5_with_DRACH, type = 'any')))==0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_cyto_confirmed_5_with_DRACH, type = 'any')))==0)) {
      heatmap_matrix_eligos <- rbind(heatmap_matrix_eligos, c(1,0,0))
    } else if ((length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_chr_ass_confirmed_5_with_DRACH, type = 'any')))==0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_nucleo_confirmed_5_with_DRACH, type = 'any')))!=0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_cyto_confirmed_5_with_DRACH, type = 'any')))!=0)) {
      heatmap_matrix_eligos <- rbind(heatmap_matrix_eligos, c(0,1,1))
    } else if ((length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_chr_ass_confirmed_5_with_DRACH, type = 'any')))==0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_nucleo_confirmed_5_with_DRACH, type = 'any')))!=0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_cyto_confirmed_5_with_DRACH, type = 'any')))==0)) {
      heatmap_matrix_eligos <- rbind(heatmap_matrix_eligos, c(0,1,0))
    } else if ((length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_chr_ass_confirmed_5_with_DRACH, type = 'any')))==0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_nucleo_confirmed_5_with_DRACH, type = 'any')))==0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_cyto_confirmed_5_with_DRACH, type = 'any')))!=0)) {
      heatmap_matrix_eligos <- rbind(heatmap_matrix_eligos, c(0,0,1))
    }
  }
  
  heatmap_matrix_eligos <- heatmap_matrix_eligos[-1,]
  colnames(heatmap_matrix_eligos) <- c('Chromatin','Nucleoplasm','Cytoplasm')
  rownames(heatmap_matrix_eligos) <- 1:nrow(heatmap_matrix_eligos)
  order <- list(c(1,0,0), c(0,1,0), c(0,0,1),c(1,1,0), c(0,1,1), c(1,0,1),c(1,1,1))
  
  order2 <- as.numeric(unlist(lapply(order, function(x) {
    rownames(heatmap_matrix_eligos[heatmap_matrix_eligos[,1]==x[1] & heatmap_matrix_eligos[,2]==x[2] & heatmap_matrix_eligos[,3]==x[3],])
  })))
  
  number <- unlist(lapply(order, function(x) {
    paste0(as.character(nrow(heatmap_matrix_eligos[heatmap_matrix_eligos[,1]==x[1] & heatmap_matrix_eligos[,2]==x[2] & heatmap_matrix_eligos[,3]==x[3],])),' - ', as.character(round(nrow(heatmap_matrix_eligos[heatmap_matrix_eligos[,1]==x[1] & heatmap_matrix_eligos[,2]==x[2] & heatmap_matrix_eligos[,3]==x[3],])*100/nrow(heatmap_matrix_eligos),2)),'%')
  }))
  
  ann_col <- data.frame(
    Fraction = factor(rep(number, c(unlist(lapply(order, function(x) {
      nrow(heatmap_matrix_eligos[heatmap_matrix_eligos[,1]==x[1] & heatmap_matrix_eligos[,2]==x[2] & heatmap_matrix_eligos[,3]==x[3],])
    })))), levels = number))
  
  rownames(ann_col) <- 1:nrow(heatmap_matrix_eligos)
  
  colors <-  brewer.pal(n = 7, name = "Dark2")
  names(colors) <- number
  mycolors <- list(Fraction = colors)
  
  heatmap_matrix_eligos2 <- heatmap_matrix_eligos[order2,]
  rownames(heatmap_matrix_eligos2) <- 1:nrow(heatmap_matrix_eligos2)
  
  pdf(paste0(path_directory,'/heatmap_DRACH_positive_ELIGOS.pdf'))
  pheatmap(angle_col = 315,heatmap_matrix_eligos2,cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = FALSE, 
           border_color = 'black', color = c(brewer.pal(8, 'Reds')[2],brewer.pal(8, 'Reds')[8]),legend_breaks = 0:1, 
           legend_labels = c('0','1'),annotation_row = ann_col,annotation_colors = mycolors)
  dev.off()
  
  # heatmap representing the spatial overlap between DRACH- hits of the three fractions
  union_hits_eligos <- c(hits_eligos_chr_ass_confirmed_5_without_DRACH, hits_eligos_nucleo_confirmed_5_without_DRACH, hits_eligos_cyto_confirmed_5_without_DRACH)
  union_hits_eligos <- reduce(union_hits_eligos, ignore.strand = FALSE)
  
  heatmap_matrix_eligos <- matrix(NA,ncol = 3)
  for (i in 1:length(union_hits_eligos)) { 
    if ((length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_chr_ass_confirmed_5_without_DRACH, type = 'any')))!=0) &
        (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_nucleo_confirmed_5_without_DRACH, type = 'any')))!=0) &
        (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_cyto_confirmed_5_without_DRACH, type = 'any')))!=0))  {
      heatmap_matrix_eligos <- rbind(heatmap_matrix_eligos, c(1,1,1))
    } else if ((length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_chr_ass_confirmed_5_without_DRACH, type = 'any')))!=0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_nucleo_confirmed_5_without_DRACH, type = 'any')))!=0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_cyto_confirmed_5_without_DRACH, type = 'any')))==0)) {
      heatmap_matrix_eligos <- rbind(heatmap_matrix_eligos, c(1,1,0))
    } else if ((length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_chr_ass_confirmed_5_without_DRACH, type = 'any')))!=0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_nucleo_confirmed_5_without_DRACH, type = 'any')))==0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_cyto_confirmed_5_without_DRACH, type = 'any')))!=0)) {
      heatmap_matrix_eligos <- rbind(heatmap_matrix_eligos, c(1,0,1))
    } else if ((length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_chr_ass_confirmed_5_without_DRACH, type = 'any')))!=0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_nucleo_confirmed_5_without_DRACH, type = 'any')))==0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_cyto_confirmed_5_without_DRACH, type = 'any')))==0)) {
      heatmap_matrix_eligos <- rbind(heatmap_matrix_eligos, c(1,0,0))
    } else if ((length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_chr_ass_confirmed_5_without_DRACH, type = 'any')))==0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_nucleo_confirmed_5_without_DRACH, type = 'any')))!=0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_cyto_confirmed_5_without_DRACH, type = 'any')))!=0)) {
      heatmap_matrix_eligos <- rbind(heatmap_matrix_eligos, c(0,1,1))
    } else if ((length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_chr_ass_confirmed_5_without_DRACH, type = 'any')))==0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_nucleo_confirmed_5_without_DRACH, type = 'any')))!=0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_cyto_confirmed_5_without_DRACH, type = 'any')))==0)) {
      heatmap_matrix_eligos <- rbind(heatmap_matrix_eligos, c(0,1,0))
    } else if ((length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_chr_ass_confirmed_5_without_DRACH, type = 'any')))==0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_nucleo_confirmed_5_without_DRACH, type = 'any')))==0) &
               (length(queryHits(findOverlaps(union_hits_eligos[i],hits_eligos_cyto_confirmed_5_without_DRACH, type = 'any')))!=0)) {
      heatmap_matrix_eligos <- rbind(heatmap_matrix_eligos, c(0,0,1))
    }
  }
  
  heatmap_matrix_eligos <- heatmap_matrix_eligos[-1,]
  colnames(heatmap_matrix_eligos) <- c('Chromatin','Nucleoplasm','Cytoplasm')
  rownames(heatmap_matrix_eligos) <- 1:nrow(heatmap_matrix_eligos)
  order <- list(c(1,0,0), c(0,1,0), c(0,0,1),c(1,1,0), c(0,1,1), c(1,0,1),c(1,1,1))
  
  order2 <- as.numeric(unlist(lapply(order, function(x) {
    rownames(heatmap_matrix_eligos[heatmap_matrix_eligos[,1]==x[1] & heatmap_matrix_eligos[,2]==x[2] & heatmap_matrix_eligos[,3]==x[3],])
  })))
  
  number <- unlist(lapply(order, function(x) {
    paste0(as.character(nrow(heatmap_matrix_eligos[heatmap_matrix_eligos[,1]==x[1] & heatmap_matrix_eligos[,2]==x[2] & heatmap_matrix_eligos[,3]==x[3],])),' - ', as.character(round(nrow(heatmap_matrix_eligos[heatmap_matrix_eligos[,1]==x[1] & heatmap_matrix_eligos[,2]==x[2] & heatmap_matrix_eligos[,3]==x[3],])*100/nrow(heatmap_matrix_eligos),2)),'%')
  }))
  
  ann_col <- data.frame(
    Fraction = factor(rep(number, c(unlist(lapply(order, function(x) {
      nrow(heatmap_matrix_eligos[heatmap_matrix_eligos[,1]==x[1] & heatmap_matrix_eligos[,2]==x[2] & heatmap_matrix_eligos[,3]==x[3],])
    })))), levels = number))
  
  rownames(ann_col) <- 1:nrow(heatmap_matrix_eligos)
  
  colors <-  brewer.pal(n = 7, name = "Dark2")
  names(colors) <- number
  mycolors <- list(Fraction = colors)
  
  heatmap_matrix_eligos2 <- heatmap_matrix_eligos[order2,]
  rownames(heatmap_matrix_eligos2) <- 1:nrow(heatmap_matrix_eligos2)
  
  pdf(paste0(path_directory,'/heatmap_DRACH_negative_ELIGOS.pdf'))
  pheatmap(angle_col = 315,heatmap_matrix_eligos2,cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = FALSE, 
           border_color = 'black', color = c(brewer.pal(8, 'Reds')[2],brewer.pal(8, 'Reds')[8]),legend_breaks = 0:1, 
           legend_labels = c('0','1'),annotation_row = ann_col,annotation_colors = mycolors)
  dev.off()
  
  return(confirmed_hits_DRACH)
}

# path_directory is the path to the directory containing three folders: /chr/, /nucleo/, /cyto/ each with the 5 tsv files 
# produced by m6Anet for each of the 5 samplings
m6Anet_results <- function(path_directory, gr_m6Anet_chr_ass, gr_m6Anet_nucleo, gr_m6Anet_cyto,
                           name_pdf_overlap_10samplings_chr_ass, name_pdf_overlap_10samplings_nucleo, name_pdf_overlap_10samplings_cyto,
                           name_pdf_histogram_chr_ass, name_pdf_histogram_nucleo, name_pdf_histogram_cyto,
                           name_pdf_overlap_3fractions) {
  
  # all the hits that map on the same chromosome, that have the same strand, same start and end positions and same gene id but 
  # different transcript id are kept only once
  gr_unique <- function(x) {
    ind_unique <- which(isUnique(paste0(seqnames(x), "_", start(x), "_", end(x), "_", strand(x),"_", mcols(x)$gene_id)))
    dup <- unique(x[-ind_unique])
    x_unique <- sort(c(x[ind_unique], dup))
    return(x_unique)
  }
  
  gr_m6anet_chr_ass <- lapply(gr_m6Anet_chr_ass, gr_unique)
  gr_m6anet_nucleo <- lapply(gr_m6Anet_nucleo, gr_unique)
  gr_m6anet_cyto <- lapply(gr_m6Anet_cyto, gr_unique)
  
  # create a heatmap that reports the overlap between the hits identified in one sampling vs the hits identified 
  # in each other sampling
  gr_m6anet_chr_ass <- unlist(gr_m6anet_chr_ass)
  pdf(file = paste0(path_directory, '/chr/', name_pdf_overlap_10samplings_chr_ass), width = 5, height = 7)
  overlapOfGRanges(gr_m6anet_chr_ass,plot = TRUE)
  dev.off()
  
  gr_m6anet_nucleo <- unlist(gr_m6anet_nucleo)
  pdf(file = paste0(path_directory, '/nucleo/', name_pdf_overlap_10samplings_nucleo), width = 5, height = 7)
  overlapOfGRanges(gr_m6anet_nucleo,plot = TRUE)
  dev.off()
  
  gr_m6anet_cyto <- unlist(gr_m6anet_cyto)
  pdf(file = paste0(path_directory, '/cyto/', name_pdf_overlap_10samplings_cyto), width = 5, height = 7)
  overlapOfGRanges(gr_m6anet_cyto,plot = TRUE)
  dev.off()
  
  m6anet_chr_ass_all_samplings <- c(gr_m6anet_chr_ass[[1]],gr_m6anet_chr_ass[[2]],gr_m6anet_chr_ass[[3]],gr_m6anet_chr_ass[[4]],gr_m6anet_chr_ass[[5]])
  
  confirmed_by_chr <- data.frame(num_samplings=seq(1,5), num_hits=rep(0,5))
  confirmed_by_5_chr <- c()
  
  # identify the hits present in at least 5 samplings (any type of overlap)
  for (r in 1:length(m6anet_chr_ass_all_samplings)) {
    gr_r_rep <- queryHits(findOverlaps(m6anet_chr_ass_all_samplings, m6anet_chr_ass_all_samplings[r], type='any')) 
    confirmed_by_chr[length(unique(m6anet_chr_ass_all_samplings[gr_r_rep]$rep)),2] <- confirmed_by_chr[length(unique(m6anet_chr_ass_all_samplings[gr_r_rep]$rep)),2] +1
    if (length(unique(m6anet_chr_ass_all_samplings[gr_r_rep]$rep)) == 5) {
      confirmed_by_5_chr <- c(confirmed_by_5_chr, r)
    }
  }
  
  # the hits confirmed in 5 samplings are resized to the original coordinates returned by m6Anet 
  hits_m6anet_chr_ass_confirmed_5 <- resize(m6anet_chr_ass_all_samplings[confirmed_by_5_chr], 2, fix = 'center')
  hits_m6anet_chr_ass_confirmed_5 <- reduce(hits_m6anet_chr_ass_confirmed_5, ignore.strand = FALSE)
  print(table(width(hits_m6anet_chr_ass_confirmed_5)))
  # resize the merged hits to have a width at least of 10 nucleotides
  for (i in 1:length(hits_m6anet_chr_ass_confirmed_5)) {
    if (width(hits_m6anet_chr_ass_confirmed_5[i])<10) {
      hits_m6anet_chr_ass_confirmed_5[i] <- resize(hits_m6anet_chr_ass_confirmed_5[i], 10, fix='center')
    }
  }
  print(table(width(hits_m6anet_chr_ass_confirmed_5)))
  
  # assign each hit to the gene on which it maps
  gene_ids <- c()
  for (i in 1:length(hits_m6anet_chr_ass_confirmed_5)) {
    over <- findOverlaps(hits_m6anet_chr_ass_confirmed_5[i],genes_txdb)
    gene_ids <- c(gene_ids, genes_txdb[subjectHits(over)[1]]$gene_id)
  }
  mcols(hits_m6anet_chr_ass_confirmed_5) <- cbind(mcols(hits_m6anet_chr_ass_confirmed_5), gene_id = gene_ids)
  
  m6anet_nucleo_all_samplings <- c(gr_m6anet_nucleo[[1]],gr_m6anet_nucleo[[2]],gr_m6anet_nucleo[[3]],gr_m6anet_nucleo[[4]],gr_m6anet_nucleo[[5]])
  
  confirmed_by_nucleo <- data.frame(num_samplings=seq(1,5), num_hits=rep(0,5))
  confirmed_by_5_nucleo <- c()
  
  for (r in 1:length(m6anet_nucleo_all_samplings)) {
    gr_r_rep <- queryHits(findOverlaps(m6anet_nucleo_all_samplings, m6anet_nucleo_all_samplings[r], type='any')) 
    confirmed_by_nucleo[length(unique(m6anet_nucleo_all_samplings[gr_r_rep]$rep)),2] <- confirmed_by_nucleo[length(unique(m6anet_nucleo_all_samplings[gr_r_rep]$rep)),2] +1
    if (length(unique(m6anet_nucleo_all_samplings[gr_r_rep]$rep)) == 5) {
      confirmed_by_5_nucleo <- c(confirmed_by_5_nucleo, r)
    }
  }
  
  hits_m6anet_nucleo_confirmed_5 <- resize(m6anet_nucleo_all_samplings[confirmed_by_5_nucleo], 2, fix = 'center')
  hits_m6anet_nucleo_confirmed_5 <- reduce(hits_m6anet_nucleo_confirmed_5, ignore.strand = FALSE)
  print(table(width(hits_m6anet_nucleo_confirmed_5)))
  for (i in 1:length(hits_m6anet_nucleo_confirmed_5)) {
    if (width(hits_m6anet_nucleo_confirmed_5[i])<10) {
      hits_m6anet_nucleo_confirmed_5[i] <- resize(hits_m6anet_nucleo_confirmed_5[i], 10, fix='center')
    }
  }
  print(table(width(hits_m6anet_nucleo_confirmed_5)))
  
  gene_ids <- c()
  for (i in 1:length(hits_m6anet_nucleo_confirmed_5)) {
    over <- findOverlaps(hits_m6anet_nucleo_confirmed_5[i],genes_txdb)
    gene_ids <- c(gene_ids, genes_txdb[subjectHits(over)[1]]$gene_id)
  }
  mcols(hits_m6anet_nucleo_confirmed_5) <- cbind(mcols(hits_m6anet_nucleo_confirmed_5), gene_id = gene_ids)
  
  m6anet_cyto_all_samplings <- c(gr_m6anet_cyto[[1]],gr_m6anet_cyto[[2]],gr_m6anet_cyto[[3]],gr_m6anet_cyto[[4]],gr_m6anet_cyto[[5]])
  
  confirmed_by_cyto <- data.frame(num_samplings=seq(1,5), num_hits=rep(0,5))
  confirmed_by_5_cyto <- c()
  
  for (r in 1:length(m6anet_cyto_all_samplings)) {
    gr_r_rep <- queryHits(findOverlaps(m6anet_cyto_all_samplings, m6anet_cyto_all_samplings[r], type='any')) 
    confirmed_by_cyto[length(unique(m6anet_cyto_all_samplings[gr_r_rep]$rep)),2] <- confirmed_by_cyto[length(unique(m6anet_cyto_all_samplings[gr_r_rep]$rep)),2] +1
    if (length(unique(m6anet_cyto_all_samplings[gr_r_rep]$rep)) == 5) {
      confirmed_by_5_cyto <- c(confirmed_by_5_cyto, r)
    }
  }
  
  hits_m6anet_cyto_confirmed_5 <- resize(m6anet_cyto_all_samplings[confirmed_by_5_cyto], 2, fix = 'center')
  hits_m6anet_cyto_confirmed_5 <- reduce(hits_m6anet_cyto_confirmed_5, ignore.strand = FALSE)
  print(table(width(hits_m6anet_cyto_confirmed_5)))
  for (i in 1:length(hits_m6anet_cyto_confirmed_5)) {
    if (width(hits_m6anet_cyto_confirmed_5[i])<10) {
      hits_m6anet_cyto_confirmed_5[i] <- resize(hits_m6anet_cyto_confirmed_5[i], 10, fix='center')
    }
  }
  print(table(width(hits_m6anet_cyto_confirmed_5)))
  
  gene_ids <- c()
  for (i in 1:length(hits_m6anet_cyto_confirmed_5)) {
    over <- findOverlaps(hits_m6anet_cyto_confirmed_5[i],genes_txdb)
    gene_ids <- c(gene_ids, genes_txdb[subjectHits(over)[1]]$gene_id)
  }
  mcols(hits_m6anet_cyto_confirmed_5) <- cbind(mcols(hits_m6anet_cyto_confirmed_5), gene_id = gene_ids)
  
  # barplot reporting the number of hits confirmed by n samplings
  p <- ggplot(confirmed_by_chr, aes(x=num_samplings, y = num_hits)) +
    geom_bar(stat="identity") +
    theme_classic(base_size = 16) +
    scale_x_continuous(breaks=seq(1,5)) +
    labs(title='Chromatin') +
    xlab('Num samplings with hit')+
    ylab('Number of hits')
  
  ggsave(paste0(path_directory, '/chr/', name_pdf_histogram_chr_ass), plot = p)
  
  p <- ggplot(confirmed_by_nucleo, aes(x=num_samplings, y = num_hits)) +
    geom_bar(stat="identity") +
    theme_classic(base_size = 16) +
    scale_x_continuous(breaks=seq(1,5)) +
    labs(title='Nucleoplasm') +
    xlab('Num samplings with hit')+
    ylab('Number of hits')
  
  ggsave(paste0(path_directory, '/nucleo/', name_pdf_histogram_nucleo), plot = p)
  
  p <- ggplot(confirmed_by_cyto, aes(x=num_samplings, y = num_hits)) +
    geom_bar(stat="identity") +
    theme_classic(base_size = 16) +
    scale_x_continuous(breaks=seq(1,5)) +
    labs(title='Cytoplasm') +
    xlab('Num samplings with hit')+
    ylab('Number of hits')
  
  ggsave(paste0(path_directory, '/cyto/', name_pdf_histogram_cyto), plot = p)
  
  mods <- matrix(0,ncol=3,nrow=3)
  colnames(mods) <- c('chr_ass','nucleo','cyto')
  rownames(mods) <- c('# hits >= 1 samplings','# hits in 5 samplings', '# hits in 5 samplings after reduce')
  
  mods[1,1] <- length(m6anet_chr_ass_all_samplings)
  mods[2,1] <- confirmed_by_chr[5,2]
  mods[3,1] <- length(hits_m6anet_chr_ass_confirmed_5)
  mods[1,2] <- length(m6anet_nucleo_all_samplings)
  mods[2,2] <- confirmed_by_nucleo[5,2]
  mods[3,2] <- length(hits_m6anet_nucleo_confirmed_5)
  mods[1,3] <- length(m6anet_cyto_all_samplings)
  mods[2,3] <- confirmed_by_cyto[5,2]
  mods[3,3] <- length(hits_m6anet_cyto_confirmed_5)
  
  write.xlsx(x = data.frame(mods),file = paste0(path_directory, 'number_hit_5_samplings.xlsx'),col.names = TRUE, row.names = TRUE)
  
  save(hits_m6anet_chr_ass_confirmed_5, file=paste0(path_directory,'/hits_m6anet_chr_ass_confirmed_5.Rda'))
  save(hits_m6anet_nucleo_confirmed_5, file=paste0(path_directory,'/hits_m6anet_nucleo_confirmed_5.Rda'))
  save(hits_m6anet_cyto_confirmed_5, file=paste0(path_directory,'/hits_m6anet_cyto_confirmed_5.Rda'))
  
  # heatmap with the pairwise overlap between m6Anet hits in the different fractions
  overlap_hits_m6anet_confirmed_5 <- list(hits_m6anet_chr_ass_confirmed_5,hits_m6anet_nucleo_confirmed_5,hits_m6anet_cyto_confirmed_5)
  names(overlap_hits_m6anet_confirmed_5) <- c('Chromatin', 'Nucleoplasm','Cytoplasm')
  
  pdf(file = paste0(path_directory, '/', name_pdf_overlap_3fractions), width = 7, height = 5)
  overlapOfGRanges(overlap_hits_m6anet_confirmed_5,plot = TRUE)
  dev.off()
  
  return(overlap_hits_m6anet_confirmed_5)
}

# path_directory is the path to the directory containing three folders: /chr/, /nucleo/, /cyto/ each with the 5 tsv files 
# produced by m6Anet for each of the 5 samplings
DRACH_overlap_m6Anet <- function(path_directory, hits_m6Anet_chr, hits_m6Anet_nucleo, hits_m6Anet_cyto) {
  
  summary_table <- matrix(ncol=3, nrow=4)
  colnames(summary_table) <- c('Total number of hits', 'DRACH+ hits', 'DRACH- hits')
  rownames(summary_table) <- c('Chromatin', 'Nucleoplasm', 'Cytoplasm', 'Overlap')
  
  summary_table[1,1] <- as.character(length(hits_m6Anet_chr))
  summary_table[2,1] <- as.character(length(hits_m6Anet_nucleo))
  summary_table[3,1] <- as.character(length(hits_m6Anet_cyto))
  
  confirmed_hits <- list(hits_m6Anet_chr, hits_m6Anet_nucleo, hits_m6Anet_cyto)
  order <- order(c(length(hits_m6Anet_chr), length(hits_m6Anet_nucleo), length(hits_m6Anet_cyto)))
  
  # compute the number of m6Anet hits present in all the fractions
  overlap_m6anet <- findOverlaps(confirmed_hits[[order[1]]],confirmed_hits[[order[2]]], type = 'any')
  overlap_m6anet_3_sets <- confirmed_hits[[order[1]]][unique(queryHits(overlap_m6anet))]
  overlap_m6anet_chr_cyto_nucleo <- findOverlaps(overlap_m6anet_3_sets,confirmed_hits[[order[3]]], type = 'any')
  hits_m6anet_confirmed_by_all <- overlap_m6anet_3_sets[unique(queryHits(overlap_m6anet_chr_cyto_nucleo))]
  summary_table[4,1] <- as.character(length(hits_m6anet_confirmed_by_all))
  
  load('/path/to/R_data/DRACH_forward_strand.Rda')
  # identify which and how many hits are DRACH+/DRACH-
  overlap_DRACH_chr <- findOverlaps(hits_m6Anet_chr, DRACH, type='any', ignore.strand=FALSE)
  hits_m6anet_chr_ass_confirmed_5_with_DRACH <- hits_m6Anet_chr[unique(queryHits(overlap_DRACH_chr))]
  save(hits_m6anet_chr_ass_confirmed_5_with_DRACH, file=paste0(path_directory,'/hits_m6anet_chr_ass_confirmed_5_with_DRACH.Rda'))
  summary_table[1,2] <- paste0(as.character(length(hits_m6anet_chr_ass_confirmed_5_with_DRACH)), ' - ', as.character(round(length(hits_m6anet_chr_ass_confirmed_5_with_DRACH)*100/length(hits_m6Anet_chr), 2)),'%')
  
  hits_m6anet_chr_ass_confirmed_5_without_DRACH <- hits_m6Anet_chr[-unique(queryHits(overlap_DRACH_chr))]
  summary_table[1,3] <- as.character(length(hits_m6anet_chr_ass_confirmed_5_without_DRACH))
  
  overlap_DRACH_nucleo <- findOverlaps(hits_m6Anet_nucleo, DRACH, type='any', ignore.strand=FALSE)
  hits_m6anet_nucleo_confirmed_5_with_DRACH <- hits_m6Anet_nucleo[unique(queryHits(overlap_DRACH_nucleo))]
  save(hits_m6anet_nucleo_confirmed_5_with_DRACH, file=paste0(path_directory,'/hits_m6anet_nucleo_confirmed_5_with_DRACH.Rda'))
  summary_table[2,2] <- paste0(as.character(length(hits_m6anet_nucleo_confirmed_5_with_DRACH)), ' - ', as.character(round(length(hits_m6anet_nucleo_confirmed_5_with_DRACH)*100/length(hits_m6Anet_nucleo), 2)),'%')
  
  hits_m6anet_nucleo_confirmed_5_without_DRACH <- hits_m6Anet_nucleo[-unique(queryHits(overlap_DRACH_nucleo))]
  summary_table[2,3] <- as.character(length(hits_m6anet_nucleo_confirmed_5_without_DRACH))
  
  overlap_DRACH_cyto <- findOverlaps(hits_m6Anet_cyto, DRACH, type='any', ignore.strand=FALSE)
  hits_m6anet_cyto_confirmed_5_with_DRACH <- hits_m6Anet_cyto[unique(queryHits(overlap_DRACH_cyto))]
  save(hits_m6anet_cyto_confirmed_5_with_DRACH, file=paste0(path_directory,'/hits_m6anet_cyto_confirmed_5_with_DRACH.Rda'))
  summary_table[3,2] <- paste0(as.character(length(hits_m6anet_cyto_confirmed_5_with_DRACH)), ' - ', as.character(round(length(hits_m6anet_cyto_confirmed_5_with_DRACH)*100/length(hits_m6Anet_cyto), 2)),'%')
  
  hits_m6anet_cyto_confirmed_5_without_DRACH <- hits_m6Anet_cyto[-unique(queryHits(overlap_DRACH_cyto))]
  summary_table[3,3] <- as.character(length(hits_m6anet_cyto_confirmed_5_without_DRACH))
  
  # compute the number of m6Anet DRACH+ hits present in all the fractions
  confirmed_hits_DRACH <- list(hits_m6anet_chr_ass_confirmed_5_with_DRACH, hits_m6anet_nucleo_confirmed_5_with_DRACH, hits_m6anet_cyto_confirmed_5_with_DRACH)
  order <- order(c(length(hits_m6anet_chr_ass_confirmed_5_with_DRACH), length(hits_m6anet_nucleo_confirmed_5_with_DRACH), length(hits_m6anet_cyto_confirmed_5_with_DRACH)))
  
  overlap_m6anet_DRACH <- findOverlaps(confirmed_hits_DRACH[[order[1]]],confirmed_hits_DRACH[[order[2]]], type = 'any')
  overlap_m6anet_3_sets_DRACH <- confirmed_hits_DRACH[[order[1]]][unique(queryHits(overlap_m6anet_DRACH))]
  overlap_m6anet_chr_cyto_nucleo_DRACH <- findOverlaps(overlap_m6anet_3_sets_DRACH,confirmed_hits_DRACH[[order[3]]], type = 'any')
  hits_m6anet_confirmed_by_all_DRACH <- overlap_m6anet_3_sets_DRACH[unique(queryHits(overlap_m6anet_chr_cyto_nucleo_DRACH))]
  summary_table[4,2] <- as.character(length(hits_m6anet_confirmed_by_all_DRACH))
  
  # compute the number of m6Anet DRACH- hits present in all the fractions
  confirmed_hits_without_DRACH <- list(hits_m6anet_chr_ass_confirmed_5_without_DRACH, hits_m6anet_nucleo_confirmed_5_without_DRACH, hits_m6anet_cyto_confirmed_5_without_DRACH)
  order <- order(c(length(hits_m6anet_chr_ass_confirmed_5_without_DRACH), length(hits_m6anet_nucleo_confirmed_5_without_DRACH), length(hits_m6anet_cyto_confirmed_5_without_DRACH)))
  
  overlap_m6anet_without_DRACH <- findOverlaps(confirmed_hits_without_DRACH[[order[1]]],confirmed_hits_without_DRACH[[order[2]]], type = 'any')
  overlap_m6anet_3_sets_without_DRACH <- confirmed_hits_without_DRACH[[order[1]]][unique(queryHits(overlap_m6anet_without_DRACH))]
  overlap_m6anet_chr_cyto_nucleo_without_DRACH <- findOverlaps(overlap_m6anet_3_sets_without_DRACH,confirmed_hits_without_DRACH[[order[3]]], type = 'any')
  hits_m6anet_confirmed_by_all_without_DRACH <- overlap_m6anet_3_sets_without_DRACH[unique(queryHits(overlap_m6anet_chr_cyto_nucleo_without_DRACH))]
  summary_table[4,3] <- as.character(length(hits_m6anet_confirmed_by_all_without_DRACH))
  
  write.xlsx(x = data.frame(summary_table),file = paste0(path_directory, 'summary_table.xlsx'),col.names = TRUE, row.names = TRUE)
  
  return(confirmed_hits_DRACH)
}

# path_directory is the path to directory containing three folders: /chr/, /nucleo/, /cyto/ each with the 5 tsv files 
# produced by m6Anet for each of the 5 samplings.
# create a folder /ELIGOS_confirmed_DRACH/ inside path_directory where you will save 
# ELIGOS DRACH+ hits confirmed by m6Anet. Create a folder /not/ inside /ELIGOS_confirmed_DRACH/ where you will save  
# ELIGOS DRACH+ hits not confirmed by m6Anet
comparison_ELIGOS_m6Anet <- function(path_directory, hits_ELIGOS_chr, hits_ELIGOS_nucleo, hits_ELIGOS_cyto, 
                                     hits_ELIGOS_chr_DRACH, hits_ELIGOS_nucleo_DRACH, hits_ELIGOS_cyto_DRACH, 
                                     hits_m6Anet_chr, hits_m6Anet_nucleo, hits_m6Anet_cyto,
                                     hits_m6Anet_chr_DRACH, hits_m6Anet_nucleo_DRACH, hits_m6Anet_cyto_DRACH, 
                                     path_pdf_overlap_ELIGOS_m6Anet_DRACH_chr, path_pdf_overlap_ELIGOS_m6Anet_DRACH_nucleo, path_pdf_overlap_ELIGOS_m6Anet_DRACH_cyto) {
  
  # represent in a heatmap, for each fraction, the overlap between ELIGOS hits confirmed by 5 samplings, m6Anet hits 
  # confirmed by 5 samplings and ELIGOS DRACH+ hits 
  overlap_eligos_all_eligos_DRACH_m6anet_all_chr_ass <- list(hits_ELIGOS_chr, hits_ELIGOS_chr_DRACH,hits_m6Anet_chr)
  names(overlap_eligos_all_eligos_DRACH_m6anet_all_chr_ass) <- c('All hits\nELIGOS','DRACH+\nELIGOS hits','All hits\nm6Anet')
  
  pdf(file = path_pdf_overlap_ELIGOS_m6Anet_DRACH_chr, width = 5, height = 7)
  overlapOfGRanges(overlap_eligos_all_eligos_DRACH_m6anet_all_chr_ass,plot = TRUE)
  dev.off()
  
  overlap_eligos_all_eligos_DRACH_m6anet_all_nucleo <- list(hits_ELIGOS_nucleo, hits_ELIGOS_nucleo_DRACH,hits_m6Anet_nucleo)
  names(overlap_eligos_all_eligos_DRACH_m6anet_all_nucleo) <- c('All hits\nELIGOS','DRACH+nELIGOS hits','All hits\nm6Anet')
  
  pdf(file = path_pdf_overlap_ELIGOS_m6Anet_DRACH_nucleo, width = 5, height = 7)
  overlapOfGRanges(overlap_eligos_all_eligos_DRACH_m6anet_all_nucleo,plot = TRUE)
  dev.off()
  
  overlap_eligos_all_eligos_DRACH_m6anet_all_cyto <- list(hits_ELIGOS_cyto, hits_ELIGOS_cyto_DRACH,hits_m6Anet_cyto)
  names(overlap_eligos_all_eligos_DRACH_m6anet_all_cyto) <- c('All hits\nELIGOS','DRACH+\nELIGOS hits','All hits\nm6Anet')
  
  pdf(file = path_pdf_overlap_ELIGOS_m6Anet_DRACH_cyto, width = 5, height = 7)
  overlapOfGRanges(overlap_eligos_all_eligos_DRACH_m6anet_all_cyto,plot = TRUE)
  dev.off()
  
  # compute the overlap between ELIGOS hits confirmed by 5 samplings and m6Anet hits confirmed by 5 samplings
  summary_table <- matrix(ncol=3, nrow=3)
  colnames(summary_table) <- c('Total number of ELIGOS hits', 'Total number of m6Anet hits', 'Overlap')
  rownames(summary_table) <- c('Chromatin', 'Nucleoplasm', 'Cytoplasm')
  
  summary_table[1,1] <- as.character(length(hits_ELIGOS_chr))
  summary_table[2,1] <- as.character(length(hits_ELIGOS_nucleo))
  summary_table[3,1] <- as.character(length(hits_ELIGOS_cyto))
  
  summary_table[1,2] <- as.character(length(hits_m6Anet_chr))
  summary_table[2,2] <- as.character(length(hits_m6Anet_nucleo))
  summary_table[3,2] <- as.character(length(hits_m6Anet_cyto))
  
  all_hits_chr <- list(hits_ELIGOS_chr, hits_m6Anet_chr)
  names(all_hits_chr) <- c('ELIGOS', 'm6Anet')
  order <- order(c(length(hits_ELIGOS_chr), length(hits_m6Anet_chr)))
  
  overlap <- findOverlaps(all_hits_chr[[order[1]]],all_hits_chr[[order[2]]], type = 'any')
  overlap_confirmed_hits_chr_ass_eligos_m6anet <- all_hits_chr[[order[1]]][unique(queryHits(overlap))]
  overlap_confirmed_hits_chr_ass_m6anet_eligos <- all_hits_chr[[order[2]]][unique(subjectHits(overlap))]
  summary_table[1,3] <- paste0(as.character(length(overlap_confirmed_hits_chr_ass_eligos_m6anet)), ' (', as.character(round(length(overlap_confirmed_hits_chr_ass_eligos_m6anet)*100/length(all_hits_chr[[order[1]]]), 2)),'%) - ', as.character(length(overlap_confirmed_hits_chr_ass_m6anet_eligos)), ' (', as.character(round(length(overlap_confirmed_hits_chr_ass_m6anet_eligos)*100/length(all_hits_chr[[order[2]]]), 2)), '%)')
  
  all_hits_nucleo <- list(hits_ELIGOS_nucleo, hits_m6Anet_nucleo)
  names(all_hits_nucleo) <- c('ELIGOS', 'm6Anet')
  order <- order(c(length(hits_ELIGOS_nucleo), length(hits_m6Anet_nucleo)))
  
  overlap <- findOverlaps(all_hits_nucleo[[order[1]]],all_hits_nucleo[[order[2]]], type = 'any')
  overlap_confirmed_hits_nucleo_eligos_m6anet <- all_hits_nucleo[[order[1]]][unique(queryHits(overlap))]
  overlap_confirmed_hits_nucleo_m6anet_eligos <- all_hits_nucleo[[order[2]]][unique(subjectHits(overlap))]
  summary_table[2,3] <- paste0(as.character(length(overlap_confirmed_hits_nucleo_eligos_m6anet)), ' (', as.character(round(length(overlap_confirmed_hits_nucleo_eligos_m6anet)*100/length(all_hits_nucleo[[order[1]]]), 2)),'%) - ', as.character(length(overlap_confirmed_hits_nucleo_m6anet_eligos)), ' (', as.character(round(length(overlap_confirmed_hits_nucleo_m6anet_eligos)*100/length(all_hits_nucleo[[order[2]]]), 2)), '%)')
  
  all_hits_cyto <- list(hits_ELIGOS_cyto, hits_m6Anet_cyto)
  names(all_hits_cyto) <- c('ELIGOS', 'm6Anet')
  order <- order(c(length(hits_ELIGOS_cyto), length(hits_m6Anet_cyto)))
  
  overlap <- findOverlaps(all_hits_cyto[[order[1]]],all_hits_cyto[[order[2]]], type = 'any')
  overlap_confirmed_hits_cyto_eligos_m6anet <- all_hits_cyto[[order[1]]][unique(queryHits(overlap))]
  overlap_confirmed_hits_cyto_m6anet_eligos <- all_hits_cyto[[order[2]]][unique(subjectHits(overlap))]
  summary_table[3,3] <- paste0(as.character(length(overlap_confirmed_hits_cyto_eligos_m6anet)), ' (', as.character(round(length(overlap_confirmed_hits_cyto_eligos_m6anet)*100/length(all_hits_cyto[[order[1]]]), 2)),'%) - ', as.character(length(overlap_confirmed_hits_cyto_m6anet_eligos)), ' (', as.character(round(length(overlap_confirmed_hits_cyto_m6anet_eligos)*100/length(all_hits_cyto[[order[2]]]), 2)), '%)')
  
  write.xlsx(x = data.frame(summary_table),file = paste0(path_directory, 'summary_table_eligos_m6anet_all.xlsx'),col.names = TRUE, row.names = TRUE)
  
  # compute the overlap between ELIGOS DRACH+ hits confirmed by 5 samplings and m6Anet DRACH+ hits 
  # confirmed by 5 samplings
  summary_table <- matrix(ncol=3, nrow=3)
  colnames(summary_table) <- c('DRACH+ ELIGOS hits', 'DRACH+ m6Anet hits', 'Overlap')
  rownames(summary_table) <- c('Chromatin', 'Nucleoplasm', 'Cytoplasm')
  
  summary_table[1,1] <- as.character(length(hits_ELIGOS_chr_DRACH))
  summary_table[2,1] <- as.character(length(hits_ELIGOS_nucleo_DRACH))
  summary_table[3,1] <- as.character(length(hits_ELIGOS_cyto_DRACH))
  
  summary_table[1,2] <- as.character(length(hits_m6Anet_chr_DRACH))
  summary_table[2,2] <- as.character(length(hits_m6Anet_nucleo_DRACH))
  summary_table[3,2] <- as.character(length(hits_m6Anet_cyto_DRACH))
  
  all_hits_chr_DRACH <- list(hits_ELIGOS_chr_DRACH, hits_m6Anet_chr_DRACH)
  names(all_hits_chr_DRACH) <- c('ELIGOS_DRACH', 'm6Anet_DRACH')
  order <- order(c(length(hits_ELIGOS_chr_DRACH), length(hits_m6Anet_chr_DRACH)))
  
  overlap <- findOverlaps(all_hits_chr_DRACH[[order[1]]],all_hits_chr_DRACH[[order[2]]], type = 'any')
  overlap_DRACH_m6anet_DRACH_eligos_chr_ass <- all_hits_chr_DRACH[[order[1]]][unique(queryHits(overlap))]
  save(overlap_DRACH_m6anet_DRACH_eligos_chr_ass, file=paste0(path_directory, '/ELIGOS_confirmed_DRACH/hits_eligos_DRACH_chr_confirmed_m6Anet.Rda'))
  eligos_DRACH_chr_non_confirmed <- all_hits_chr_DRACH[[order[1]]][-unique(queryHits(overlap))]
  save(eligos_DRACH_chr_non_confirmed, file=paste0(path_directory, '/ELIGOS_confirmed_DRACH/not/hits_eligos_DRACH_chr_not_confirmed.Rda'))
  overlap_hits_DRACH_eligos_DRACH_m6anet_chr_ass <- all_hits_chr_DRACH[[order[2]]][unique(subjectHits(overlap))]
  summary_table[1,3] <- paste0(as.character(length(overlap_DRACH_m6anet_DRACH_eligos_chr_ass)), ' (', as.character(round(length(overlap_DRACH_m6anet_DRACH_eligos_chr_ass)*100/length(all_hits_chr_DRACH[[order[1]]]), 2)),'%) - ', as.character(length(overlap_hits_DRACH_eligos_DRACH_m6anet_chr_ass)), ' (', as.character(round(length(overlap_hits_DRACH_eligos_DRACH_m6anet_chr_ass)*100/length(all_hits_chr_DRACH[[order[2]]]), 2)), '%)')
  
  all_hits_nucleo_DRACH <- list(hits_ELIGOS_nucleo_DRACH, hits_m6Anet_nucleo_DRACH)
  names(all_hits_nucleo_DRACH) <- c('ELIGOS_DRACH', 'm6Anet_DRACH')
  order <- order(c(length(hits_ELIGOS_nucleo_DRACH), length(hits_m6Anet_nucleo_DRACH)))
  
  overlap <- findOverlaps(all_hits_nucleo_DRACH[[order[1]]],all_hits_nucleo_DRACH[[order[2]]], type = 'any')
  overlap_DRACH_m6anet_DRACH_eligos_nucleo <- all_hits_nucleo_DRACH[[order[1]]][unique(queryHits(overlap))]
  save(overlap_DRACH_m6anet_DRACH_eligos_nucleo, file=paste0(path_directory, '/ELIGOS_confirmed_DRACH/hits_eligos_DRACH_nucleo_confirmed_m6Anet.Rda'))
  eligos_DRACH_nucleo_non_confirmed <- all_hits_nucleo_DRACH[[order[1]]][-unique(queryHits(overlap))]
  save(eligos_DRACH_nucleo_non_confirmed, file=paste0(path_directory, '/ELIGOS_confirmed_DRACH/not/hits_eligos_DRACH_nucleo_not_confirmed.Rda'))
  overlap_hits_DRACH_eligos_DRACH_m6anet_nucleo <- all_hits_nucleo_DRACH[[order[2]]][unique(subjectHits(overlap))]
  summary_table[2,3] <- paste0(as.character(length(overlap_DRACH_m6anet_DRACH_eligos_nucleo)), ' (', as.character(round(length(overlap_DRACH_m6anet_DRACH_eligos_nucleo)*100/length(all_hits_nucleo_DRACH[[order[1]]]), 2)),'%) - ', as.character(length(overlap_hits_DRACH_eligos_DRACH_m6anet_nucleo)), ' (', as.character(round(length(overlap_hits_DRACH_eligos_DRACH_m6anet_nucleo)*100/length(all_hits_nucleo_DRACH[[order[2]]]), 2)), '%)')
  
  all_hits_cyto_DRACH <- list(hits_ELIGOS_cyto_DRACH, hits_m6Anet_cyto_DRACH)
  names(all_hits_cyto_DRACH) <- c('ELIGOS_DRACH', 'm6Anet_DRACH')
  order <- order(c(length(hits_ELIGOS_cyto_DRACH), length(hits_m6Anet_cyto_DRACH)))
  
  overlap <- findOverlaps(all_hits_cyto_DRACH[[order[1]]],all_hits_cyto_DRACH[[order[2]]], type = 'any')
  overlap_DRACH_m6anet_DRACH_eligos_cyto <- all_hits_cyto_DRACH[[order[1]]][unique(queryHits(overlap))]
  save(overlap_DRACH_m6anet_DRACH_eligos_cyto, file=paste0(path_directory, '/ELIGOS_confirmed_DRACH/hits_eligos_DRACH_cyto_confirmed_m6Anet.Rda'))
  eligos_DRACH_cyto_non_confirmed <- all_hits_cyto_DRACH[[order[1]]][-unique(queryHits(overlap))]
  save(eligos_DRACH_cyto_non_confirmed, file=paste0(path_directory, '/ELIGOS_confirmed_DRACH/not/hits_eligos_DRACH_cyto_not_confirmed.Rda'))
  overlap_hits_DRACH_eligos_DRACH_m6anet_cyto <- all_hits_cyto_DRACH[[order[2]]][unique(subjectHits(overlap))]
  summary_table[3,3] <- paste0(as.character(length(overlap_DRACH_m6anet_DRACH_eligos_cyto)), ' (', as.character(round(length(overlap_DRACH_m6anet_DRACH_eligos_cyto)*100/length(all_hits_cyto_DRACH[[order[1]]]), 2)),'%) - ', as.character(length(overlap_hits_DRACH_eligos_DRACH_m6anet_cyto)), ' (', as.character(round(length(overlap_hits_DRACH_eligos_DRACH_m6anet_cyto)*100/length(all_hits_cyto_DRACH[[order[2]]]), 2)), '%)')
  
  write.xlsx(x = data.frame(summary_table),file = paste0(path_directory, 'summary_table_eligos_m6anet_DRACH.xlsx'),col.names = TRUE, row.names = TRUE)
}


# ELIGOS ANALYSIS ON NASCENT+PRE-EXISTING READS (pvalue<0.05, ad.pval<0.05, OR>1)
confirmed_hits_ELIGOS_min05_min05_mag1 <- ELIGOS_results(path_directory = '/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/',
                                                         name_pdf_overlap_10samplings_cyto = 'overlap_cyto_replicates_4sU_10nt.pdf', name_pdf_overlap_10samplings_chr_ass = 'overlap_chr_ass_replicates_4sU_10nt.pdf', name_pdf_overlap_10samplings_nucleo = 'overlap_nucleo_replicates_4sU_10nt.pdf',
                                                         name_pdf_histogram_cyto = 'number_rep_with_hit_cyto_4sU.pdf', name_pdf_histogram_chr_ass = 'number_rep_with_hit_chr_ass_4sU.pdf', name_pdf_histogram_nucleo = 'number_rep_with_hit_nucleo_4sU.pdf',
                                                         name_pdf_overlap_3fractions = 'overlap_hits_eligos_confirmed_5_4sU_10nt.pdf', p=0.05,ap=0.05,OR=1)

confirmed_hits_DRACH_ELIGOS_min05_min05_mag1 <- DRACH_overlap_ELIGOS(path_directory = '/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/',
                                                                     hits_ELIGOS_chr = confirmed_hits_ELIGOS_min05_min05_mag1[[1]], hits_ELIGOS_nucleo = confirmed_hits_ELIGOS_min05_min05_mag1[[2]], hits_ELIGOS_cyto = confirmed_hits_ELIGOS_min05_min05_mag1[[3]])

# m6Anet ANALYSIS ON NASCENT+PRE-EXISTING READS (prob.modification>0.75)
# load the tsv file with the output of m6Anet for the three fractions for all the 5 samplings
m6Anet_fractions_4sU_chr_ass <- list.files(path = '/path/to/fractions_m6anet_4sU_library_gene_subsampling_prob0.75/chr', pattern = 'tsv', full.names = TRUE)
m6Anet_fractions_4sU_nucleo <- list.files(path = '/path/to/fractions_m6anet_4sU_library_gene_subsampling_prob0.75/nucleo', pattern = 'tsv', full.names = TRUE)
m6Anet_fractions_4sU_cyto <- list.files(path = '/path/to/fractions_m6anet_4sU_library_gene_subsampling_prob0.75/cyto', pattern = 'tsv', full.names = TRUE)

# load the vector with the transcript names as names and the names of the gene from which 
# they are transcribed as values
load('/path/to/R_data/tx_gene.Rda')
tx_txdb <- GenomicFeatures::transcripts(txdb)

m6Anet_output <- function(dir) {
  list_gr <- list()
  
  for (file in dir) {
    t <- read.table(file, header=TRUE)
    
    gr <- GRanges(seqnames = Rle(t$Chr),
                  ranges = IRanges(start = t$Start, end=t$End),
                  strand = Rle(t$Strand),
                  prob_mod = t$Prob_mod,
                  gene_id = unname(tx_gene[t$Transcript_id]),
                  tx_id = t$Transcript_id,
                  rep = unlist(strsplit(strsplit(file, split = '/')[[1]][8], "_")[[1]][5]))
    
    # extract the analysed sites with a probability of modification of at least 0.75
    gr <- gr[which(gr$prob_mod >= 0.75)]
    gr <- resize(gr, 6, fix = 'center')
    
    list_gr <- c(list_gr, gr)
  }
  
  return(list_gr)
}

gr_m6Anet_chr_ass_4sU <- m6Anet_output(m6Anet_fractions_4sU_chr_ass)
gr_m6Anet_nucleo_4sU <- m6Anet_output(m6Anet_fractions_4sU_nucleo)
gr_m6Anet_cyto_4sU <- m6Anet_output(m6Anet_fractions_4sU_cyto)

# remove the hits mapping on transcripts that have been assigned to the wrong chromosome by m6Anet
gr_m6Anet_chr_ass_4sU <- lapply(gr_m6Anet_chr_ass_4sU, function(x) {
  rem <- c()
  for (i in 1:length(x)) {
    tx <- x[i]$tx_id
    chr <- as.character(seqnames(x[i]))
    if (chr != as.character(seqnames(tx_txdb[tx_txdb$tx_name == tx]))) {
      rem <- c(rem, i)
    }
  }
  if (length(rem) != 0) {
    return(x[-rem])
  } else {
    return(x)
  }
})

gr_m6Anet_nucleo_4sU <- lapply(gr_m6Anet_nucleo_4sU, function(x) {
  rem <- c()
  for (i in 1:length(x)) {
    tx <- x[i]$tx_id
    chr <- as.character(seqnames(x[i]))
    if (chr != as.character(seqnames(tx_txdb[tx_txdb$tx_name == tx]))) {
      rem <- c(rem, i)
    }
  }
  if (length(rem) != 0) {
    return(x[-rem])
  } else {
    return(x)
  }
})

gr_m6Anet_cyto_4sU <- lapply(gr_m6Anet_cyto_4sU, function(x) {
  rem <- c()
  for (i in 1:length(x)) {
    tx <- x[i]$tx_id
    chr <- as.character(seqnames(x[i]))
    if (chr != as.character(seqnames(tx_txdb[tx_txdb$tx_name == tx]))) {
      rem <- c(rem, i)
    }
  }
  if (length(rem) != 0) {
    return(x[-rem])
  } else {
    return(x)
  }
})

confirmed_hits_m6Anet_prob0.75 <- m6Anet_results(path_directory = '/path/to/fractions_m6anet_4sU_library_gene_subsampling_prob0.75/',
                                                 gr_m6Anet_chr_ass = gr_m6Anet_chr_ass_4sU, gr_m6Anet_nucleo = gr_m6Anet_nucleo_4sU, gr_m6Anet_cyto = gr_m6Anet_cyto_4sU,
                                                 name_pdf_overlap_10samplings_chr_ass = 'overlap_chr_ass_replicates_4sU_10nt.pdf', name_pdf_overlap_10samplings_nucleo = 'overlap_nucleo_replicates_4sU_10nt.pdf', name_pdf_overlap_10samplings_cyto = 'overlap_cyto_replicates_4sU_10nt.pdf',
                                                 name_pdf_histogram_chr_ass = 'number_rep_with_hit_chr_ass_4sU.pdf', name_pdf_histogram_nucleo = 'number_rep_with_hit_nucleo_4sU.pdf', name_pdf_histogram_cyto = 'number_rep_with_hit_cyto_4sU.pdf',
                                                 name_pdf_overlap_3fractions = 'overlap_hits_m6anet_confirmed_5_4sU_10nt.pdf')

confirmed_hits_DRACH_m6Anet_prob0.75 <- DRACH_overlap_m6Anet(path_directory = '/path/to/fractions_m6anet_4sU_library_gene_subsampling_prob0.75/',
                                                             hits_m6Anet_chr = confirmed_hits_m6Anet_prob0.75[[1]], hits_m6Anet_nucleo = confirmed_hits_m6Anet_prob0.75[[2]], hits_m6Anet_cyto = confirmed_hits_m6Anet_prob0.75[[3]])

ELIGOS_vs_m6Anet <- comparison_ELIGOS_m6Anet(path_directory = '/path/to/fractions_m6anet_4sU_library_gene_subsampling_prob0.75/',hits_ELIGOS_chr = confirmed_hits_ELIGOS_min05_min05_mag1[[1]], hits_ELIGOS_nucleo = confirmed_hits_ELIGOS_min05_min05_mag1[[2]], hits_ELIGOS_cyto = confirmed_hits_ELIGOS_min05_min05_mag1[[3]],
                                             hits_ELIGOS_chr_DRACH = confirmed_hits_DRACH_ELIGOS_min05_min05_mag1[[1]], hits_ELIGOS_nucleo_DRACH = confirmed_hits_DRACH_ELIGOS_min05_min05_mag1[[2]], hits_ELIGOS_cyto_DRACH = confirmed_hits_DRACH_ELIGOS_min05_min05_mag1[[3]],
                                             hits_m6Anet_chr = confirmed_hits_m6Anet_prob0.75[[1]], hits_m6Anet_nucleo = confirmed_hits_m6Anet_prob0.75[[2]], hits_m6Anet_cyto = confirmed_hits_m6Anet_prob0.75[[3]],
                                             hits_m6Anet_chr_DRACH = confirmed_hits_DRACH_m6Anet_prob0.75[[1]], hits_m6Anet_nucleo_DRACH = confirmed_hits_DRACH_m6Anet_prob0.75[[2]], hits_m6Anet_cyto_DRACH = confirmed_hits_DRACH_m6Anet_prob0.75[[3]],
                                             path_pdf_overlap_ELIGOS_m6Anet_DRACH_chr = '/path/to/fractions_m6anet_4sU_library_gene_subsampling_prob0.75/overlap_all_hits_eligos_hits_eligos_DRACH_all_hits_m6anet_chr_ass_10nt.pdf',
                                             path_pdf_overlap_ELIGOS_m6Anet_DRACH_nucleo = '/path/to/fractions_m6anet_4sU_library_gene_subsampling_prob0.75/overlap_all_hits_eligos_hits_eligos_DRACH_all_hits_m6anet_nucleo_10nt.pdf',
                                             path_pdf_overlap_ELIGOS_m6Anet_DRACH_cyto = '/path/to/fractions_m6anet_4sU_library_gene_subsampling_prob0.75/overlap_all_hits_eligos_hits_eligos_DRACH_all_hits_m6anet_cyto_10nt.pdf')

######################
# m6Anet ANALYSIS ON NASCENT+PRE-EXISTING READS (prob.modification>0.9)
# load the tsv file with the output of m6Anet for the three fractions for all the 5 samplings
m6Anet_fractions_4sU_chr_ass <- list.files(path = '/path/to/fractions_m6anet_4sU_library_gene_subsampling_prob0.9/chr', pattern = 'tsv', full.names = TRUE)
m6Anet_fractions_4sU_nucleo <- list.files(path = '/path/to/fractions_m6anet_4sU_library_gene_subsampling_prob0.9/nucleo', pattern = 'tsv', full.names = TRUE)
m6Anet_fractions_4sU_cyto <- list.files(path = '/path/to/fractions_m6anet_4sU_library_gene_subsampling_prob0.9/cyto', pattern = 'tsv', full.names = TRUE)

# load the vector with the transcript names as names and the names of the gene from which 
# they are transcribed as values
load('/path/to/R_data/tx_gene.Rda')
tx_txdb <- GenomicFeatures::transcripts(txdb)

m6Anet_output <- function(dir) {
  list_gr <- list()
  
  for (file in dir) {
    t <- read.table(file, header=TRUE)
    
    gr <- GRanges(seqnames = Rle(t$Chr),
                  ranges = IRanges(start = t$Start, end=t$End),
                  strand = Rle(t$Strand),
                  prob_mod = t$Prob_mod,
                  gene_id = unname(tx_gene[t$Transcript_id]),
                  tx_id = t$Transcript_id,
                  rep = unlist(strsplit(strsplit(file, split = '/')[[1]][8], "_")[[1]][5]))
    
    # extract the analysed sites with a probability of modification of at least 0.9
    gr <- gr[which(gr$prob_mod >= 0.9)]
    gr <- resize(gr, 6, fix = 'center')
    
    list_gr <- c(list_gr, gr)
  }
  
  return(list_gr)
}

gr_m6Anet_chr_ass_4sU <- m6Anet_output(m6Anet_fractions_4sU_chr_ass)
gr_m6Anet_nucleo_4sU <- m6Anet_output(m6Anet_fractions_4sU_nucleo)
gr_m6Anet_cyto_4sU <- m6Anet_output(m6Anet_fractions_4sU_cyto)

# remove the hits mapping on transcripts that have been assigned to the wrong chromosome by m6Anet
gr_m6Anet_chr_ass_4sU <- lapply(gr_m6Anet_chr_ass_4sU, function(x) {
  rem <- c()
  for (i in 1:length(x)) {
    tx <- x[i]$tx_id
    chr <- as.character(seqnames(x[i]))
    if (chr != as.character(seqnames(tx_txdb[tx_txdb$tx_name == tx]))) {
      rem <- c(rem, i)
    }
  }
  if (length(rem) != 0) {
    return(x[-rem])
  } else {
    return(x)
  }
})

gr_m6Anet_nucleo_4sU <- lapply(gr_m6Anet_nucleo_4sU, function(x) {
  rem <- c()
  for (i in 1:length(x)) {
    tx <- x[i]$tx_id
    chr <- as.character(seqnames(x[i]))
    if (chr != as.character(seqnames(tx_txdb[tx_txdb$tx_name == tx]))) {
      rem <- c(rem, i)
    }
  }
  if (length(rem) != 0) {
    return(x[-rem])
  } else {
    return(x)
  }
})

gr_m6Anet_cyto_4sU <- lapply(gr_m6Anet_cyto_4sU, function(x) {
  rem <- c()
  for (i in 1:length(x)) {
    tx <- x[i]$tx_id
    chr <- as.character(seqnames(x[i]))
    if (chr != as.character(seqnames(tx_txdb[tx_txdb$tx_name == tx]))) {
      rem <- c(rem, i)
    }
  }
  if (length(rem) != 0) {
    return(x[-rem])
  } else {
    return(x)
  }
})

confirmed_hits_m6Anet_prob0.9 <- m6Anet_results(path_directory = '/path/to/fractions_m6anet_4sU_library_gene_subsampling_prob0.9/',
                                                gr_m6Anet_chr_ass = gr_m6Anet_chr_ass_4sU, gr_m6Anet_nucleo = gr_m6Anet_nucleo_4sU, gr_m6Anet_cyto = gr_m6Anet_cyto_4sU,
                                                name_pdf_overlap_10samplings_chr_ass = 'overlap_chr_ass_replicates_4sU_10nt.pdf', name_pdf_overlap_10samplings_nucleo = 'overlap_nucleo_replicates_4sU_10nt.pdf', name_pdf_overlap_10samplings_cyto = 'overlap_cyto_replicates_4sU_10nt.pdf',
                                                name_pdf_histogram_chr_ass = 'number_rep_with_hit_chr_ass_4sU.pdf', name_pdf_histogram_nucleo = 'number_rep_with_hit_nucleo_4sU.pdf', name_pdf_histogram_cyto = 'number_rep_with_hit_cyto_4sU.pdf',
                                                name_pdf_overlap_3fractions = 'overlap_hits_m6anet_confirmed_5_4sU_10nt.pdf')

confirmed_hits_DRACH_m6Anet_prob0.9 <- DRACH_overlap_m6Anet(path_directory = '/path/to/fractions_m6anet_4sU_library_gene_subsampling_prob0.9/',
                                                            hits_m6Anet_chr = confirmed_hits_m6Anet_prob0.9[[1]], hits_m6Anet_nucleo = confirmed_hits_m6Anet_prob0.9[[2]], hits_m6Anet_cyto = confirmed_hits_m6Anet_prob0.9[[3]])

ELIGOS_vs_m6Anet <- comparison_ELIGOS_m6Anet(path_directory = '/path/to/fractions_m6anet_4sU_library_gene_subsampling_prob0.9/',hits_ELIGOS_chr = confirmed_hits_ELIGOS_min05_min05_mag1[[1]], hits_ELIGOS_nucleo = confirmed_hits_ELIGOS_min05_min05_mag1[[2]], hits_ELIGOS_cyto = confirmed_hits_ELIGOS_min05_min05_mag1[[3]],
                                             hits_ELIGOS_chr_DRACH = confirmed_hits_DRACH_ELIGOS_min05_min05_mag1[[1]], hits_ELIGOS_nucleo_DRACH = confirmed_hits_DRACH_ELIGOS_min05_min05_mag1[[2]], hits_ELIGOS_cyto_DRACH = confirmed_hits_DRACH_ELIGOS_min05_min05_mag1[[3]],
                                             hits_m6Anet_chr = confirmed_hits_m6Anet_prob0.9[[1]], hits_m6Anet_nucleo = confirmed_hits_m6Anet_prob0.9[[2]], hits_m6Anet_cyto = confirmed_hits_m6Anet_prob0.9[[3]],
                                             hits_m6Anet_chr_DRACH = confirmed_hits_DRACH_m6Anet_prob0.9[[1]], hits_m6Anet_nucleo_DRACH = confirmed_hits_DRACH_m6Anet_prob0.9[[2]], hits_m6Anet_cyto_DRACH = confirmed_hits_DRACH_m6Anet_prob0.9[[3]],
                                             path_pdf_overlap_ELIGOS_m6Anet_DRACH_chr = '/path/to/fractions_m6anet_4sU_library_gene_subsampling_prob0.9/overlap_all_hits_eligos_hits_eligos_DRACH_all_hits_m6anet_chr_ass_10nt.pdf',
                                             path_pdf_overlap_ELIGOS_m6Anet_DRACH_nucleo = '/path/to/fractions_m6anet_4sU_library_gene_subsampling_prob0.9/overlap_all_hits_eligos_hits_eligos_DRACH_all_hits_m6anet_nucleo_10nt.pdf',
                                             path_pdf_overlap_ELIGOS_m6Anet_DRACH_cyto = '/path/to/fractions_m6anet_4sU_library_gene_subsampling_prob0.9/overlap_all_hits_eligos_hits_eligos_DRACH_all_hits_m6anet_cyto_10nt.pdf')

#######################
# ELIGOS ANALYSIS ON NASCENT+PRE-EXISTING READS from replicate 1 (pvalue<0.05, ad.pval<0.05, OR>1)

confirmed_hits_ELIGOS_rep1_min05_min05_mag1 <- ELIGOS_results(path_directory = '/path/to/fractions_eligos_4sU_library_gene_subsampling_rep1_min05_min05_mag1/',
                                                              name_pdf_overlap_10samplings_cyto = 'overlap_cyto_replicates_4sU_10nt_rep1.pdf', name_pdf_overlap_10samplings_chr_ass = 'overlap_chr_ass_replicates_4sU_10nt_rep1.pdf', name_pdf_overlap_10samplings_nucleo = 'overlap_nucleo_replicates_4sU_10nt_rep1.pdf',
                                                              name_pdf_histogram_cyto = 'number_rep_with_hit_cyto_4sU_rep1.pdf', name_pdf_histogram_chr_ass = 'number_rep_with_hit_chr_ass_4sU_rep1.pdf', name_pdf_histogram_nucleo = 'number_rep_with_hit_nucleo_4sU_rep1.pdf',
                                                              name_pdf_overlap_3fractions = 'overlap_hits_eligos_confirmed_5_4sU_10nt_rep1.pdf', p=0.05, ap=0.05, OR=1)

confirmed_hits_DRACH_ELIGOS_rep1_min05_min05_mag1 <- DRACH_overlap_ELIGOS(path_directory = '/path/to/fractions_eligos_4sU_library_gene_subsampling_rep1_min05_min05_mag1/',
                                                                          hits_ELIGOS_chr = confirmed_hits_ELIGOS_rep1_min05_min05_mag1[[1]], hits_ELIGOS_nucleo = confirmed_hits_ELIGOS_rep1_min05_min05_mag1[[2]], hits_ELIGOS_cyto = confirmed_hits_ELIGOS_rep1_min05_min05_mag1[[3]])

#######################
# ELIGOS ANALYSIS ON NASCENT+PRE-EXISTING READS from replicate 2 (pvalue<0.05, ad.pval<0.05, OR>1)

confirmed_hits_ELIGOS_rep2_min05_min05_mag1 <- ELIGOS_results(path_directory = '/path/to/fractions_eligos_4sU_library_gene_subsampling_rep2_min05_min05_mag1/',
                                                              name_pdf_overlap_10samplings_cyto = 'overlap_cyto_replicates_4sU_10nt_rep2.pdf', name_pdf_overlap_10samplings_chr_ass = 'overlap_chr_ass_replicates_4sU_10nt_rep2.pdf', name_pdf_overlap_10samplings_nucleo = 'overlap_nucleo_replicates_4sU_10nt_rep2.pdf',
                                                              name_pdf_histogram_cyto = 'number_rep_with_hit_cyto_4sU_rep2.pdf', name_pdf_histogram_chr_ass = 'number_rep_with_hit_chr_ass_4sU_rep2.pdf', name_pdf_histogram_nucleo = 'number_rep_with_hit_nucleo_4sU_rep2.pdf',
                                                              name_pdf_overlap_3fractions = 'overlap_hits_eligos_confirmed_5_4sU_10nt_rep2.pdf', p=0.05, ap=0.05, OR=1)

confirmed_hits_DRACH_ELIGOS_rep2_min05_min05_mag1 <- DRACH_overlap_ELIGOS(path_directory = '/path/to/fractions_eligos_4sU_library_gene_subsampling_rep2_min05_min05_mag1/',
                                                                          hits_ELIGOS_chr = confirmed_hits_ELIGOS_rep2_min05_min05_mag1[[1]], hits_ELIGOS_nucleo = confirmed_hits_ELIGOS_rep2_min05_min05_mag1[[2]], hits_ELIGOS_cyto = confirmed_hits_ELIGOS_rep2_min05_min05_mag1[[3]])


########################
# ELIGOS ANALYSIS ON NASCENT READS (pvalue<0.05, ad.pval<0.05, OR>1)

confirmed_hits_ELIGOS_nascent_min05_min05_mag1 <- ELIGOS_results(path_directory = '/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/',
                                                                 name_pdf_overlap_10samplings_cyto = 'overlap_cyto_replicates_4sU_10nt.pdf', name_pdf_overlap_10samplings_chr_ass = 'overlap_chr_ass_replicates_4sU_10nt.pdf', name_pdf_overlap_10samplings_nucleo = 'overlap_nucleo_replicates_4sU_10nt.pdf',
                                                                 name_pdf_histogram_cyto = 'number_rep_with_hit_cyto_4sU.pdf', name_pdf_histogram_chr_ass = 'number_rep_with_hit_chr_ass_4sU.pdf', name_pdf_histogram_nucleo = 'number_rep_with_hit_nucleo_4sU.pdf',
                                                                 name_pdf_overlap_3fractions = 'overlap_hits_eligos_confirmed_5_4sU_10nt_nascent.pdf', p=0.05, ap=0.05, OR=1)

confirmed_hits_DRACH_ELIGOS_nascent_min05_min05_mag1 <- DRACH_overlap_ELIGOS(path_directory = '/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/',
                                                                             hits_ELIGOS_chr = confirmed_hits_ELIGOS_nascent_min05_min05_mag1[[1]], hits_ELIGOS_nucleo = confirmed_hits_ELIGOS_nascent_min05_min05_mag1[[2]], hits_ELIGOS_cyto = confirmed_hits_ELIGOS_nascent_min05_min05_mag1[[3]])

# identify which ELIGOS hits (DRACH+ and DRACH-) from nascent reads are only in Chromatin fraction
# create a folder /only_chr/ inside the directory containing the results of the analysis on nascent reads in which you will save
# the RData objects with ELIGOS hits only present in chromatin fraction.
# inside /only_chr/ create /with_DRACH/ to save ELIGOS DRACH+ hits only present in chromatin fraction and 
# /without_DRACH/ to save ELIGOS DRACH- hits only present in chromatin fraction
only_chr <- confirmed_hits_DRACH_ELIGOS_nascent_min05_min05_mag1[[1]][-unique(subjectHits(findOverlaps(confirmed_hits_DRACH_ELIGOS_nascent_min05_min05_mag1[[2]],confirmed_hits_DRACH_ELIGOS_nascent_min05_min05_mag1[[1]], type='any')))]
only_chr <- only_chr[-unique(queryHits(findOverlaps(only_chr, confirmed_hits_DRACH_ELIGOS_nascent_min05_min05_mag1[[3]],type='any')))]
save(only_chr, file='/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/only_chr/with_DRACH/hits_eligos_only_chr_with_DRACH.Rda')

load('/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/without_DRACH/hits_eligos_chr_ass_confirmed_5_without_DRACH.Rda')
load('/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/without_DRACH/hits_eligos_nucleo_confirmed_5_without_DRACH.Rda')
load('/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/without_DRACH/hits_eligos_cyto_confirmed_5_without_DRACH.Rda')
only_chr <- hits_eligos_chr_ass_confirmed_5_without_DRACH[-unique(subjectHits(findOverlaps(hits_eligos_nucleo_confirmed_5_without_DRACH,hits_eligos_chr_ass_confirmed_5_without_DRACH, type='any')))]
only_chr <- only_chr[-unique(queryHits(findOverlaps(only_chr, hits_eligos_cyto_confirmed_5_without_DRACH, type='any')))]
save(only_chr, file='/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/only_chr/without_DRACH/hits_eligos_only_chr_without_DRACH.Rda')

######################## 
# ELIGOS ANALYSIS ON NASCENT+PRE-EXISTING READS with library-level subsampling of nascent ELIGOS (pvalue<0.05, ad.pval<0.05, or>1)

confirmed_hits_ELIGOS_2_min05_min05_mag1 <- ELIGOS_results(path_directory = '/path/to/fractions_eligos_4sU_library_gene_subsampling_total_510645_min05_min05_mag1/',
                                                           name_pdf_overlap_10samplings_cyto = 'overlap_cyto_replicates_4sU_10nt.pdf', name_pdf_overlap_10samplings_chr_ass = 'overlap_chr_ass_replicates_4sU_10nt.pdf', name_pdf_overlap_10samplings_nucleo = 'overlap_nucleo_replicates_4sU_10nt.pdf',
                                                           name_pdf_histogram_cyto = 'number_rep_with_hit_cyto_4sU.pdf', name_pdf_histogram_chr_ass = 'number_rep_with_hit_chr_ass_4sU.pdf', name_pdf_histogram_nucleo = 'number_rep_with_hit_nucleo_4sU.pdf',
                                                           name_pdf_overlap_3fractions = 'overlap_hits_eligos_confirmed_5_4sU_10nt.pdf', p=0.05, ap=0.05, OR=1)

confirmed_hits_DRACH_ELIGOS_2_min05_min05_mag1 <- DRACH_overlap_ELIGOS(path_directory = '/path/to/fractions_eligos_4sU_library_gene_subsampling_total_510645_min05_min05_mag1/',
                                                                       hits_ELIGOS_chr = confirmed_hits_ELIGOS_2_min05_min05_mag1[[1]], hits_ELIGOS_nucleo = confirmed_hits_ELIGOS_2_min05_min05_mag1[[2]], hits_ELIGOS_cyto = confirmed_hits_ELIGOS_2_min05_min05_mag1[[3]])
