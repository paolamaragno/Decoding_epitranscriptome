# Plot the distribution of the median and max probability of modification of m6Anet analysed sites
# overlapping with each ELIGOS DRACH+ hit. ELIGOS DRACH+ hits are defined in different ways:
# p-value ≤ 0.05, adjusted p-value ≤ 0.05, odds-ratio ≥ 1;
# p-value ≤ 0.05, adjusted p-value ≤ 0.05, odds-ratio < 1;
# p-value ≤ 0.05, adjusted p-value > 0.05, odds-ratio ≥ 1;
# p-value ≤ 0.05, adjusted p-value > 0.05, odds-ratio < 1.


library('GenomicFeatures')
library('GenomicRanges')

gtf_file <- "/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/references/Homo_sapiens.GRCh38.104.gtf"
txdb <- makeTxDbFromGFF(gtf_file)
genes_txdb <- GenomicFeatures::genes(txdb)

# path_directory is the directory containing three folders: /chr/, /nucleo/, /cyto/ each with the 10 txt files 
# produced by ELIGOS for each of the 10 samplings
ELIGOS_results_min05_min_05_min1 <- function(path_directory, p, ap, OR) {
  
  # load the txt files with ELIGOS output for all the 10 samplings of each fraction
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
      gr <- gr[(gr$pval <= p) & (gr$pvalAdj <= ap) & (gr$oddR < OR)]
      
      list_gr <- c(list_gr, gr)
    }
    return(list_gr)
  }
  
  gr_eligos_chr_ass <- eligos_output(eligos_chr_ass, 2)
  gr_eligos_nucleo <- eligos_output(eligos_nucleo, 2)
  gr_eligos_cyto <- eligos_output(eligos_cyto, 2)
  
  gr <- list(gr_eligos_chr_ass,gr_eligos_nucleo,gr_eligos_cyto)
  return(gr)
}

# path_directory is the directory containing three folders: /chr/, /nucleo/, /cyto/ each with the 10 txt files 
# produced by ELIGOS for each of the 10 samplings
ELIGOS_results_min05_mag05_mag1 <- function(path_directory, p, ap, OR) {
  
  # load the txt files with ELIGOS output for all the 10 samplings of each fraction
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
      gr <- gr[(gr$pval <= p) & (gr$pvalAdj > ap) & (gr$oddR >= OR)]
      
      list_gr <- c(list_gr, gr)
    }
    return(list_gr)
  }
  gr_eligos_chr_ass <- eligos_output(eligos_chr_ass, 2)
  gr_eligos_nucleo <- eligos_output(eligos_nucleo, 2)
  gr_eligos_cyto <- eligos_output(eligos_cyto, 2)
  
  gr <- list(gr_eligos_chr_ass,gr_eligos_nucleo,gr_eligos_cyto)
  return(gr)
}

# path_directory is the directory containing three folders: /chr/, /nucleo/, /cyto/ each with the 10 txt files 
# produced by ELIGOS for each of the 10 samplings
ELIGOS_results_min05_mag05_min1 <- function(path_directory, p, ap, OR) {
  
  # load the txt files with ELIGOS output for all the 10 samplings of each fraction
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
      gr <- gr[(gr$pval <= p) & (gr$pvalAdj > ap) & (gr$oddR < OR)]
      
      list_gr <- c(list_gr, gr)
    }
    return(list_gr)
  }
  
  gr_eligos_chr_ass <- eligos_output(eligos_chr_ass, 2)
  gr_eligos_nucleo <- eligos_output(eligos_nucleo, 2)
  gr_eligos_cyto <- eligos_output(eligos_cyto, 2)
  
  gr <- list(gr_eligos_chr_ass,gr_eligos_nucleo,gr_eligos_cyto)
  return(gr)
}

# create a folder /hits_ELIGOS/ inside path_directory in which you will save the hits of each fraction 
# confirmed by at least 5 samplings
ELIGOS_results <- function(path_directory, gr) {
  
  gr_eligos_chr_ass <- gr[[1]]
  # assign each hit to the gene on which it maps (I did this step again since there are errors in ELIGOS results)
  for (i in 1:length(gr_eligos_chr_ass)) {
    gene_ids <- c()
    for (j in 1:length(gr_eligos_chr_ass[[i]])){
      over <- findOverlaps(gr_eligos_chr_ass[[i]][j],genes_txdb)
      gene_ids <- c(gene_ids, genes_txdb[subjectHits(over)[1]]$gene_id)
    }
    mcols(gr_eligos_chr_ass[[i]]) <- cbind(mcols(gr_eligos_chr_ass[[i]]), gene_id = gene_ids)
  }
  
  gr_eligos_nucleo <- gr[[2]]
  for (i in 1:length(gr_eligos_nucleo)) {
    gene_ids <- c()
    for (j in 1:length(gr_eligos_nucleo[[i]])){
      over <- findOverlaps(gr_eligos_nucleo[[i]][j],genes_txdb)
      gene_ids <- c(gene_ids, genes_txdb[subjectHits(over)[1]]$gene_id)
    }
    mcols(gr_eligos_nucleo[[i]]) <- cbind(mcols(gr_eligos_nucleo[[i]]), gene_id = gene_ids)
  }
  
  gr_eligos_cyto <- gr[[3]]
  for (i in 1:length(gr_eligos_cyto)) {
    gene_ids <- c()
    for (j in 1:length(gr_eligos_cyto[[i]])){
      over <- findOverlaps(gr_eligos_cyto[[i]][j],genes_txdb)
      gene_ids <- c(gene_ids, genes_txdb[subjectHits(over)[1]]$gene_id)
    }
    mcols(gr_eligos_cyto[[i]]) <- cbind(mcols(gr_eligos_cyto[[i]]), gene_id = gene_ids)
  }
  
  gr_eligos_chr_ass <- unlist(gr_eligos_chr_ass)
  gr_eligos_nucleo <- unlist(gr_eligos_nucleo)
  gr_eligos_cyto <- unlist(gr_eligos_cyto)
  
  eligos_chr_ass_all_replicate_10nt <- c(gr_eligos_chr_ass[[1]],gr_eligos_chr_ass[[2]],gr_eligos_chr_ass[[3]],gr_eligos_chr_ass[[4]],gr_eligos_chr_ass[[5]],
                                         gr_eligos_chr_ass[[6]],gr_eligos_chr_ass[[7]],gr_eligos_chr_ass[[8]],gr_eligos_chr_ass[[9]],gr_eligos_chr_ass[[10]])
  
  confirmed_by_chr <- data.frame(num_samplings=seq(1,10), num_hits=rep(0,10))
  confirmed_by_at_least_5_chr <- c()
  
  # identify the hits present in at least 5 samplings (any type of overlap)
  for (r in 1:length(eligos_chr_ass_all_replicate_10nt)) {
    gr_r_rep <- queryHits(findOverlaps(eligos_chr_ass_all_replicate_10nt, eligos_chr_ass_all_replicate_10nt[r], type='any')) 
    confirmed_by_chr[length(unique(eligos_chr_ass_all_replicate_10nt[gr_r_rep]$rep)),2] <- confirmed_by_chr[length(unique(eligos_chr_ass_all_replicate_10nt[gr_r_rep]$rep)),2] +1
    if (length(unique(eligos_chr_ass_all_replicate_10nt[gr_r_rep]$rep)) >= 5) {
      confirmed_by_at_least_5_chr <- c(confirmed_by_at_least_5_chr, r)
    }
  }
  
  # the hits confirmed in at least 5 samplings are resized to the original coordinates returned by ELIGOS 
  hits_eligos_chr_ass_confirmed_5 <- resize(eligos_chr_ass_all_replicate_10nt[confirmed_by_at_least_5_chr], 2, fix = 'center')
  hits_eligos_chr_ass_confirmed_5 <- reduce(hits_eligos_chr_ass_confirmed_5, ignore.strand = FALSE)
  print(table(width(hits_eligos_chr_ass_confirmed_5)))
  # resize the merged hits to have a width at least of 10 nucleotides
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
  
  eligos_nucleo_all_replicate_10nt <- c(gr_eligos_nucleo[[1]],gr_eligos_nucleo[[2]],gr_eligos_nucleo[[3]],gr_eligos_nucleo[[4]],gr_eligos_nucleo[[5]],
                                        gr_eligos_nucleo[[6]],gr_eligos_nucleo[[7]],gr_eligos_nucleo[[8]],gr_eligos_nucleo[[9]],gr_eligos_nucleo[[10]])
  
  confirmed_by_nucleo <- data.frame(num_samplings=seq(1,10), num_hits=rep(0,10))
  confirmed_by_at_least_5_nucleo <- c()
  
  for (r in 1:length(eligos_nucleo_all_replicate_10nt)) {
    gr_r_rep <- queryHits(findOverlaps(eligos_nucleo_all_replicate_10nt, eligos_nucleo_all_replicate_10nt[r], type='any')) 
    confirmed_by_nucleo[length(unique(eligos_nucleo_all_replicate_10nt[gr_r_rep]$rep)),2] <- confirmed_by_nucleo[length(unique(eligos_nucleo_all_replicate_10nt[gr_r_rep]$rep)),2] +1
    if (length(unique(eligos_nucleo_all_replicate_10nt[gr_r_rep]$rep)) >= 5) {
      confirmed_by_at_least_5_nucleo <- c(confirmed_by_at_least_5_nucleo, r)
    }
  }
  
  hits_eligos_nucleo_confirmed_5 <- resize(eligos_nucleo_all_replicate_10nt[confirmed_by_at_least_5_nucleo], 2, fix = 'center')
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
  
  eligos_cyto_all_replicate_10nt <- c(gr_eligos_cyto[[1]],gr_eligos_cyto[[2]],gr_eligos_cyto[[3]],gr_eligos_cyto[[4]],gr_eligos_cyto[[5]],
                                      gr_eligos_cyto[[6]],gr_eligos_cyto[[7]],gr_eligos_cyto[[8]],gr_eligos_cyto[[9]],gr_eligos_cyto[[10]])
  
  confirmed_by_cyto <- data.frame(num_samplings=seq(1,10), num_hits=rep(0,10))
  confirmed_by_at_least_5_cyto <- c()
  
  for (r in 1:length(eligos_cyto_all_replicate_10nt)) {
    gr_r_rep <- queryHits(findOverlaps(eligos_cyto_all_replicate_10nt, eligos_cyto_all_replicate_10nt[r], type='any')) 
    confirmed_by_cyto[length(unique(eligos_cyto_all_replicate_10nt[gr_r_rep]$rep)),2] <- confirmed_by_cyto[length(unique(eligos_cyto_all_replicate_10nt[gr_r_rep]$rep)),2] +1
    if (length(unique(eligos_cyto_all_replicate_10nt[gr_r_rep]$rep)) >= 5) {
      confirmed_by_at_least_5_cyto <- c(confirmed_by_at_least_5_cyto, r)
    }
  }
  
  hits_eligos_cyto_confirmed_5 <- resize(eligos_cyto_all_replicate_10nt[confirmed_by_at_least_5_cyto], 2, fix = 'center')
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
  
  save(hits_eligos_chr_ass_confirmed_5, file=paste0(path_directory,'/hits_ELIGOS/hits_eligos_chr_ass_confirmed_5.Rda'))
  save(hits_eligos_nucleo_confirmed_5, file=paste0(path_directory,'/hits_ELIGOS/hits_eligos_nucleo_confirmed_5.Rda'))
  save(hits_eligos_cyto_confirmed_5, file=paste0(path_directory,'/hits_ELIGOS/hits_eligos_cyto_confirmed_5.Rda'))
  
  hits_eligos_confirmed_5 <- list(hits_eligos_chr_ass_confirmed_5,hits_eligos_nucleo_confirmed_5,hits_eligos_cyto_confirmed_5)
  return(hits_eligos_confirmed_5)
}

DRACH_overlap_ELIGOS <- function(path_directory, hits_ELIGOS_chr, hits_ELIGOS_nucleo, hits_ELIGOS_cyto) {
  
  load('/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/R_data/DRACH_forward_strand.Rda')
  overlap_DRACH_chr <- findOverlaps(hits_ELIGOS_chr, DRACH, type='any', ignore.strand=FALSE)
  hits_eligos_chr_ass_confirmed_5_with_DRACH <- hits_ELIGOS_chr[unique(queryHits(overlap_DRACH_chr))]
  save(hits_eligos_chr_ass_confirmed_5_with_DRACH, file=paste0(path_directory,'/hits_ELIGOS/hits_eligos_chr_ass_confirmed_5_with_DRACH.Rda'))
  
  overlap_DRACH_nucleo <- findOverlaps(hits_ELIGOS_nucleo, DRACH, type='any', ignore.strand=FALSE)
  hits_eligos_nucleo_confirmed_5_with_DRACH <- hits_ELIGOS_nucleo[unique(queryHits(overlap_DRACH_nucleo))]
  save(hits_eligos_nucleo_confirmed_5_with_DRACH, file=paste0(path_directory,'/hits_ELIGOS/hits_eligos_nucleo_confirmed_5_with_DRACH.Rda'))
  
  overlap_DRACH_cyto <- findOverlaps(hits_ELIGOS_cyto, DRACH, type='any', ignore.strand=FALSE)
  hits_eligos_cyto_confirmed_5_with_DRACH <- hits_ELIGOS_cyto[unique(queryHits(overlap_DRACH_cyto))]
  save(hits_eligos_cyto_confirmed_5_with_DRACH, file=paste0(path_directory,'/hits_ELIGOS/hits_eligos_cyto_confirmed_5_with_DRACH.Rda'))
}

####################

gr_ELIGOS_min05_min05_min1 <- ELIGOS_results_min05_min_05_min1(path_directory = '/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_eligos_4sU_library_gene_subsampling_min05_min05_min1/', p=0.05,ap=0.05,OR=1)
confirmed_hits_ELIGOS_min05_min05_min1 <- ELIGOS_results(path_directory = '/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_eligos_4sU_library_gene_subsampling_min05_min05_min1/',gr_ELIGOS_min05_min05_min1)
confirmed_hits_DRACH_ELIGOS_min05_min05_min1 <- DRACH_overlap_ELIGOS(path_directory = '/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_eligos_4sU_library_gene_subsampling_min05_min05_min1/',
                                                                     hits_ELIGOS_chr = confirmed_hits_ELIGOS_min05_min05_min1[[1]], hits_ELIGOS_nucleo = confirmed_hits_ELIGOS_min05_min05_min1[[2]], hits_ELIGOS_cyto = confirmed_hits_ELIGOS_min05_min05_min1[[3]])

gr_ELIGOS_min05_mag05_mag1 <- ELIGOS_results_min05_mag05_mag1(path_directory = '/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_eligos_4sU_library_gene_subsampling_min05_mag05_mag1/', p=0.05,ap=0.05,OR=1)
confirmed_hits_ELIGOS_min05_mag05_mag1 <- ELIGOS_results(path_directory = '/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_eligos_4sU_library_gene_subsampling_min05_mag05_mag1/', gr_ELIGOS_min05_mag05_mag1)
confirmed_hits_DRACH_ELIGOS_min05_mag05_mag1 <- DRACH_overlap_ELIGOS(path_directory = '/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_eligos_4sU_library_gene_subsampling_min05_mag05_mag1/',
                                                                     hits_ELIGOS_chr = confirmed_hits_ELIGOS_min05_mag05_mag1[[1]], hits_ELIGOS_nucleo = confirmed_hits_ELIGOS_min05_mag05_mag1[[2]], hits_ELIGOS_cyto = confirmed_hits_ELIGOS_min05_mag05_mag1[[3]])

gr_ELIGOS_min05_mag05_min1 <- ELIGOS_results_min05_mag05_min1(path_directory = '/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_eligos_4sU_library_gene_subsampling_min05_mag05_min1/', p=0.05,ap=0.05,OR=1)
confirmed_hits_ELIGOS_min05_mag05_min1 <- ELIGOS_results(path_directory = '/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_eligos_4sU_library_gene_subsampling_min05_mag05_min1/',gr_ELIGOS_min05_mag05_min1)
confirmed_hits_DRACH_ELIGOS_min05_mag05_min1 <- DRACH_overlap_ELIGOS(path_directory = '/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_eligos_4sU_library_gene_subsampling_min05_mag05_min1/',
                                                                     hits_ELIGOS_chr = confirmed_hits_ELIGOS_min05_mag05_min1[[1]], hits_ELIGOS_nucleo = confirmed_hits_ELIGOS_min05_mag05_min1[[2]], hits_ELIGOS_cyto = confirmed_hits_ELIGOS_min05_mag05_min1[[3]])

##############

# load the tsv file with the output of m6Anet for the three fractions for all the 10 samplings
m6Anet_fractions_4sU_chr_ass <- list.files(path = '/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_m6anet_4sU_library_gene_subsampling_prob0.9/chr', pattern = 'tsv', full.names = TRUE)
m6Anet_fractions_4sU_nucleo <- list.files(path = '/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_m6anet_4sU_library_gene_subsampling_prob0.9/nucleo', pattern = 'tsv', full.names = TRUE)
m6Anet_fractions_4sU_cyto <- list.files(path = '/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_m6anet_4sU_library_gene_subsampling_prob0.9/cyto', pattern = 'tsv', full.names = TRUE)

# load the vector with the transcript names as names and the names of the gene from which 
# they are transcribed as values
load('/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/R_data/tx_gene.Rda')
gtf_file <- "/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/references/Homo_sapiens.GRCh38.104.gtf"
txdb <- makeTxDbFromGFF(gtf_file)
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
    
    list_gr <- c(list_gr, gr)
  }
  return(list_gr)
}

gr_m6Anet_chr_ass_4sU <- m6Anet_output(m6Anet_fractions_4sU_chr_ass)
gr_m6Anet_chr_ass_4sU <- unlist(as(gr_m6Anet_chr_ass_4sU, 'GRangesList'))
gr_m6Anet_nucleo_4sU <- m6Anet_output(m6Anet_fractions_4sU_nucleo)
gr_m6Anet_nucleo_4sU <- unlist(as(gr_m6Anet_nucleo_4sU, 'GRangesList'))
gr_m6Anet_cyto_4sU <- m6Anet_output(m6Anet_fractions_4sU_cyto)
gr_m6Anet_cyto_4sU <- unlist(as(gr_m6Anet_cyto_4sU, 'GRangesList'))

# remove the hits mapping on transcripts that have been assigned to the wrong chromosome by m6Anet
rem <- c()
for (i in 1:length(gr_m6Anet_chr_ass_4sU)) {
  tx <- gr_m6Anet_chr_ass_4sU[i]$tx_id
  chr <- as.character(seqnames(gr_m6Anet_chr_ass_4sU[i]))
  if (chr != as.character(seqnames(tx_txdb[tx_txdb$tx_name == tx]))) {
    rem <- c(rem, i)
  }
}
gr_m6Anet_chr_ass_4sU <- gr_m6Anet_chr_ass_4sU[-rem]

rem <- c()
for (i in 1:length(gr_m6Anet_nucleo_4sU)) {
  tx <- gr_m6Anet_nucleo_4sU[i]$tx_id
  chr <- as.character(seqnames(gr_m6Anet_nucleo_4sU[i]))
  if (chr != as.character(seqnames(tx_txdb[tx_txdb$tx_name == tx]))) {
    rem <- c(rem, i)
  }
}
gr_m6Anet_nucleo_4sU <- gr_m6Anet_nucleo_4sU[-rem]

rem <- c()
for (i in 1:length(gr_m6Anet_cyto_4sU)) {
  tx <- gr_m6Anet_cyto_4sU[i]$tx_id
  chr <- as.character(seqnames(gr_m6Anet_cyto_4sU[i]))
  if (chr != as.character(seqnames(tx_txdb[tx_txdb$tx_name == tx]))) {
    rem <- c(rem, i)
  }
}
gr_m6Anet_cyto_4sU <- gr_m6Anet_cyto_4sU[-rem]

# all the hits that map on the same chromosome, that have the same strand, same start and end positions and same gene id but 
# different transcript id are kept only once
gr_unique <- function(x) {
  ind_unique <- which(isUnique(paste0(seqnames(x), "_", start(x), "_", end(x), "_", strand(x),"_", mcols(x)$gene_id)))
  dup <- unique(x[-ind_unique])
  x_unique <- sort(c(x[ind_unique], dup))
  return(x_unique)
}

gr_m6Anet_chr_ass_4sU <- gr_unique(gr_m6Anet_chr_ass_4sU)
gr_m6Anet_nucleo_4sU <- gr_unique(gr_m6Anet_nucleo_4sU)
gr_m6Anet_cyto_4sU <- gr_unique(gr_m6Anet_cyto_4sU)

# path_directory is the directory containing three folders: /chr/, /nucleo/, /cyto/ each with the 10 txt files 
# produced by ELIGOS for each of the 10 samplings
# create a folder /plot_prob_distribution_m6Anet/ inside path_directory where you will save the histrograms
# with the distribution of the median and maximum probability of modification of m6Anet analysed sites 
# overlapping with each ELIGOS DRACH+ hit
prob_distribution <- function(path_directory) {
  load(paste0(path_directory, '/hits_ELIGOS/hits_eligos_chr_ass_confirmed_5_with_DRACH.Rda'))
  load(paste0(path_directory, '/hits_ELIGOS/hits_eligos_nucleo_confirmed_5_with_DRACH.Rda'))
  load(paste0(path_directory, '/hits_ELIGOS/hits_eligos_cyto_confirmed_5_with_DRACH.Rda'))
  
  prob_median_chr <- c()
  prob_max_chr <- c()
  
  # identify which sites analysed by m6Anet overlap with each ELIGOS DRACH+ hit and compute the median and maximum probability of modification
  # of this set of m6Anet sites
  lapply(seq_along(hits_eligos_chr_ass_confirmed_5_with_DRACH), function(i,x) {
    m6Anet_over_prob <- gr_m6Anet_chr_ass_4sU[unique(subjectHits(findOverlaps(x[i],gr_m6Anet_chr_ass_4sU)))]$prob_mod
    if (length(m6Anet_over_prob) != 0) {
      prob_median_chr <<- c(prob_median_chr, median(m6Anet_over_prob))
      prob_max_chr <<- c(prob_max_chr, max(m6Anet_over_prob))
    }
  }, x = hits_eligos_chr_ass_confirmed_5_with_DRACH)
  
  jpeg(paste0(path_directory, '/plot_prob_distribution_m6Anet/median_chr.jpeg'))
  hist(prob_median_chr, breaks=100, main = 'Median probability of m6Anet hits overlapping\nwith ELIGOS DRACH+ hits - C.Associated',
       xlab = 'Median probability')
  dev.off()
  
  jpeg(paste0(path_directory, '/plot_prob_distribution_m6Anet/max_chr.jpeg'))
  hist(prob_max_chr, breaks=100, main = 'Max probability of m6Anet hits overlapping\nwith ELIGOS DRACH+ hits - C.Associated',
       xlab = 'Max probability')
  dev.off()
  
  prob_median_nucleo <- c()
  prob_max_nucleo <- c()
  
  lapply(seq_along(hits_eligos_nucleo_confirmed_5_with_DRACH), function(i,x) {
    m6Anet_over_prob <- gr_m6Anet_nucleo_4sU[unique(subjectHits(findOverlaps(x[i],gr_m6Anet_nucleo_4sU)))]$prob_mod
    if (length(m6Anet_over_prob) != 0) {
      prob_median_nucleo <<- c(prob_median_nucleo, median(m6Anet_over_prob))
      prob_max_nucleo <<- c(prob_max_nucleo, max(m6Anet_over_prob))
    }
  }, x = hits_eligos_nucleo_confirmed_5_with_DRACH)
  
  jpeg(paste0(path_directory,'/plot_prob_distribution_m6Anet/median_nucleo.jpeg'))
  hist(prob_median_nucleo, breaks=100, main = 'Median probability of m6Anet hits overlapping\nwith ELIGOS DRACH+ hits - Nucleoplasmic',
       xlab = 'Median probability')
  dev.off()
  
  jpeg(paste0(path_directory,'/plot_prob_distribution_m6Anet/max_nucleo.jpeg'))
  hist(prob_max_nucleo, breaks=100, main = 'Max probability of m6Anet hits overlapping\nwith ELIGOS DRACH+ hits - Nucleoplasmic',
       xlab = 'Max probability')
  dev.off()
  
  prob_median_cyto <- c()
  prob_max_cyto <- c()
  
  lapply(seq_along(hits_eligos_cyto_confirmed_5_with_DRACH), function(i,x) {
    m6Anet_over_prob <- gr_m6Anet_cyto_4sU[unique(subjectHits(findOverlaps(x[i],gr_m6Anet_cyto_4sU)))]$prob_mod
    if (length(m6Anet_over_prob) != 0) {
      prob_median_cyto <<- c(prob_median_cyto, median(m6Anet_over_prob))
      prob_max_cyto <<- c(prob_max_cyto, max(m6Anet_over_prob))
    }
  }, x = hits_eligos_cyto_confirmed_5_with_DRACH)
  
  jpeg(paste0(path_directory,'/plot_prob_distribution_m6Anet/median_cyto.jpeg'))
  hist(prob_median_cyto, breaks=100, main = 'Median probability of m6Anet hits overlapping\nwith ELIGOS DRACH+ hits - Cytoplasmic',
       xlab = 'Median probability')
  dev.off()
  
  jpeg(paste0(path_directory,'/plot_prob_distribution_m6Anet/max_cyto.jpeg'))
  hist(prob_max_cyto, breaks=100, main = 'Max probability of m6Anet hits overlapping\nwith ELIGOS DRACH+ hits - Cytoplasmic',
       xlab = 'Max probability')
  dev.off()
}

prob_distribution('/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/')
prob_distribution('/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_eligos_4sU_library_gene_subsampling_min05_min05_min1/')
prob_distribution('/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_eligos_4sU_library_gene_subsampling_min05_mag05_mag1/')
prob_distribution('/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_eligos_4sU_library_gene_subsampling_min05_mag05_min1/')

