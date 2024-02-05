# perform the overlap between ELIGOS DRACH+ hits not confirmed by m6Anet and the coordinates of the m6A marks
# and of the the binding sites of specific effectors associated with m6A from RMBase3 and RMVar. 
# This allows to compute how many ELIGOS DRACH+ hits non confirmed by m6Anet are m6A sites according to the databases

library('GenomicRanges')
library('GenomicFeatures')
library('xlsx')

# overlap between ELIGOS DRACH+ hits (not confirmed by m6Anet) and m6a marks from the two databases
# directory_hits_not_confirmed is the path to the directory containing ELIGOS DRACH+ hits not confirmed by m6Anet
overlap_marks <- function(hits_chr_grange,hits_nucleo_grange,hits_cyto_grange,directory_hits_not_confirmed) {
  
  # initiate a matrix that will report, for each fraction, the number (and percentage) of ELIGOS hits overlapping with m6A sites.
  # The matrix will also report the number of ELIGOS hits containing m6A and shared by all the fractions
  number_hits_per_mod <- matrix(0,nrow=4, ncol=2)
  rownames(number_hits_per_mod) <- c('Chromatin Associated', 'Nucleoplasmic', 'Cytoplasmic', 'Overlap')
  colnames(number_hits_per_mod) <- c('m6A', 'Tot DRACH+ hits')
  
  number_hits_per_mod[1,2] <- as.character(length(hits_chr_grange))
  number_hits_per_mod[2,2] <- as.character(length(hits_nucleo_grange))
  number_hits_per_mod[3,2] <- as.character(length(hits_cyto_grange))
  
  # compute the number of ELIGOS DRACH+ hits (not confirmed by m6Anet) shared by all the fractions
  all_hits <- list(hits_chr_grange, hits_nucleo_grange, hits_cyto_grange)
  order <- order(c(length(hits_chr_grange), length(hits_nucleo_grange), length(hits_cyto_grange)))
  
  overlap <- findOverlaps(all_hits[[order[1]]],all_hits[[order[2]]], type = 'any')
  overlap2 <- all_hits[[order[1]]][unique(queryHits(overlap))]
  overlap_all <- findOverlaps(overlap2,all_hits[[order[3]]], type = 'any')
  number_hits_per_mod[4,2] <- length(overlap2[unique(queryHits(overlap_all))])
  
  mod_file <- read.table('/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/mod/m6A_total.bed')
  
  mod_grange <- GRanges(seqnames = mod_file$chr,
                        ranges = IRanges(start = mod_file$start, end=mod_file$end),
                        strand = Rle(mod_file$strand),
                        mod_type =  mod_file$mod_type)
  
  # overlap between ELIGOS DRACH+ hits of a fraction and the coordinates of m6A from the two databases
  suppressWarnings(overlap_chr <- findOverlaps(hits_chr_grange,mod_grange))
  hits_chr_with_mod <- c()
  if (length(overlap_chr) != 0) {
    hits_chr_with_mod <- hits_chr_grange[unique(queryHits(overlap_chr))]
    number_hits_per_mod[1,'m6A'] <- paste0(as.character(length(unique(queryHits(overlap_chr)))), ' - ', as.character(round(length(unique(queryHits(overlap_chr)))*100/as.numeric(number_hits_per_mod[1,2]),2)),'%')
    # update the metadata field of ELIGOS DRACH+ hits overlapping with m6A
    for (i in unique(queryHits(overlap_chr))) {
      if (hits_chr_grange[i]$mod_type == 'Unknown') {
        hits_chr_grange[i]$mod_type <- 'm6A'
      } 
    }
  }
  
  suppressWarnings(overlap_nucleo <- findOverlaps(hits_nucleo_grange,mod_grange))  
  hits_nucleo_with_mod <- c()
  if (length(overlap_nucleo) != 0) {
    hits_nucleo_with_mod <- hits_nucleo_grange[unique(queryHits(overlap_nucleo))]
    number_hits_per_mod[2,'m6A'] <- paste0(as.character(length(unique(queryHits(overlap_nucleo)))), ' - ', as.character(round(length(unique(queryHits(overlap_nucleo)))*100/as.numeric(number_hits_per_mod[2,2]),2)),'%')
    # update the metadata field of ELIGOS hits overlapping with m6A
    for (i in unique(queryHits(overlap_nucleo))) {
      if (hits_nucleo_grange[i]$mod_type == 'Unknown') {
        hits_nucleo_grange[i]$mod_type <- 'm6A'
      } 
    }
  } 
  
  suppressWarnings(overlap_cyto <- findOverlaps(hits_cyto_grange,mod_grange))
  hits_cyto_with_mod <- c()
  if (length(overlap_cyto) != 0) {
    hits_cyto_with_mod <- hits_cyto_grange[unique(queryHits(overlap_cyto))]
    number_hits_per_mod[3,'m6A'] <- paste0(as.character(length(unique(queryHits(overlap_cyto)))), ' - ', as.character(round(length(unique(queryHits(overlap_cyto)))*100/as.numeric(number_hits_per_mod[3,2]),2)),'%')
    # update the metadata field of ELIGOS hits overlapping with m6A
    for (i in unique(queryHits(overlap_cyto))) {
      if (hits_cyto_grange[i]$mod_type == 'Unknown') {
        hits_cyto_grange[i]$mod_type <- 'm6A'
      } 
    }
  }
  
  # compute how many ELIGOS DRACH+ hits overlapping with m6A are shared by all the fractions
  if (length(hits_chr_with_mod) !=0 & length(hits_nucleo_with_mod) != 0 & length(hits_cyto_with_mod) != 0) {
    hits_with_mod <- list(hits_chr_with_mod, hits_nucleo_with_mod, hits_cyto_with_mod)
    order <- order(c(length(hits_chr_with_mod), length(hits_nucleo_with_mod), length(hits_cyto_with_mod)))
    
    overlap <- findOverlaps(hits_with_mod[[order[1]]],hits_with_mod[[order[2]]], type = 'any')
    overlap2 <- hits_with_mod[[order[1]]][unique(queryHits(overlap))]
    overlap_all <- findOverlaps(overlap2,hits_with_mod[[order[3]]], type = 'any')
    number_hits_per_mod[4,'m6A'] <- length(overlap2[unique(queryHits(overlap_all))]) 
  } 
  
  write.xlsx(x = data.frame(number_hits_per_mod),file = paste0(directory_hits_not_confirmed, '/ELIGOS_DRACH_not_confirmed_m6A_mark.xlsx'),col.names = TRUE, row.names = TRUE)
  
  unique_mod_types <- function(hits) {
    hits$mod_type <- unlist(lapply(seq_along(hits), function(i,x) {
      paste(unique(unlist(strsplit(x[i]$mod_type, split =';'))), collapse = ';')
    }, x=hits))
    return(hits)
  }
  
  hits_chr_grange <- unique_mod_types(hits_chr_grange)
  hits_nucleo_grange <- unique_mod_types(hits_nucleo_grange)
  hits_cyto_grange <- unique_mod_types(hits_cyto_grange)
  
  l <- list(hits_chr_grange,hits_nucleo_grange,hits_cyto_grange)
  names(l) <- c('hits_chr_grange','hits_nucleo_grange','hits_cyto_grange')
  return(l)
}

# overlap between ELIGOS DRACH+ hits (not confirmed by m6Anet) and the binding sites of specific effectors
# associated with m6A from the two databases
overlap_binding_sites <- function(hits, n) {
  
  over <-findOverlaps(hits,all_sites_RMBase_RMvar, ignore.strand=TRUE, type='any')
  hits_overlapping_specific_effector_m6A <- c()
  # report the index of ELIGOS DRACH+ hits overlapping with specific effectors associated with m6A 
  for (i in 1:length(over)) {
    if (all_sites_RMBase_RMvar[subjectHits(over[i])]$RBP_name %in% specific_effectors_only_m6A) {
      hits_overlapping_specific_effector_m6A <- c(hits_overlapping_specific_effector_m6A, queryHits(over[i]))
    }
  }
  number_hits_per_mod_RBP_DRACH[n,1] <<- paste0(as.character(length(unique(hits_overlapping_specific_effector_m6A))), ' - ', as.character(round(length(unique(hits_overlapping_specific_effector_m6A))/length(hits)*100,2)),'%')
  
  # update the field "mod_type" of ELIGOS DRACH+ hits that didn't overlap with m6A mark from the databases
  for (i in 1:length(hits)) {
    if (hits[i]$mod_type == 'Unknown') {
      if (i %in% queryHits(over)) {
        binding_sites <- unique(all_sites_RMBase_RMvar[subjectHits(over[queryHits(over)==i])])
        RBP_specific <- unique(binding_sites$RBP_name[which(binding_sites$RBP_name %in% specific_effectors_only_m6A)])
        hits[i]$effector <- paste(RBP_specific, collapse = ';')
        hits[i]$mod_type <- paste(unique(enzyme_type_mod_unique[enzyme_type_mod_unique$Name_effector %in% RBP_specific,3]), collapse=';')
        hits[i]$type_effector <- paste(unique(enzyme_type_mod_unique[enzyme_type_mod_unique$Name_effector %in% RBP_specific,2]), collapse=';')
      }}}
  
  hits$mod_type <- unlist(lapply(seq_along(hits), function(i,x) {
    paste(unique(unlist(strsplit(x[i]$mod_type, split =';'))), collapse = ';')
  }, x=hits))
  
  l <- list(unique(hits_overlapping_specific_effector_m6A), hits)
  names(l) <- c('specific m6A', 'hits')
  return(l)
}

# combine the information from the overlap with m6A marks and the overlap with the binding sites of specific
# effectors associated with m6A from the databases.
# directory_hits_not_confirmed is the path to the directory containing ELIGOS DRACH+ hits not confirmed by m6Anet
final_summary <- function(directory_hits_not_confirmed) {
  
  load(paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_chr_not_confirmed_m6A_mark_RBPs.Rda'))
  load(paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_nucleo_not_confirmed_m6A_mark_RBPs.Rda'))
  load(paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_cyto_not_confirmed_m6A_mark_RBPs.Rda'))
  
  # initiate a matrix reporting the number (and percentage) of ELIGOS DRACH+ hits annotated with m6A.
  # Report also the number of ELIGOS DRACH+ hits annotated to m6A and shared by all the fractions
  number_hits_DRACH <- matrix(0,nrow=4, ncol=2)
  rownames(number_hits_DRACH) <- c('Chromatin Associated', 'Nucleoplasmic', 'Cytoplasmic','Overlap')
  colnames(number_hits_DRACH) <- c('m6A', 'Tot number of DRACH+ hits')
  
  number_hits_DRACH[1,'Tot number of DRACH+ hits'] <- length(eligos_DRACH_chr_non_confirmed_with_bindings)
  number_hits_DRACH[2,'Tot number of DRACH+ hits'] <- length(eligos_DRACH_nucleo_non_confirmed_with_bindings)
  number_hits_DRACH[3,'Tot number of DRACH+ hits'] <- length(eligos_DRACH_cyto_non_confirmed_with_bindings)
  
  count_mods_DRACH <- function(hits_chr, hits_nucleo, hits_cyto) {
    mod_types_chr <- unlist(strsplit(hits_chr$mod_type, split=';'))
    mod_types_nucleo <- unlist(strsplit(hits_nucleo$mod_type, split=';'))
    mod_types_cyto <- unlist(strsplit(hits_cyto$mod_type, split=';'))
    # for each fraction, compute the number and percentage of ELIGOS DRACH+ hits annotated to m6A
    number_hits_DRACH[1,'m6A'] <<- paste0(as.character(length(mod_types_chr[mod_types_chr=='m6A'])),' - ', as.character(round(length(mod_types_chr[mod_types_chr=='m6A'])*100/as.numeric(number_hits_DRACH[1,'Tot number of DRACH+ hits']),2)),'%')
    number_hits_DRACH[2,'m6A'] <<- paste0(as.character(length(mod_types_nucleo[mod_types_nucleo=='m6A'])),' - ', as.character(round(length(mod_types_nucleo[mod_types_nucleo=='m6A'])*100/as.numeric(number_hits_DRACH[2,'Tot number of DRACH+ hits']),2)),'%')
    number_hits_DRACH[3,'m6A'] <<- paste0(as.character(length(mod_types_cyto[mod_types_cyto=='m6A'])),' - ', as.character(round(length(mod_types_cyto[mod_types_cyto=='m6A'])*100/as.numeric(number_hits_DRACH[3,'Tot number of DRACH+ hits']),2)),'%')
    
    # compute how many ELIGOS DRACH+ hits annotated to m6A are shared by all the fractions
    hits_mod_chr <- hits_chr[unlist(lapply(seq_along(hits_chr),function(i) {'m6A' %in% unlist(strsplit(hits_chr[i]$mod_type, split=';'))}))] 
    hits_mod_nucleo <- hits_nucleo[unlist(lapply(seq_along(hits_nucleo),function(i) {'m6A' %in% unlist(strsplit(hits_nucleo[i]$mod_type, split=';'))}))] 
    hits_mod_cyto <-  hits_cyto[unlist(lapply(seq_along(hits_cyto),function(i) {'m6A' %in% unlist(strsplit(hits_cyto[i]$mod_type, split=';'))}))] 
    hits_mod <- list(hits_mod_chr, hits_mod_nucleo, hits_mod_cyto)
    order <- order(c(length(hits_mod_chr), length(hits_mod_nucleo), length(hits_mod_cyto)))
    
    overlap <- findOverlaps(hits_mod[[order[1]]],hits_mod[[order[2]]], type = 'any')
    overlap2 <- hits_mod[[order[1]]][unique(queryHits(overlap))]
    overlap_all <- findOverlaps(overlap2,hits_mod[[order[3]]], type = 'any')
    number_hits_DRACH[4,'m6A'] <<- length(overlap2[unique(queryHits(overlap_all))])
  }
  
  count_mods_DRACH(eligos_DRACH_chr_non_confirmed_with_bindings,eligos_DRACH_nucleo_non_confirmed_with_bindings,eligos_DRACH_cyto_non_confirmed_with_bindings)
  
  write.xlsx(x = data.frame(number_hits_DRACH),file = paste0(directory_hits_not_confirmed, '/ELIGOS_DRACH_not_confirmed_summary.xlsx'),col.names = TRUE, row.names = TRUE)
}

############
# directory_eligos_DRACH is the path to the directory containing ELIGOS DRACH+ hits not confirmed by m6Anet (p>0.75)
directory_hits_not_confirmed <-  "/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_m6anet_4sU_library_gene_subsampling_prob0.75/ELIGOS_confirmed_DRACH/not/"

load(paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_chr_not_confirmed.Rda'))
load(paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_nucleo_not_confirmed.Rda'))
load(paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_cyto_not_confirmed.Rda'))

# add a metadata column to each GRanges object (initiated with Unknown) in which 'm6A' will be 
# reported only for the hits overlapping with any m6A mark from the databases
mcols(eligos_DRACH_chr_non_confirmed) <- cbind(mcols(eligos_DRACH_chr_non_confirmed), mod_type= 'Unknown')
mcols(eligos_DRACH_nucleo_non_confirmed) <- cbind(mcols(eligos_DRACH_nucleo_non_confirmed), mod_type= 'Unknown')
mcols(eligos_DRACH_cyto_non_confirmed) <- cbind(mcols(eligos_DRACH_cyto_non_confirmed), mod_type= 'Unknown')

hits_DRACH <- overlap_marks(eligos_DRACH_chr_non_confirmed,eligos_DRACH_nucleo_non_confirmed,eligos_DRACH_cyto_non_confirmed,directory_hits_not_confirmed)
eligos_DRACH_chr_non_confirmed <- hits_DRACH[[1]]
eligos_DRACH_nucleo_non_confirmed <- hits_DRACH[[2]]
eligos_DRACH_cyto_non_confirmed <- hits_DRACH[[3]]

save(eligos_DRACH_chr_non_confirmed, file =paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_chr_not_confirmed_m6A_mark.Rda'))
save(eligos_DRACH_nucleo_non_confirmed, file =paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_nucleo_not_confirmed_m6A_mark.Rda'))
save(eligos_DRACH_cyto_non_confirmed, file =paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_cyto_not_confirmed_m6A_mark.Rda'))

# initiate a matrix to report, for each fraction, the percentage of ELIGOS DRACH+ hits (not confirmed by m6Anet) overlapping with specific effectors
# associated with m6A and the number of ELIGOS DRACH+ hits (not confirmed by m6Anet) overlapping with specific effectors associated with m6A and
# shared by all the fractions
number_hits_per_mod_RBP_DRACH <- matrix(0,nrow=4, ncol=2)
rownames(number_hits_per_mod_RBP_DRACH) <- c('Chromatin Associated', 'Nucleoplasmic', 'Cytoplasmic','Overlap')
colnames(number_hits_per_mod_RBP_DRACH) <- c('% of hits overlapping with\nspecific effectors m6A',
                                             'Tot number of DRACH+ hits')
number_hits_per_mod_RBP_DRACH[1,2] <- length(eligos_DRACH_chr_non_confirmed)
number_hits_per_mod_RBP_DRACH[2,2] <-length(eligos_DRACH_nucleo_non_confirmed)
number_hits_per_mod_RBP_DRACH[3,2] <-length(eligos_DRACH_cyto_non_confirmed)

# add a metadata column to each GRanges object in which the names of the specific effectors associated with m6A 
# with which each hit overlaps will be reported as well as the function of these effectors (both initiated with Unknown)
mcols(eligos_DRACH_chr_non_confirmed) <- cbind(mcols(eligos_DRACH_chr_non_confirmed), effector = 'Unknown', type_effector = 'Unknown')
mcols(eligos_DRACH_nucleo_non_confirmed) <- cbind(mcols(eligos_DRACH_nucleo_non_confirmed), effector = 'Unknown', type_effector = 'Unknown')
mcols(eligos_DRACH_cyto_non_confirmed) <- cbind(mcols(eligos_DRACH_cyto_non_confirmed), effector = 'Unknown', type_effector = 'Unknown')

load('/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/RBPs/specific_effectors_only_m6A.Rda')
load('/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/RBPs/all_sites_RMBase_RMvar.Rda')
load('/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/RBPs/enzyme_type_mod_unique.RDa')

results_chr <- overlap_binding_sites(eligos_DRACH_chr_non_confirmed,1)
chr_specific_m6A <- eligos_DRACH_chr_non_confirmed[results_chr[[1]]]
eligos_DRACH_chr_non_confirmed_with_bindings <- results_chr[[2]]

results_nucleo <- overlap_binding_sites(eligos_DRACH_nucleo_non_confirmed,2)
nucleo_specific_m6A <- eligos_DRACH_nucleo_non_confirmed[results_nucleo[[1]]]
eligos_DRACH_nucleo_non_confirmed_with_bindings <- results_nucleo[[2]]

results_cyto <- overlap_binding_sites(eligos_DRACH_cyto_non_confirmed,3)
cyto_specific_m6A <- eligos_DRACH_cyto_non_confirmed[results_cyto[[1]]]
eligos_DRACH_cyto_non_confirmed_with_bindings <- results_cyto[[2]]

# compute the number of ELIGOS DRACH+ hits (not confirmed by m6Anet) overlapping with specific effectors associated
# with m6A and shared by all the fractions 
all_hits_specific_m6A <- list(chr_specific_m6A,nucleo_specific_m6A,cyto_specific_m6A)
order <- order(c(length(chr_specific_m6A), length(nucleo_specific_m6A), length(cyto_specific_m6A)))
overlap <- findOverlaps(all_hits_specific_m6A[[order[1]]],all_hits_specific_m6A[[order[2]]], type = 'any')
overlap2 <- all_hits_specific_m6A[[order[1]]][unique(queryHits(overlap))]
overlap_all <- findOverlaps(overlap2,all_hits_specific_m6A[[order[3]]], type = 'any')
number_hits_per_mod_RBP_DRACH[4,1] <- length(overlap2[unique(queryHits(overlap_all))])

write.xlsx(x = data.frame(number_hits_per_mod_RBP_DRACH),file = paste0(directory_hits_not_confirmed, 'ELIGOS_DRACH_not_confirmed_m6A_mark_sp_effectors.xlsx'),col.names = TRUE, row.names=TRUE)

save(eligos_DRACH_chr_non_confirmed_with_bindings, file =paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_chr_not_confirmed_m6A_mark_RBPs.Rda'))
save(eligos_DRACH_nucleo_non_confirmed_with_bindings, file =paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_nucleo_not_confirmed_m6A_mark_RBPs.Rda'))
save(eligos_DRACH_cyto_non_confirmed_with_bindings, file =paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_cyto_not_confirmed_m6A_mark_RBPs.RDa'))

final_summary(directory_hits_not_confirmed)

#############
# directory_eligos_DRACH is the path to the directory containing ELIGOS DRACH+ hits not confirmed by m6Anet (p>0.9)
directory_hits_not_confirmed <-  "/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_m6anet_4sU_library_gene_subsampling_prob0.9/ELIGOS_confirmed_DRACH/not/"

load(paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_chr_not_confirmed.Rda'))
load(paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_nucleo_not_confirmed.Rda'))
load(paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_cyto_not_confirmed.Rda'))

# add a metadata column to each GRanges object (initiated with Unknown) in which 'm6A' will be 
# reported only for the hits overlapping with any m6A mark from the databases
mcols(eligos_DRACH_chr_non_confirmed) <- cbind(mcols(eligos_DRACH_chr_non_confirmed), mod_type= 'Unknown')
mcols(eligos_DRACH_nucleo_non_confirmed) <- cbind(mcols(eligos_DRACH_nucleo_non_confirmed), mod_type= 'Unknown')
mcols(eligos_DRACH_cyto_non_confirmed) <- cbind(mcols(eligos_DRACH_cyto_non_confirmed), mod_type= 'Unknown')

hits_DRACH <- overlap_marks(eligos_DRACH_chr_non_confirmed,eligos_DRACH_nucleo_non_confirmed,eligos_DRACH_cyto_non_confirmed,directory_hits_not_confirmed)
eligos_DRACH_chr_non_confirmed <- hits_DRACH[[1]]
eligos_DRACH_nucleo_non_confirmed <- hits_DRACH[[2]]
eligos_DRACH_cyto_non_confirmed <- hits_DRACH[[3]]

save(eligos_DRACH_chr_non_confirmed, file =paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_chr_not_confirmed_m6A_mark.Rda'))
save(eligos_DRACH_nucleo_non_confirmed, file =paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_nucleo_not_confirmed_m6A_mark.Rda'))
save(eligos_DRACH_cyto_non_confirmed, file =paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_cyto_not_confirmed_m6A_mark.Rda'))

# initiate a matrix to report, for each fraction, the percentage of ELIGOS DRACH+ hits (not confirmed by m6Anet) overlapping with specific effectors
# associated with m6A and the number of ELIGOS DRACH+ hits (not confirmed by m6Anet) overlapping with specific effectors associated with m6A and
# shared by all the fractions
number_hits_per_mod_RBP_DRACH <- matrix(0,nrow=4, ncol=2)
rownames(number_hits_per_mod_RBP_DRACH) <- c('Chromatin Associated', 'Nucleoplasmic', 'Cytoplasmic','Overlap')
colnames(number_hits_per_mod_RBP_DRACH) <- c('% of hits overlapping with\nspecific effectors m6A',
                                             'Tot number of DRACH+ hits')
number_hits_per_mod_RBP_DRACH[1,2] <- length(eligos_DRACH_chr_non_confirmed)
number_hits_per_mod_RBP_DRACH[2,2] <-length(eligos_DRACH_nucleo_non_confirmed)
number_hits_per_mod_RBP_DRACH[3,2] <-length(eligos_DRACH_cyto_non_confirmed)

# add a metadata column to each GRanges object in which the names of the specific effectors associated with m6A 
# with which each hit overlaps will be reported as well as the function of these effectors (both initiated with Unknown)
mcols(eligos_DRACH_chr_non_confirmed) <- cbind(mcols(eligos_DRACH_chr_non_confirmed), effector = 'Unknown', type_effector = 'Unknown')
mcols(eligos_DRACH_nucleo_non_confirmed) <- cbind(mcols(eligos_DRACH_nucleo_non_confirmed), effector = 'Unknown', type_effector = 'Unknown')
mcols(eligos_DRACH_cyto_non_confirmed) <- cbind(mcols(eligos_DRACH_cyto_non_confirmed), effector = 'Unknown', type_effector = 'Unknown')

load('/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/RBPs/specific_effectors_only_m6A.Rda')
load('/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/RBPs/all_sites_RMBase_RMvar.Rda')
load('/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/RBPs/enzyme_type_mod_unique.RDa')

results_chr <- overlap_binding_sites(eligos_DRACH_chr_non_confirmed,1)
chr_specific_m6A <- eligos_DRACH_chr_non_confirmed[results_chr[[1]]]
eligos_DRACH_chr_non_confirmed_with_bindings <- results_chr[[2]]

results_nucleo <- overlap_binding_sites(eligos_DRACH_nucleo_non_confirmed,2)
nucleo_specific_m6A <- eligos_DRACH_nucleo_non_confirmed[results_nucleo[[1]]]
eligos_DRACH_nucleo_non_confirmed_with_bindings <- results_nucleo[[2]]

results_cyto <- overlap_binding_sites(eligos_DRACH_cyto_non_confirmed,3)
cyto_specific_m6A <- eligos_DRACH_cyto_non_confirmed[results_cyto[[1]]]
eligos_DRACH_cyto_non_confirmed_with_bindings <- results_cyto[[2]]

# compute the number of ELIGOS DRACH+ hits (not confirmed by m6Anet) overlapping with specific effectors associated
# with m6A and shared by all the fractions 
all_hits_specific_m6A <- list(chr_specific_m6A,nucleo_specific_m6A,cyto_specific_m6A)
order <- order(c(length(chr_specific_m6A), length(nucleo_specific_m6A), length(cyto_specific_m6A)))
overlap <- findOverlaps(all_hits_specific_m6A[[order[1]]],all_hits_specific_m6A[[order[2]]], type = 'any')
overlap2 <- all_hits_specific_m6A[[order[1]]][unique(queryHits(overlap))]
overlap_all <- findOverlaps(overlap2,all_hits_specific_m6A[[order[3]]], type = 'any')
number_hits_per_mod_RBP_DRACH[4,1] <- length(overlap2[unique(queryHits(overlap_all))])

write.xlsx(x = data.frame(number_hits_per_mod_RBP_DRACH),file = paste0(directory_hits_not_confirmed, 'ELIGOS_DRACH_not_confirmed_m6A_mark_sp_effectors.xlsx'),col.names = TRUE, row.names=TRUE)

save(eligos_DRACH_chr_non_confirmed_with_bindings, file =paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_chr_not_confirmed_m6A_mark_RBPs.Rda'))
save(eligos_DRACH_nucleo_non_confirmed_with_bindings, file =paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_nucleo_not_confirmed_m6A_mark_RBPs.Rda'))
save(eligos_DRACH_cyto_non_confirmed_with_bindings, file =paste0(directory_hits_not_confirmed,'hits_eligos_DRACH_cyto_not_confirmed_m6A_mark_RBPs.Rda'))

final_summary(directory_hits_not_confirmed)

