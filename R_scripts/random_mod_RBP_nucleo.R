# overlap between each of the 1,000 sets of randomly generated DRACH+ and DRACH- hits (nucleoplasm fraction) 
# and the coordinates of effectors' binding sites from RMBase3 and RMVar

library('GenomicRanges')
library('GenomicFeatures')

enzyme_type_mod_unique <- read.table('/home/pmaragno/random/RBPs/enzyme_type_mod_unique.RDa')
load('/home/pmaragno/random/RBPs/specific_effectors.Rda')
load('/home/pmaragno/random/RBPs/specific_effectors_only_m6A.Rda')
load('/home/pmaragno/random/RBPs/specific_effectors_non_m6A.Rda')
load('/home/pmaragno/random/RBPs/ambiguous_effectors.Rda')
load('/home/pmaragno/random/RBPs/all_sites_RMBase_RMvar.Rda')

compute_overlapping_hits_binding_sites <- function(hits) {
  # add to the i-th GRanges of random hits a metadata column in which the names of the effectors with which
  # each hit overlaps will be reported as well as the function of these effectors (both initiated with Unknown)
  mcols(hits) <- cbind(mcols(hits), effector = 'Unknown', type_effector = 'Unknown')
  
  # overlap between the i-th GRanges of random hits and all the effectors' binding sites
  over <- findOverlaps(hits,all_sites_RMBase_RMvar, ignore.strand=TRUE, type='any')
  hits_overlapping_m6A <- unique(queryHits(over[all_sites_RMBase_RMvar[subjectHits(over)]$RBP_name %in% specific_effectors_only_m6A]))
  hits_overlapping_non_m6A <- unique(queryHits(over[all_sites_RMBase_RMvar[subjectHits(over)]$RBP_name %in% specific_effectors_non_m6A]))
  # report the percentage of random hits of the i-th GRanges overlapping with specific effectors associated with m6A/
  # specific effectors not associated with m6A/ambiguous effectors/whatever kind of effector
  m6A <<- c(m6A, round(length(unique(hits_overlapping_m6A))/length(hits)*100,2))
  non_m6A <<- c(non_m6A, round(length(unique(hits_overlapping_non_m6A))/length(hits)*100,2))
  specific <<- c(specific, round(length(unique(c(hits_overlapping_m6A,hits_overlapping_non_m6A)))/length(hits)*100,2))
  total <<- c(total, round(length(unique(queryHits(over)))/length(hits)*100,2))
  
  # update the field "mod_type" of the i-th GRanges of random hits, only for the hits not overlapping with any RNA mark from 
  # the databases
  lapply(as.list(unique(queryHits(over))), function(y){
    if (hits[y]$mod_type == 'Unknown') {
      binding_sites <- unique(all_sites_RMBase_RMvar[subjectHits(over[queryHits(over)==y])])
      RBP_specific <- unique(binding_sites$RBP_name[which(binding_sites$RBP_name %in% specific_effectors)])
      RBP_ambiguous <- unique(binding_sites$RBP_name[which(!binding_sites$RBP_name %in% specific_effectors)])
      if (length(RBP_specific) == 0) {
        # if the hit overlaps only with one or more ambiguous effectors, the names of these effectors are added in the metadata separated by ";"
        # and the modification type of the hit is annotated as Ambiguous
        hits[y]$effector <<- paste(unique(RBP_ambiguous), collapse = ';')
        hits[y]$mod_type <<- 'Ambiguous'
        hits[y]$type_effector <<- 'Ambiguous'
      } else if (length(RBP_specific) == length(unique(binding_sites$RBP_name))) {
        # if the hit overlaps only with one or more specific effectors (associated either with m6A or with another RNA mark), 
        # the names of these effectors are added in the metadata separated by ";" and the modification type of the hit is annotated with the 
        # names of the RNA marks related to these specific effectors separated by ";"
        hits[y]$effector <<- paste(unique(RBP_specific), collapse = ';')
        hits[y]$mod_type <<- paste(unique(enzyme_type_mod_unique[enzyme_type_mod_unique$Name_effector %in% RBP_specific,3]), collapse=';')
        hits[y]$type_effector <<- paste(unique(enzyme_type_mod_unique[enzyme_type_mod_unique$Name_effector %in% RBP_specific,2]), collapse=';')
      } else {
        # if the hit overlaps both with one or more specific effectors and with one or more ambiguous effectors, the names of these effectors 
        # are added in the metadata separated by ";" and the modification type of the hit is annotated with the names of the RNA marks related 
        # to the specific effectors with which the hit overlaps separated by ";"  also including the modification type Ambiguous
        hits[y]$effector <<- paste(c(RBP_specific,RBP_ambiguous), collapse = ';')
        hits[y]$mod_type <<- paste(c(unique(enzyme_type_mod_unique[enzyme_type_mod_unique$Name_effector %in% RBP_specific,3]), 'Ambiguous'), collapse=';')
        hits[y]$type_effector <<- paste(c(unique(enzyme_type_mod_unique[enzyme_type_mod_unique$Name_effector %in% RBP_specific,2]), 'Ambiguous'), collapse=';')
      }
    }
  })
  
  lapply(seq_along(1:length(hits)), function(i){
    hits[i]$mod_type <<- paste(unique(unlist(strsplit(hits[i]$mod_type, split =';'))), collapse = ';')
  })
  return(hits)
}

load('/home/pmaragno/random/nucleo_mod_type_DRACH.Rda')
m6A <- c()
non_m6A <- c()
specific <- c()
total<- c()
# overlap between each of the 1,000 sets of randomly generated DRACH+ hits and the binding sites of effectors associated with RNA mods
# from the two databases
hits_nucleo_mod_RBP_DRACH <- lapply(nucleo_mod_type_DRACH[[12]], compute_overlapping_hits_binding_sites)
hits_nucleo_mod_RBP_DRACH <- c(list(m6A), list(non_m6A), list(specific), list(total), hits_nucleo_mod_RBP_DRACH)
names(hits_nucleo_mod_RBP_DRACH) <- c('m6A', 'non_m6A', 'specific', 'total', 'hits')
save(hits_nucleo_mod_RBP_DRACH, file='/home/pmaragno/random/hits_nucleo_mod_RBP_DRACH.Rda')

load('/home/pmaragno/random/nucleo_mod_type_without_DRACH.Rda')
m6A <- c()
non_m6A <- c()
specific <- c()
total<- c()
# overlap between each of the 1,000 sets of randomly generated DRACH- hits and the binding sites of effectors associated with RNA mods
# from the two databases
hits_nucleo_mod_RBP_without_DRACH <- lapply(nucleo_mod_type_without_DRACH[[12]], compute_overlapping_hits_binding_sites)
hits_nucleo_mod_RBP_without_DRACH <- c(list(m6A), list(non_m6A), list(specific), list(total), hits_nucleo_mod_RBP_without_DRACH)
names(hits_nucleo_mod_RBP_without_DRACH) <- c('m6A', 'non_m6A', 'specific', 'total', 'hits')
save(hits_nucleo_mod_RBP_without_DRACH, file='/home/pmaragno/random/hits_nucleo_mod_RBP_without_DRACH.Rda')

