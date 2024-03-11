library('GenomicRanges')
library('GenomicFeatures')
library('xlsx')

# overlap between ELIGOS hits (only in Chromatin Associated nascent RNAs) and RNA marks, of each RNA mod type, from the two databases.
# directory_hits_only_chr is the path to the directory containing ELIGOS hits only in chromatin associated nascent RNAs.
# directory_mod is the path to the directory containing, for each RNA mod type, the bed file with
# the coordinates of all the RNA marks of that type from RMBase3+RMVar.
# m6A is TRUE when analysing DRACH+ hits, otherwise FALSE
overlap_marks <- function(hits_chr_grange,m6A,directory_hits_only_chr,directory_mod) {
  
  # initiate a matrix that will report the number (and percentage) of ELIGOS hits only present in chromatin associated nascent 
  # RNAs containing each RNA mod type and the overall number (and percentage) of ELIGOS hits only present in chromatin associated 
  # nascent RNAs annotated at least with one RNA modification type
  number_hits_per_mod <- matrix(0,nrow=2, ncol=12)
  rownames(number_hits_per_mod) <- c('Chromatin', 'Number of annotated marks')
  if (m6A) {
    colnames(number_hits_per_mod) <- c('m6A','Y','m1A','m5C', 'm7G','A-I', 'Nm', 'm6Am', 'm5U', 'Others', 'DRACH+ hits annotated', 'Tot DRACH+ hits')
  } else {
    colnames(number_hits_per_mod) <- c('m6A','Y','m1A','m5C', 'm7G','A-I', 'Nm', 'm6Am', 'm5U', 'Others', 'DRACH- hits annotated', 'Tot DRACH- hits')
  }
  
  number_hits_per_mod[1,12] <- as.character(length(hits_chr_grange))
  
  # initiate the vector that will report the index of ELIGOS hits only present in chromatin associated nascent RNAs overlapping with each RNA mod type
  query_matching_mark_chr <- c()
  
  mod <- list.files(path = directory_mod, pattern = "total.bed", full.names = TRUE)
  
  for (f in mod) {
    
    mod_file <- read.table(f)
    
    if (f == "/path/to/mod/other_RNA_mods_total.bed") {
      mod_grange <- GRanges(seqnames = gsub(pattern = 'M', replacement = 'MT', x = gsub(pattern='chr', replacement = '', x=mod_file[,1])),
                            ranges = IRanges(start = mod_file[,2], end=mod_file[,3]),
                            strand = Rle(mod_file[,6]),
                            mod_type = mod_file[,7])
      number_hits_per_mod[2,'Others'] <- as.character(length(mod_grange))
    } else {
      mod_grange <- GRanges(seqnames = gsub(pattern = 'M', replacement = 'MT', x = gsub(pattern='chr', replacement = '', x=mod_file$chr)),
                            ranges = IRanges(start = mod_file$start, end=mod_file$end),
                            strand = Rle(mod_file$strand),
                            mod_type = mod_file$mod_type)
      number_hits_per_mod[2,unique(mod_grange$mod_type)] <- as.character(length(mod_grange))
    }
    
    # overlap between ELIGOS hits and the coordinates of a specific RNA mark from the two databases
    suppressWarnings(overlap_chr <- findOverlaps(hits_chr_grange,mod_grange))
    hits_chr_with_mod <- c()
    if (length(overlap_chr) != 0) {
      # save inside the vector the index of ELIGOS hits overlapping with this RNA mod type
      query_matching_mark_chr <- c(query_matching_mark_chr, unique(queryHits(overlap_chr)))
      hits_chr_with_mod <- hits_chr_grange[unique(queryHits(overlap_chr))]
      # report in the matrix the number (and percentage) of hits that can be annotated with this RNA mod type
      if (f == "/path/to/mod/other_RNA_mods_total.bed") {
        number_hits_per_mod[1,'Others'] <- paste0(as.character(length(unique(queryHits(overlap_chr)))), ' - ', as.character(round(length(unique(queryHits(overlap_chr)))*100/as.numeric(number_hits_per_mod[1,12]),2)),'%')
      } else {
        number_hits_per_mod[1,unique(mod_grange$mod_type)] <- paste0(as.character(length(unique(queryHits(overlap_chr)))), ' - ', as.character(round(length(unique(queryHits(overlap_chr)))*100/as.numeric(number_hits_per_mod[1,12]),2)),'%')
      }
      # update the metadata field of ELIGOS hits overlapping with this modification type
      for (i in unique(queryHits(overlap_chr))) {
        if (hits_chr_grange[i]$mod_type == 'Unknown') {
          hits_chr_grange[i]$mod_type <- paste(unique(mod_grange[subjectHits(overlap_chr[queryHits(overlap_chr) == i])]$mod_type), collapse = ';')
        } else {
          hits_chr_grange[i]$mod_type <- paste(c(hits_chr_grange[i]$mod_type,unique(mod_grange[subjectHits(overlap_chr[queryHits(overlap_chr) == i])]$mod_type)), collapse = ';')
        }
      }
    } 
  }
  
  # compute the number (and percentage) of ELIGOS hits annotated at least with one RNA mod type
  number_hits_per_mod[1,11] <- paste0(as.character(length(unique(query_matching_mark_chr))), ' - ', as.character(round(length(unique(query_matching_mark_chr))*100/as.numeric(number_hits_per_mod[1,12]),2)), '%')
  
  if (m6A) {
    write.xlsx(x = data.frame(number_hits_per_mod),file = paste0(directory_hits_only_chr, '/with_DRACH/number_hit_per_mod_DRACH_bothBD_only_chr.xlsx'),col.names = TRUE, row.names =TRUE)
  } else {
    write.xlsx(x = data.frame(number_hits_per_mod),file = paste0(directory_hits_only_chr, '/without_DRACH/number_hit_per_mod_non_DRACH_bothBD_only_chr.xlsx'),col.names = TRUE, row.names =TRUE)
  }
  
  unique_mod_types <- function(hits) {
    hits$mod_type <- unlist(lapply(seq_along(hits), function(i,x) {
      paste(unique(unlist(strsplit(x[i]$mod_type, split =';'))), collapse = ';')
    }, x=hits))
    return(hits)
  }
  
  hits_chr_grange <- unique_mod_types(hits_chr_grange)
  
  return(hits_chr_grange)
}

# directory_hits_only_chr is the path to the directory containing ELIGOS hits that are only in nascent chromatin associated RNAs.
# directory_mod is the path to the directory containing, for each RNA mod type, the bed file with
# the coordinates of all the RNA marks of that type from RMBase3+RMVar
ELIGOS_hits_DB_marks <- function(directory_mod, directory_hits_only_chr) {
  
  # load the RData containing ELIGOS hits only present in Chromatin Associated nascent RNAs, DRACH+ and DRACH- separately
  load(paste0(directory_hits_only_chr,'/without_DRACH/hits_eligos_only_chr_without_DRACH.Rda'))
  only_chr_without_DRACH <- only_chr
  load(paste0(directory_hits_only_chr,'/with_DRACH/hits_eligos_only_chr_with_DRACH.Rda'))
  only_chr_with_DRACH <- only_chr
  
  # add a metadata column to each GRanges object in which the RNA marks with which each hit overlaps will be reported (initiated 
  # with "Unknown")
  mcols(only_chr_without_DRACH) <- cbind(mcols(only_chr_without_DRACH), mod_type= 'Unknown')
  mcols(only_chr_with_DRACH) <- cbind(mcols(only_chr_with_DRACH), mod_type= 'Unknown')
  
  hits_DRACH <- overlap_marks(only_chr_with_DRACH,TRUE,directory_hits_only_chr,directory_mod)
  hits_non_DRACH <- overlap_marks(only_chr_without_DRACH,FALSE,directory_hits_only_chr,directory_mod)
  
  save(hits_DRACH, file=paste0(directory_hits_only_chr,'/with_DRACH/hits_eligos_only_chr_with_DRACH_mod_type.Rda'))
  save(hits_non_DRACH, file =paste0(directory_hits_only_chr,'/without_DRACH/hits_eligos_only_chr_without_DRACH_mod_type.Rda'))
}

ELIGOS_hits_DB_marks(directory_mod = '/path/to/mod', 
                     directory_hits_only_chr = '/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/only_chr')

################

load('/path/to/RBPs/specific_effectors.Rda')
load('/path/to/RBPs/specific_effectors_only_m6A.Rda')
load('/path/to/RBPs/specific_effectors_non_m6A.Rda')
load('/path/to/RBPs/enzyme_type_mod_unique.RDa')
load('/path/to/RBPs/all_sites_RMBase_RMvar.Rda')

# overlap between ELIGOS hits (only in Chromatin Associated nascent RNAs) and effectors' binding sites from the two databases.
# m6A is TRUE when analysing DRACH+ hits, otherwise FALSE
overlap_binding_sites <- function(hits, n, m, p, q, r, m6A) {
  
  over <-findOverlaps(hits,all_sites_RMBase_RMvar, ignore.strand=TRUE, type='any')
  hits_overlapping_specific_effector_m6A <- c()
  hits_overlapping_specific_effector_non_m6A <- c()
  # report the index of ELIGOS hits overlapping with specific effectors associated with m6A or specific effectors not associated
  # with m6A
  for (i in 1:length(over)) {
    if (all_sites_RMBase_RMvar[subjectHits(over[i])]$RBP_name %in% specific_effectors_only_m6A) {
      hits_overlapping_specific_effector_m6A <- c(hits_overlapping_specific_effector_m6A, queryHits(over[i]))
    }
    if (all_sites_RMBase_RMvar[subjectHits(over[i])]$RBP_name %in% specific_effectors_non_m6A) {
      hits_overlapping_specific_effector_non_m6A <- c(hits_overlapping_specific_effector_non_m6A, queryHits(over[i]))
    }
  }
  if (m6A) {
    number_hits_per_mod_RBP_DRACH[n,m] <<-paste0(as.character(length(unique(hits_overlapping_specific_effector_non_m6A))), ' - ', as.character(round(length(unique(hits_overlapping_specific_effector_non_m6A))/length(hits)*100,2)),'%')
    number_hits_per_mod_RBP_DRACH[n,p] <<- paste0(as.character(length(unique(hits_overlapping_specific_effector_m6A))),' - ', as.character(round(length(unique(hits_overlapping_specific_effector_m6A))/length(hits)*100,2)),'%')
    number_hits_per_mod_RBP_DRACH[n,q] <<-paste0(as.character(length(unique(c(hits_overlapping_specific_effector_non_m6A,hits_overlapping_specific_effector_m6A)))),' - ', as.character(round(length(unique(c(hits_overlapping_specific_effector_non_m6A,hits_overlapping_specific_effector_m6A)))/length(hits)*100,2)),'%')
    number_hits_per_mod_RBP_DRACH[n,r] <<- paste0(as.character(length(unique(queryHits(over)))), ' - ', as.character(round(length(unique(queryHits(over)))/length(hits)*100,2)),'%')
  } else {
    number_hits_per_mod_RBP_non_DRACH[n,m] <<-paste0(as.character(length(unique(hits_overlapping_specific_effector_non_m6A))), ' - ',as.character(round(length(unique(hits_overlapping_specific_effector_non_m6A))/length(hits)*100,2)),'%')
    number_hits_per_mod_RBP_non_DRACH[n,p] <<- paste0(as.character(length(unique(hits_overlapping_specific_effector_m6A))),' - ', as.character(round(length(unique(hits_overlapping_specific_effector_m6A))/length(hits)*100,2)),'%')
    number_hits_per_mod_RBP_non_DRACH[n,q] <<-paste0(as.character(length(unique(c(hits_overlapping_specific_effector_non_m6A,hits_overlapping_specific_effector_m6A)))),' - ',as.character(round(length(unique(c(hits_overlapping_specific_effector_non_m6A,hits_overlapping_specific_effector_m6A)))/length(hits)*100,2)),'%')
    number_hits_per_mod_RBP_non_DRACH[n,r] <<- paste0(as.character(length(unique(queryHits(over)))), ' - ',as.character(round(length(unique(queryHits(over)))/length(hits)*100,2)),'%')
  }
  
  # update the field "mod_type" of ELIGOS hits that didn't overlap with any RNA mark from the databases
  for (i in 1:length(hits)) {
    if (hits[i]$mod_type == 'Unknown') {
      if (i %in% queryHits(over)) {
        binding_sites <- unique(all_sites_RMBase_RMvar[subjectHits(over[queryHits(over)==i])])
        RBP_specific <- unique(binding_sites$RBP_name[which(binding_sites$RBP_name %in% specific_effectors)])
        RBP_ambiguous <- unique(binding_sites$RBP_name[which(!binding_sites$RBP_name %in% specific_effectors)])
        if (length(RBP_specific) == 0) {
          # if the hit overlaps only with one or more ambiguous effectors, the names of these effectors are added in the metadata separated by “;” 
          # and the modification type of the hit is annotated as “Ambiguous”
          hits[i]$effector <- paste(RBP_ambiguous, collapse = ';')
          hits[i]$mod_type <- 'Ambiguous'
          hits[i]$type_effector <- 'Ambiguous'
        } else if (length(RBP_specific) == length(unique(binding_sites$RBP_name))) {
          # if the hit overlaps only with one or more specific effectors (associated either with m6A or with another RNA mark), 
          # the names of these effectors are added in the metadata separated by “;” and the modification type of the hit is annotated with the 
          # names of the RNA marks related to these specific effectors separated by “;”
          hits[i]$effector <- paste(RBP_specific, collapse = ';')
          hits[i]$mod_type <- paste(unique(enzyme_type_mod_unique[enzyme_type_mod_unique$Name_effector %in% RBP_specific,3]), collapse=';')
          hits[i]$type_effector <- paste(unique(enzyme_type_mod_unique[enzyme_type_mod_unique$Name_effector %in% RBP_specific,2]), collapse=';')
        } else {
          # if the hit overlaps both with one or more specific effectors and with one or more ambiguous effectors, the names of these effectors 
          # are added in the metadata separated by “;” and the modification type of the hit is annotated with the names of the RNA marks related 
          # to the specific effectors with which the hit overlaps separated by “;” also including the modification type “Ambiguous”
          hits[i]$effector <- paste(c(RBP_specific,RBP_ambiguous), collapse = ';')
          hits[i]$mod_type <- paste(c(unique(enzyme_type_mod_unique[enzyme_type_mod_unique$Name_effector %in% RBP_specific,3]), 'Ambiguous'), collapse=';')
          hits[i]$type_effector <- paste(c(unique(enzyme_type_mod_unique[enzyme_type_mod_unique$Name_effector %in% RBP_specific,2]), 'Ambiguous'), collapse=';')
        }
      }
    }}
  
  hits$mod_type <- unlist(lapply(seq_along(hits), function(i,x) {
    paste(unique(unlist(strsplit(x[i]$mod_type, split =';'))), collapse = ';')
  }, x=hits))
  
  l <- list(unique(hits_overlapping_specific_effector_non_m6A), unique(hits_overlapping_specific_effector_m6A),unique(c(hits_overlapping_specific_effector_non_m6A,hits_overlapping_specific_effector_m6A)),unique(queryHits(over)), hits)
  names(l) <- c('specific non m6A', 'specific m6A','specific', 'total', 'hits')
  return(l)
}

# directory_hits_only_chr is the path to the directory containing ELIGOS hits that are only in nascent chromatin associated RNAs
ELIGOS_hits_DB_effectors <- function(directory_hits_only_chr) {
  
  load(paste0(directory_hits_only_chr,'/without_DRACH/hits_eligos_only_chr_without_DRACH_mod_type.Rda'))
  load(paste0(directory_hits_only_chr,'/with_DRACH/hits_eligos_only_chr_with_DRACH_mod_type.Rda'))
  
  # initiate a matrix to report the percentage of ELIGOS DRACH- hits overlapping with each category of effectors and
  # the percentage of ELIGOS DRACH- hits overlapping at least with one category
  number_hits_per_mod_RBP_non_DRACH <<- matrix(0,nrow=2, ncol=5)
  rownames(number_hits_per_mod_RBP_non_DRACH) <<- c('Chromatin', 'Number of RBPs annotated')
  colnames(number_hits_per_mod_RBP_non_DRACH) <<- c('% of hits overlapping with\nspecific effectors non m6A', 
                                                    '% of hits overlapping with\nspecific effectors m6A',
                                                    '% of hits overlapping with\nspecific effectors',
                                                    '% of hits overlapping with\ntotal effectors',
                                                    'Tot number of DRACH- hits')
  
  # the same for ELIGOS DRACH+ hits 
  number_hits_per_mod_RBP_DRACH <<- matrix(0,nrow=2, ncol=5)
  rownames(number_hits_per_mod_RBP_DRACH) <<- c('Chromatin', 'Number of RBPs annotated')
  colnames(number_hits_per_mod_RBP_DRACH) <<- c('% of hits overlapping with\nspecific effectors non m6A', 
                                                '% of hits overlapping with\nspecific effectors m6A',
                                                '% of hits overlapping with\nspecific effectors',
                                                '% of hits overlapping with\ntotal effectors',
                                                'Tot number of DRACH+ hits')
  
  number_hits_per_mod_RBP_non_DRACH[1,5]<<- length(hits_non_DRACH)
  number_hits_per_mod_RBP_non_DRACH[2,1] <<-length(specific_effectors_non_m6A)
  number_hits_per_mod_RBP_non_DRACH[2,2] <<-length(specific_effectors_only_m6A)
  number_hits_per_mod_RBP_non_DRACH[2,3] <<-length(specific_effectors)
  number_hits_per_mod_RBP_non_DRACH[2,4] <<-length(unique(enzyme_type_mod_unique[,1]))
  number_hits_per_mod_RBP_DRACH[1,5] <<- length(hits_DRACH)
  number_hits_per_mod_RBP_DRACH[2,1] <<-length(specific_effectors_non_m6A)
  number_hits_per_mod_RBP_DRACH[2,2] <<-length(specific_effectors_only_m6A)
  number_hits_per_mod_RBP_DRACH[2,3] <<-length(specific_effectors)
  number_hits_per_mod_RBP_DRACH[2,4] <<-length(unique(enzyme_type_mod_unique[,1]))
  
  # add a metadata column to each GRanges object in which the names of the effectors with hit overlaps will be reported as well as the 
  # function of these effectors (both initiated with "Unknown")
  mcols(hits_non_DRACH) <- cbind(mcols(hits_non_DRACH), effector = 'Unknown', type_effector = 'Unknown')
  mcols(hits_DRACH) <- cbind(mcols(hits_DRACH), effector = 'Unknown', type_effector = 'Unknown')
  
  results_chr <- overlap_binding_sites(hits_non_DRACH,1,1,2,3,4, FALSE)
  chr_specific_non_m6A <- hits_non_DRACH[results_chr[[1]]]
  chr_specific_m6A <- hits_non_DRACH[results_chr[[2]]]
  chr_specific <- hits_non_DRACH[results_chr[[3]]]
  chr_total <- hits_non_DRACH[results_chr[[4]]]
  hits_non_DRACH_with_bindings <- results_chr[[5]]
  
  results_chr <- overlap_binding_sites(hits_DRACH,1,1,2,3,4, TRUE)
  chr_specific_non_m6A <- hits_DRACH[results_chr[[1]]]
  chr_specific_m6A <- hits_DRACH[results_chr[[2]]]
  chr_specific <- hits_DRACH[results_chr[[3]]]
  chr_total <- hits_DRACH[results_chr[[4]]]
  hits_DRACH_with_bindings <- results_chr[[5]]
  
  write.xlsx(x = data.frame(number_hits_per_mod_RBP_non_DRACH),file = paste0(directory_hits_only_chr, '/without_DRACH/hits_overlapping_RBPs_without_DRACH_only_chr.xlsx'),col.names = TRUE, row.names = TRUE)
  write.xlsx(x = data.frame(number_hits_per_mod_RBP_DRACH),file = paste0(directory_hits_only_chr, '/with_DRACH/hits_overlapping_RBPs_with_DRACH_only_chr.xlsx'),col.names = TRUE, row.names = TRUE)
  
  save(hits_non_DRACH_with_bindings, file =paste0(directory_hits_only_chr,'/without_DRACH/hits_eligos_only_chr_without_DRACH_mod_type_RBP.Rda'))
  save(hits_DRACH_with_bindings, file=paste0(directory_hits_only_chr,'/with_DRACH/hits_eligos_only_chr_with_DRACH_mod_type_RBP.Rda'))
}

ELIGOS_hits_DB_effectors(directory_hits_only_chr = '/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/only_chr/')

################

# combine the information from the overlap with the RNA marks and the overlap with the binding sites of the
# effectors from the public databases.
# directory_hits_only_chr is the path to the directory containing ELIGOS hits thare only in nascent chromatin associated RNAs
final_summary <- function(all_mods, directory_hits_only_chr) {
  
  load(paste0(directory_hits_only_chr,'/without_DRACH/hits_eligos_only_chr_without_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits_only_chr,'/with_DRACH/hits_eligos_only_chr_with_DRACH_mod_type_RBP.Rda'))
  
  # initiate a matrix reporting, for each RNA mod type, the number (and percentage) of ELIGOS DRACH- hits 
  # annotated to that RNA mod
  number_hits_non_DRACH <- matrix(0,nrow=1, ncol=length(all_mods)+2)
  rownames(number_hits_non_DRACH) <- c('Chromatin')
  colnames(number_hits_non_DRACH) <- c(all_mods, 'Tot number of DRACH- annotated hits', 'Tot number of DRACH- hits')
  number_hits_non_DRACH[1,'Tot number of DRACH- hits'] <- length(hits_non_DRACH_with_bindings)
  
  # initiate a matrix reporting, for each RNA mod type, the number (and percentage) of ELIGOS DRACH+ hits 
  # annotated to that RNA mod
  number_hits_DRACH <- matrix(0,nrow=1, ncol=length(all_mods)+2)
  rownames(number_hits_DRACH) <- c('Chromatin')
  colnames(number_hits_DRACH) <- c(all_mods, 'Tot number of DRACH+ annotated hits', 'Tot number of DRACH+ hits')
  number_hits_DRACH[1,'Tot number of DRACH+ hits'] <- length(hits_DRACH_with_bindings)
  
  count_mods_non_DRACH <- function(hits_chr) {
    mod_types_chr <- unlist(strsplit(hits_chr$mod_type, split=';'))
    # iterate over the modification types
    for (mod in colnames(number_hits_non_DRACH)[1:length(all_mods)]) {
      # compute the number and percentage of ELIGOS DRACH- hits annotated to that RNA mod type
      number_hits_non_DRACH[1,mod] <<- paste0(as.character(length(mod_types_chr[mod_types_chr==mod])),' - ', as.character(round(length(mod_types_chr[mod_types_chr==mod])*100/as.numeric(number_hits_non_DRACH[1,'Tot number of DRACH- hits']),2)),'%')
    }
  }
  
  count_mods_DRACH <- function(hits_chr) {
    mod_types_chr <- unlist(strsplit(hits_chr$mod_type, split=';'))
    # iterate over the modification types
    for (mod in colnames(number_hits_DRACH)[1:length(all_mods)]) {
      # compute the number and percentage of ELIGOS DRACH+ hits annotated to that RNA mod type
      number_hits_DRACH[1,mod] <<- paste0(as.character(length(mod_types_chr[mod_types_chr==mod])),' - ', as.character(round(length(mod_types_chr[mod_types_chr==mod])*100/as.numeric(number_hits_DRACH[1,'Tot number of DRACH+ hits']),2)),'%')
    }
  }
  
  count_mods_non_DRACH(hits_non_DRACH_with_bindings)
  count_mods_DRACH(hits_DRACH_with_bindings)
  
  # remove ELIGOS DRACH-/DRACH+ hits that don't overlap with any RNA mark/effectors' binding site (Unknown) 
  # and those that overlap only with one or more ambiguous effectors' binding sites (Ambiguous)
  hits_non_DRACH_with_bindings_marks_specific <- hits_non_DRACH_with_bindings[!hits_non_DRACH_with_bindings$mod_type %in% c('Unknown', 'Ambiguous')]
  hits_DRACH_with_bindings_marks_specific <- hits_DRACH_with_bindings[!hits_DRACH_with_bindings$mod_type  %in% c('Unknown', 'Ambiguous')]
  
  # compute how many ELIGOS DRACH-/DRACH+ hits are annotated at least with one RNA mark or overlap with the 
  # binding site of at least one specific effector
  number_hits_non_DRACH[1,'Tot number of DRACH- annotated hits'] <- paste0(as.character(length(hits_non_DRACH_with_bindings_marks_specific)), ' - ', as.character(round(length(hits_non_DRACH_with_bindings_marks_specific)*100/as.numeric(number_hits_non_DRACH[1,'Tot number of DRACH- hits']),2)),'%')
  number_hits_DRACH[1,'Tot number of DRACH+ annotated hits'] <- paste0(as.character(length(hits_DRACH_with_bindings_marks_specific)), ' - ', as.character(round(length(hits_DRACH_with_bindings_marks_specific)*100/as.numeric(number_hits_DRACH[1,'Tot number of DRACH+ hits']),2)),'%')
  
  write.xlsx(x = data.frame(number_hits_non_DRACH),file = paste0(directory_hits_only_chr, '/without_DRACH/hits_without_DRACH_only_chr_summary_DBs.xlsx'),col.names = TRUE, row.names = TRUE)
  write.xlsx(x = data.frame(number_hits_DRACH),file = paste0(directory_hits_only_chr, '/with_DRACH/hits_with_DRACH_only_chr_summary_DBs.xlsx'),col.names = TRUE, row.names = TRUE)
}

directory_hits_only_chr = '/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/only_chr/'

load(paste0(directory_hits_only_chr,'/without_DRACH/hits_eligos_only_chr_without_DRACH_mod_type_RBP.Rda'))
load(paste0(directory_hits_only_chr,'/with_DRACH/hits_eligos_only_chr_with_DRACH_mod_type_RBP.Rda'))

# identify which are the RNA mod types with which ELIGOS hits are annotated after the overlap with the RNA marks and 
# the effectors from the databases
all_mods <- c(unique(hits_DRACH_with_bindings$mod_type),
              unique(hits_non_DRACH_with_bindings$mod_type))
all_mods <- unique(all_mods)
all_mods <- unique(unlist(strsplit(all_mods, split=';')))
all_mods <- c('m6A', 'm1A', 'm5C', 'm7G', 'A-I', 'Nm', 'Ambiguous')

final_summary(all_mods,directory_hits_only_chr)

