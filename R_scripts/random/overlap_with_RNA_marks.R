library('GenomicRanges')
library('GenomicFeatures')

mod_position <- function(directory_mod, hits) {
  
  # initiate a vector for each RNA mod type in which, for each of the 1,000 random sets, the corresponding percentage of sequences
  # overlapping with at least one RNA mark of that mod type will be saved
  m6A <- c()
  Y <- c()
  m1A <- c()
  m5C <- c()
  m7G <- c()
  AI <- c()
  Nm <- c()
  m6Am <- c()
  m5U <- c()
  Others <- c()
  # initiate the vector that will report, for each of the 1,000 sets, the index of the random sequences overlapping with 
  # each RNA mod type
  Characterized_hits <- list()
  
  # add to each of the 1,000 GRanges a metadata column in which the RNA marks with which each sequence overlaps will be reported (initiated 
  # with "Unknown")
  hits <- lapply(hits, function(x) {
    mcols(x) <- cbind(mcols(x), mod_type = 'Unknown')
    return(x)
  })
  
  mod <- list.files(path = directory_mod, pattern = "total.bed", full.names = TRUE)
  
  # iterate over the RNA modification types
  lapply(as.list(mod), function(f) {
    
    mod_file <- read.table(f)
    
    if (f == "/path/to/mod/other_RNA_mods_total.bed") {
      mod_grange <- GRanges(seqnames = gsub(pattern = 'M', replacement = 'MT', x = gsub(pattern='chr', replacement = '', x=mod_file$V1)),
                            ranges = IRanges(start = mod_file$V2, end=mod_file$V3),
                            strand = Rle(mod_file$V6),
                            mod_type = mod_file$V7)
     } else {
      mod_grange <- GRanges(seqnames = gsub(pattern = 'M', replacement = 'MT', x = gsub(pattern='chr', replacement = '', x=mod_file$chr)),
                            ranges = IRanges(start = mod_file$start, end=mod_file$end),
                            strand = Rle(mod_file$strand),
                            mod_type = mod_file$mod_type)
     }
    
    # iterate over the 1,000 sets of random sequences
    lapply(seq_along(1:1000), function(i) {
      # overlap between the i-th set of random sequences and the coordinates of a specific RNA mark from the two databases
      suppressWarnings(overlap <- findOverlaps(hits[[i]], mod_grange))
      # save inside the vector the index of the random sequences of the i-th set overlapping with this RNA mod type
      if (f == mod[1]) {
        if (length(overlap) == 0) {
          Characterized_hits <<- c(Characterized_hits, list(c(0)))
        } else {
          Characterized_hits <<- c(Characterized_hits, list(unique(queryHits(overlap))))
        }
      } else {
        Characterized_hits[[i]] <<- c(Characterized_hits[[i]], unique(queryHits(overlap)))
      }
      
      # save the percentage of random sequences overlapping with the marks of each RNA mod type in the corresponding vector
      if (length(overlap) != 0) {
        if (f == "/path/to/mod/other_RNA_mods_total.bed") {
          Others <<- c(Others, length(unique(queryHits(overlap))) * 100 / length(hits[[i]]))
        } else if (f == "/path/to/mod/A-I_total.bed") {
          AI <<- c(AI, length(unique(queryHits(overlap))) * 100 / length(hits[[i]]))
        } else if (f == "/path/to/mod/m1A_total.bed") {
          m1A <<- c(m1A, length(unique(queryHits(overlap))) * 100 / length(hits[[i]]))
        } else if (f == "/path/to/mod/m5C_total.bed") {
          m5C <<- c(m5C, length(unique(queryHits(overlap))) * 100 / length(hits[[i]]))
        } else if (f == "/path/to/mod/m5U_total.bed") {
          m5U <<- c(m5U, length(unique(queryHits(overlap))) * 100 / length(hits[[i]]))
        } else if (f == "/path/to/mod/m6A_total.bed") {
          m6A <<- c(m6A, length(unique(queryHits(overlap))) * 100 / length(hits[[i]]))
        } else if (f == "/path/to/mod/m6Am_total.bed") {
          m6Am <<- c(m6Am, length(unique(queryHits(overlap))) * 100 / length(hits[[i]]))
        } else if (f == "/path/to/mod/m7G_total.bed") {
          m7G <<- c(m7G, length(unique(queryHits(overlap))) * 100 / length(hits[[i]]))
        } else if (f == "/path/to/mod/Nm_total.bed") {
          Nm <<- c(Nm, length(unique(queryHits(overlap))) * 100 / length(hits[[i]]))
        } else {
          Y <<- c(Y, length(unique(queryHits(overlap))) * 100 / length(hits[[i]]))
        }} else {
          # if the i-th set of random hits doesn't have any sequence overlapping with the coordinates of that specific RNA mark, 
          # add 0 in the corresponding vector
          if (f == "/path/to/mod/other_RNA_mods_total.bed") {
            Others <<- c(Others, 0)
          } else if (f == "/path/to/mod/A-I_total.bed") {
            AI <<- c(AI, 0)
          } else if (f == "/path/to/mod/m1A_total.bed") {
            m1A <<- c(m1A,0)
          } else if (f == "/path/to/mod/m5C_total.bed") {
            m5C <<- c(m5C, 0)
          } else if (f == "/path/to/mod/m5U_total.bed") {
            m5U <<- c(m5U, 0)
          } else if (f == "/path/to/mod/m6A_total.bed") {
            m6A <<- c(m6A, 0)
          } else if (f == "/path/to/mod/m6Am_total.bed") {
            m6Am <<- c(m6Am, 0)
          } else if (f == "/path/to/mod/m7G_total.bed") {
            m7G <<- c(m7G, 0)
          } else if (f == "/path/to/mod/Nm_total.bed") {
            Nm <<- c(Nm, 0)
          } else {
            Y <<- c(Y, 0)
          }
        }
      
      # update the field "mod_type" of the i-th set of random sequences
      lapply(unique(queryHits(overlap)), function(o) {
        if (hits[[i]][o]$mod_type == 'Unknown') {
          hits[[i]][o]$mod_type <<- paste(unique(mod_grange[subjectHits(overlap[queryHits(overlap) == o])]$mod_type), collapse = ';')
        } else {
          hits[[i]][o]$mod_type <<- paste(unique(c(hits[[i]][o]$mod_type, unique(mod_grange[subjectHits(overlap[queryHits(overlap) == o])]$mod_type))), collapse = ';')
        }
      })
    })
  })
  
  # for each of the 1,000 random sets, compute the percentage of sequences annotated at least with one RNA mod type
  Characterized_hits <- lapply(Characterized_hits , function(x) {
    length(unique(x)[unique(x) != 0]) * 100 / length(hits[[1]])
  })
  
  l <- list(m6A,Y,m1A,m5C,m7G,AI,Nm,m6Am,m5U,Others,Characterized_hits,hits)
  names(l) <- c('m6A','Y','m1A','m5C','m7G','A-I','Nm','m6Am','m5U','Others','Characterized_hits','hits'
  )
  
  return(l)
}

load('/path/to/random_hits_chr_DRACH.Rda')
load('/path/to/random_hits_chr_without_DRACH.Rda')

# compute the minimum number mn of random DRACH+ sequences from chromatin fraction across the 1,000 sets
min_DRACH <- min(unlist(lapply(random_hits_chr_DRACH, length)))
# for each of the 1,000 sets of random DRACH+ sequences from chromatin fraction extract the first mn random sequences
random_hits_chr_DRACH_reduced <- lapply(random_hits_chr_DRACH, function(x) {
  x[1:min_DRACH]
})
save(random_hits_chr_DRACH_reduced, file='/path/to/random_hits_chr_DRACH_reduced.Rda')

# overlap between each of the 1,000 sets of randomly generated DRACH+ sequences and RNA marks, of each RNA mod type, from the two databases
chr_ass_mod_type_DRACH <- mod_position('/path/to/mod',random_hits_chr_DRACH_reduced)  
save(chr_ass_mod_type_DRACH, file = '/path/to/chr_ass_mod_type_DRACH.Rda')

# compute the minimum number mn of random DRACH- sequences from chromatin fraction across the 1,000 sets
min_without_DRACH <- min(unlist(lapply(random_hits_chr_without_DRACH, length)))
# for each of the 1,000 sets of random DRACH- sequences from chromatin fraction extract the first mn random sequences
random_hits_chr_without_DRACH_reduced <- lapply(random_hits_chr_without_DRACH, function(x) {
  x[1:min_without_DRACH]
})
save(random_hits_chr_without_DRACH_reduced, file='/path/to/random_hits_chr_without_DRACH_reduced.Rda')

# overlap between each of the 1,000 sets of randomly generated DRACH- sequences and RNA marks, of each RNA mod type, from the two databases
chr_ass_mod_type_without_DRACH <- mod_position('/path/to/mod',random_hits_chr_without_DRACH_reduced)  
save(chr_ass_mod_type_without_DRACH, file='/path/to/chr_ass_mod_type_without_DRACH.Rda')

load('/path/to/random_hits_nucleo_DRACH.Rda')
load('/path/to/random_hits_nucleo_without_DRACH.Rda')

# compute the minimum number mn of random DRACH+ sequences from nucleoplasm fraction across the 1,000 sets
min_DRACH <- min(unlist(lapply(random_hits_nucleo_DRACH, length)))
# for each of the 1,000 sets of random DRACH+ sequences from nucleoplasm fraction extract the first mn random sequences
random_hits_nucleo_DRACH_reduced <- lapply(random_hits_nucleo_DRACH, function(x) {
  x[1:min_DRACH]
})
save(random_hits_nucleo_DRACH_reduced, file='/path/to/random_hits_nucleo_DRACH_reduced.Rda')

# overlap between each of the 1,000 sets of randomly generated DRACH+ sequences and RNA marks, of each RNA mod type, from the two databases
nucleo_mod_type_DRACH <- mod_position('/path/to/mod',random_hits_nucleo_DRACH_reduced)  
save(nucleo_mod_type_DRACH, file = '/path/to/nucleo_mod_type_DRACH.Rda')

# compute the minimum number mn of random DRACH- sequences from nucleoplasm fraction across the 1,000 sets
min_without_DRACH <- min(unlist(lapply(random_hits_nucleo_without_DRACH, length)))
# for each of the 1,000 sets of random DRACH- sequences from nucleoplasm fraction extract the first mn random sequences
random_hits_nucleo_without_DRACH_reduced <- lapply(random_hits_nucleo_without_DRACH, function(x) {
  x[1:min_without_DRACH]
})
save(random_hits_nucleo_without_DRACH_reduced, file='/path/to/random_hits_nucleo_without_DRACH_reduced.Rda')

# overlap between each of the 1,000 sets of randomly generated DRACH- sequences and RNA marks, of each RNA mod type, from the two databases
nucleo_mod_type_without_DRACH <- mod_position('/path/to/mod',random_hits_nucleo_without_DRACH_reduced)  
save(nucleo_mod_type_without_DRACH, file='/path/to/nucleo_mod_type_without_DRACH.Rda')

load('/path/to/random_hits_cyto_DRACH.Rda')
load('/path/to/random_hits_cyto_without_DRACH.Rda')

# compute the minimum number mn of random DRACH+ sequences from cytoplasm fraction across the 1,000 sets
min_DRACH <- min(unlist(lapply(random_hits_cyto_DRACH, length)))
# for each of the 1,000 sets of random DRACH+ sequences from cytoplasm fraction extract the first mn random sequences
random_hits_cyto_DRACH_reduced <- lapply(random_hits_cyto_DRACH, function(x) {
  x[1:min_DRACH]
})
save(random_hits_cyto_DRACH_reduced, file='/path/to/random_hits_cyto_DRACH_reduced.Rda')

# overlap between each of the 1,000 sets of randomly generated DRACH+ sequences and RNA marks, of each RNA mod type, from the two databases
cyto_mod_type_DRACH <- mod_position('/path/to/mod',random_hits_cyto_DRACH_reduced)
save(cyto_mod_type_DRACH, file = '/path/to/cyto_mod_type_DRACH.Rda')

# compute the minimum number mn of random DRACH- sequences from cytoplasm fraction across the 1,000 sets
min_without_DRACH <- min(unlist(lapply(random_hits_cyto_without_DRACH, length)))
# for each of the 1,000 sets of random DRACH- sequences from cytoplasm fraction extract the first mn random sequences
random_hits_cyto_without_DRACH_reduced <- lapply(random_hits_cyto_without_DRACH, function(x) {
  x[1:min_without_DRACH]
})
save(random_hits_cyto_without_DRACH_reduced, file='/path/to/random_hits_cyto_without_DRACH_reduced.Rda')

# overlap between each of the 1,000 sets of randomly generated DRACH- sequences and RNA marks, of each RNA mod type, from the two databases
cyto_mod_type_without_DRACH <- mod_position('/path/to/mod',random_hits_cyto_without_DRACH_reduced)
save(cyto_mod_type_without_DRACH, file='/path/to/cyto_mod_type_without_DRACH.Rda')


