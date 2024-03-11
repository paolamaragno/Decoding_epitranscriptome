library('GenomicRanges')
library('GenomicFeatures')

# load, for each fraction, ELIGOS DRACH+ and DRACH- hits 
load('/path/to/hits_eligos_chr_ass_confirmed_5_without_DRACH.Rda')
load('/path/to/hits_eligos_cyto_confirmed_5_without_DRACH.Rda')
load('/path/to/hits_eligos_nucleo_confirmed_5_without_DRACH.Rda')
load('/path/to/hits_eligos_chr_ass_confirmed_5_with_DRACH.Rda')
load('/path/to/hits_eligos_nucleo_confirmed_5_with_DRACH.Rda')
load('/path/to/hits_eligos_cyto_confirmed_5_with_DRACH.Rda')

# compute the overall number of ELIGOS hits (both DRACH+ and DRACH-) from each fraction
tot_hits_chr <- length(hits_eligos_chr_ass_confirmed_5_without_DRACH) +length(hits_eligos_chr_ass_confirmed_5_with_DRACH)
tot_hits_nucleo <- length(hits_eligos_nucleo_confirmed_5_without_DRACH) +length(hits_eligos_nucleo_confirmed_5_with_DRACH)
tot_hits_cyto <- length(hits_eligos_cyto_confirmed_5_without_DRACH) +length(hits_eligos_cyto_confirmed_5_with_DRACH)

# identify, for each fraction, the genes on which at least one ELIGOS hit maps
genes_chr <- unique(c(unique(hits_eligos_chr_ass_confirmed_5_with_DRACH$gene_id), unique(hits_eligos_chr_ass_confirmed_5_without_DRACH$gene_id)))
genes_nucleo <- unique(c(unique(hits_eligos_nucleo_confirmed_5_with_DRACH$gene_id),unique(hits_eligos_nucleo_confirmed_5_without_DRACH$gene_id)))
genes_cyto <- unique(c(unique(hits_eligos_cyto_confirmed_5_with_DRACH$gene_id), unique(hits_eligos_cyto_confirmed_5_without_DRACH$gene_id)))

load('/path/to/exon_coord_red.RDa')
# for each fraction identify the coordinates of the exons of the genes on which at least one ELIGOS hit maps
coordinates_exons_genes_hits_chr <- exon_coord_red[genes_chr]
# remove the exons with a width lower than 10 nucleotides
coordinates_exons_genes_hits_chr <- lapply(coordinates_exons_genes_hits_chr, function(x) {
  x[width(x) >= 10] 
})
coordinates_exons_genes_hits_chr <- unlist(as(coordinates_exons_genes_hits_chr, 'GRangesList'))
# create a vector reporting, for each exon of each selected gene, the chr_start_end_strand
coordinates_exons_genes_vector_chr <- paste0(seqnames(coordinates_exons_genes_hits_chr), "_", start(coordinates_exons_genes_hits_chr), "_", end(coordinates_exons_genes_hits_chr), '_', strand(coordinates_exons_genes_hits_chr))

coordinates_exons_genes_hits_nucleo <- exon_coord_red[genes_nucleo]
coordinates_exons_genes_hits_nucleo <- lapply(coordinates_exons_genes_hits_nucleo, function(x) {
  x[width(x) >= 10] 
})
coordinates_exons_genes_hits_nucleo <- unlist(as(coordinates_exons_genes_hits_nucleo, 'GRangesList'))
coordinates_exons_genes_vector_nucleo <- paste0(seqnames(coordinates_exons_genes_hits_nucleo), "_", start(coordinates_exons_genes_hits_nucleo), "_", end(coordinates_exons_genes_hits_nucleo), '_', strand(coordinates_exons_genes_hits_nucleo))

coordinates_exons_genes_hits_cyto <- exon_coord_red[genes_cyto]
coordinates_exons_genes_hits_cyto <- lapply(coordinates_exons_genes_hits_cyto, function(x) {
  x[width(x) >= 10] 
})
coordinates_exons_genes_hits_cyto <- unlist(as(coordinates_exons_genes_hits_cyto, 'GRangesList'))
coordinates_exons_genes_vector_cyto <- paste0(seqnames(coordinates_exons_genes_hits_cyto), "_", start(coordinates_exons_genes_hits_cyto), "_", end(coordinates_exons_genes_hits_cyto), '_', strand(coordinates_exons_genes_hits_cyto))

# function to generate, for each fraction, for 1,000 times a number of random sequences equal to
# the overall number of ELIGOS hits (both DRACH+ and DRACH-) from the same fraction
generate_random_hits <- function(number_hits, fraction) {
  
  random_hits_DRACH <- c()
  random_hits_without_DRACH <- c()
  
  if (fraction == 'chr') {
    coordinates_exons_genes_vector <- coordinates_exons_genes_vector_chr
  } else if (fraction == 'nucleo') {
    coordinates_exons_genes_vector <- coordinates_exons_genes_vector_nucleo
  } else {
    coordinates_exons_genes_vector <- coordinates_exons_genes_vector_cyto
  }
  
  # generate a number of random sequences of 10 nucleotides equal to the overall number of ELIGOS hits
  for (i in 1:1000){
    set.seed(i)
    # randomly sample a number of exons - of the genes on which at least one ELIGOS hit maps in that fraction -
    # equal to the overall number of ELIGOS hits (both DRACH+ and DRACH-) from the same fraction
    random_coordinates <- sample(coordinates_exons_genes_vector,number_hits, replace = TRUE)
    hits <- unlist(lapply(as.list(random_coordinates), function(x){
      chr <- unlist(strsplit(x, split='_'))[1]
      start <- as.numeric(unlist(strsplit(x, split='_'))[2])
      end <- as.numeric(unlist(strsplit(x, split='_'))[3])
      strand <- unlist(strsplit(x, split='_'))[4]
      # extract a random position between the starting position of the exon and the final position - 9
      start_random_hit <- sample(start:(end-9), 1)
      return(paste0(chr,'_', as.character(start_random_hit), '_', as.character(start_random_hit+9), '_', strand))
    }))
    
    gr <- GRanges(seqnames = unlist(lapply(as.list(hits), function(x) unlist(strsplit(x, split='_'))[1])),
                  ranges = IRanges(as.numeric(unlist(lapply(as.list(hits), function(x) unlist(strsplit(x, split='_'))[2]))),
                                   as.numeric(unlist(lapply(as.list(hits), function(x) unlist(strsplit(x, split='_'))[3])))),
                  strand = unlist(lapply(as.list(hits), function(x) unlist(strsplit(x, split='_'))[4])))
    
    # identify which random hits contain the DRACH motif and which not
    overlap_DRACH <- findOverlaps(gr, DRACH, type='any', ignore.strand=FALSE)
    random_hits_DRACH <- c(random_hits_DRACH, gr[unique(queryHits(overlap_DRACH))])
    random_hits_without_DRACH <- c(random_hits_without_DRACH, gr[-unique(queryHits(overlap_DRACH))])
  }
  
  l <- list(random_hits_DRACH,random_hits_without_DRACH)
  names(l) <- c('random_hits_DRACH','random_hits_without_DRACH')
  return(l)
  
}

load('/path/to/DRACH_forward_strand.Rda')

# for each fraction generate for 1,000 times a number of random sequences equal to
# the overall number of ELIGOS hits from the same fraction 
random_hits_chr <- generate_random_hits(tot_hits_chr, 'chr')
random_hits_chr_DRACH <- random_hits_chr[[1]]
random_hits_chr_without_DRACH <- random_hits_chr[[2]]
save(random_hits_chr_DRACH, file='/path/to/random_hits_chr_DRACH.Rda')
save(random_hits_chr_without_DRACH, file='/path/to/random_hits_chr_without_DRACH.Rda')

random_hits_nucleo <- generate_random_hits(tot_hits_nucleo, 'nucleo')
random_hits_nucleo_DRACH <- random_hits_nucleo[[1]]
random_hits_nucleo_without_DRACH <- random_hits_nucleo[[2]]
save(random_hits_nucleo_DRACH, file='/path/to/random_hits_nucleo_DRACH.Rda')
save(random_hits_nucleo_without_DRACH, file='/path/to/random_hits_nucleo_without_DRACH.Rda')

random_hits_cyto <- generate_random_hits(tot_hits_cyto, 'cyto')
random_hits_cyto_DRACH <- random_hits_cyto[[1]]
random_hits_cyto_without_DRACH <- random_hits_cyto[[2]]
save(random_hits_cyto_DRACH, file='/path/to/random_hits_cyto_DRACH.Rda')
save(random_hits_cyto_without_DRACH, file='/path/to/random_hits_cyto_without_DRACH.Rda')
