library('Biostrings')
library('GenomicAlignments')
library('GenomicRanges')
library('GenomicFeatures')

# path to the genome reference file in fasta format
genome <- readDNAStringSet('/path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa', format="fasta")
names(genome) <- unlist(lapply(strsplit(names(genome), split=' '), function(x) x[1]))
# look for DRACH motif on the genome positive strand
DRACH_forward <- GRanges(vmatchPattern(pattern = "DRACH", subject = genome, fixed = 'subject'))
strand(DRACH_forward) <- '+'
# look for the reverse complement of the DRACH motif on the genome positive strand (corresponding to DRACH
# sites on the genome negative strand)
DRACH_reverse <- GRanges(vmatchPattern(pattern = reverseComplement(DNAString('DRACH')), subject = genome, fixed = 'subject'))
strand(DRACH_reverse) <- '-'
DRACH <- c(DRACH_forward,DRACH_reverse)
save(DRACH, file='/path/to/R_data/DRACH_forward_strand.Rda')

##########

# load the tsv file with the output of m6Anet for the three fractions for all the 10 samplings
m6Anet_fractions_4sU_chr_ass <- list.files(path = '/path/to/fractions_m6anet_4sU_library_gene_subsampling_prob0.9/chr', pattern = 'tsv', full.names = TRUE)
m6Anet_fractions_4sU_nucleo <- list.files(path = '/path/to/fractions_analysis_Paola_SUM159/fractions_m6anet_4sU_library_gene_subsampling_prob0.9/nucleo', pattern = 'tsv', full.names = TRUE)
m6Anet_fractions_4sU_cyto <- list.files(path = '/path/to/fractions_analysis_Paola_SUM159/fractions_m6anet_4sU_library_gene_subsampling_prob0.9/cyto', pattern = 'tsv', full.names = TRUE)

# load the vector with the transcript names as names and the names of the gene from which 
# they are transcribed as values
load('/path/to/R_data/tx_gene.Rda')

# extract, from each tsv file, the sites that have been analysed by m6Anet
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
                  rep = unlist(strsplit(strsplit(file, split = '/')[[1]][7], "_")[[1]][5]))
    
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

gtf_file <- "/path/to/Homo_sapiens.GRCh38.104.gtf"
txdb <- makeTxDbFromGFF(gtf_file)
tx_txdb <- GenomicFeatures::transcripts(txdb)

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

# compute the percentage of m6Anet analysed sites containing the DRACH motif
chr_ass_m6Anet_DRACH <- findOverlaps(gr_m6Anet_chr_ass_4sU, DRACH, type='any', ignore.strand=FALSE)
length(unique(queryHits(chr_ass_m6Anet_DRACH)))/length(gr_m6Anet_chr_ass_4sU)

nucleo_m6Anet_DRACH <- findOverlaps(gr_m6Anet_nucleo_4sU, DRACH, type='any')
length(unique(queryHits(nucleo_m6Anet_DRACH)))/length(gr_m6Anet_nucleo_4sU)

cyto_m6Anet_DRACH <- findOverlaps(gr_m6Anet_cyto_4sU, DRACH, type='any')
length(unique(queryHits(cyto_m6Anet_DRACH)))/length(gr_m6Anet_cyto_4sU)


