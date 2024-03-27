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
