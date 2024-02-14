library('GenomicFeatures')
library('GenomicRanges')

# download the full tables that you can find at these links: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6801 and 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM967817

# SNPs SUM159 on hg19
snps_SUM <- read.table('GSM967817-19203.txt', fill=TRUE, header = TRUE)[,c(1,2)]
snps_SUM <- snps_SUM[snps_SUM$VALUE != 'NoCall',]

coordinates_snps_SUM <- read.table('GPL6801-4019.txt', fill=TRUE, header = TRUE, sep = '\t')
coordinates_snps_SUM <- coordinates_snps_SUM[coordinates_snps_SUM$Chromosome != '---',]
snps_SUM_final <- coordinates_snps_SUM[coordinates_snps_SUM$ID %in% snps_SUM$ID_REF, c(2,7,7,9,10)]

snps_SUM_final <- data.frame(chromosome = snps_SUM_final$Chromosome, start=snps_SUM_final$Physical.Position, stop=snps_SUM_final$Physical.Position.1, alleleA = snps_SUM_final$Allele.A, alleleB=snps_SUM_final$Allele.B)
write.table(snps_SUM_final, '/path/to/output/SNPs_SUM.bed', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


