# from the alignment of all the reads of each fraction on Human genome only the library-level
# subsampling is performed. Then, for the genes with coverage at least 20X in all the fractions, 
# the Spearman correlation between the number of reads mapping on them in each pairwise comparison between fractions
# is computed. The hierarchical clustering is done between the different fractions based on the number
# of reads mapping, in each fraction, on the genes with coverage > 20X

library('GenomicFeatures')
library('GenomicAlignments')

# path to GTF annotation file
gtf_file <- "/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/references/Homo_sapiens.GRCh38.104.gtf"
txdb <- makeTxDbFromGFF(gtf_file)
genes_txdb <- GenomicFeatures::genes(txdb)

# import, for each fraction, the filtered bam file with the reads of all the replicates of that fraction
bam_chr <- readGAlignments('/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/chromatin_associated_PASS_rep1_rep2_rep3_rep4_4sU_spliced_mapped_to_Homo_sapiens.GRCh38.dna.primary_assembly_filtered.bam', use.names=TRUE)
bam_nucleo <- readGAlignments('/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/nucleoplasmic_PASS_rep1_rep2_4sU_spliced_mapped_to_Homo_sapiens.GRCh38.dna.primary_assembly_filtered.bam', use.names=TRUE)
bam_cyto <- readGAlignments('/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/cytoplasmic_PASS_rep1_rep2_4sU_spliced_mapped_to_Homo_sapiens.GRCh38.dna.primary_assembly_filtered.bam', use.names=TRUE)

# overlap between the coordinates of each mapping read and the coordinates of all the genes
chr_ass_genes_over <- findOverlaps(bam_chr,genes_txdb)
nucleo_genes_over <- findOverlaps(bam_nucleo,genes_txdb)
cyto_genes_over <- findOverlaps(bam_cyto,genes_txdb)

# remove all the reads that map on multiple genes
bam_chr <- bam_chr[queryHits(chr_ass_genes_over)][isUnique(queryHits(chr_ass_genes_over))]
bam_nucleo <- bam_nucleo[queryHits(nucleo_genes_over)][isUnique(queryHits(nucleo_genes_over))]
bam_cyto <- bam_cyto[queryHits(cyto_genes_over)][isUnique(queryHits(cyto_genes_over))]

# count the minimum number of mapping reads across the three fractions
threshold_library_subsampling <- min(length(bam_chr), length(bam_nucleo), length(bam_cyto))

# for each fraction perform the subsampling at library level
set.seed(1)
bam_chr_after_library_sub <- bam_chr[sample(names(bam_chr), threshold_library_subsampling),]
set.seed(1)
bam_nucleo_after_library_sub <- bam_nucleo[sample(names(bam_nucleo), threshold_library_subsampling),]
set.seed(1)
bam_cyto_after_library_sub <- bam_cyto[sample(names(bam_cyto), threshold_library_subsampling),]

# overlap between the coordinates of each selected mapping read and the coordinates of all the genes
chr_ass_genes_over2 <- findOverlaps(bam_chr_after_library_sub,genes_txdb)
nucleo_genes_over2 <- findOverlaps(bam_nucleo_after_library_sub,genes_txdb)
cyto_genes_over2 <- findOverlaps(bam_cyto_after_library_sub,genes_txdb)

# extract the names of the genes on which the reads overlap
gene_names_chr_ass <- names(genes_txdb[subjectHits(chr_ass_genes_over2)])
gene_names_nucleo <- names(genes_txdb[subjectHits(nucleo_genes_over2)])
gene_names_cyto <- names(genes_txdb[subjectHits(cyto_genes_over2)])

# split the reads on the base of the gene on which they map
chr_ass_genes_over_split <- split(queryHits(chr_ass_genes_over2), gene_names_chr_ass)
nucleo_genes_over_split <- split(queryHits(nucleo_genes_over2), gene_names_nucleo)
cyto_genes_over_split <- split(queryHits(cyto_genes_over2), gene_names_cyto)

# count the number of reads on each gene
chr_ass_number_reads_per_gene <- unlist(lapply(chr_ass_genes_over_split, length))
nucleo_number_reads_per_gene <- unlist(lapply(nucleo_genes_over_split, length))
cyto_number_reads_per_gene <- unlist(lapply(cyto_genes_over_split, length))

# identify all the genes
all_genes <- union(names(chr_ass_number_reads_per_gene), names(cyto_number_reads_per_gene))
all_genes <- union(all_genes, names(nucleo_number_reads_per_gene))

# create a matrix with a row per gene and three columns, one for each fraction (initialize to 0)
genes_counts <- matrix(data = 0, nrow = length(all_genes), ncol=3)
rownames(genes_counts) <- all_genes
colnames(genes_counts) <- c("Chromatin\nAssociated", "Nucleoplasmic", "Cytoplasmic")

# update the number of reads mapping on each gene in each fraction (letting 0 when there are no reads mapping on the gene)
genes_counts[names(chr_ass_number_reads_per_gene), "Chromatin\nAssociated"] <- unname(chr_ass_number_reads_per_gene)
genes_counts[names(nucleo_number_reads_per_gene), "Nucleoplasmic"] <- unname(nucleo_number_reads_per_gene)
genes_counts[names(cyto_number_reads_per_gene), "Cytoplasmic"] <- unname(cyto_number_reads_per_gene)

# keep the genes with at least 20 mapping reads in all the fractions
genes_counts_20 <- genes_counts[which(genes_counts[,1]>=20 & genes_counts[,2]>=20 & genes_counts[,3]>=20),]

# compute Spearman correlation between the number of reads per gene considering each pair of fractions
cor_chr_cyto <- cor(genes_counts_20[,1],genes_counts_20[,3], method = "s")
cor_chr_nucleo <- cor(genes_counts_20[,1],genes_counts_20[,2], method = "s")
cor_nucleo_cyto <- cor(genes_counts_20[,2],genes_counts_20[,3], method = "s")

# scatter plot of all pairwise comparisons between the number of reads per gene
pdf('/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/sp_cor_pairs_of_fractions.pdf')
par(mfrow = c(2,2))

plot(log = 'xy', x=genes_counts_20[,3], y=genes_counts_20[,1], xlab = c('Number of reads from cytoplasmic RNAs'), ylab = c('Number of reads from chr. associated RNAs'), main = paste0("Cytoplasmic VS Chromatin Associated\n", "rSp. = ", sprintf("%.2f", cor_chr_cyto)), cex=0.1, cex.main=1)
abline(0, 1, col = 2)

plot(log = 'xy', x=genes_counts_20[,3],y=genes_counts_20[,2], xlab = c('Number of reads from cytoplasmic RNAs'), ylab = c('Number of reads from nucleoplasmic RNAs'), main = paste0("Cytoplasmic VS Nucleoplasmic\n", "rSp. = ", sprintf("%.2f", cor_nucleo_cyto)), cex=0.1, cex.main=1)
abline(0, 1, col = 2)

plot(log = 'xy',x=genes_counts_20[,1],y=genes_counts_20[,2], xlab = c('Number of reads from chr. associated RNAs'), ylab = c('Number of reads from nucleoplasmic RNAs'), main = paste0("Chromatin Associated VS Nucleoplasmic\n", "rSp. = ", sprintf("%.2f", cor_chr_nucleo)), cex=0.1, cex.main=1)
abline(0, 1, col = 2)

dev.off()

# hierarchical clustering 
pdf('/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/Hclust_number_reads_per_gene_after_library_sub.pdf')
plot(hclust(as.dist(1 - cor(genes_counts_20, use = "na.or.complete", method = "spearman"))), xlab = "", ylab = "1-Spearman_correlation", main = "Hierarchical clustering between the number of reads per gene\nin each fraction after the library-level subsampling")
dev.off()
