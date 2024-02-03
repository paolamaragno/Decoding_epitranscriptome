library('GenomicFeatures')
library('GenomicRanges')

# function that acts on the single GRanges relative to a gene reporting all its cds (from all its transcripts) 
# and produces a new GRanges that, for each transcript of that gene, reports its last cds 
f <- function(x) {
  tx_unique <- unique(names(x))
  gr_new <- GRanges()
  
  if (unique(strand(x)) == '+') {
    for (tx in tx_unique) {
      x_tx <- x[which(names(x) == tx)]
      x_tx <- x_tx[which(end(ranges(x_tx)) == max(end(ranges(x_tx))))]
      gr_new <- append(gr_new, x_tx)
    }
  } else {
    for (tx in tx_unique) {
      x_tx <- x[which(names(x) == tx)]
      x_tx <- x_tx[which(start(ranges(x_tx)) == min(start(ranges(x_tx))))]
      gr_new <- append(gr_new, x_tx)
    }
  }
  return(gr_new)
}

# function applied on each GRanges relative to a single gene and reporting the last cds of each single transcript annotated to 
# that gene. The stop codon is extracted as the 200 nucleotides around the end of the last cds
f_resize <- function(x) {
  if (unique(strand(x))=='+') {
    for (i in 1:length(x)) {
      start(x[i]) <- end(x[i]) - 100
      end(x[i]) <- end(x[i]) + 100
    } 
  } else {
    for (i in 1:length(x)) {
      start(x[i]) <- start(x[i]) - 100
      end(x[i]) <- start(x[i]) + 200
    }
  }
  return(x)
}

########################
gtf_file <- "/path/to/Homo_sapiens.GRCh38.104.gtf"
txdb <- makeTxDbFromGFF(gtf_file)

genes_txdb <- GenomicFeatures::genes(txdb)
tx <- transcriptsBy(txdb, by = "gene")

# create a vector with the transcript names as names and the names of the gene from which they are transcribed as values
tmp <- lapply(seq_along(tx), function(y, n, i) {
  cbind(y[[i]]$tx_name, n[i]) 
}, y = tx, n = names(tx))
tx_gene <- do.call(rbind, tmp)
rownames(tx_gene) <- tx_gene[, 1]
tx_gene <- tx_gene[, 2]
save(tx_gene, file="path/to/output_dir/tx_gene.RDa")

# EXONS
exon_coord <- exonsBy(txdb, by = 'gene')
# merge the overlapping ranges 
exon_coord_red <- lapply(exon_coord, reduce)
save(exon_coord_red, file="/path/to/output_dir/exon_coord_red.RDa")

# create a vector with the exon names as names and the names of the gene to which they belong as values
tmp <- lapply(seq_along(exon_coord), function(y, n, i) {
  cbind(y[[i]]$exon_name, n[i]) 
}, y = exon_coord, n = names(exon_coord))
exon_gene <- do.call(rbind, tmp)
rownames(exon_gene) <- exon_gene[, 1]
exon_gene <- exon_gene[, 2]
save(exon_gene, file="/path/to/output_dir/exon_gene.RDa")

# 5'UTR
five_UTR_coord <- fiveUTRsByTranscript(txdb)
# create a vector with the names of the gene to which each 5'UTR belongs, 
# keeping also the repeated genes in case more consecutive 5'UTRs belong to the same gene
five_UTR_coord_genes_names <- lapply(five_UTR_coord, function(x) exon_gene[mcols(x)$exon_name])
# split the vector of 5'UTR names on the base of the the vector of gene names
five_UTR_coord_bygene <- split(unlist(five_UTR_coord), unname(unlist(five_UTR_coord_genes_names)))
# merge the overlapping ranges
five_UTR_coord_bygene_red <- lapply(five_UTR_coord_bygene, reduce)
five_UTR_coord_bygene_red <- lapply(five_UTR_coord_bygene_red, function(x) {
  mcols(x) <- cbind(mcols(x), feature='5UTR')
  x
})
save(five_UTR_coord_bygene_red, file="/path/to/output_dir/five_UTR_coord_bygene_red.RDa")

# 3'UTR
three_UTR_coord <- threeUTRsByTranscript(txdb)
# create a vector with the names of the gene to which each 3'UTR belongs, 
# keeping also the repeated genes in case more consecutive 3'UTRs belong to the same gene
three_UTR_coord_genes_names <- lapply(three_UTR_coord, function(x) exon_gene[mcols(x)$exon_name])
# split the vector of 3'UTR names on the base of the the vector of gene names
three_UTR_coord_bygene <- split(unlist(three_UTR_coord), unname(unlist(three_UTR_coord_genes_names)))
# merge the overlapping ranges
three_UTR_coord_bygene_red <- lapply(three_UTR_coord_bygene, reduce)
three_UTR_coord_bygene_red <- lapply(three_UTR_coord_bygene_red, function(x) {
  mcols(x) <- cbind(mcols(x), feature='3UTR')
  x
})
save(three_UTR_coord_bygene_red, file="/path/to/output_dir/three_UTR_coord_bygene_red.RDa")

# CODING EXONS
# coding exons are defined as the exon regions that are neither 5'UTR nor 3'UTR
coding_exons <- lapply(seq_along(exon_coord_red), function(y, n, i, f, t) {
  gene <- n[i]
  if (gene %in% names(f)) {
    coding <- setdiff(y[[i]], f[[gene]])
  } else {coding <- y[[i]]}
  if (gene %in% names(t)) {
    coding2 <- setdiff(coding, t[[gene]])
  } else {coding2 <- coding}
  return(coding2)
}, y = exon_coord_red, n = names(exon_coord_red), f = five_UTR_coord_bygene_red, t = three_UTR_coord_bygene_red)

names(coding_exons) <- names(exon_coord_red)

for (i in 1:length(coding_exons)) {
  if (length(coding_exons[[i]])==0) {
    print(i)
  }
}
coding_exons <- coding_exons[-10597]
coding_exons <- lapply(coding_exons, function(x) {
  mcols(x) <- cbind(mcols(x), feature='coding exon')
  x
})
save(coding_exons, file="/path/to/output_dir/coding_exons.RDa")

# STOP CODON
cds_coord <- cdsBy(txdb, by = 'tx',use.names=TRUE)

# create a vector with the names of the gene to which each coding region belongs, 
# keeping also the repeated genes in case more consecutive coding regions belong to the same gene
tmp <- lapply(seq_along(cds_coord), function(y, n, i) {
  rep(tx_gene[n[i]], length(y[[i]]))
}, y = cds_coord, n = names(cds_coord))
gene_name_for_split <- unlist(tmp)

# split the vector of coding regions' names on the base of the the vector of gene names
cds_coord_bygene <- split(unlist(cds_coord), unname(gene_name_for_split))
# extract for each transcript annotated to each gene the most 3' coding region
cds_coord_bygene_last_cds_per_tx <- lapply(cds_coord_bygene, f)
# extract for each transcript annotated to each gene the 200 nucleotides around the end of its most 3' coding region
stop_codon_coord <- lapply(cds_coord_bygene_last_cds_per_tx, f_resize) 
# merge overlapping ranges
stop_codon_coord_red <- lapply(stop_codon_coord, reduce)
stop_codon_coord_red <- lapply(stop_codon_coord_red, function(x) {
  mcols(x) <- cbind(mcols(x), feature='stop codon')
  x
})
save(stop_codon_coord_red, file="/path/to/output_dir/stop_codon_coord_red.RDa")

# INTRONS
# introns are defined as the transcript regions that are never exons (coding exons/5'UTR/3'UTR)

# for each gene extract the start coordinate of its transcript that is the most 5' and 
# end coordinate of its transcript that is the most 3'
tx_start_end <- lapply(tx, function(x) {
  start <- min(start(x))
  end <- max(end(x))
  GRanges(seqnames = seqnames(x)[1], ranges = IRanges(start=start, end = end), strand = strand(x)[1])
})

introns <- lapply(seq_along(tx_start_end), function(y, n, i, f) {
  gene <- n[i]
  tx_overall_region <- y[[gene]]
  exons <- f[[gene]]
  return(setdiff(tx_overall_region, exons))
}, y = tx_start_end, n = names(tx_start_end), f = exon_coord_red)

names(introns) <- names(tx_start_end)
# merge the overlapping ranges 
introns_coord_red <- lapply(introns, reduce)

# remove the genes without introns
not_null <- lapply(introns_coord_red, function(x) {length(x) != 0})
introns_coord_red <- introns_coord_red[unname(unlist(not_null))]
introns_coord_red <- lapply(introns_coord_red, function(x) {
  mcols(x) <- cbind(mcols(x), feature='intron')
  x
})
save(introns_coord_red, file="/path/to/output_dir/introns_coord_red.RDa")

all_genes_annotated <- c(names(five_UTR_coord_bygene_red), names(three_UTR_coord_bygene_red), names(stop_codon_coord_red),
                         names(coding_exons), names(introns_coord_red))
all_genes_annotated <- unique(all_genes_annotated)

# for each annotated gene create a new GRanges object containing the coordinates of its 5'UTRs, 3'UTRs, stop codon, introns and exons
protein_coding_genes_5UTR_3UTR_introns_exons_stop <- c()
for (gene in all_genes_annotated) {
  gr <- c()
  if (gene %in% names(five_UTR_coord_bygene_red)) {
    five_utr <- five_UTR_coord_bygene_red[gene][[1]]
    gr <- c(gr,five_utr)
  } 
  if (gene %in% names(coding_exons)) {
    exons <- coding_exons[gene][[1]]
    gr <- c(gr, exons)
  }
  if (gene %in% names(introns_coord_red)) {
    introns <- introns_coord_red[gene][[1]]
    gr <- c(gr, introns)
  }
  if (gene %in% names(stop_codon_coord_red)) {
    stop_cod <- stop_codon_coord_red[gene][[1]]
    gr <- c(gr,stop_cod)
  }
  if (gene %in% names(three_UTR_coord_bygene_red)) {
    three_utr <- three_UTR_coord_bygene_red[gene][[1]]
    gr <- c(gr,three_utr)
  }
  if (length(gr) > 0) {
    gr <- unlist(as(gr, "GRangesList"))
    protein_coding_genes_5UTR_3UTR_introns_exons_stop <- c(protein_coding_genes_5UTR_3UTR_introns_exons_stop, gr)
  }
}

# assign to each GRanges of the list the name of the corresponding gene
names(protein_coding_genes_5UTR_3UTR_introns_exons_stop) <- all_genes_annotated

# convert into GRangesList
protein_coding_genes_5UTR_3UTR_introns_exons_stop <- as(protein_coding_genes_5UTR_3UTR_introns_exons_stop,'GRangesList')

save(protein_coding_genes_5UTR_3UTR_introns_exons_stop, file="/path/to/output_dir/protein_coding_genes_5UTR_3UTR_introns_exons_stop.RDa")
