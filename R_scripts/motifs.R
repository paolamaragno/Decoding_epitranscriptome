library('readxl')
library('xlsx')
library('ggplot2')
library('GenomicRanges')
library('pheatmap')

load("/path/to/R_data/protein_coding_genes_5UTR_3UTR_introns_exons_stop.Rda")

# path_directory is the path to the directory containing three folders: /chr/, /nucleo/, /cyto/ each with the 5 txt files 
# produced by ELIGOS for each of the 5 samplings
# /xstreme/ is a folder in path_directory with the results of XSTREME analysis
motif_analysis <- function(path_directory) {
  
  load(paste0(path_directory,'/without_DRACH/hits_eligos_chr_ass_confirmed_5_without_DRACH.Rda'))
  load(paste0(path_directory,'/without_DRACH/hits_eligos_cyto_confirmed_5_without_DRACH.Rda'))
  load(paste0(path_directory,'/without_DRACH/hits_eligos_nucleo_confirmed_5_without_DRACH.Rda'))
  load(paste0(path_directory,'/without_DRACH/union_hits_without_DRACH.Rda'))
  
  # load the excel file with the information regarding all the selected motifs
  motifs <- read_excel(paste0(path_directory,'/xstreme/motifs.xlsx'), col_names = c('RANK','SEED','CLUSTER','SOURCE', 'ID',	'ALT_ID',	'CONSENSUS',	'WIDTH',	'SITES',	'SEA_PVALUE',	'EVALUE',	'EVALUE_ACC',	'SIM_SOURCE',	'SIM_MOTIF', 'MOTIF_URL'))
  
  # create a data frame reporting, for each fraction, the number of ELIGOS DRACH- hits containing each of the selected motifs
  seq <- rep(motifs$CONSENSUS, each=3) 
  Fraction <- factor(rep(c("Chromatin", "Nucleoplasm", "Cytoplasm"), length(motifs$CONSENSUS)), levels = c("Chromatin", "Nucleoplasm", "Cytoplasm"))
  sites <- rep(0,length(motifs$CONSENSUS)*3)
  stats <- data.frame(seq,Fraction,sites)
  
  significance <- c()
  
  # iterate over the selected motifs
  for (i in 1:length(motifs$CONSENSUS)) {
    motif <- read_excel(paste0(path_directory,'/xstreme/', motifs$CONSENSUS[i], '.xlsx'), col_names = FALSE)
    # identify which are the 20-nucleotides regions that contain that motif according to XSTREME
    seqs_motif <- union_hits_without_DRACH[unique(as.numeric(unlist(motif[,1])))]
    
    # identify which are ELIGOS DRACH- hits overlapping with the 20-nucleotides regions containing that motif
    seqs_chr_motif <- hits_eligos_chr_ass_confirmed_5_without_DRACH[unique(queryHits(findOverlaps(hits_eligos_chr_ass_confirmed_5_without_DRACH,seqs_motif,type='any')))]
    stats[which(stats$seq == motifs$CONSENSUS[i] & stats$Fraction =='Chromatin'),3]  <- length(seqs_chr_motif)
    mcols(seqs_chr_motif) <- cbind(mcols(seqs_chr_motif),consensus = motifs$CONSENSUS[i], fraction='chr') 
    save(seqs_chr_motif, file=paste0(path_directory,'/xstreme/',motifs$CONSENSUS[i],'_chr.Rda'))
    
    seqs_nucleo_motif <- hits_eligos_nucleo_confirmed_5_without_DRACH[unique(queryHits(findOverlaps(hits_eligos_nucleo_confirmed_5_without_DRACH,seqs_motif,type='any')))]
    stats[which(stats$seq == motifs$CONSENSUS[i] & stats$Fraction =='Nucleoplasm'),3]  <- length(seqs_nucleo_motif)
    mcols(seqs_nucleo_motif) <- cbind(mcols(seqs_nucleo_motif),consensus = motifs$CONSENSUS[i], fraction='nucleo') 
    save(seqs_nucleo_motif, file=paste0(path_directory,'/xstreme/',motifs$CONSENSUS[i],'_nucleo.Rda'))
    
    seqs_cyto_motif <- hits_eligos_cyto_confirmed_5_without_DRACH[unique(queryHits(findOverlaps(hits_eligos_cyto_confirmed_5_without_DRACH,seqs_motif,type='any')))]
    stats[which(stats$seq == motifs$CONSENSUS[i] & stats$Fraction =='Cytoplasm'),3] <- length(seqs_cyto_motif)
    mcols(seqs_cyto_motif) <- cbind(mcols(seqs_cyto_motif),consensus = motifs$CONSENSUS[i], fraction='cyto') 
    save(seqs_cyto_motif, file=paste0(path_directory,'/xstreme/',motifs$CONSENSUS[i],'_cyto.Rda'))
    
    if (length(seqs_chr_motif) == max(length(seqs_chr_motif), length(seqs_nucleo_motif), length(seqs_cyto_motif))) {
      significance <- c(significance,paste0('e', unlist(strsplit(motifs[motifs$CONSENSUS == motifs$CONSENSUS[i],]$EVALUE, split='e'))[2]),'','')
    } else if (length(seqs_nucleo_motif) == max(length(seqs_chr_motif), length(seqs_nucleo_motif), length(seqs_cyto_motif))) {
      significance <- c(significance,'',paste0('e', unlist(strsplit(motifs[motifs$CONSENSUS == motifs$CONSENSUS[i],]$EVALUE, split='e'))[2]),'')
    } else {
      significance <- c(significance,'','',paste0('e', unlist(strsplit(motifs[motifs$CONSENSUS == motifs$CONSENSUS[i],]$EVALUE, split='e'))[2]))
    }
    }
  
  stats <- cbind(stats, significance)
  
  order <- lapply(unique(stats$seq), function(x) median(stats$sites[stats$seq == x]))
  names(order) <- unique(stats$seq)
  order <- sort(unlist(order), decreasing = TRUE)
  
  # create a barplot that reports the number of ELIGOS DRACH- hits that contain each motif in each fraction
  p <- ggplot(data = stats, aes(fill=Fraction, y=sites, x=seq)) + 
    geom_bar(position="dodge", stat="identity") +
    theme_classic() +
    xlab('')+
    ylab('Number of DRACH- hits with the motif') +
    scale_x_discrete(limits=names(order)) +
    theme(axis.text.x = element_text(angle=45, hjust = 1)) +
    geom_text(aes(label = significance), vjust = -1)
  
  ggsave(paste0(path_directory, '/xstreme/number_sites.pdf'), plot = p, width = 13, height = 7)
  
  # function to represent the distribution in the different gene parts of ELIGOS DRACH- hits containing each motif
  heatmap_gene <- function(R_objects) {
    
    # create a unique GRanges with all the ELIGOS DRACH- hits containing at least one motif
    all_hits <- GRanges()
    for (obj in R_objects) {
      n <- load(obj) 
      hits <- get(n)
      all_hits <- c(all_hits, hits)
    }
    
    genes <- intersect(names(protein_coding_genes_5UTR_3UTR_introns_exons_stop),unique(all_hits$gene_id))
      
    # initiate a matrix that, for each motif, will report the number of ELIGOS DRACH- hits containing that
    # motif and mapping on each gene part
    m <- matrix(0, ncol = 5, nrow = length(R_objects)) 
    colnames(m) <- c('5UTR', 'coding exon', 'intron', 'stop codon', '3UTR')
    rownames(m) <- unique(all_hits$consensus) 
      
    # iterate over the selected motifs
    for (consensus in unique(all_hits$consensus)) {
      # identify which hits contain that motif
      hits_consensus <- all_hits[all_hits$consensus == consensus]
      counts <- matrix(0, ncol = 5, nrow = 1) 
      colnames(counts) <- c('5UTR', 'coding exon', 'intron', 'stop codon', '3UTR')
      
      # iterate over the hits that contain a specific motif and identify in which gene part they fall
      for (hit in 1:length(hits_consensus)) {
        gene <- hits_consensus[hit]$gene_id
        parts_gene <- protein_coding_genes_5UTR_3UTR_introns_exons_stop[gene][[1]]
        over <- findOverlaps(hits_consensus[hit], parts_gene)
        
        if (length(subjectHits(over)) != 0) {
          which_features <- unique(parts_gene[subjectHits(over)]$feature)   
          for (colname in colnames(counts)) {
            if (colname %in% which_features) {
              counts[1,colname] <- counts[1,colname] +1
            } 
          }
        }
      }
      
      m[consensus,] <- counts[1,] 
    }
    
    motif_evalue <- motifs$EVALUE
    names(motif_evalue) <- motifs$CONSENSUS
    
    if (unique(all_hits$fraction) == 'chr') {
      fraction = c('Chromatin')
    } else if (unique(all_hits$fraction) == 'nucleo') {
      fraction = c('Nucleoplasm')
    } else {
      fraction = c('Cytoplasm')
    }
    
    ann <- data.frame(motif_evalue = motif_evalue[rownames(m)])
    rownames(ann) <- rownames(m)
    print(ann)
    
    if (fraction == 'Chromatin') {
      pdf(paste0(path_directory, 'xstreme/motifs_distribution_on_gene_', fraction, '.pdf'), height = 7, width = 6)
      pheatmap(m, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, 
             annotation_row = ann, annotation_legend = FALSE, treeheight_row = 0, main = fraction, fontsize = 10, color = colorRampPalette(c("white", "red"))(50), 
               legend_breaks = round(seq(0,max(m),length.out=13),0), legend_labels = as.character(round(seq(0,max(m),length.out=13),0)), fontsize_row = 10, fontsize_col = 10,display_numbers = m)
    } else {
      pdf(paste0(path_directory, 'xstreme/motifs_distribution_on_gene_', fraction, '.pdf'), height = 7, width = 6)
      pheatmap(m, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, 
               treeheight_row = 0, main = fraction, fontsize = 10, color = colorRampPalette(c("white", "red"))(50), 
               legend_breaks = round(seq(0,max(m),length.out=13),0), legend_labels = as.character(round(seq(0,max(m),length.out=13),0)), fontsize_row = 10, fontsize_col = 10,display_numbers = m)
    }
    dev.off()
  } 
  
  R_objects <- list.files(path = paste0(path_directory,'xstreme'), pattern = 'Rda', full.names = TRUE)
  R_objects_chr <- R_objects[unlist(lapply(strsplit(R_objects,'_'), function(x) 'chr.Rda' %in% x))]
  heatmap_gene(R_objects_chr)
  R_objects_nucleo <- R_objects[unlist(lapply(strsplit(R_objects,'_'), function(x) 'nucleo.Rda' %in% x))]
  heatmap_gene(R_objects_nucleo)
  R_objects_cyto <- R_objects[unlist(lapply(strsplit(R_objects,'_'), function(x) 'cyto.Rda' %in% x))]
  heatmap_gene(R_objects_cyto)
}

#################

# in the file /sea_out/sequences.tsv, for each motif identified as enriched by SEA, there is the id of the input sequences
# in which that motif is present. Save all the columns but the first in a xlsx file (as CGGAGGR.xlsx) 

t <- read_xlsx('/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/xstreme/CGGAGGR.xlsx',col_names = FALSE)

t <- t[-which(t$...3 =='fp'),]
t <- t[-which(t$...4 ==1),]
t <- t[,1]
write.xlsx(x = data.frame(t),file = '/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/xstreme/CGGAGGR.xlsx',col.names = FALSE )

t <- read_xlsx('/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/xstreme/GAAGGAA.xlsx',col_names = FALSE)

t <- t[-which(t$...3 =='fp'),]
t <- t[-which(t$...4 ==1),]
t <- t[,1]
write.xlsx(x = data.frame(t),file = '/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/xstreme/GAAGGAA.xlsx',col.names = FALSE )

t <- read_xlsx('/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/xstreme/AGAGARR.xlsx',col_names = FALSE)

t <- t[-which(t$...3 =='fp'),]
t <- t[-which(t$...4 ==1),]
t <- t[,1]
write.xlsx(x = data.frame(t),file = '/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/xstreme/AGAGARR.xlsx',col.names = FALSE )

t <- read_xlsx('/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/xstreme/CCUYCCC.xlsx',col_names = FALSE)

t <- t[-which(t$...3 =='fp'),]
t <- t[-which(t$...4 ==1),]
t <- t[,1]
write.xlsx(x = data.frame(t),file = '/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/xstreme/CCUYCCC.xlsx',col.names = FALSE )

t <- read_xlsx('/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/xstreme/UGGGGAU.xlsx',col_names = FALSE)

t <- t[-which(t$...4 =='fp'),]
t <- t[-which(t$...5 ==1),]
t <- t[,2]
write.xlsx(x = data.frame(t),file = '/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/xstreme/UGGGGAU.xlsx',col.names = FALSE )

t <- read_xlsx('/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/xstreme/DUAGGGA.xlsx',col_names = FALSE)

t <- t[-which(t$...3 =='fp'),]
t <- t[-which(t$...4 ==1),]
t <- t[,1]
write.xlsx(x = data.frame(t),file = '/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/xstreme/DUAGGGA.xlsx',col.names = FALSE )

t <- read_xlsx('/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/xstreme/AGAAGAN.xlsx',col_names = FALSE)

t <- t[-which(t$...3 =='fp'),]
t <- t[-which(t$...4 ==1),]
t <- t[,1]
write.xlsx(x = data.frame(t),file = '/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/xstreme/AGAAGAN.xlsx',col.names = FALSE )

t <- read_xlsx('/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/xstreme/GCUGGMC.xlsx',col_names = FALSE)

t <- t[-which(t$...3 =='fp'),]
t <- t[-which(t$...4 ==1),]
t <- t[,1]
write.xlsx(x = data.frame(t),file = '/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/xstreme/GCUGGMC.xlsx',col.names = FALSE )

t <- read_xlsx('/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/xstreme/AGGAGCA.xlsx',col_names = FALSE)

t <- t[-which(t$...3 =='fp'),]
t <- t[-which(t$...4 ==1),]
t <- t[,1]
write.xlsx(x = data.frame(t),file = '/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/xstreme/AGGAGCA.xlsx',col.names = FALSE )

motif_analysis('/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/')

