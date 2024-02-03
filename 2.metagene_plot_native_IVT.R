# generate a metagene plot showing the coverage profile of the reads from the RNAs of one of the
# three fractions of SUM159 cells and the coverage profile of the reads from IVT RNAs

library(Rsamtools)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
library(Guitar)
library(ggplot2)

read_bam_file <- function(file_path) {
  galign <- readGAlignments(file_path)
  granges <- as(galign, "GRanges")
  return(granges)
}

# Function to create a list of GRanges objects from multiple bam files
create_granges_list <- function(file_paths) {
  granges_list <- GRangesList()
  for (file_path in file_paths) {
    granges <- read_bam_file(file_path)
    granges_list[[file_path]] <- granges
  }
  return(granges_list)
}

# path to the GTF annotation file
gtf_file<-"/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/references/Homo_sapiens.GRCh38.104.gtf"
txdb_from_gtf <- makeTxDbFromGFF(file=gtf_file)

# Specify the file paths of the bam files
bam_files <- c("/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/coverage_native_IVT/chromatin_associated_PASS_rep1_rep2_rep3_rep4_4sU_spliced_mapped_to_Homo_sapiens.GRCh38.dna.primary_assembly_filtered.bam",
               "/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/coverage_native_IVT/IVT_gridION_PASS_all_reads_filtered.bam")

# Specify the file labels
bam_labels <- c("Chromatin SUM159", "IVT SUM159")

# Create the list of GRanges objects
granges_list <- create_granges_list(bam_files)

metagene_plot_of_coverage <- GuitarPlot(stGRangeLists = granges_list,
                                        stGroupName = bam_labels,
                                        txTxdb = txdb_from_gtf,
                                        #pltTxType="mrna" allows GuitarPlot to map reads onto 5'UTR,CSD,and 3'UTR
                                        #pltTxType="ncrna" allows GuitarPlot to map reads onto ncRNA
                                        pltTxType = "mrna",
                                        #Change to TRUE if plotting "ncrna"
                                        headOrtail = TRUE,
                                        enableCI = FALSE,
                                        legend()) 


pdf("/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/coverage_native_IVT/coverage_native_IVT_chr.pdf")
metagene_plot_of_coverage + 
  scale_fill_manual(values=rep(c("transparent"),length(bam_labels))) +
  theme_void() +
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) 
#then, you can add title, specify color using just ggplot2 functions
dev.off()


# Specify the file paths of the bam files
bam_files <- c("/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/coverage_native_IVT/nucleoplasmic_PASS_rep1_rep2_4sU_spliced_mapped_to_Homo_sapiens.GRCh38.dna.primary_assembly_filtered.bam",
               "/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/coverage_native_IVT/IVT_gridION_PASS_all_reads_filtered.bam")

# Specify the file labels
bam_labels <- c("Nucleoplasm SUM159", "IVT SUM159")

# Create the list of GRanges objects
granges_list <- create_granges_list(bam_files)

metagene_plot_of_coverage <- GuitarPlot(stGRangeLists = granges_list,
                                        stGroupName = bam_labels,
                                        txTxdb = txdb_from_gtf,
                                        #pltTxType="mrna" allows GuitarPlot to map reads onto 5'UTR,CSD,and 3'UTR
                                        #pltTxType="ncrna" allows GuitarPlot to map reads onto ncRNA
                                        pltTxType = "mrna",
                                        #Change to TRUE if plotting "ncrna"
                                        headOrtail = TRUE,
                                        enableCI = FALSE) 


pdf("/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/coverage_native_IVT/coverage_native_IVT_nucleo.pdf")
metagene_plot_of_coverage + 
  scale_fill_manual(values=rep(c("transparent"),length(bam_labels))) +
  theme_void()+
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())
#then, you can add title, specify color using just ggplot2 functions
dev.off()


# Specify the file paths of the bam files
bam_files <- c("/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/coverage_native_IVT/cytoplasmic_PASS_rep1_rep2_4sU_spliced_mapped_to_Homo_sapiens.GRCh38.dna.primary_assembly_filtered.bam",
               "/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/coverage_native_IVT/IVT_gridION_PASS_all_reads_filtered.bam")

# Specify the file labels
bam_labels <- c("Cytoplasm SUM159", "IVT SUM159")

# Create the list of GRanges objects
granges_list <- create_granges_list(bam_files)

metagene_plot_of_coverage <- GuitarPlot(stGRangeLists = granges_list,
                                        stGroupName = bam_labels,
                                        txTxdb = txdb_from_gtf,
                                        #pltTxType="mrna" allows GuitarPlot to map reads onto 5'UTR,CSD,and 3'UTR
                                        #pltTxType="ncrna" allows GuitarPlot to map reads onto ncRNA
                                        pltTxType = "mrna",
                                        #Change to TRUE if plotting "ncrna"
                                        headOrtail = TRUE,
                                        enableCI = FALSE) 


pdf("/Users/paolamarango/Desktop/fractions_analysis_Paola_SUM159/coverage_native_IVT/coverage_native_IVT_cyto.pdf")
metagene_plot_of_coverage + 
  scale_fill_manual(values=rep(c("transparent"),length(bam_labels))) +
  theme_void()+
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())
#then, you can add title, specify color using just ggplot2 functions
dev.off()
