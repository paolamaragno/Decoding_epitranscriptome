library('Rsamtools')
library('GenomicRanges')
library('GenomicFeatures')
library('GenomicAlignments')
library('Guitar')
library('ggplot2')

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
gtf_file<-"/path/to/Homo_sapiens.GRCh38.104.gtf"
txdb_from_gtf <- makeTxDbFromGFF(file=gtf_file)

# Specify the file paths of the bam files
bam_files <- c("/path/to/native_filtered.bam",
               "/path/to/IVT_filtered.bam")

# Specify the file labels
bam_labels <- c("Native SUM159", "IVT K562 - SUM159")

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


pdf("/path/to/coverage_native_IVT.pdf")
metagene_plot_of_coverage + 
  scale_fill_manual(values=rep(c("transparent"),length(bam_labels))) +
  theme_void() +
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) 
#then, you can add title, specify color using just ggplot2 functions
dev.off()
