#!/usr/bin/env Rscript
### load input variables ###
args = commandArgs(trailingOnly=TRUE)

for(v in args)
{
	vTmp <- strsplit(v,"=")[[1]]
	assign(vTmp[[1]],vTmp[[2]])
}

### environment set-up ###
source(paste0(scriptsFolder,"nanoID.environmentSetUp.R"))

# Bam names
#
load(file=paste0(resultsFolder,"nanoid/bam.RData"))
bamNames <- names(bam)
rm(bam)
gc()

### create sequencing summary matrices ###
sequencing.files = list.files(dirname(FAST5paths),pattern = "sequencing_summary",full.names = TRUE,recursive = TRUE)

sequencing.summary <- read.delim(sequencing.files[[1]],sep="\t",row.names = NULL,header = TRUE,stringsAsFactors = FALSE)

if(length(sequencing.files>1))
{
	for(sequencing.file in sequencing.files[-1])
	{
		sequencing.summary <- rbind(sequencing.summary,read.delim(sequencing.file,sep="\t",row.names = NULL,header = TRUE,stringsAsFactors = FALSE))
	}
}

rownames(sequencing.summary) = sequencing.summary[,"read_id"]


sequencing.summary <- sequencing.summary[bamNames,]

save(sequencing.summary,file=paste0(resultsFolder,"nanoid/sequencing.summary.RData"))