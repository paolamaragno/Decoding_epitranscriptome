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

### load objects ###
load(paste0(resultsFolder,"nanoid/fast5.files.RData"))

load(file=paste0(resultsFolder,"nanoid/bam.RData"))
bamNames <- names(bam)
rm(bam)
gc()

### for any read we associate the fast5 file and the nucleotide sequence predicted by guppy ###

read.to.fast5.name = c()
read.sequence.list = list()
saveList(object = read.sequence.list, file = paste0(resultsFolder,"nanoid/read.sequence.list.llo"), append = FALSE, compress = FALSE)

for(fast5.file in fast5.files){
	## extraction of fast5 reads names
	fast5ls <- h5ls(fast5.file)

	### For each fast5 I take the most recent Basecall_1D slot
	slotName <- sort(unique(fast5ls$name[grep("Basecall_1D",fast5ls$name)]),decreasing=TRUE)[[1]]

	fast5reads <- unique(fast5ls$name[grep("read_",fast5ls$name)])

	fast5reads <- fast5reads[fast5reads %in% paste0("read_",bamNames)]

	if(length(fast5reads)==0){print("No reads of interest!");next}
	
	print(paste0(length(fast5reads)," reads of interest in fast5: ",fast5.file))
	system("grep MemFree /proc/meminfo")

	## instead of a for we use a mclapply to parallelize
	tmp <- unlist(mclapply(fast5reads, function(fast5read)
	{
		h5 <- h5read(fast5.file,paste0(fast5read,"/Analyses/",slotName,"/BaseCalled_template/Fastq"))

		fast5read <- strsplit(fast5read,"_")[[1]][[2]]

		out <- gsub("U","T",strsplit(h5,split = "\n")[[1]][2])
		names(out) <- fast5.file

		out
	},mc.cores=mc.cores))

 	namesTmp <- names(tmp)
 	names(namesTmp) <- sapply(strsplit(fast5reads,"_"),function(i)i[[2]])
 	names(tmp) <- sapply(strsplit(fast5reads,"_"),function(i)i[[2]])

  	read.to.fast5.name = c(read.to.fast5.name, namesTmp)
  	read.sequence.list = c(read.sequence.list, tmp)

  	saveList(object = as.list(tmp), file = paste0(resultsFolder,"nanoid/read.sequence.list.llo"), append = TRUE, compress = FALSE)

  	#@ rm(fast5ls)
  	#@ rm(fast5reads)
  	#@ rm(namesTmp)
  	#@ rm(tmp)
  	#@ gc()
}

### results save
# save(read.to.fast5.name,file=paste0(resultsFolder,"nanoid/read.to.fast5.name.RData"))

# read.sequence.list <- readList(file = paste0(resultsFolder,"nanoid/read.sequence.list.llo"))
names(read.sequence.list) <- names(read.to.fast5.name)

# MF I decide to keep only the aligned reads
#
read.to.fast5.name <- read.to.fast5.name[bamNames]
read.sequence.list <- read.sequence.list[bamNames]

save(read.to.fast5.name,file=paste0(resultsFolder,"nanoid/read.to.fast5.name.RData"))
save(read.sequence.list,file=paste0(resultsFolder,"nanoid/read.sequence.list.RData"))
saveList(read.sequence.list,paste0(resultsFolder,"nanoid/read.sequence.list.llo"))
