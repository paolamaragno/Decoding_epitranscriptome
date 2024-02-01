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

### custom function ###
rna.to.dna = function(h5,suffix = ""){
	  fastq = strsplit(h5,split = "\n")
  read.name = strsplit(fastq[[1]][1],split = " ")
    read.name[[1]][1] = paste0(read.name[[1]][1],suffix)
    fastq[[1]][1] = paste(c(unlist(read.name),""),collapse = " ")
      fastq[[1]][2] = gsub("U","T",fastq[[1]][2])
      paste(c(unlist(fastq),""),collapse = "\n")
}

### create fast5.files vector ###
fast5.files = list.files(FAST5paths,pattern = ".fast5$",full.names = TRUE,recursive = TRUE)
save(fast5.files,file=paste0(resultsFolder,"nanoid/fast5.files.RData"))

## for any multi-reads fast5 we create a fastq

registerDoParallel(cores = mc.cores)

foreach(fast5.file=fast5.files) %dopar% {
	print(fast5.file)

	nameTmp <- gsub("\\.fast5","",basename(fast5.file))
	
	# if(file.exists(paste0(resultsFolder,"FASTQ_DNA/",nameTmp,".fastq")))
	# {
	# 	print("Exists")
	# }else{
		## extraction of fast5 reads names
		fast5ls <- h5ls(fast5.file)

		### For each fast5 I take the most recent Basecall_1D slot
		slotName <- sort(unique(fast5ls$name[grep("Basecall_1D",fast5ls$name)]),decreasing=TRUE)[[1]]

		fast5reads <- unique(fast5ls$name[grep("read_",fast5ls$name,ignore.case=TRUE)])

		sink(paste0(resultsFolder,"FASTQ_DNA/",nameTmp,".fastq"))
			
		for(fast5read in fast5reads)
		{
			cat(rna.to.dna(h5read(fast5.file,paste0(fast5read,"/Analyses/",slotName,"/BaseCalled_template/Fastq")),suffix = ""))
		}
		sink()
	# }
}
