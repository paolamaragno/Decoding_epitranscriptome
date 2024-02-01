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

### alignment.reconstruction.function ###
alignment.reconstruction.function = function(read.name){
  fast5.name = read.to.fast5.name[read.name]
  read = read.sequence.list[[read.name]]
  which.strand = as.character(strand(bam[read.name]))
  if (which.strand == "-"){read = as.character(reverseComplement(as(read,'DNAString')))}
  reference = aligned.sequence.list[[read.name]]
  cigar = cigar.list[[read.name]]
  exploded.cigar = as.vector(Rle(names(cigar),cigar))
  alignment.length = length(exploded.cigar)
  aligned.read = rep("-",alignment.length)
  aligned.read[which(exploded.cigar %in% c("M","I","S","H"))] = strsplit(read,split = "")[[1]]
  aligned.reference = rep("-",alignment.length)
  aligned.reference[which(exploded.cigar %in% c("M","D","N"))] = strsplit(reference,split = "")[[1]]
  alignment = rbind(aligned.read,aligned.reference)
  colnames(alignment) = exploded.cigar
  alignment.uncut = alignment
  alignment = alignment[,which(colnames(alignment) != "H"),drop = FALSE] # cigar H is for hard-clipping
  alignment = alignment[,which(colnames(alignment) != "S"),drop = FALSE] # cigar S is for soft-clipping
  alignment = alignment[,which(colnames(alignment) != "N"),drop = FALSE] # cigar N is for splicing
  alignment = alignment[,which(alignment["aligned.reference",] != "-"),drop = FALSE] # "-" is for insertions in the reference
  
  if (which.strand == "-"){
    complement.map = c("A" = "T","C" = "G","G" = "C","T" = "A","N" = "N","-" = "-")
    
    alignment["aligned.read",] = complement.map[alignment["aligned.read",]]
    alignment["aligned.reference",] = complement.map[alignment["aligned.reference",]]
  }
  
  return(alignment)
}

### alignment.reconstruction ###
bam = get(load(paste0(resultsFolder,"nanoid/bam.RData")))
cigar.list = get(load(paste0(resultsFolder,"nanoid/cigar.list.RData")))
read.to.fast5.name = get(load(paste0(resultsFolder,"nanoid/read.to.fast5.name.RData")))

bamNames <- names(bam)
index.list = split(1:length(bamNames), ceiling(1:length(bamNames)/100))

if(file.exists(where=paste0(resultsFolder,"nanoid/alignment.reconstruction.list.llo")))
{
   start = (getListLength(paste0(resultsFolder,"nanoid/alignment.reconstruction.list.llo"))/100 + 1)
   print(paste0("First index: ",start))
}else{
  alignment.reconstruction = list()
  saveList(object = alignment.reconstruction, file = paste0(resultsFolder,"nanoid/alignment.reconstruction.list.llo"), append = FALSE, compress = FALSE)
  start = 1
}

registerDoParallel(cores = mc.cores)

for (j in start:length(index.list)){
  if(j%%50==0){print(j);system("grep MemFree /proc/meminfo")}

  aligned.sequence.list <- readList(paste0(resultsFolder,"nanoid/aligned.sequence.list.llo"),index.list[[j]])
  names(aligned.sequence.list) <- bamNames[index.list[[j]]]

  read.sequence.list <- readList(paste0(resultsFolder,"nanoid/read.sequence.list.llo"),index.list[[j]])
  names(read.sequence.list) <- bamNames[index.list[[j]]]

  alignment.reconstruction.tmp = foreach(n = bamNames[index.list[[j]]],.noexport = setdiff(ls(),c("bam","read.to.fast5.name","read.sequence.list","sequencing.summary"))) %dopar% alignment.reconstruction.function(n)
  saveList(object = alignment.reconstruction.tmp, file = paste0(resultsFolder,"nanoid/alignment.reconstruction.list.llo"), append = TRUE, compress = FALSE)

  #@ rm(alignment.reconstruction.tmp)
  #@ rm(aligned.sequence.list)
  #@ rm(read.sequence.list)
  #@ gc()
}
alignment.reconstruction <- readList(file = paste0(resultsFolder,"nanoid/alignment.reconstruction.list.llo"))
names(alignment.reconstruction) = bamNames

save(alignment.reconstruction,file=paste0(resultsFolder,"nanoid/alignment.reconstruction.RData"))
