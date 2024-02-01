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

load(paste0(nanoIDSupplementals,"five.mers.minus.RData"))
load(paste0(nanoIDSupplementals,"cigar.five.mers.RData"))
load(paste0(nanoIDSupplementals,"cigar.five.mers.S.containing.RData"))
load(paste0(nanoIDSupplementals,"cigar.five.mers.H.S.N.I.centered.RData"))

### five.mer.alignment.reconstruction.function ###
five.mer.alignment.reconstruction.function = function(read.name){
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
  
  if (which.strand == "-"){
    complement.map = c("A" = "T","C" = "G","G" = "C","T" = "A","N" = "N","-" = "-")
    
    alignment["aligned.read",] = complement.map[alignment["aligned.read",]]
    alignment["aligned.reference",] = complement.map[alignment["aligned.reference",]]
  }
  
  aligned.read.sequence.five.mers = as.character(five.mers.minus[as.numeric(stats::filter(as.numeric(c("A" = 1,"C" = 2,"G" = 3,"T" = 4,"-" = 5)[alignment["aligned.read",which(!colnames(alignment) %in% c("H","S"))]]),5^(0:8)[1:5],sides = 1)-(sum(5^(0:8)[1:5])-1))[-c(1:4)]])
  aligned.reference.sequence.five.mers = as.character(five.mers.minus[as.numeric(stats::filter(as.numeric(c("A" = 1,"C" = 2,"G" = 3,"T" = 4,"-" = 5)[alignment["aligned.reference",which(!colnames(alignment) %in% c("H","S"))]]),5^(0:8)[1:5],sides = 1)-(sum(5^(0:8)[1:5])-1))[-c(1:4)]])
  
  exploded.cigar.sequence.five.mers = as.character(cigar.five.mers[as.numeric(stats::filter(as.numeric(c("M" = 1,"I" = 2,"D" = 3,"H" = 4,"S" = 5,"N" = 6)[exploded.cigar]),6^(0:8)[1:5],sides = 1)-(sum(6^(0:8)[1:5])-1))[-c(1:4)]])
  
  five.mers.alignment = rbind(exploded.cigar.sequence.five.mers,
                              "aligned.read.sequence.five.mers" = NA,
                              "aligned.reference.sequence.five.mers" = NA)
  five.mers.alignment["aligned.read.sequence.five.mers",!(exploded.cigar.sequence.five.mers %in% c(cigar.five.mers.S.containing))] = aligned.read.sequence.five.mers
  five.mers.alignment["aligned.reference.sequence.five.mers",!(exploded.cigar.sequence.five.mers %in% c(cigar.five.mers.S.containing))] = aligned.reference.sequence.five.mers
  five.mers.alignment = five.mers.alignment[,!(exploded.cigar.sequence.five.mers %in% cigar.five.mers.H.S.N.I.centered)]
  
  return(five.mers.alignment)
}

### five.mer.alignment.reconstruction ###
bam = get(load(paste0(resultsFolder,"nanoid/bam.RData")))
cigar.list = get(load(paste0(resultsFolder,"nanoid/cigar.list.RData")))
read.to.fast5.name = get(load(paste0(resultsFolder,"nanoid/read.to.fast5.name.RData")))

bamNames <- names(bam)

index.list = split(1:length(bamNames), ceiling(1:length(bamNames)/100))

if(file.exists(where=paste0(resultsFolder,"nanoid/five.mer.alignment.reconstruction.list.llo")))
{
   start = getListLength(paste0(resultsFolder,"nanoid/five.mer.alignment.reconstruction.list.llo"))/100 + 1
   print(paste0("First index: ",start))
}else{
  five.mer.alignment.reconstruction.list = list()
  saveList(object = five.mer.alignment.reconstruction.list, file = paste0(resultsFolder,"nanoid/five.mer.alignment.reconstruction.list.llo"), append = FALSE, compress = FALSE)
  start = 1
}

registerDoParallel(cores = mc.cores)
for (j in start:length(index.list)){

  if(j%%50==0){print(j);system("grep MemFree /proc/meminfo")}

  aligned.sequence.list <- readList(paste0(resultsFolder,"nanoid/aligned.sequence.list.llo"),index.list[[j]])
  names(aligned.sequence.list) <- bamNames[index.list[[j]]]

  read.sequence.list <- readList(paste0(resultsFolder,"nanoid/read.sequence.list.llo"),index.list[[j]])
  names(read.sequence.list) <- bamNames[index.list[[j]]]

  five.mer.alignment.reconstruction.list.tmp = foreach(n = bamNames[index.list[[j]]],.noexport = setdiff(ls(),c("bam","aligned.sequence.list","cigar.list","read.sequence.list","read.to.fast5.name"))) %dopar% five.mer.alignment.reconstruction.function(n)
	
  saveList(object = five.mer.alignment.reconstruction.list.tmp, file = paste0(resultsFolder,"nanoid/five.mer.alignment.reconstruction.list.llo"), append = TRUE, compress = FALSE)

  #@ rm(five.mer.alignment.reconstruction.list.tmp)
  #@ rm(aligned.sequence.list)
  #@ rm(read.sequence.list)
  #@ gc()
}

five.mer.alignment.reconstruction.list = readList(paste0(resultsFolder,"nanoid/five.mer.alignment.reconstruction.list.llo"))
names(five.mer.alignment.reconstruction.list) <- bamNames

save(five.mer.alignment.reconstruction.list,file=file.path(paste0(resultsFolder,"nanoid/five.mer.alignment.reconstruction.RData")))