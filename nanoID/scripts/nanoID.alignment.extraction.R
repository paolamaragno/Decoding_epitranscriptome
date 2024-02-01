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

### extract aligned sequences ###
reference.genome = readDNAStringSet(filepath = genomeFasta)
names(reference.genome) <- sapply(strsplit(names(reference.genome)," "),"[[",1)

bam = readGAlignments(bamPath,param = NULL,use.names = TRUE)
bam = bam[isUnique(names(bam))]

aligned.sequence.list = as.list(reference.genome[as(bam,"GRanges")])
aligned.sequence.list = lapply(aligned.sequence.list,as.character)
names(aligned.sequence.list) = names(bam)

cigar.strings = cigar(bam)
cigar.ops = explodeCigarOps(cigar.strings)
cigar.list = explodeCigarOpLengths(cigar.strings)
cigar.list = lapply(1:length(cigar.strings),function(x){vec = cigar.list[[x]];names(vec) = cigar.ops[[x]];vec})
names(cigar.list) = names(bam)

save(bam,file=paste0(resultsFolder,"nanoid/bam.RData"))
save(aligned.sequence.list,file=paste0(resultsFolder,"nanoid/aligned.sequence.list.RData"))
save(cigar.list,file=paste0(resultsFolder,"nanoid/cigar.list.RData"))

saveList(aligned.sequence.list,paste0(resultsFolder,"nanoid/aligned.sequence.list.llo"))
