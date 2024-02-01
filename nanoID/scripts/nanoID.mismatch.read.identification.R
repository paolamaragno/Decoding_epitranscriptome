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
objects = list.files(nanoIDSupplementals)
for(object in objects){load(paste0(nanoIDSupplementals,object))}

### build.read.based.mismatch.identification ###
build.read.based.mismatch.identification = function(read.name){
  alignment = alignment.reconstruction[[read.name]]
  
  read.T.positions = which(alignment["aligned.reference",] == "T")
  read.non.T.positions = which(alignment["aligned.reference",] %in% c("A","C","G"))
  read.T.containing = intersect(1:ncol(alignment),unique(c(read.T.positions-2,read.T.positions-1,read.T.positions,read.T.positions+1,read.T.positions+2)))
  read.non.T.containing = setdiff(1:ncol(alignment),unique(c(read.T.positions-2,read.T.positions-1,read.T.positions,read.T.positions+1,read.T.positions+2)))
  
  read.A.positions = which(alignment["aligned.reference",] == "A")
  read.C.positions = which(alignment["aligned.reference",] == "C")
  read.G.positions = which(alignment["aligned.reference",] == "G")
  
  alignment.T.positions = alignment[,read.T.positions,drop = FALSE]
  alignment.non.T.positions = alignment[,read.non.T.positions,drop = FALSE]
  alignment.T.containing = alignment[,read.T.containing,drop = FALSE]
  alignment.non.T.containing = alignment[,read.non.T.containing,drop = FALSE]
  
  alignment.A.positions = alignment[,read.A.positions,drop = FALSE]
  alignment.C.positions = alignment[,read.C.positions,drop = FALSE]
  alignment.G.positions = alignment[,read.G.positions,drop = FALSE]
  
  aligned.read.sequence.five.mers = as.character(five.mers.minus[as.numeric(stats::filter(as.numeric(c("A" = 1,"C" = 2,"G" = 3,"T" = 4,"-" = 5)[alignment["aligned.read",which(!colnames(alignment) %in% c("H","S"))]]),5^(0:8)[1:5],sides = 1)-(sum(5^(0:8)[1:5])-1))[-c(1:4)]])
  aligned.reference.sequence.five.mers = as.character(five.mers.minus[as.numeric(stats::filter(as.numeric(c("A" = 1,"C" = 2,"G" = 3,"T" = 4,"-" = 5)[alignment["aligned.reference",which(!colnames(alignment) %in% c("H","S"))]]),5^(0:8)[1:5],sides = 1)-(sum(5^(0:8)[1:5])-1))[-c(1:4)]])
  
  tryCatch({
    five.mers.alignment = five.mer.alignment.reconstruction[[read.name]]
    five.mers.alignment = five.mers.alignment[,which(five.mers.alignment["aligned.read.sequence.five.mers",] != "NA" & five.mers.alignment["aligned.reference.sequence.five.mers",] != "NA")]
    five.mers.alignment = five.mers.alignment[,which(five.mers.alignment["exploded.cigar.sequence.five.mers",] == "MMMMM")]
    
    five.mers.alignment.mismatches = as.numeric(five.mers.alignment["aligned.read.sequence.five.mers",] != five.mers.alignment["aligned.reference.sequence.five.mers",])
    names(five.mers.alignment.mismatches) = five.mers.alignment["aligned.read.sequence.five.mers",]
  },error = function(x){five.mers.alignment.mismatches = 0})
  
  results = c(
    "Length" = length(alignment["aligned.read",]),
    "A" = length(alignment["aligned.read",which(alignment["aligned.reference",] == "A")]),
    "C" = length(alignment["aligned.read",which(alignment["aligned.reference",] == "C")]),
    "G" = length(alignment["aligned.read",which(alignment["aligned.reference",] == "G")]),
    "T" = length(alignment.T.positions["aligned.read",]),
    "TtoA" = sum(alignment.T.positions["aligned.read",] == "A"),
    "TtoC" = sum(alignment.T.positions["aligned.read",] == "C"),
    "TtoG" = sum(alignment.T.positions["aligned.read",] == "G"),
    "Tto-" = sum(alignment.T.positions["aligned.read",] == "-"),
    "TtoACG" = sum(alignment.T.positions["aligned.read",] %in% c("A","C","G")),
    "TtoACG-" = sum(alignment.T.positions["aligned.read",] %in% c("A","C","G","-")),
    "binom.test.TtoACG" = tryCatch({binom.test(sum(alignment.T.positions["aligned.read",] %in% c("A","C","G")),length(alignment.T.positions["aligned.read",]),p = mismatch.probability.TtoACG,alternative = "greater")$p.value},error = function(x){1}),
    "prop.test.TtoACG" = tryCatch({prop.test(sum(alignment.T.positions["aligned.read",] %in% c("A","C","G")),length(alignment.T.positions["aligned.read",]),p = mismatch.probability.TtoACG,alternative = "greater")$p.value},error = function(x){1}),
    "estimate.TtoACG" = sum(alignment.T.positions["aligned.read",] %in% c("A","C","G"))/length(alignment.T.positions["aligned.read",]),
    
    "fisher.test.TtoACG" = fisher.test(matrix(c(sum(alignment.T.positions["aligned.read",] %in% c("A","C","G")),
                                                sum(alignment.T.positions["aligned.read",] %in% c("T")),
                                                sum(alignment.non.T.positions["aligned.reference",which(alignment.non.T.positions["aligned.read",] != "-")] != alignment.non.T.positions["aligned.read",which(alignment.non.T.positions["aligned.read",] != "-")]),
                                                sum(alignment.non.T.positions["aligned.reference",which(alignment.non.T.positions["aligned.read",] != "-")] == alignment.non.T.positions["aligned.read",which(alignment.non.T.positions["aligned.read",] != "-")])),nrow = 2,dimnames = list(c("D","notD"),c("A","notA"))),alternative = "greater")$p.value,
    
    "binom.test.TtoACG-" = tryCatch({binom.test(sum(alignment.T.positions["aligned.read",] %in% c("A","C","G","-")),length(alignment.T.positions["aligned.read",]),p = mismatch.probability.TtoACG.minus,alternative = "greater")$p.value},error = function(x){1}),
    "prop.test.TtoACG-" = tryCatch({prop.test(sum(alignment.T.positions["aligned.read",] %in% c("A","C","G","-")),length(alignment.T.positions["aligned.read",]),p = mismatch.probability.TtoACG.minus,alternative = "greater")$p.value},error = function(x){1}),
    "estimate.TtoACG-" = sum(alignment.T.positions["aligned.read",] %in% c("A","C","G","-"))/length(alignment.T.positions["aligned.read",]),
    "estimate.Tto-" = sum(alignment.T.positions["aligned.read",] %in% c("-"))/length(alignment.T.positions["aligned.read",]),
    
    "fisher.test.TtoACG-" = fisher.test(matrix(c(sum(alignment.T.positions["aligned.read",] %in% c("A","C","G","-")),
                                                 sum(alignment.T.positions["aligned.read",] %in% c("T")),
                                                 sum(alignment.non.T.positions["aligned.reference",] != alignment.non.T.positions["aligned.read",]),
                                                 sum(alignment.non.T.positions["aligned.reference",] == alignment.non.T.positions["aligned.read",])),nrow = 2,dimnames = list(c("D","notD"),c("A","notA"))),alternative = "greater")$p.value,
    
    "AtoC" = sum(alignment.A.positions["aligned.read",] == "C"),
    "AtoG" = sum(alignment.A.positions["aligned.read",] == "G"),
    "AtoT" = sum(alignment.A.positions["aligned.read",] == "T"),
    "Ato-" = sum(alignment.A.positions["aligned.read",] == "-"),
    "AtoTCG" = sum(alignment["aligned.read",which(alignment["aligned.reference",] == "A")] %in% c("T","C","G")),
    "AtoTCG-" = sum(alignment["aligned.read",which(alignment["aligned.reference",] == "A")] %in% c("T","C","G","-")),
    "CtoA" = sum(alignment.C.positions["aligned.read",] == "A"),
    "CtoG" = sum(alignment.C.positions["aligned.read",] == "G"),
    "CtoT" = sum(alignment.C.positions["aligned.read",] == "T"),
    "Cto-" = sum(alignment.C.positions["aligned.read",] == "-"),
    "CtoATG" = sum(alignment["aligned.read",which(alignment["aligned.reference",] == "C")] %in% c("A","T","G")),
    "CtoATG-" = sum(alignment["aligned.read",which(alignment["aligned.reference",] == "C")] %in% c("A","T","G","-")),
    "GtoA" = sum(alignment.G.positions["aligned.read",] == "A"),
    "GtoC" = sum(alignment.G.positions["aligned.read",] == "C"),
    "GtoT" = sum(alignment.G.positions["aligned.read",] == "T"),
    "Gto-" = sum(alignment.G.positions["aligned.read",] == "-"),
    "GtoATC" = sum(alignment["aligned.read",which(alignment["aligned.reference",] == "G")] %in% c("A","T","C")),
    "GtoATC-" = sum(alignment["aligned.read",which(alignment["aligned.reference",] == "G")] %in% c("A","T","C","-")),
    "estimate.AtoTCG" = sum(alignment["aligned.read",which(alignment["aligned.reference",] == "A")] %in% c("T","C","G"))/length(alignment["aligned.read",which(alignment["aligned.reference",] == "A")]),
    "estimate.AtoTCG-" = sum(alignment["aligned.read",which(alignment["aligned.reference",] == "A")] %in% c("T","C","G","-"))/length(alignment["aligned.read",which(alignment["aligned.reference",] == "A")]),
    "estimate.CtoATG" = sum(alignment["aligned.read",which(alignment["aligned.reference",] == "C")] %in% c("A","T","G"))/length(alignment["aligned.read",which(alignment["aligned.reference",] == "C")]),
    "estimate.CtoATG-" = sum(alignment["aligned.read",which(alignment["aligned.reference",] == "C")] %in% c("A","T","G","-"))/length(alignment["aligned.read",which(alignment["aligned.reference",] == "C")]),
    "estimate.GtoATC" = sum(alignment["aligned.read",which(alignment["aligned.reference",] == "G")] %in% c("A","T","C"))/length(alignment["aligned.read",which(alignment["aligned.reference",] == "G")]),
    "estimate.GtoATC-" = sum(alignment["aligned.read",which(alignment["aligned.reference",] == "G")] %in% c("A","T","C","-"))/length(alignment["aligned.read",which(alignment["aligned.reference",] == "G")]),
    "estimate.mismatches" = sum(alignment["aligned.reference",which(alignment["aligned.reference",] != "-" & alignment["aligned.read",] != "-")] != alignment["aligned.read",which(alignment["aligned.reference",] != "-" & alignment["aligned.read",] != "-")])/length(alignment["aligned.read",]),
    "estimate.mismatches.-" = sum(alignment["aligned.reference",which(alignment["aligned.reference",] != "-")] != alignment["aligned.read",which(alignment["aligned.reference",] != "-")])/length(alignment["aligned.read",]),
    "estimate.mismatches.T.containing" = sum(alignment.T.containing["aligned.reference",which(alignment.T.containing["aligned.reference",] != "-" & alignment.T.containing["aligned.read",] != "-")] != alignment.T.containing["aligned.read",which(alignment.T.containing["aligned.reference",] != "-" & alignment.T.containing["aligned.read",] != "-")])/length(alignment.T.containing["aligned.read",]),
    "estimate.mismatches.T.containing.-" = sum(alignment.T.containing["aligned.reference",which(alignment.T.containing["aligned.reference",] != "-")] != alignment.T.containing["aligned.read",which(alignment.T.containing["aligned.reference",] != "-")])/length(alignment.T.containing["aligned.read",]),
    "estimate.mismatches.non.T.containing" = sum(alignment.non.T.containing["aligned.reference",which(alignment.non.T.containing["aligned.reference",] != "-" & alignment.non.T.containing["aligned.read",] != "-")] != alignment.non.T.containing["aligned.read",which(alignment.non.T.containing["aligned.reference",] != "-" & alignment.non.T.containing["aligned.read",] != "-")])/length(alignment.non.T.containing["aligned.read",]),
    "estimate.mismatches.non.T.containing.-" = sum(alignment.non.T.containing["aligned.reference",which(alignment.non.T.containing["aligned.reference",] != "-")] != alignment.non.T.containing["aligned.read",which(alignment.non.T.containing["aligned.reference",] != "-")])/length(alignment.non.T.containing["aligned.read",]),
    
    "binom.test.T.containing" = tryCatch({binom.test(sum(alignment.T.containing["aligned.reference",which(alignment.T.containing["aligned.reference",] != "-" & alignment.T.containing["aligned.read",] != "-")] != alignment.T.containing["aligned.read",which(alignment.T.containing["aligned.reference",] != "-" & alignment.T.containing["aligned.read",] != "-")]),length(alignment.T.containing["aligned.read",]),p = mismatch.probability.T.containing,alternative = "greater")$p.value},error = function(x){1}),
    "prop.test.T.containing" = tryCatch({prop.test(sum(alignment.T.containing["aligned.reference",which(alignment.T.containing["aligned.reference",] != "-" & alignment.T.containing["aligned.read",] != "-")] != alignment.T.containing["aligned.read",which(alignment.T.containing["aligned.reference",] != "-" & alignment.T.containing["aligned.read",] != "-")]),length(alignment.T.containing["aligned.read",]),p = mismatch.probability.T.containing,alternative = "greater")$p.value},error = function(x){1}),
    
    "fisher.test.T.containing" = fisher.test(matrix(c(sum(alignment.T.containing["aligned.reference",which(alignment.T.containing["aligned.read",] != "-")] != alignment.T.containing["aligned.read",which(alignment.T.containing["aligned.read",] != "-")]),
                                                      sum(alignment.T.containing["aligned.reference",which(alignment.T.containing["aligned.read",] != "-")] == alignment.T.containing["aligned.read",which(alignment.T.containing["aligned.read",] != "-")]),
                                                      sum(alignment.non.T.containing["aligned.reference",which(alignment.non.T.containing["aligned.read",] != "-")] != alignment.non.T.containing["aligned.read",which(alignment.non.T.containing["aligned.read",] != "-")]),
                                                      sum(alignment.non.T.containing["aligned.reference",which(alignment.non.T.containing["aligned.read",] != "-")] == alignment.non.T.containing["aligned.read",which(alignment.non.T.containing["aligned.read",] != "-")])),nrow = 2,dimnames = list(c("D","notD"),c("A","notA"))),alternative = "greater")$p.value,
    
    "binom.test.T.containing.-" = tryCatch({binom.test(sum(alignment.T.containing["aligned.reference",which(alignment.T.containing["aligned.reference",] != "-")] != alignment.T.containing["aligned.read",which(alignment.T.containing["aligned.reference",] != "-")]),length(alignment.T.containing["aligned.read",]),p = mismatch.probability.T.containing.minus,alternative = "greater")$p.value},error = function(x){1}),
    "prop.test.T.containing.-" = tryCatch({prop.test(sum(alignment.T.containing["aligned.reference",which(alignment.T.containing["aligned.reference",] != "-")] != alignment.T.containing["aligned.read",which(alignment.T.containing["aligned.reference",] != "-")]),length(alignment.T.containing["aligned.read",]),p = mismatch.probability.T.containing.minus,alternative = "greater")$p.value},error = function(x){1}),
    
    "fisher.test.T.containing.-" = fisher.test(matrix(c(sum(alignment.T.containing["aligned.reference",] != alignment.T.containing["aligned.read",]),
                                                        sum(alignment.T.containing["aligned.reference",] == alignment.T.containing["aligned.read",]),
                                                        sum(alignment.non.T.containing["aligned.reference",] != alignment.non.T.containing["aligned.read",]),
                                                        sum(alignment.non.T.containing["aligned.reference",] == alignment.non.T.containing["aligned.read",])),nrow = 2,dimnames = list(c("D","notD"),c("A","notA"))),alternative = "greater")$p.value)
  
  three.mers.T.containing = sapply(five.mers.three.mers.centered.T.containing.list,function(x){sum(five.mers.alignment.mismatches[x],na.rm = TRUE)/sum(names(five.mers.alignment.mismatches) %in% x)})
  names(three.mers.T.containing) = paste0(names(five.mers.three.mers.centered.T.containing.list),".T.containing")
  
  three.mers.non.T.containing = sapply(five.mers.three.mers.centered.non.T.containing.list,function(x){sum(five.mers.alignment.mismatches[x],na.rm = TRUE)/sum(names(five.mers.alignment.mismatches) %in% x)})
  names(three.mers.non.T.containing) = paste0(names(five.mers.three.mers.centered.non.T.containing.list),".non.T.containing")
  
  results = c(results,three.mers.T.containing,three.mers.non.T.containing)

  return(results)
}

### read based mismatch identification ###
bam = get(load(paste0(resultsFolder,"nanoid/bam.RData")))
bamNames <- names(bam)
rm(bam)
gc()

index.list = split(1:length(bamNames), ceiling(1:length(bamNames)/100))

# alignment.reconstruction = get(load(paste0(resultsFolder,"nanoid/alignment.reconstruction.RData")))
# five.mer.alignment.reconstruction = get(load(paste0(resultsFolder,"nanoid/five.mer.alignment.reconstruction.RData")))

start = 1

if(file.exists(where=paste0(resultsFolder,"nanoid/read.based.mismatch.identification.csv")))
{
  nrowMat <- system(paste0("cat ",paste0(resultsFolder,"nanoid/read.based.mismatch.identification.csv")," | wc -l"), intern=TRUE)
  nrowMat <- (as.numeric(nrowMat)-1)

  if(nrowMat%in%sapply(index.list,max))
  {
    start = (which(nrowMat==sapply(index.list,max))+1)
  }
}

if(start<=length(index.list))
{
  registerDoParallel(cores = mc.cores)

  for (j in start:length(index.list)){
    if(j%%50==0){print(j);system("grep MemFree /proc/meminfo")}

    alignment.reconstruction = readList(paste0(resultsFolder,"nanoid/alignment.reconstruction.list.llo"),index=index.list[[j]])
    names(alignment.reconstruction) <- bamNames[index.list[[j]]]

    five.mer.alignment.reconstruction = readList(paste0(resultsFolder,"nanoid/five.mer.alignment.reconstruction.list.llo"),index=index.list[[j]])
    names(five.mer.alignment.reconstruction) <- bamNames[index.list[[j]]]

    read.based.mismatch.identification.list.tmp = foreach(n = bamNames[index.list[[j]]],.noexport = setdiff(ls(),c("alignment.reconstruction"))) %dopar% build.read.based.mismatch.identification(n)

    read.based.mismatch.identification.list.tmp <- do.call("rbind", read.based.mismatch.identification.list.tmp)
    rownames(read.based.mismatch.identification.list.tmp) <- bamNames[index.list[[j]]]

    fwrite(x=read.based.mismatch.identification.list.tmp
          ,file=paste0(resultsFolder,"nanoid/read.based.mismatch.identification.csv")
          ,append=TRUE)

    #@ rm(read.based.mismatch.identification.list.tmp)
    #@ rm(alignment.reconstruction)
    #@ rm(five.mer.alignment.reconstruction)
    #@ gc()
  }
}

if(!file.exists(where=paste0(resultsFolder,"nanoid/rownames.csv")))
{
  fwrite(x=as.matrix(bamNames,ncol=1),file=paste0(resultsFolder,"nanoid/rownames.csv"),row.names=TRUE,col.names=TRUE,quote=FALSE)
  system(paste0("cut --complement -f1 -d ',' ",resultsFolder,"nanoid/rownames.csv > ",resultsFolder,"nanoid/rownamesB.csv"))
  system(paste0("mv ",resultsFolder,"nanoid/rownamesB.csv ",resultsFolder,"nanoid/rownames.csv"))
}
