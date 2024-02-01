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

### To remove the info related to the Move slot ### 
if(moveSlot=="true"){moveSlot="Move"}else{moveSlot="Foe"}

### load objects ###
objects = list.files(nanoIDSupplementals)
for(object in objects){load(paste0(nanoIDSupplementals,object))}

### build.raw.signal.list ###
build.raw.signal.list = function(read.name){
  if(moveSlot=="Move") # In this way, if the moveSlot is "foe", the production of the NaN matrix is faster.
  {
    try({
      fast5.file = read.to.fast5.name[read.name]
      
      raw.signal = h5read(fast5.file,paste0("read_",read.name))
      slotName = sort(names(raw.signal$Analyses)[grep("Basecall_1D",names(raw.signal$Analyses))],decreasing=TRUE)[[1]]
      raw.signal = raw.signal$Raw$Signal
      # raw.signal = h5read(fast5.file,paste0("read_",read.name,"/Raw/Signal"))

      move = h5read(fast5.file,paste0("read_",read.name,"/Analyses/",slotName,"/BaseCalled_template/",moveSlot))
      move.rle = Rle(move)
      
      read.sequence = rev(strsplit(read.sequence.list[[read.name]],split = "")[[1]])
      
      event.repeats = move[move == 1]
      event.repeats[cumsum(runLength(move.rle)[runValue(move.rle) == 1])] = runLength(move.rle)[runValue(move.rle) == 0] + 1
      if (move[length(move)] == 1){event.repeats[length(event.repeats)] = 1}
      
      stride = 10
      
      num_events_template = sequencing.summary[strsplit(read.name,split = "\\.")[[1]][1],"num_events_template"]
      num_events = sequencing.summary[strsplit(read.name,split = "\\.")[[1]][1],"num_events"]
      
      anchor = (num_events - num_events_template)*stride
      anchored.raw.signal = raw.signal[anchor:length(raw.signal)]
      
      read.sequence.five.mers = as.character(five.mers.minus[as.numeric(stats::filter(as.numeric(c("A" = 1,"C" = 2,"G" = 3,"T" = 4,"-" = 5)[read.sequence]),5^(0:8)[1:5],sides = 1)-(sum(5^(0:8)[1:5])-1))[-c(1:4)]])
      read.sequence.five.mers = c(NA,NA,read.sequence.five.mers,NA,NA)
      
      events.to.raw = cut(1:length(anchored.raw.signal),(0:num_events_template)*stride,right = FALSE,labels = FALSE)
      event.raw.signal = lapply(tapply(1:length(anchored.raw.signal),INDEX = events.to.raw,identity),function(x){anchored.raw.signal[x]})
      raw.means = sapply(event.raw.signal,mean)
      
      cumsum.event.repeats = cumsum(event.repeats)
      
      reevent = cbind(c(1,cumsum.event.repeats[-length(cumsum.event.repeats)] + 1),cumsum.event.repeats)
      colnames(reevent) = c("start","end")

      reevents.to.raw = cut(1:length(anchored.raw.signal),c(((0:num_events_template)*stride)[reevent[,"start"]],num_events_template*stride+1),right = FALSE,labels = FALSE)
      reevent.raw.signal = lapply(tapply(1:length(anchored.raw.signal),INDEX = reevents.to.raw,identity),function(x){anchored.raw.signal[x]})
      reevent.raw.means = sapply(reevent.raw.signal,mean)
      
      events = cbind("move" = move,
                     "read.sequence" = rep(read.sequence,times = event.repeats),
                     "read.sequence.five.mers" = rep(read.sequence.five.mers,times = event.repeats))
      
      means = cbind("raw.means" = raw.means,
                    "reevent.raw.means" = rep(reevent.raw.means,times = event.repeats))
      
      events.move = events[move == 1,]
      means.move = means[move == 1,]
      
      events.no.move = events[move == 0,]
      means.no.move = means[move == 0,]
      
      results = c(read.name,
                  apply(means.move[events.move[,"read.sequence.five.mers"] %in% five.mers.T.centered,"reevent.raw.means",drop = FALSE],2,mean),
                  apply(means.move[events.move[,"read.sequence.five.mers"] %in% five.mers.T.containing,"reevent.raw.means",drop = FALSE],2,mean),
                  apply(means.move[events.move[,"read.sequence.five.mers"] %in% five.mers.non.T.containing,"reevent.raw.means",drop = FALSE],2,mean),
                  
                  apply(means.move[events.move[,"read.sequence.five.mers"] %in% setdiff(five.mers.T.containing,five.mers.T.centered),"reevent.raw.means",drop = FALSE],2,mean),
                  
                  apply(means.move[events.move[,"read.sequence.five.mers"] %in% five.mers.T.leading,"reevent.raw.means",drop = FALSE],2,mean),
                  apply(means.move[events.move[,"read.sequence.five.mers"] %in% five.mers.T.lagging,"reevent.raw.means",drop = FALSE],2,mean),
                  
                  apply(means.move[events.move[,"read.sequence.five.mers"] %in% five.mers.T.centered.leading,"reevent.raw.means",drop = FALSE],2,mean),
                  apply(means.move[events.move[,"read.sequence.five.mers"] %in% five.mers.T.centered.lagging,"reevent.raw.means",drop = FALSE],2,mean),
                  
                  apply(means.move[events.move[,"read.sequence.five.mers"] %in% intersect(five.mers.A.centered,five.mers.T.containing),"reevent.raw.means",drop = FALSE],2,mean),
                  apply(means.move[events.move[,"read.sequence.five.mers"] %in% intersect(five.mers.C.centered,five.mers.T.containing),"reevent.raw.means",drop = FALSE],2,mean),
                  apply(means.move[events.move[,"read.sequence.five.mers"] %in% intersect(five.mers.G.centered,five.mers.T.containing),"reevent.raw.means",drop = FALSE],2,mean),
                  
                  apply(means.move[events.move[,"read.sequence.five.mers"] %in% intersect(five.mers.A.centered,five.mers.non.T.containing),"reevent.raw.means",drop = FALSE],2,mean),
                  apply(means.move[events.move[,"read.sequence.five.mers"] %in% intersect(five.mers.C.centered,five.mers.non.T.containing),"reevent.raw.means",drop = FALSE],2,mean),
                  apply(means.move[events.move[,"read.sequence.five.mers"] %in% intersect(five.mers.G.centered,five.mers.non.T.containing),"reevent.raw.means",drop = FALSE],2,mean),
                  
                  sapply(five.mers.three.mers.centered.list,function(x){apply(means.move[events.move[,"read.sequence.five.mers"] %in% x,"reevent.raw.means",drop = FALSE],2,mean)}),

                  sapply(five.mers.three.mers.centered.T.containing.list,function(x){apply(means.move[events.move[,"read.sequence.five.mers"] %in% x,"reevent.raw.means",drop = FALSE],2,mean)}),
                  sapply(five.mers.three.mers.centered.non.T.containing.list,function(x){apply(means.move[events.move[,"read.sequence.five.mers"] %in% x,"reevent.raw.means",drop = FALSE],2,mean)})
      )
      names(results) = NULL
      return(results)
    },silent = TRUE)
  }
  if(!exists("results")) # Once I got an error, event.repeats longer than reevent.raw.means...
  {
    c(read.name,rep("NaN",169))
  }
}

col.names = c("read.name",
              "five.mers.T.centered","five.mers.T.containing","five.mers.non.T.containing",
              "five.mers.non.T.centered",
              "five.mers.T.leading","five.mers.T.lagging","five.mers.T.centered.leading","five.mers.T.centered.lagging",
              "five.mers.A.centered.T.containing","five.mers.C.centered.T.containing","five.mers.G.centered.T.containing",
              "five.mers.A.centered.non.T.containing","five.mers.C.centered.non.T.containing","five.mers.G.centered.non.T.containing",
              mkAllStrings(c("A","C","G","T"),3),
              paste0(names(five.mers.three.mers.centered.T.containing.list),".T.containing"),
              paste0(names(five.mers.three.mers.centered.non.T.containing.list),".non.T.containing"))

z.score.rescaling = function(x,y,z){((z - mean(y,na.rm = TRUE))*(sd(x,na.rm = TRUE)/sd(y,na.rm = TRUE))) + mean(x,na.rm = TRUE)}

### raw signal of kmers ###
bam = get(load(paste0(resultsFolder,"nanoid/bam.RData")))
read.to.fast5.name = get(load(paste0(resultsFolder,"nanoid/read.to.fast5.name.RData")))

sequencing.summary = get(load(paste0(resultsFolder,"nanoid/sequencing.summary.RData")))

bamNames <- names(bam)
index.list = split(1:length(bamNames), ceiling(1:length(bamNames)/100))
#@ read.sequence.list = get(load(paste0(resultsFolder,"nanoid/read.sequence.list.RData")))

### GARR optimized code

start = 1

if(file.exists(where=paste0(resultsFolder,"nanoid/raw.signal.csv")))
{
  nrowMat <- system(paste0("cat ",paste0(resultsFolder,"nanoid/raw.signal.csv")," | wc -l"), intern=TRUE)
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
    if(j%%25==0){print(j);system("grep MemFree /proc/meminfo")}

    read.sequence.list = readList(paste0(resultsFolder,"nanoid/read.sequence.list.llo"),index=index.list[[j]])
    names(read.sequence.list) <- bamNames[index.list[[j]]]

    raw.signal.tmp = foreach(n = bamNames[index.list[[j]]],.noexport = setdiff(ls(),c("bam","read.to.fast5.name","read.sequence.list","sequencing.summary"))) %dopar% build.raw.signal.list(n)
    raw.signal.tmp <- do.call("rbind", raw.signal.tmp)
    
    rownames(raw.signal.tmp) <- bamNames[index.list[[j]]]
 
    colnames(raw.signal.tmp) = col.names
 
    fwrite(x=raw.signal.tmp
           ,file=paste0(resultsFolder,"nanoid/raw.signal.csv")
           ,append=TRUE)
    
    #@ rm(raw.signal.tmp)
    #@ rm(raw.signal.list)
    #@ gc()
  }
}
