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

### build.trace.add.on.list ###
build.trace.add.on.list = function(read.name){
  if(moveSlot=="Move") # In this way, if the moveSlot is "foe", the production of the NaN matrix is faster.
  {
    try({
      fast5.file = read.to.fast5.name[read.name]
      read.sequence = rev(strsplit(read.sequence.list[[read.name]],split = "")[[1]])

      analyses = h5read(fast5.file,paste0("read_",read.name,"/Analyses"))
      slotName = sort(names(analyses)[grep("Basecall_1D",names(analyses))],decreasing=TRUE)[[1]]

      trace = h5read(fast5.file,paste0("read_",read.name,"/Analyses/",slotName,"/BaseCalled_template/Trace"))
      rownames(trace) = c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T")
      move = h5read(fast5.file,paste0("read_",read.name,"/Analyses/",slotName,"/BaseCalled_template/",moveSlot))
      move.rle = Rle(move)
      event.repeats = move[move == 1]
      event.repeats[cumsum(runLength(move.rle)[runValue(move.rle) == 1])] = runLength(move.rle)[runValue(move.rle) == 0] + 1
      if (move[length(move)] == 1){event.repeats[length(event.repeats)] = 1}
      colnames(trace) = rep(read.sequence,times = event.repeats)
      
      move.rle.A = Rle(move[colnames(trace) == "A"])
      move.rle.C = Rle(move[colnames(trace) == "C"])
      move.rle.G = Rle(move[colnames(trace) == "G"])
      move.rle.T = Rle(move[colnames(trace) == "T"])
      
      results = c(read.name,
                  
                  mean(runLength(move.rle)[runValue(move.rle) == 1]),
                  quantile(runLength(move.rle)[runValue(move.rle) == 1],seq(0,1,length.out = 100)),
                  
                  mean(runLength(move.rle)[runValue(move.rle) == 0]),
                  quantile(runLength(move.rle)[runValue(move.rle) == 0],seq(0,1,length.out = 100)),
                  
                  tryCatch({mean(runLength(move.rle.A)[runValue(move.rle.A) == 1])},error = function(x){0}),
                  tryCatch({quantile(runLength(move.rle.A)[runValue(move.rle.A) == 1],seq(0,1,length.out = 100))},error = function(x){0}),
                  tryCatch({mean(runLength(move.rle.A)[runValue(move.rle.A) == 0])},error = function(x){0}),
                  tryCatch({quantile(runLength(move.rle.A)[runValue(move.rle.A) == 0],seq(0,1,length.out = 100))},error = function(x){0}),
                  
                  tryCatch({mean(runLength(move.rle.C)[runValue(move.rle.C) == 1])},error = function(x){0}),
                  tryCatch({quantile(runLength(move.rle.C)[runValue(move.rle.C) == 1],seq(0,1,length.out = 100))},error = function(x){0}),
                  tryCatch({mean(runLength(move.rle.C)[runValue(move.rle.C) == 0])},error = function(x){0}),
                  tryCatch({quantile(runLength(move.rle.C)[runValue(move.rle.C) == 0],seq(0,1,length.out = 100))},error = function(x){0}),
                  
                  tryCatch({mean(runLength(move.rle.G)[runValue(move.rle.G) == 1])},error = function(x){0}),
                  tryCatch({quantile(runLength(move.rle.G)[runValue(move.rle.G) == 1],seq(0,1,length.out = 100))},error = function(x){0}),
                  tryCatch({mean(runLength(move.rle.G)[runValue(move.rle.G) == 0])},error = function(x){0}),
                  tryCatch({quantile(runLength(move.rle.G)[runValue(move.rle.G) == 0],seq(0,1,length.out = 100))},error = function(x){0}),
                  
                  tryCatch({mean(runLength(move.rle.T)[runValue(move.rle.T) == 1])},error = function(x){0}),
                  tryCatch({quantile(runLength(move.rle.T)[runValue(move.rle.T) == 1],seq(0,1,length.out = 100))},error = function(x){0}),
                  tryCatch({mean(runLength(move.rle.T)[runValue(move.rle.T) == 0])},error = function(x){0}),
                  tryCatch({quantile(runLength(move.rle.T)[runValue(move.rle.T) == 0],seq(0,1,length.out = 100))},error = function(x){0})
      )
      names(results) = NULL
      return(results)
    },silent = TRUE)
  }
  if(!exists("results")) # Once I got an error, event.repeats longer than reevent.raw.means...
  {
    c(read.name,rep("NaN",1010))
  }
}

col.names = c("read.name",
              paste0("1s","_mean"),
              paste0("1s_",1:100),
              paste0("0s","_mean"),
              paste0("0s_",1:100),
              paste0("1s","_A_mean"),
              paste0("1s_A_",1:100),
              paste0("0s","_A_mean"),
              paste0("0s_A_",1:100),
              paste0("1s","_C_mean"),
              paste0("1s_C_",1:100),
              paste0("0s","_C_mean"),
              paste0("0s_C_",1:100),
              paste0("1s","_G_mean"),
              paste0("1s_G_",1:100),
              paste0("0s","_G_mean"),
              paste0("0s_G_",1:100),
              paste0("1s","_T_mean"),
              paste0("1s_T_",1:100),
              paste0("0s","_T_mean"),
              paste0("0s_T_",1:100))

### probability of traces ###
bam = get(load(paste0(resultsFolder,"nanoid/bam.RData")))
bamNames <- names(bam)
gc()

read.to.fast5.name = get(load(paste0(resultsFolder,"nanoid/read.to.fast5.name.RData")))
#@ read.sequence.list = get(load(paste0(resultsFolder,"nanoid/read.sequence.list.RData")))

index.list = split(1:length(bamNames), ceiling(1:length(bamNames)/100))

start = 1

if(file.exists(where=paste0(resultsFolder,"nanoid/traces.add.on.csv")))
{
  nrowMat <- system(paste0("cat ",paste0(resultsFolder,"nanoid/traces.add.on.csv")," | wc -l"), intern=TRUE)
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
  
    if(j%%50==0){print(j);gc();system("grep MemFree /proc/meminfo")}
  
    read.sequence.list <- readList(paste0(resultsFolder,"nanoid/read.sequence.list.llo"),index.list[[j]])
    names(read.sequence.list) <- bamNames[index.list[[j]]]
  
    trace.add.on.list.tmp = foreach(n = bamNames[index.list[[j]]],.noexport = setdiff(ls(),c("bam","read.to.fast5.name","read.sequence.list"))) %dopar% build.trace.add.on.list(n)
    trace.add.on.list.tmp <- do.call("rbind", trace.add.on.list.tmp)
    
    rownames(trace.add.on.list.tmp) <- bamNames[index.list[[j]]]
 
    colnames(trace.add.on.list.tmp) = col.names
 
    fwrite(x=trace.add.on.list.tmp
           ,file=paste0(resultsFolder,"nanoid/traces.add.on.csv")
           ,append=TRUE)
    
    #@ rm(read.sequence.list)
    #@ rm(trace.add.on.list.tmp)
    #@ gc()
  }
}
