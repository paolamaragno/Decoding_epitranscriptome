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

### build.trace.list ###
build.trace.list = function(read.name){
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
      trace.move = trace[,move > 0]
      trace.no.move = trace[,move == 0]
      
      results = c(read.name,
                  
                  tryCatch({apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "A",drop = FALSE],1,mean)},error = function(x){0}),
                  tryCatch({apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "A",drop = FALSE],1,mean)},error = function(x){0}),
                  tryCatch({apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "A",drop = FALSE],1,mean)},error = function(x){0}),
                  
                  tryCatch({as.vector(t(apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "A",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                  tryCatch({as.vector(t(apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "A",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                  tryCatch({as.vector(t(apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "A",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                  
                  tryCatch({apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "C",drop = FALSE],1,mean)},error = function(x){0}),
                  tryCatch({apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "C",drop = FALSE],1,mean)},error = function(x){0}),
                  tryCatch({apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "C",drop = FALSE],1,mean)},error = function(x){0}),
                  
                  tryCatch({as.vector(t(apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "C",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                  tryCatch({as.vector(t(apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "C",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                  tryCatch({as.vector(t(apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "C",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                  
                  tryCatch({apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "G",drop = FALSE],1,mean)},error = function(x){0}),
                  tryCatch({apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "G",drop = FALSE],1,mean)},error = function(x){0}),
                  tryCatch({apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "G",drop = FALSE],1,mean)},error = function(x){0}),
                  
                  tryCatch({as.vector(t(apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "G",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                  tryCatch({as.vector(t(apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "G",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                  tryCatch({as.vector(t(apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "G",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                  
                  tryCatch({apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "T",drop = FALSE],1,mean)},error = function(x){0}),
                  tryCatch({apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "T",drop = FALSE],1,mean)},error = function(x){0}),
                  tryCatch({apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "T",drop = FALSE],1,mean)},error = function(x){0}),
                  
                  tryCatch({as.vector(t(apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "T",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                  tryCatch({as.vector(t(apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "T",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                  tryCatch({as.vector(t(apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "T",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                  
                  tryCatch({apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) %in% c("A","C","G"),drop = FALSE],1,mean)},error = function(x){0}),
                  tryCatch({apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) %in% c("A","C","G"),drop = FALSE],1,mean)},error = function(x){0}),
                  tryCatch({apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) %in% c("A","C","G"),drop = FALSE],1,mean)},error = function(x){0}),
                  
                  tryCatch({as.vector(t(apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) %in% c("A","C","G"),drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                  tryCatch({as.vector(t(apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) %in% c("A","C","G"),drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                  tryCatch({as.vector(t(apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) %in% c("A","C","G"),drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)})
      )
      names(results) = NULL
      return(results)
    },silent = TRUE)
  }
  if(!exists("results")) # Once I got an error, event.repeats longer than reevent.raw.means...
  {
    c(read.name,rep("NaN",1320))
  }
}

col.names = c("read.name",
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_A_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_A_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_A_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_A_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_A_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_A_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_C_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_C_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_C_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_C_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_C_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_C_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_G_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_G_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_G_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_G_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_G_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_G_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_T_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_T_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_T_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_T_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_T_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_T_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_ACG_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_ACG_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_ACG_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_ACG_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_ACG_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_ACG_"),seq(0,1,length.out = 10),paste0)))

### probability of traces ###
bam = get(load(paste0(resultsFolder,"nanoid/bam.RData")))
bamNames <- names(bam)
read.to.fast5.name = get(load(paste0(resultsFolder,"nanoid/read.to.fast5.name.RData")))
#@ read.sequence.list = get(load(paste0(resultsFolder,"nanoid/read.sequence.list.RData")))

index.list = split(1:length(bamNames), ceiling(1:length(bamNames)/100))

start = 1

if(file.exists(where=paste0(resultsFolder,"nanoid/traces.csv")))
{
  nrowMat <- system(paste0("cat ",paste0(resultsFolder,"nanoid/traces.csv")," | wc -l"), intern=TRUE)
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
  
   read.sequence.list <- readList(paste0(resultsFolder,"nanoid/read.sequence.list.llo"),index.list[[j]])
   names(read.sequence.list) <- bamNames[index.list[[j]]]
  
   trace.list.tmp = foreach(n = bamNames[index.list[[j]]],.noexport = setdiff(ls(),c("bam","read.to.fast5.name","read.sequence.list"))) %dopar% build.trace.list(n)
   trace.list.tmp <- do.call("rbind", trace.list.tmp)
   
   rownames(trace.list.tmp) <- bamNames[index.list[[j]]]

   colnames(trace.list.tmp) = col.names

   fwrite(x=trace.list.tmp
          ,file=paste0(resultsFolder,"nanoid/traces.csv")
          ,append=TRUE)

   #@ rm(trace.list.tmp)
   #@ rm(read.sequence.list)
   #@ gc()
  }
}
