#!/usr/bin/env Rscript
### load input variables ###
args = commandArgs(trailingOnly=TRUE)

for(v in args)
{
	vTmp <- strsplit(v,"=")[[1]]
	assign(vTmp[[1]],vTmp[[2]])
}

source(paste0(scriptsFolder,"nanoID.environmentSetUp.R"))

system(paste0("cut --complement -f1 -d ',' ",resultsFolder,"nanoid/traces.csv > ",resultsFolder,"nanoid/tracesB.csv"))

system(paste0("cut --complement -f1 -d ',' ",resultsFolder,"nanoid/traces.add.on.csv > ",resultsFolder,"nanoid/traces.add.onB.csv"))

system(paste0("cut --complement -f1 -d ',' ",resultsFolder,"nanoid/raw.signal.csv > ",resultsFolder,"nanoid/raw.signalB.csv"))

system(paste0("cut --complement -f1 -d ',' ",resultsFolder,"nanoid/raw.signal.five.mers.csv > ",resultsFolder,"nanoid/raw.signal.five.mersB.csv"))

system(paste0("cut --complement -f1 -d ',' ",resultsFolder,"nanoid/raw.signal.five.mers.add.on.csv > ",resultsFolder,"nanoid/raw.signal.five.mers.add.onB.csv"))

system(paste0("paste -d ',' ",resultsFolder,"nanoid/rownames.csv ",resultsFolder,"nanoid/read.based.mismatch.identification.csv > ",resultsFolder,"nanoid/read.based.mismatch.identificationB.csv"))
system(paste0("mv ",resultsFolder,"nanoid/read.based.mismatch.identificationB.csv ",resultsFolder,"nanoid/read.based.mismatch.identification.csv"))

system(paste0("paste ",resultsFolder,"nanoid/read.based.mismatch.identification.csv "
					  ,resultsFolder,"nanoid/tracesB.csv "
					  ,resultsFolder,"nanoid/traces.add.onB.csv "
					  ,resultsFolder,"nanoid/raw.signalB.csv "
					  ,resultsFolder,"nanoid/raw.signal.five.mersB.csv "
					  ,resultsFolder,"nanoid/raw.signal.five.mers.add.onB.csv -d ',' > "
					  ,resultsFolder,"nanoid/read.based.parameters.csv"))

system(paste0("rm ",resultsFolder,"nanoid/tracesB.csv "
					  ,resultsFolder,"nanoid/traces.add.onB.csv "
					  ,resultsFolder,"nanoid/raw.signalB.csv "
					  ,resultsFolder,"nanoid/raw.signal.five.mersB.csv "
					  ,resultsFolder,"nanoid/raw.signal.five.mers.add.onB.csv"))

read.based.parameters <- as.matrix(fread(paste0(resultsFolder,"nanoid/read.based.parameters.csv"),sep=",",header=TRUE,nrows=10))
rownames(read.based.parameters) <- read.based.parameters[,1]
read.based.parameters <- read.based.parameters[,-1]

colnamesTmp = c("read.name",paste0("P",substr(rep("0000",ncol(read.based.parameters)),1,4-nchar(1:ncol(read.based.parameters))),1:ncol(read.based.parameters),"-",colnames(read.based.parameters)))

fwrite(x=as.list(colnamesTmp),file=paste0(resultsFolder,"nanoid/colnames.csv"),row.names=TRUE,col.names=TRUE,quote=FALSE)
system(paste0("cut --complement -f1 -d ',' ",resultsFolder,"nanoid/colnames.csv > ",resultsFolder,"nanoid/colnamesB.csv"))
system(paste0("mv ",resultsFolder,"nanoid/colnamesB.csv ",resultsFolder,"nanoid/colnames.csv"))

system(paste0("tail -n+2 ",resultsFolder,"nanoid/read.based.parameters.csv > ",resultsFolder,"nanoid/read.based.parametersB.csv"))
system(paste0("cat ",resultsFolder,"nanoid/colnames.csv ",resultsFolder,"nanoid/read.based.parametersB.csv > ",resultsFolder,"read.based.parameters.csv"))
system(paste0("rm ",resultsFolder,"nanoid/read.based.parametersB.csv"))
