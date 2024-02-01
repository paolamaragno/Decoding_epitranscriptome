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

### Model loading
## All available models
train.models <- list.files(trainedModelFolder)
train.models <- train.models[grepl("h5$",train.models)]

### Features data
read.based.parameters.files <- list.files(paste0(resultsFolder,"nanoid"))[grep("read.based.parameters",list.files(paste0(resultsFolder,"nanoid")))]

for(train.model in train.models)
{
	modificationProbabilities <- sapply(read.based.parameters.files,function(read.based.parameters.file)
	{
		nrowMat <- system(paste0("cat ",resultsFolder,"nanoid/",read.based.parameters.file," | wc -l"), intern=TRUE)
		nrowMat <- (as.numeric(nrowMat)-1)

		index.list = split(1:nrowMat, ceiling(1:nrowMat/10000))
		index.list = append(0,index.list)

		probaTmp = tryCatch(mclapply(seq_along(index.list)[-1],function(i)
		{
			print(i)
			read.based.parameters = fread(file=paste0(resultsFolder,"nanoid/",read.based.parameters.file),header=TRUE,nrow=10000,skip=max(index.list[[(i-1)]]))
			rownames(read.based.parameters) <- unlist(read.based.parameters[,1])
			read.based.parameters <- read.based.parameters[,-1]
			read.based.parameters[is.na(read.based.parameters)] = 0
			train.model <- load_model_hdf5(paste0(trainedModelFolder,"/",train.model))
			p=predict_proba(train.model,as.matrix(read.based.parameters))
			rm(read.based.parameters)
			gc()
			p
		},mc.cores=mc.cores),error=function(e)return(rep(NaN,nrowMat)))
		probaTmp <- unlist(probaTmp)

		rownamesMat <- as.character(unlist(read.table(paste0(resultsFolder,"nanoid/rownames.csv"),header=TRUE)))

		names(probaTmp) <- rownamesMat
		probaTmp
	})

	if(is.list(modificationProbabilities))
	{
		modificationProbabilities <- unlist(unname(modificationProbabilities))
	}else{
		modificationProbabilities <- modificationProbabilities[,1]
	}

	save(file=paste0(resultsFolder,"nanoid/nanoID.modification.probabilities_",train.model,".RData"),x=modificationProbabilities)
}