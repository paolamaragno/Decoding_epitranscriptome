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

### Bed file for counts estimation
bedFile <- read.table(genesBed,header=FALSE)
colnames(bedFile) <- c("chr","start","end","id","score","strand")
bedFile <- makeGRangesFromDataFrame(as.data.frame(bedFile),keep.extra.columns=TRUE)
names(bedFile) <- bedFile$id

### Features data
read.based.parameters.files <- list.files(path=resultsDir,pattern=,"nanoID.modification.probabilities",recursive=TRUE)

## I remove training conditions
read.based.parameters.files <- read.based.parameters.files[!(grepl(paste0("^",unlabeled_time,"/"),read.based.parameters.files)|grepl(paste0("^",fullylabeled_time,"/"),read.based.parameters.files))]

if(length(read.based.parameters.files)==0)
{
	print("Only conditions involved in training")
}else{

	## Models to loop for kinetic rates estimation
	train.models <- list.files(trainedModelFolder)
	train.models <- train.models[grepl("h5$",train.models)]
	
	train.models <- sapply(train.models,function(model)
	{
		read.based.parameters.files[grepl(paste0("nanoID.modification.probabilities_",model,".RData$"),read.based.parameters.files)]
	})

	for(j in colnames(train.models))
	{
		read.based.parameters.files <- train.models[,j]
		
		### Nascent and PreExisting counts
		countMatrices <- lapply(dirname(read.based.parameters.files),function(i)
		{
			alignmentsTmp <- get(load(paste0(resultsDir,i,"/bam.RData")))
			overlapsTmp <- findOverlaps(alignmentsTmp,bedFile)
			
			overlapsTmp <- overlapsTmp[isUnique(overlapsTmp@from),]

			alignmentsTmp <- cbind("readId"=names(alignmentsTmp)[overlapsTmp@from]
								 , "tuId"=names(bedFile)[overlapsTmp@to])

			classificationsTmp <- get(load(paste0(resultsDir,i,"/nanoID.modification.probabilities_",j,".RData")))

			commonReadsTmp <- intersect(alignmentsTmp[,"readId"],names(classificationsTmp))

			alignmentsTmp <- alignmentsTmp[alignmentsTmp[,"readId"]%in%commonReadsTmp,]


			classificationsTmp <- classificationsTmp[alignmentsTmp[,"readId"]]

			preexistingCountsTmp <- table(alignmentsTmp[classificationsTmp<0.5,"tuId"])
			nascentCountsTmp <- table(alignmentsTmp[classificationsTmp>0.5,"tuId"])

			allTuTmp <- union(names(preexistingCountsTmp),names(nascentCountsTmp))

			countsTmp <- cbind(preexistingCountsTmp[allTuTmp],nascentCountsTmp[allTuTmp])
			rownames(countsTmp) <- allTuTmp
			colnames(countsTmp) <- c("PreExisting","Nascent")
			countsTmp[!is.finite(countsTmp)] <- 0
			countsTmp
		})

		# countMatricesNorm <- lapply(countMatrices,function(m){normTmp <- apply(m,2,sum);t(apply(m,1,function(i)i/normTmp*10^6))})
		countMatricesNorm <- lapply(countMatrices,function(m){(10^6*m)/sum(m)})

		### Calculate synthesis rate and stability
		## Cell cycle
		alpha = log(2)/as.numeric(CCL)

		## Labeling times
		labelingTimes = as.numeric(sapply(strsplit(read.based.parameters.files,"/"),"[[",1))

		outMatrices <- lapply(seq_along(countMatrices),function(i)
		{
			m <- countMatrices[[i]]
			labelingTimes <- labelingTimes[[i]]	
			## Decay rate
			decay.rate.single.moleculeTmp = - alpha - (1/labelingTimes)*log(1 - m[,"Nascent"]/apply(m,1,sum))
			decay.rate.single.moleculeTmp[decay.rate.single.moleculeTmp <= 0] = NA
			
			## Synthesis rate
			synthesis.rate.single.moleculeTmp = apply(m,1,sum)*(alpha + decay.rate.single.moleculeTmp)

			## Half life
			half.lives.single.moleculeTmp = log(2)/decay.rate.single.moleculeTmp

			cbind(m,"Synthesis"=synthesis.rate.single.moleculeTmp,"Degradation"=decay.rate.single.moleculeTmp,"HalfLife"=half.lives.single.moleculeTmp)
		})

		outMatricesNorm <- lapply(seq_along(countMatricesNorm),function(i)
		{
			m <- countMatricesNorm[[i]]
			labelingTimes <- labelingTimes[[i]]	
			## Decay rate
			decay.rate.single.moleculeTmp = - alpha - (1/labelingTimes)*log(1 - m[,"Nascent"]/apply(m,1,sum))
			decay.rate.single.moleculeTmp[decay.rate.single.moleculeTmp <= 0] = NA
			
			## Synthesis rate
			synthesis.rate.single.moleculeTmp = apply(m,1,sum)*(alpha + decay.rate.single.moleculeTmp)

			## Half life
			half.lives.single.moleculeTmp = log(2)/decay.rate.single.moleculeTmp

			cbind(m,"Synthesis"=synthesis.rate.single.moleculeTmp,"Degradation"=decay.rate.single.moleculeTmp,"HalfLife"=half.lives.single.moleculeTmp)
		})

		names(outMatrices) <- names(outMatricesNorm) <- sapply(read.based.parameters.files,function(k)basename(dirname(dirname(k))))

		## Output
		save(outMatrices,file=paste0(resultsDir,"/kineticRates/","kineticRates_",j,".RData"))
		save(outMatricesNorm,file=paste0(resultsDir,"/kineticRates/","kineticRatesNorm_",j,".RData"))
	}
}
