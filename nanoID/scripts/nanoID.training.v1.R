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

### Features data
## Unlabeled
read.based.parameters.files.unlabeled <- list.files(paste0(resultsDir,"/",unlabeled_time),recursive=TRUE)
read.based.parameters.files.unlabeled <- read.based.parameters.files.unlabeled[grepl("nanoid/read.based.parameters.csv",read.based.parameters.files.unlabeled)]

## Fully labeled
read.based.parameters.files.fullylabeled <- list.files(paste0(resultsDir,"/",fullylabeled_time),recursive=TRUE)
read.based.parameters.files.fullylabeled <- read.based.parameters.files.fullylabeled[grepl("nanoid/read.based.parameters.csv",read.based.parameters.files.fullylabeled)]

### Parameters
nFeatures = 4694
epochs = 40
shuffle = TRUE

### Training data
## Unlabeled
read.based.parameters.unlabeled <- lapply(read.based.parameters.files.unlabeled,function(i)
{
	read.based.parameters.unlabeled.tmp <- fread(paste0(resultsDir,"/",unlabeled_time,"/",i))
	reads.tmp <- unlist(read.based.parameters.unlabeled.tmp[,1])
	read.based.parameters.unlabeled.tmp <- read.based.parameters.unlabeled.tmp[,-1]
	read.based.parameters.unlabeled.tmp[is.na(read.based.parameters.unlabeled.tmp)] = 0
	read.based.parameters.unlabeled.tmp <- as.matrix(read.based.parameters.unlabeled.tmp)
	rownames(read.based.parameters.unlabeled.tmp) <- reads.tmp
	read.based.parameters.unlabeled.tmp
})
read.based.parameters.unlabeled <- do.call("rbind", read.based.parameters.unlabeled)

## Fullylabeled
read.based.parameters.fullylabeled <- lapply(read.based.parameters.files.fullylabeled,function(i)
{
	read.based.parameters.fullylabeled.tmp <- fread(paste0(resultsDir,"/",fullylabeled_time,"/",i))
	reads.tmp <- unlist(read.based.parameters.fullylabeled.tmp[,1])
	read.based.parameters.fullylabeled.tmp <- read.based.parameters.fullylabeled.tmp[,-1]
	read.based.parameters.fullylabeled.tmp[is.na(read.based.parameters.fullylabeled.tmp)] = 0
	read.based.parameters.fullylabeled.tmp <- as.matrix(read.based.parameters.fullylabeled.tmp)
	rownames(read.based.parameters.fullylabeled.tmp) <- reads.tmp
	read.based.parameters.fullylabeled.tmp
})
read.based.parameters.fullylabeled <- do.call("rbind", read.based.parameters.fullylabeled)

### Training set sizes
sizes <- round(seq(0,min(nrow(read.based.parameters.unlabeled),nrow(read.based.parameters.fullylabeled)),length.out=(as.numeric(Nsamplings)+1)))
sizes <- unique(c(as.numeric(unlist(strsplit(trainingSizes,"_"))),sizes))
sizes <- sort(sizes,decreasing=FALSE) # I start from the smallest
sizes <- sizes[-1] # I remove the 0

### Seeds
seeds <- 1:as.numeric(Nseeds)

### Training
for(seedTmp in seeds)
{
	for(sizeTmp in sizes)
	{
		set.seed(seedTmp)
		labeledTmp <- read.based.parameters.fullylabeled[sample(1:nrow(read.based.parameters.fullylabeled),size=sizeTmp),]
		
		set.seed(seedTmp)
		unlabeledTmp <- read.based.parameters.unlabeled[sample(1:nrow(read.based.parameters.unlabeled),size=sizeTmp),]

		### Train and test with 5-fold cross validation
		y_all = c(rep(1,nrow(labeledTmp)),rep(0,nrow(unlabeledTmp)))
		x_all = rbind(labeledTmp,unlabeledTmp)
		
		rm(labeledTmp);gc()
		rm(unlabeledTmp);gc()

		### Model definition
		train.model <- keras_model_sequential()
		train.model %>% 
		layer_batch_normalization(input_shape = nFeatures) %>%
		layer_dense(units = 512, activation = "relu", input_shape = nFeatures) %>%
		layer_batch_normalization() %>%
		layer_dropout(rate = 0.25) %>%
		layer_dense(units = 256, activation = "relu") %>%
		layer_batch_normalization() %>%
		layer_dropout(rate = 0.25) %>%
		layer_dense(units = 128, activation = "relu") %>%
		layer_batch_normalization() %>%
		layer_dropout(rate = 0.25) %>%
		layer_dense(units = 64, activation = "relu") %>%
		layer_batch_normalization() %>%
		layer_dropout(rate = 0.25) %>%
		layer_dense(units = 32, activation = "relu") %>%
		layer_batch_normalization() %>%
		layer_dropout(rate = 0.25) %>%
		layer_dense(units = 16, activation = "relu") %>%
		layer_batch_normalization() %>%
		layer_dropout(rate = 0.25) %>%
		layer_dense(units = 8, activation = "relu") %>%
		layer_batch_normalization() %>%
		layer_dropout(rate = 0.25) %>%
		layer_dense(units = 1, activation = "sigmoid")
	  
		train.model %>% summary()
	  
		train.model %>% compile(
			#optimizer = 'adam',
			#optimizer = optimizer_adam(),
			optimizer = optimizer_rmsprop(),
			#optimizer = optimizer_adadelta(),
			#optimizer = optimizer_adagrad(),
			#optimizer = optimizer_adamax(),
			#optimizer = optimizer_nadam(),
			#optimizer = optimizer_sgd(),
			loss = 'binary_crossentropy',
			metrics = list('accuracy')
		)

		set.seed(seedTmp)
		y_train_tmp_0 = sample(which(y_all==0),round(0.7*length(which(y_all==0))))
		y_train_tmp_1 = sample(which(y_all==1),round(0.7*length(which(y_all==1))))

		y_test_tmp_0 = setdiff(which(y_all==0),y_train_tmp_0)
		y_test_tmp_1 = setdiff(which(y_all==1),y_train_tmp_1)

		y_train_tmp = c(rep(1,length(y_train_tmp_1)),rep(0,length(y_train_tmp_0)))
		x_train_tmp = rbind(x_all[y_train_tmp_1,],x_all[y_train_tmp_0,])

		train.model %>% fit(
			x = x_train_tmp,
			y = y_train_tmp,
			epochs = epochs,
			batch_size = 32,
			validation_split = 0.2,
			shuffle = shuffle,
			verbose = 1
		)

		pred.table.train = table("model"=predict_classes(train.model,x_train_tmp),"true"=y_train_tmp)
		accuracy.train = (sum(diag(pred.table.train)))/sum(pred.table.train)

		rm(x_train_tmp);gc()

		y_test_tmp = c(rep(1,length(y_test_tmp_1)),rep(0,length(y_test_tmp_0)))
		x_test_tmp = rbind(x_all[y_test_tmp_1,],x_all[y_test_tmp_0,])

		test_reads_tmp = rownames(x_test_tmp)

		pred.table.test = table("model"=predict_classes(train.model,x_test_tmp),"true"=y_test_tmp)
		accuracy.test = (sum(diag(pred.table.test)))/sum(pred.table.test)

		rm(x_test_tmp);gc()

		save_model_hdf5(train.model,paste0(trainedModelFolder,"/train.model","_",length(y_all),"_",seedTmp,".h5"),overwrite = TRUE,include_optimizer = TRUE)
		saveRDS(list("train_reads_unlabeled" = train_reads_tmp[y_train_tmp==0],
					 "train_reads_labeled" = train_reads_tmp[y_train_tmp==1],
					 "test_reads_unlabeled" = test_reads_tmp[y_test_tmp==0],
					 "test_reads_labeled" = test_reads_tmp[y_test_tmp==1],
					 "accuracy.train" = accuracy.train,
					 "accuracy.test" = accuracy.test)
			   ,paste0(trainedModelFolder,"/train.test.reads","_",length(y_all),"_",seedTmp,".rds"))
	}
}
