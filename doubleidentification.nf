params.python_path='/home/hyx/system/anaconda3/bin/python'
params.type = 'double'
params.file_input = '/home/hyx/nf/data/doublet/HM_Mix_6K.rds'
params.file_output = '/home/hyx/nf/output/doublet/6k.csv'
params.file_label = '/home/hyx/nf/data/emptyDroplet/pbmc3k.csv'
params.fdr = 0.005

// Channel DRead_data
// Channel DPreprocess_data

process Read() {
	debug true
	tag 'Read'

	input:
	val pp

	output:
	stdout

	script:
	"""
	#!/usr/bin/env Rscript
	rm(list=ls())
	setwd('/home/hyx/nf/script')
	#Dataset <- readRDS("/home/hyx/nf/data/doublet/HM_Mix_12K.rds")
	Dataset <- readRDS("$pp")
	source('Lib.R')
	library(Seurat)
	H <- CreateSeuratObject(Dataset[[1]])
	H\$label <- Dataset[[2]]
	save(H, file='/home/hyx/nf/temp/H_0.RData')
	"""
}

process Preprocessing {
	debug true
	tag 'Preprocessing'
	input:
	val x

	output:
	stdout
	script:
	'''
	#!/usr/bin/env Rscript
	rm(list=ls())
	setwd('/home/hyx/nf/script')
	#Dataset <- readRDS("/home/hyx/nf/data/doublet/HM_Mix_12K.rds")
    #load(file="/home/hyx/nf/temp/tttt.RData")
	load(file='/home/hyx/nf/temp/H_0.RData')
	source('Lib.R')
	library(Seurat)
	#H <- CreateSeuratObject(Dataset[[1]])
	#H$label <- Dataset[[2]]
	H <- NormalizeData(H)
	H <- FindVariableFeatures(H, selection.method = "vst", nfeatures = 2000)
	H <- ScaleData(H)
	H <- RunPCA(H)
	H <- RunUMAP(H, dims = 1:10)
	#A <- gsub('doublet','Doublet',H@meta.data$label)
	#A <- gsub('singlet','Singlet',A)
	save(H,file="/home/hyx/nf/temp/H_1.RData")
	'''
}
// process EmptyDrops{
// 	debug true
// 	tag 'EmptyDrops'

// 	input:
// 	val x

// 	output:
// 	stdout

// 	script:
// 		"""
// 		#! $params.python_path


// 		"""
// }
process NeighborsDefine {
	debug true
	tag 'NeighborsDefine'

	input:
	val x

	output:
	stdout

	script:
	'''
	#!/usr/bin/env Rscript
	setwd('/home/hyx/nf/script')
	source('Lib.R')
	library(Seurat)	
	#从目录导入H
	load(file="/home/hyx/nf/temp/H_1.RData")
	A <- gsub('doublet','Doublet',H@meta.data$label)
	A <- gsub('singlet','Singlet',A)
	save(A,file="/home/hyx/nf/temp/A.RData")
	'''
}
process DoubletSim {
	debug true
	tag 'DoubletSim'

	input:
	val x

	output:
	stdout

	script:
	'''
	#!/usr/bin/env Rscript
	setwd('/home/hyx/nf/script')
	source('Lib.R')
	library(Seurat)	

	load(file="/home/hyx/nf/temp/H_1.RData")	
	load(file="/home/hyx/nf/temp/A.RData")
	sweep.res.list_H <- paramSweep_v3(H, PCs = 1:10, sct = FALSE)
	gt.calls <- A
	sweep.stats_H <- summarizeSweep(sweep.res.list_H, GT = TRUE, GT.calls = gt.calls)
	bcmvn_H <- find.pK(sweep.stats_H)	

	#H -> H_2.RData
	save(H,file="/home/hyx/nf/temp/H_2.RData")
	#print('stage_3')
	'''


}
process DoubletProb{
	debug true
	tag 'DoubletProb'

	input:
	val x

	output:
	stdout

	script:
	'''
	#!/usr/bin/env Rscript
	setwd('/home/hyx/nf/script')
	source('Lib.R')
	library(Seurat)	

	load(file="/home/hyx/nf/temp/H_2.RData")
	#load(file="/home/hyx/nf/temp/A.RData")
	#Est
	true_double <- H@meta.data$label[H@meta.data$label=='doublet']
	B <- length(true_double)/length(H@meta.data$label)
	nExp_poi <- round(B*nrow(H@meta.data))

	#Run
	H <- doubletFinder_v3(H, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

	save(H,file="/home/hyx/nf/temp/H_3.RData")
	#print('stage_4')

	'''
}
process Validation {
	debug true
	tag 'Validation'

	input:
	val params.file_output

	output:
	stdout

	script:
	"""
	#!/usr/bin/env Rscript
	setwd('/home/hyx/nf/script')
	source('Lib.R')
	library(Seurat)

	load(file="/home/hyx/nf/temp/H_3.RData")
	load(file="/home/hyx/nf/temp/A.RData")

	true_value <- A
	predict_value <- H@meta.data[[length(H@meta.data)]]
	#简单的建立字典的方法
	y = c("Singlet"=0, "Doublet"=1)
	r = y[as.character(H@meta.data[[length(H@meta.data)]])]
	names(r) = NULL
	H@meta.data\$DF <- r
	r1 = y[as.character(A)]
	names(r1)=NULL
	H\$label_true <- r1


	true_value <- H@meta.data\$label_true
	predict_value <- H@meta.data\$DF#分别将16和17列取出来
	write.csv(predict_value,"$params.file_output")

	retrieved=sum(predict_value)
	precision=sum(true_value & predict_value)/retrieved#0.814
	recall=sum(predict_value & true_value)/sum(true_value)#0.970
	F1_score = precision*recall

	library(modEvA)
	auc=AUC(obs=H@meta.data\$label_true,pred=H@meta.data[[5]],curve = "ROC", simplif=TRUE, main = "ROC curve")
	#auc
	aupr=AUC(obs=H@meta.data\$label_true,pred=H@meta.data[[5]],curve = "PR", simplif=TRUE, main = "PR curve")
	#aupr

	paste("precision: ",precision)
	paste("recall: ",recall)
	paste("F1_score: ",F1_score)
	paste("auc: ",auc)
	paste("aupr: ",aupr)


	"""
}
process ERead{
	debug true
	tag 'ERead'

	input:
	val ppp

	output:
	stdout

	script:
		"""
		#!/usr/bin/env Rscript
		rm(list=ls())
		setwd('/home/hyx/nf/script')
		library(DropletUtils)
		fname <- file.path('./$ppp')
		sce <- read10xCounts(fname, col.names=TRUE)
		save(sce,file="/home/hyx/nf/temp/e_1.RData")
		"""
}

process EPreprocessing {
	debug true
	tag 'EPreprocessing'

	input:
	val x

	output:
	stdout

	script:

	'''
	#!/usr/bin/env Rscript
	#rm(list=ls())
	setwd('/home/hyx/nf/script')
    library(DropletUtils)
	set.seed(100)
	'''

}
process EmptyDrops{
	debug true
	tag 'EmptyDrops'

	input:
	val x

	output:
	stdout

	script:
	"""
	#!/usr/bin/env Rscript
	library(DropletUtils)
	load(file="/home/hyx/nf/temp/e_1.RData")
	e.out <- emptyDrops(counts(sce))
	save(e.out,file="/home/hyx/nf/temp/e_2.RData")
	save(sce,file="/home/hyx/nf/temp/e_3.RData")

	"""
}

process Filtering {
	debug true
	tag 'Filtering'

	input:
	val x

	output:
	stdout

	script:
	"""
	#!/usr/bin/env Rscript
    library(DropletUtils)
	load(file="/home/hyx/nf/temp/e_2.RData")
    load(file="/home/hyx/nf/temp/e_3.RData")
    sce\$empty <- ifelse(is.na(e.out\$FDR) | e.out\$FDR <= $params.fdr, "False", "True")
	save(sce,file="/home/hyx/nf/temp/e_4.RData")

	"""
}

process EValidation {
	debug true
	tag 'EValidation'

	input:
	val x

	output:
	stdout

	script:
	"""
	#!/usr/bin/env Rscript
	setwd('/home/hyx/nf/script')
    labels <- read.csv("$params.file_label", row.names = 1)
    load(file="/home/hyx/nf/temp/e_4.RData")
    true_value <- ifelse(labels\$simulated == "False", 0, 1)
    predict_value <- ifelse(sce\$empty == "False", 0, 1)


	source('/home/hyx/nf/script/Lib.R')
	library(Seurat)
	
    write.csv(predict_value,'$params.file_output')

	retrieved=sum(predict_value)
	precision=sum(true_value & predict_value)/retrieved#0.814
	recall=sum(predict_value & true_value)/sum(true_value)#0.970
	F1_score = precision*recall

	library(modEvA)
	auc=AUC(obs=true_value,pred=predict_value,curve = "ROC", simplif=TRUE, main = "ROC curve")
	#auc
	aupr=AUC(obs=true_value,pred=predict_value,curve = "PR", simplif=TRUE, main = "PR curve")
	#aupr

	paste("precision: ",precision)
	paste("recall: ",recall)
	paste("F1_score: ",F1_score)
	paste("auc: ",auc)
	paste("aupr: ",aupr)


	"""
}






workflow{
	// wk = Channel.of("")
	// def p = Channel.fromPath(params.file_input)
	// def o = Channel.fromPath(params.file_output)
    if(params.type=='double')
        Read(params.file_input)|Preprocessing|NeighborsDefine|DoubletSim|DoubletProb|Validation
    else if(params.type=='empty')
        // ERead(params.file_input)|EPreprocessing|EmptyDrops|Filter|EValidation
		ERead(params.file_input)|EPreprocessing|EmptyDrops|Filtering|EValidation
		// Filtering(1)|EValidation

}