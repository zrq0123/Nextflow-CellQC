#install.packages('Seurat')
paramSweep_v3 <- function(seu, PCs=1:10, sct = FALSE, num.cores=1) {
  require(Seurat); require(fields);
  ## Set pN-pK param sweep ranges
  pK <- c(0.0005, 0.001, 0.005, seq(0.01,0.3,by=0.01))
  pN <- seq(0.05,0.3,by=0.05)
  
  ## Remove pK values with too few cells
  min.cells <- round(nrow(seu@meta.data)/(1-0.05) - nrow(seu@meta.data))
  pK.test <- round(pK*min.cells)
  pK <- pK[which(pK.test >= 1)]
  
  ## Extract pre-processing parameters from original data analysis workflow
  orig.commands <- seu@commands
  
  ## Down-sample cells to 20000 (when applicable) for computational effiency
  if (nrow(seu@meta.data) > 20000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 20000, replace=FALSE)]
    data <- seu@assays$RNA$counts[ , real.cells]
    n.real.cells <- ncol(data)
  }
  
  if (nrow(seu@meta.data) <= 20000){
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA$counts
    n.real.cells <- ncol(data)
  }
  
  ## Iterate through pN, computing pANN vectors at varying pK
  #no_cores <- detectCores()-1
  if(num.cores>1){
    require(parallel)
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)),
                        FUN = parallel_paramSweep_v3,
                        n.real.cells,
                        real.cells,
                        pK,
                        pN,
                        data,
                        orig.commands,
                        PCs,
                        sct,mc.cores=num.cores)
    stopCluster(cl)
  }else{
    output2 <- lapply(as.list(1:length(pN)),
                      FUN = parallel_paramSweep_v3,
                      n.real.cells,
                      real.cells,
                      pK,
                      pN,
                      data,
                      orig.commands,
                      PCs,
                      sct)
  }
  
  ## Write parallelized output into list
  sweep.res.list <- list()
  list.ind <- 0
  for(i in 1:length(output2)){
    for(j in 1:length(output2[[i]])){
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }
  
  ## Assign names to list of results
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_" ))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
  
}


summarizeSweep <- function(sweep.list, GT = FALSE, GT.calls = NULL) {
  require(KernSmooth); require(ROCR)
  ## Set pN-pK param sweep ranges
  name.vec <- names(sweep.list)
  name.vec <- unlist(strsplit(name.vec, split="pN_"))
  name.vec <- name.vec[seq(2, length(name.vec), by=2)]
  name.vec <- unlist(strsplit(name.vec, split="_pK_"))
  pN <- as.numeric(unique(name.vec[seq(1, length(name.vec), by=2)]))
  pK <- as.numeric(unique(name.vec[seq(2, length(name.vec), by=2)]))
  
  ## Initialize data structure w/ or w/o AUC column, depending on whether ground-truth doublet classifications are available
  if (GT == TRUE) {
    sweep.stats <- as.data.frame(matrix(0L, nrow=length(sweep.list), ncol=4))
    colnames(sweep.stats) <- c("pN","pK","AUC","BCreal")
    sweep.stats$pN <- factor(rep(pN, each=length(pK), levels = pN))
    sweep.stats$pK <- factor(rep(pK, length(pN),levels = pK))
  }
  
  if (GT == FALSE) {
    sweep.stats <- as.data.frame(matrix(0L, nrow=length(sweep.list), ncol=3))
    colnames(sweep.stats) <- c("pN","pK","BCreal")
    sweep.stats$pN <- factor(rep(pN, each=length(pK), levels = pN))
    sweep.stats$pK <- factor(rep(pK, length(pN),levels = pK))
  }
  
  ## Perform pN-pK parameter sweep summary
  for (i in 1:length(sweep.list)) {
    res.temp <- sweep.list[[i]]
    
    ## Use gaussian kernel density estimation of pANN vector to compute bimodality coefficient
    gkde <- approxfun(bkde(res.temp$pANN, kernel="normal"))
    x <- seq(from=min(res.temp$pANN), to=max(res.temp$pANN), length.out=nrow(res.temp))
    sweep.stats$BCreal[i] <- bimodality_coefficient(gkde(x))
    
    if (GT == FALSE) { next }
    
    ## If ground-truth doublet classifications are available, perform ROC analysis on logistic
    ## regression model trained using pANN vector
    meta <- as.data.frame(matrix(0L, nrow=nrow(res.temp), ncol=2))
    meta[,1] <- GT.calls
    meta[,2] <- res.temp$pANN
    train.ind <- sample(1:nrow(meta), round(nrow(meta)/2), replace=FALSE)
    test.ind <- (1:nrow(meta))[-train.ind]
    colnames(meta) <- c("SinDub","pANN")
    meta$SinDub <- factor(meta$SinDub, levels = c("Doublet","Singlet"))
    model.lm <- glm(SinDub ~ pANN, family="binomial"(link='logit'), data=meta, subset=train.ind)
    prob <- predict(model.lm, newdata=meta[test.ind, ], type="response")
    ROCpred <- ROCR::prediction(predictions=prob, labels=meta$SinDub[test.ind])
    perf.auc <- ROCR::performance(ROCpred, measure="auc")
    sweep.stats$AUC[i] <- perf.auc@y.values[[1]]
  }
  
  return(sweep.stats)
  
}

find.pK <- function(sweep.stats) {
  
  ## Implementation for data without ground-truth doublet classifications 
  '%ni%' <- Negate('%in%')
  if ("AUC" %ni% colnames(sweep.stats) == TRUE) {
    ## Initialize data structure for results storage
    bc.mvn <- as.data.frame(matrix(0L, nrow=length(unique(sweep.stats$pK)), ncol=5))
    colnames(bc.mvn) <- c("ParamID","pK","MeanBC","VarBC","BCmetric")
    bc.mvn$pK <- unique(sweep.stats$pK)
    bc.mvn$ParamID <- 1:nrow(bc.mvn)
    
    ## Compute bimodality coefficient mean, variance, and BCmvn across pN-pK sweep results
    x <- 0
    for (i in unique(bc.mvn$pK)) {
      x <- x + 1
      ind <- which(sweep.stats$pK == i)
      bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
      bc.mvn$VarBC[x] <- sd(sweep.stats[ind, "BCreal"])^2
      bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"])/(sd(sweep.stats[ind, "BCreal"])^2)
    }
    
    ## Plot for visual validation of BCmvn distribution
    par(mar=rep(1,4))
    x <- plot(x=bc.mvn$ParamID, y=bc.mvn$BCmetric, pch=16, col="#41b6c4", cex=0.75)
    x <- lines(x=bc.mvn$ParamID, y=bc.mvn$BCmetric, col="#41b6c4")
    print(x)
    
    return(bc.mvn)
    
  }
  
  ## Implementation for data with ground-truth doublet classifications (e.g., MULTI-seq, CellHashing, Demuxlet, etc.)
  if ("AUC" %in% colnames(sweep.stats) == TRUE) {
    ## Initialize data structure for results storage
    bc.mvn <- as.data.frame(matrix(0L, nrow=length(unique(sweep.stats$pK)), ncol=6))
    colnames(bc.mvn) <- c("ParamID","pK","MeanAUC","MeanBC","VarBC","BCmetric")
    bc.mvn$pK <- unique(sweep.stats$pK)
    bc.mvn$ParamID <- 1:nrow(bc.mvn)
    
    ## Compute bimodality coefficient mean, variance, and BCmvn across pN-pK sweep results
    x <- 0
    for (i in unique(bc.mvn$pK)) {
      x <- x + 1
      ind <- which(sweep.stats$pK == i)
      bc.mvn$MeanAUC[x] <- mean(sweep.stats[ind, "AUC"])
      bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
      bc.mvn$VarBC[x] <- sd(sweep.stats[ind, "BCreal"])^2
      bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"])/(sd(sweep.stats[ind, "BCreal"])^2)
    }
    
    ## Plot for visual validation of BCmvn distribution
    par(mar=rep(1,4))
    x <- plot(x=bc.mvn$ParamID, y=bc.mvn$MeanAUC, pch=18, col="black", cex=0.75,xlab=NA, ylab = NA)
    x <- lines(x=bc.mvn$ParamID, y=bc.mvn$MeanAUC, col="black", lty=2)
    par(new=TRUE)
    x <- plot(x=bc.mvn$ParamID, y=bc.mvn$BCmetric, pch=16, col="#41b6c4", cex=0.75)
    axis(side=4)
    x <- lines(x=bc.mvn$ParamID, y=bc.mvn$BCmetric, col="#41b6c4")
    print(x)
    
    return(bc.mvn)
    
  }
}

doubletFinder_v3 <- function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, sct = FALSE, annotations = NULL) {
  require(Seurat); require(fields); require(KernSmooth)
  
  ## Generate new list of doublet classificatons from existing pANN vector to save time
  if (reuse.pANN != FALSE ) {
    pANN.old <- seu@meta.data[ , reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    return(seu)
  }
  
  if (reuse.pANN == FALSE) {
    ## Make merged real-artifical data
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA$counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    ##print(paste("Creating",n_doublets,"artificial doublets...",sep=" "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    # Keep track of the types of the simulated doublets
    if(!is.null(annotations)){
      stopifnot(typeof(annotations)=="character")
      stopifnot(length(annotations)==length(Cells(seu)))
      stopifnot(!any(is.na(annotations)))
      annotations <- factor(annotations)
      names(annotations) <- Cells(seu)
      doublet_types1 <- annotations[real.cells1]
      doublet_types2 <- annotations[real.cells2]
    }
    ## Store important pre-processing information
    orig.commands <- seu@commands
    
    ## Pre-process Seurat object
    if (sct == FALSE) {
      ##print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      
      ##print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets,
                                     normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
                                     scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
                                     margin = orig.commands$NormalizeData.RNA@params$margin)
      
      ##print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(seu_wdoublets,
                                            selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
                                            loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
                                            clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
                                            mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
                                            dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
                                            num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
                                            binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
                                            nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
                                            mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
                                            dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
      
      ##print("Scaling data...")
      seu_wdoublets <- ScaleData(seu_wdoublets,
                                 features = orig.commands$ScaleData.RNA$features,
                                 model.use = orig.commands$ScaleData.RNA$model.use,
                                 do.scale = orig.commands$ScaleData.RNA$do.scale,
                                 do.center = orig.commands$ScaleData.RNA$do.center,
                                 scale.max = orig.commands$ScaleData.RNA$scale.max,
                                 block.size = orig.commands$ScaleData.RNA$block.size,
                                 min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
      
      ##print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets,
                              features = orig.commands$ScaleData.RNA$features,
                              npcs = length(PCs),
                              rev.pca =  orig.commands$RunPCA.RNA$rev.pca,
                              weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
                              verbose=FALSE)
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc() # Free up memory
    }
    
    if (sct == TRUE) {
      require(sctransform)
      ##print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      
      ##print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)
      
      ##print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc()
    }
    
    ## Compute PC distance matrix
    ##print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)
    
    ## Compute pANN
    ##print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    if(!is.null(annotations)){
      neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = length(levels(doublet_types1))))
    }
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
      if(!is.null(annotations)){
        for(ct in unique(annotations)){
          neighbor_types[i,] <- 
            table( doublet_types1[neighbors - n_real.cells] ) +
            table( doublet_types2[neighbors - n_real.cells] )
          neighbor_types[i,] <- neighbor_types[i,] / sum( neighbor_types[i,] )
        }
      }
    }
    
    ##print("Classifying doublets..")
    classifications <- rep("Singlet",n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN",pN,pK,nExp,sep="_")] <- pANN[rownames(seu@meta.data), 1]
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    if(!is.null(annotations)){
      colnames(neighbor_types) = levels(doublet_types1)
      for(ct in levels(doublet_types1)){
        seu@meta.data[, paste("DF.doublet.contributors",pN,pK,nExp,ct,sep="_")] <- neighbor_types[,ct] 
      }
    }
    return(seu)
  }
}
parallel_paramSweep_v3 <- function(n, n.real.cells, real.cells, pK, pN, data, orig.commands, PCs, sct)  {

  sweep.res.list = list()
  list.ind = 0

  ## Make merged real-artifical data
  ##print(paste("Creating artificial doublets for pN = ", pN[n]*100,"%",sep=""))
  n_doublets <- round(n.real.cells/(1 - pN[n]) - n.real.cells)
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  doublets <- (data[, real.cells1] + data[, real.cells2])/2
  colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
  data_wdoublets <- cbind(data, doublets)

  ## Pre-process Seurat object
  if (sct == FALSE) {
    ##print("Creating Seurat object...")
    seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)

    ##print("Normalizing Seurat object...")
    seu_wdoublets <- NormalizeData(seu_wdoublets,
                                   normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
                                   scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
                                   margin = orig.commands$NormalizeData.RNA@params$margin)

    ##print("Finding variable genes...")
    seu_wdoublets <- FindVariableFeatures(seu_wdoublets,
                                          selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
                                          loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
                                          clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
                                          mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
                                          dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
                                          num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
                                          binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
                                          nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
                                          mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
                                          dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)

    ##print("Scaling data...")
    seu_wdoublets <- ScaleData(seu_wdoublets,
                               features = orig.commands$ScaleData.RNA$features,
                               model.use = orig.commands$ScaleData.RNA$model.use,
                               do.scale = orig.commands$ScaleData.RNA$do.scale,
                               do.center = orig.commands$ScaleData.RNA$do.center,
                               scale.max = orig.commands$ScaleData.RNA$scale.max,
                               block.size = orig.commands$ScaleData.RNA$block.size,
                               min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)

    ##print("Running PCA...")
    seu_wdoublets <- RunPCA(seu_wdoublets,
                            features = orig.commands$ScaleData.RNA$features,
                            npcs = length(PCs),
                            rev.pca =  orig.commands$RunPCA.RNA$rev.pca,
                            weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
                            verbose=FALSE)
  }

  if (sct == TRUE) {
    require(sctransform)
    ##print("Creating Seurat object...")
    seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)

    ##print("Running SCTransform...")
    seu_wdoublets <- SCTransform(seu_wdoublets)

    ##print("Running PCA...")
    seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
  }

  ## Compute PC distance matrix
  ##print("Calculating PC distance matrix...")
  nCells <- nrow(seu_wdoublets@meta.data)
  pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
  rm(seu_wdoublets)
  gc()
  dist.mat <- fields::rdist(pca.coord)[,1:n.real.cells]

  ## Pre-order PC distance matrix prior to iterating across pK for pANN computations
  ##print("Defining neighborhoods...")
  for (i in 1:n.real.cells) {
    dist.mat[,i] <- order(dist.mat[,i])
  }

  ## Trim PC distance matrix for faster manipulations
  ind <- round(nCells * max(pK))+5
  dist.mat <- dist.mat[1:ind, ]

  ## Compute pANN across pK sweep
  ##print("Computing pANN across all pK...")
  for (k in 1:length(pK)) {
    ##print(paste("pK = ", pK[k], "...", sep = ""))
    pk.temp <- round(nCells * pK[k])
    pANN <- as.data.frame(matrix(0L, nrow = n.real.cells, ncol = 1))
    colnames(pANN) <- "pANN"
    rownames(pANN) <- real.cells
    list.ind <- list.ind + 1

    for (i in 1:n.real.cells) {
      neighbors <- dist.mat[2:(pk.temp + 1),i]
      pANN$pANN[i] <- length(which(neighbors > n.real.cells))/pk.temp
    }

    sweep.res.list[[list.ind]] <- pANN

  }

  return(sweep.res.list)
}
bimodality_coefficient <- function(x) {
  G <- skewness(x)
  sample.excess.kurtosis <- kurtosis(x)
  K <- sample.excess.kurtosis
  n <- length(x)
  B <- ((G^2)+1)/(K+((3*((n-1)^2))/((n-2)*(n-3))))
  return(B)
}
doubletFinder <- function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE) {
  require(Seurat); require(fields); require(KernSmooth)

  ## Generate new list of doublet classificatons from existing pANN vector to save time
  if (reuse.pANN != FALSE ) {
    pANN.old <- seu@meta.data[ , reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[,ncol(seu@meta.data)+1] <- classifications
    colnames(seu@meta.data)[ncol(seu@meta.data)] <- paste("DF.classifications",pN,pK,nExp,sep="_")
    return(seu)
  }

  if (reuse.pANN == FALSE) {
    ## Make merged real-artifical data
    real.cells <- colnames(seu@data)
    data <- seu@raw.data[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    ##print(paste("Creating",n_doublets,"artificial doublets...",sep=" "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)

    ## Pre-process Seurat object
    ##print("Creating Seurat object...")
    seu_wdoublets <- CreateSeuratObject(raw.data = data_wdoublets)

    ##print("Normalizing Seurat object...")
    seu_wdoublets <- NormalizeData(seu_wdoublets,
                                   normalization.method = seu@calc.params$NormalizeData$normalization.method,
                                   scale.factor = seu@calc.params$NormalizeData$scale.factor)

    ##print("Finding variable genes...")
    seu_wdoublets <- FindVariableGenes(seu_wdoublets,
                                       do.plot = FALSE, x.low.cutoff = seu@calc.params$FindVariableGenes$x.low.cutoff,
                                       x.high.cutoff = seu@calc.params$FindVariableGenes$x.high.cutoff,
                                       y.high.cutoff = seu@calc.params$FindVariableGenes$y.high.cutoff,
                                       y.cutoff = seu@calc.params$FindVariableGenes$y.cutoff)

    ##print("Scaling data...")
    seu_wdoublets <- ScaleData(seu_wdoublets, display.progress = TRUE)

    ##print("Running PCA...")
    seu_wdoublets <- RunPCA(seu_wdoublets,
                            pc.genes = seu_wdoublets@var.genes,
                            pcs.print = 0,
                            pcs.compute = length(PCs))
    cell.names <- seu_wdoublets@cell.names
    nCells <- length(cell.names)

    ## Compute PC distance matrix
    #print("Calculating PC distance matrix")
    pca.coord <- seu_wdoublets@dr$pca@cell.embeddings[, PCs]
    rm(seu_wdoublets) ## Free-up memory
    gc() ## Free-up memory
    dist.mat <- fields::rdist(pca.coord)
    colnames(dist.mat) <- cell.names
    rownames(dist.mat) <- cell.names

    ## Compute pANN
    ##print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      neighbor.names <- rownames(dist.mat)[neighbors]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
    }

    ##print("Classifying doublets..")
    classifications <- rep("Singlet",n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[,ncol(seu@meta.data)+1] <- pANN[seu@cell.names, 1]
    colnames(seu@meta.data)[ncol(seu@meta.data)] <- paste("pANN",pN,pK,nExp,sep="_")
    seu@meta.data[,ncol(seu@meta.data)+1] <- classifications
    colnames(seu@meta.data)[ncol(seu@meta.data)] <- paste("DF.classifications",pN,pK,nExp,sep="_")
    return(seu)
  }
}
kurtosis <- function (x) {
  n <- length(x)
  K <- (1/n)*sum((x-mean(x))^4)/(((1/n)*sum((x-mean(x))^2))^2)-3
  K <- ((n - 1)*((n+1)*K-3*(n-1))/((n-2)*(n-3)))+3
  return(K)
}
modelHomotypic <- function(annotations) {
  anno.freq <- table(annotations)/length(annotations)
  x <- sum(anno.freq^2)
  return(x)
}
paramSweep <- function(seu, PCs) {
  require(Seurat); require(fields)
  ## Set pN-pK param sweep ranges
  pK <- c(0.0005, 0.001, 0.005, seq(0.01,0.3,by=0.01))
  pN <- seq(0.05,0.3,by=0.05)

  sweep.res.list <- list()
  list.ind <- 0

  ## Down-sample cells to 10000 (when applicable) for computational effiency
  if (ncol(seu@data) > 10000) {
    real.cells <- colnames(seu@data)[sample(1:ncol(seu@data), 10000, replace=FALSE)]
    data <- seu@raw.data[ , real.cells]
    n.real.cells <- ncol(data)
  }

  if (ncol(seu@data) <= 10000){
    real.cells <- colnames(seu@data)
    data <- seu@raw.data[, real.cells]
    n.real.cells <- ncol(data)
  }

  ## Iterate through pN, computing pANN vectors at varying pK
  for (n in 1:length(pN)) {

    ## Make merged real-artifical data
    ##print(paste("Creating artificial doublets for pN = ", pN[n]*100,"%",sep=""))
    n_doublets <- round(n.real.cells/(1 - pN[n]) - n.real.cells)
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)

    ## Pre-process Seurat object
    ##print("Creating Seurat object...")
    seu_wdoublets <- CreateSeuratObject(raw.data = data_wdoublets)

    ##print("Normalizing Seurat object...")
    seu_wdoublets <- NormalizeData(seu_wdoublets,
                                   normalization.method = seu@calc.params$NormalizeData$normalization.method,
                                   scale.factor = seu@calc.params$NormalizeData$scale.factor)

    ##print("Finding variable genes...")
    seu_wdoublets <- FindVariableGenes(seu_wdoublets,
                                       do.plot = FALSE,
                                       x.low.cutoff = seu@calc.params$FindVariableGenes$x.low.cutoff,
                                       x.high.cutoff = seu@calc.params$FindVariableGenes$x.high.cutoff,
                                       y.high.cutoff = seu@calc.params$FindVariableGenes$y.high.cutoff,
                                       y.cutoff = seu@calc.params$FindVariableGenes$y.cutoff)

    ##print("Scaling data...")
    seu_wdoublets <- ScaleData(seu_wdoublets, display.progress = TRUE)

    ##print("Running PCA...")
    seu_wdoublets <- RunPCA(seu_wdoublets,
                            pc.genes = seu_wdoublets@var.genes,
                            pcs.print = 0,
                            pcs.compute = length(PCs))

    ## Compute PC distance matrix
    ##print("Calculating PC distance matrix...")
    nCells <- length(seu_wdoublets@cell.names)
    pca.coord <- seu_wdoublets@dr$pca@cell.embeddings[, PCs]
    rm(seu_wdoublets)
    gc()
    dist.mat <- fields::rdist(pca.coord)[,1:n.real.cells]

    ## Pre-order PC distance matrix prior to iterating across pK for pANN computations
    ##print("Defining neighborhoods...")
    for (i in 1:n.real.cells) {
      dist.mat[,i] <- order(dist.mat[,i])
    }

    ## Trim PC distance matrix for faster manipulations
    ind <- round(nCells * max(pK))+5
    dist.mat <- dist.mat[1:ind, ]

    ## Compute pANN across pK sweep
    ##print("Computing pANN across all pK...")
    for (k in 1:length(pK)) {
      ##print(paste("pK = ", pK[k], "...", sep = ""))
      pk.temp <- round(nCells * pK[k])
      pANN <- as.data.frame(matrix(0L, nrow = n.real.cells, ncol = 1))
      colnames(pANN) <- "pANN"
      rownames(pANN) <- real.cells
      list.ind <- list.ind + 1

      for (i in 1:n.real.cells) {
        neighbors <- dist.mat[2:(pk.temp + 1),i]
        pANN$pANN[i] <- length(which(neighbors > n.real.cells))/pk.temp
      }

      sweep.res.list[[list.ind]] <- pANN

    }
  }

  ## Assign names to list of results
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_" ))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)

}

skewness <- function(x) {
  n <- length(x)
  S <- (1/n)*sum((x-mean(x))^3)/(((1/n)*sum((x-mean(x))^2))^1.5)
  S <- S*(sqrt(n*(n-1)))/(n-2)
  return(S)
}