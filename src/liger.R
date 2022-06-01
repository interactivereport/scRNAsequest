PKGloading <- function(){
  require(Seurat)
  require(liger)
  require(rhdf5)
  require(Matrix)
  options(stringsAsFactors = FALSE)
}

getobs <- function(strH5ad){
  message("\tobtainning obs ...")
  obs <- h5read(strH5ad,"obs")
  meta <- do.call(cbind.data.frame, obs[grep("^_",names(obs),invert=T)])
  dimnames(meta) <- list(obs[["_index"]],grep("^_",names(obs),invert=T,value=T))
  for(one in names(obs[["__categories"]])){
    meta[,one] <- obs[["__categories"]][[one]][1+meta[,one]]
  }
  return(meta)
}
getX <- function(strH5ad){
  message("\tobtainning X ...")
  X <- h5read(strH5ad,"X")
  gID <- h5read(strH5ad,"var/_index")
  cID <- h5read(strH5ad,"obs/_index")
  if((max(X$indices)+1)==length(gID)){ # CSR sparse matrix
    M <- sparseMatrix(i=X$indices+1,p=X$indptr,x=as.numeric(X$data),
                      dims=c(length(gID),length(cID)),
                      dimnames=list(gID,cID))
  }else if((max(X1$indices)+1)==length(cID)){#CSC sparse matrix
    M <- sparseMatrix(j=X$indices+1,p=X$indptr,x=as.numeric(X$data),
                      dims=c(length(gID),length(cID)),
                      dimnames=list(gID,cID))
  }
  return(M)
}

runLiger <- function(strH5ad,strHVG,batchKey,strOut){
  X <- getX(strH5ad)
  meta <- getobs(strH5ad)
  Xlist <- list()
  for(one in unique(meta[,batchKey])) Xlist[[one]] <- X[,(1:nrow(meta))[meta[,batchKey]==one]]
  ## create liger obj and norm
  message("\tcreating liger object")
  adatalg <- createLiger(raw.data = Xlist, take.gene.union = FALSE, remove.missing = FALSE, make.sparse = TRUE)
  adatalg <- normalize(adatalg)
  ## get hvg from scanpy (all data har workflows should have the same hvg)
  hvg <- data.table::fread(strHVG)
  adatalg@var.genes <- unlist(hvg[,1],use.names=F)[unlist(hvg[,2])]
  ## scale
  adatalg <- scaleNotCenter(adatalg)
  ### Factorization, first step is to determine k
  ## running suggestK on multiple cores can greatly decrease the runtime
  #k.suggest <- suggestK(adatalg, num.cores = 50, gen.new = T, plot.log2 = T, nrep = 3)
  ### this takes a long time
  #
  #k.suggest
  #
  message("\toptimizing ALS ...")
  adatalg <- optimizeALS(adatalg, k=50)
  adatalg <- quantileAlignSNF(adatalg)
  # using H.norm to simulate PCA
  layout <- setNames(as.data.frame(adatalg@H.norm),paste0("pca_",1:ncol(adatalg@H.norm)))
  
  message("\trunning UMAP ...")
  #adatalg <- runUMAP(adatalg) # this create a tmp file in the conda env root folder, so permission denied
  UMAP <- uwot::umap(adatalg@H.norm,
                     n_components = 2, metric = "euclidean",
                     n_neighbors = 10, min_dist = 0.1)
  rownames(UMAP) <- rownames(adatalg@H.norm)
  message(tempdir())
  #saveRDS(adatalg, file = "working_data/concat_liger_runumap.Rdata")
  layout <- cbind(layout,
                  setNames(as.data.frame(UMAP),paste0("umap_",1:ncol(UMAP))))

  message("\trunning tSNE ...")
  adatalg <- runTSNE(adatalg)
  layout <- cbind(layout,
                  setNames(as.data.frame(adatalg@tsne.coords),paste0("tsne_",1:ncol(adatalg@tsne.coords))))
  
  layout <- cbind(layout,liger_cluster=adatalg@clusters)
  layout <- cbind(cID=rownames(layout),layout)
  data.table::fwrite(layout,strOut)
}
main <- function(){
  suppressMessages(suppressWarnings(PKGloading()))
  batchKey <- "library_id" #"batch"#
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)<3) stop("Path to 3 files are required!")
  strH5ad <- args[1]
  strHVG <- args[2]
  if(!file.exists(strH5ad) || !file.exists(strHVG)) stop("either raw h5ad or HVG file is missing!")
  strOut <- args[3]
  runLiger(strH5ad,strHVG,batchKey,strOut)
}


main()