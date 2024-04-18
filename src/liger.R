PKGloading <- function(){
  #require(Seurat)
  require(rliger)
  require(dplyr)
  require(stringr)
  options(stringsAsFactors = FALSE)
}
#source("../../src/readH5ad.R")
#strH5ad <- "raw/TST11837_oyoung_single_raw_postfilter.h5ad"
#strHVG <- "Liger/TST11837_oyoung_single_hvg.csv"
#batchKey <- "library_id"
#clusterMethod <- "Leiden"
#clusterResolution <- "0.8"

runLiger <- function(strH5ad,strHVG,batchKey,
                     clusterMethod,clusterResolution,
                     strOut){
  X <- getX(strH5ad)
  meta <- getobs(strH5ad)
  Xlist <- list()
  for(one in unique(meta[,batchKey])) Xlist[[one]] <- X[,(1:nrow(meta))[meta[,batchKey]==one]]
  ## create liger obj and norm
  message("\tcreating liger object")
  adatalg <- suppressMessages(createLiger(rawData = Xlist,removeMissing=FALSE,addPrefix=F))
  # normalization and scale
  message("\tnormalizing & scaling")
  adatalg <- suppressMessages(normalize(adatalg))
  ## get hvg from scanpy (all data har workflows should have the same hvg)
  hvg <- data.table::fread(strHVG)
  varFeatures(adatalg) <- unlist(hvg[,1],use.names=F)[unlist(hvg[,2])]
  adatalg <- suppressMessages(scaleNotCenter(adatalg))

  message("\tIntegration with Joint Matrix Factorization  ...")
  adatalg <- suppressMessages(runIntegration(adatalg, k=50))
  adatalg <- suppressMessages(quantileNorm(adatalg))
  
  # clustering
  message("\tclustering: ",clusterMethod,"(",clusterResolution,")")
  adatalg <- suppressMessages(runCluster(adatalg,resolution=as.numeric(clusterResolution),
                        method=tolower(clusterMethod)))
  # UMAP
  message("\tUMAP")
  adatalg <- suppressMessages(runUMAP(adatalg,nNeighbors=30,minDist=0.3))
  
  # save Cluster, embedding (PCA/UMAP)
  message("\tsaving")
  X <- dplyr::bind_cols(data.frame(adatalg@cellMeta[,paste0(tolower(clusterMethod),"_cluster"),drop=F]) %>% 
                          dplyr::rename_at(1,~paste0("Liger_",.)),
                        data.frame(adatalg@H.norm) %>% 
                          dplyr::rename_at(vars(contains("Factor_")),list(~str_replace(.,"Factor_","pca_"))),
                        setNames(data.frame(adatalg@dimReds$UMAP),tolower(colnames(adatalg@dimReds$UMAP))))
  saveRDS(X,strOut)
}
main <- function(){
  suppressMessages(suppressWarnings(PKGloading()))
  batchKey <- "library_id" #"batch"#
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)<5) stop("Path to 3 files and two cluster parameters are required!")
  strH5ad <- args[1]
  strHVG <- args[2]
  if(!file.exists(strH5ad) || !file.exists(strHVG)) stop("either raw h5ad or HVG file is missing!")
  strOut <- args[3]
  clusterMethod <- args[4]
  clusterResolution <- args[5]

  source(paste0(dirname(gsub("--file=","",grep("file=",commandArgs(),value=T))),"/readH5ad.R"))
  runLiger(strH5ad,strHVG,batchKey,clusterMethod,clusterResolution,strOut)
}


main()