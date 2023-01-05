# Adding meta and embedding from h5ad into seurat object
PKGloading <- function(){
  require(Seurat)
  require(SeuratObject)
  require(SeuratDisk)
  require(dplyr)
}

seuratObj <- function(strRDS,strH5ad){
  D <- readRDS(strRDS)
  # add additional the cell meta information
  message("Obtain and append cell meta information")
  meta <- getobs(strH5ad)
  meta <- meta[,!colnames(meta)%in%colnames(D@meta.data)]
  meta <- merge(D@meta.data,meta,by=0,all.x=T)
  rownames(meta) <- unlist(meta[,1])
  D@meta.data <- meta[,-1]
  # add embedding
  message("Obtain and append cell integration embedding")
  for(one in getobsmKey(strH5ad)){
    keyname <- gsub("X_","",one)
    message("\t",keyname)
    key <- gsub("[[:punct:]]","",keyname)
    X <- getobsm(strH5ad,one)
    colnames(X) <- paste(key,1:ncol(X),sep="_")
    D[[keyname]] <- CreateDimReducObject(embeddings=X,key=paste0(key,"_"),assay="RNA")
  }
  message("Saving h5seurat")
  suppressMessages(SeuratDisk::SaveH5Seurat(D,filename=gsub("h5ad$","h5seurat",strH5ad),overwrite=T))
}

main <- function(){
  # functions of accessing h5ad
  a <- commandArgs()
  strPath <- gsub("^--file=","",grep("^--file",a,value=T)[1])
  scRNAseq_DE_path <- dirname(normalizePath(strPath))
  source(paste0(scRNAseq_DE_path,"/readH5ad.R"))
  
  suppressMessages(suppressWarnings(PKGloading()))
  args = commandArgs(trailingOnly=TRUE)
  seuratObj(args[1],args[2])
}


main()