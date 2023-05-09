PKGloading <- function(){
  require(Seurat)
  require(SeuratObject)
  require(SeuratDisk)
  require(sctransform)
  require(rhdf5)
  require(Matrix)
  #require(future) #no effect on SCT
  #plan("multiprocess", workers = 8)
  options(future.globals.maxsize=3145728000)
}

processH5ad <- function(strH5ad,batch,strOut,expScale,bPrepSCT){
  if(is.null(bPrepSCT)) bPrepSCT <- F
  X <- getX(strH5ad)
  gID <- setNames(rownames(X),gsub("_","-",rownames(X)))
  rownames(X) <- names(gID) #gsub("_","-",rownames(X))
  D <- CreateSeuratObject(counts=X,
                          project="SCT",
                          meta.data=getobs(strH5ad))
  if(is.null(expScale) || expScale==0){
    message("\n\n===== SCT normalization is selected =====")
    Dmedian <- NA
    if(!bPrepSCT){
      message("\t***median UMI is used to normalize across samples by scale_factor from vst***")
      expScale <- Dmedian <- median(colSums(D@assays$RNA@counts))#min(colSums(D@assays$RNA@counts))#
    }
    Dlist <- SplitObject(D,split.by=batch)
    Dlist <- sapply(Dlist,function(one,medianUMI){
      message("\t\tSCT ",one@meta.data[1,batch]," ...")
      #after checking/testing, the above would return proper SCT nomralized data
      #/edgehpc/dept/compbio/users/zouyang/process/PRJNA544731/src/SCT_scale_batch.ipynb
      #https://github.com/satijalab/sctransform/issues/128
      return(suppressMessages(suppressWarnings(
        SCTransform(one,method = 'glmGamPoi',
                    new.assay.name="SCT",
                    return.only.var.genes = FALSE,
                    scale_factor=medianUMI,
                    verbose = FALSE)
      )))
    },Dmedian)
    if(length(Dlist)==1){
      D <- Dlist[[1]]
    }else{
      #D <- merge(Dlist[[1]], y=Dlist[-1])
      # according to the test the blow save almost half of memory comparing the above
      D <- Dlist[[1]]
      Dlist[[1]] <- NULL
      for(i in 1:length(Dlist)){
        D <- merge(D,Dlist[[1]])
        Dlist[[1]] <- NULL
      }
      if(!is.null(bPrepSCT) && bPrepSCT){
        message("\t***PrepSCTFindMarkers***\n\t\tMight take a while ...")
        D <- PrepSCTFindMarkers(D)
        expScale <- D@assays$SCT@SCTModel.list[[1]]@median_umi
      }
    }
  }else{
    message("\n\n===== LogNormal normalization is selected =====")
    if(expScale<=100)
      expScale <- round(quantile(D@meta.data$nCount_RNA,expScale/100)/1e3)*1e3
    message("\tScale: ",expScale)
    D <- NormalizeData(D,assay="RNA",normalization.method = "LogNormalize", scale.factor = expScale)
  }
  saveX(D,strOut,gID,expScale)
}
saveX <- function(D,strH5,gID,expScale){
  message("\tsaving expression: ",strH5)
  saveRDS(D,paste0(strH5,".rds"))
  if("SCT"%in%names(D)){
    X <- Matrix::t(Matrix::Matrix(D@assays$SCT@data,sparse=T))
  }else{
    X <- Matrix::t(Matrix::Matrix(D@assays$RNA@data,sparse=T))
  }
  suppressMessages(suppressWarnings({
    a <- file.remove(strH5)
    h5write(X@x,file=strH5,name="data")
    h5write(X@i,file=strH5,name="indices")
  }))
  h5write(X@p,file=strH5,name="indptr")
  h5write(dim(X),file=strH5,name="shape")
  #h5write(rownames(X),file=strH5,name="row_names")
  #h5write(colnames(X),file=strH5,name="col_names")
  h5closeAll()
  
  #strF <- gsub("h5$","info.h5",strH5)
  #a <- file.remove(strF)
  #h5write(rownames(X),file=strF,name="cID")
  #h5write(gID[colnames(X)],file=strF,name="gID")
  #h5write(expScale,file=strF,name="scaleFactor")
  #h5closeAll()
  
  cat(paste(rownames(X),collapse="\n"),sep="",file=gsub("h5$","cID",strH5))
  cat(paste(gID[colnames(X)],collapse="\n"),sep="",file=gsub("h5$","gID",strH5))
  cat(expScale,file=gsub("h5$","scaleF",strH5))
}
main <- function(){
  suppressMessages(suppressWarnings(PKGloading()))
  batchKey="library_id" #"batch"
  args = commandArgs(trailingOnly=TRUE)
  if(length(args)<3) stop("Path to h5ad file, the output file and config file are required!")
  strH5ad <- args[1]
  if(!file.exists(strH5ad)) stop(paste0("H5ad file (",strH5ad,") does not exist!"))
  strOut <- args[2]
  strConfig <- args[3]
  if(!file.exists(strConfig)){
    stop(paste0("Config file (",strConfig,") does not exist!"))
  }else{
    config <- yaml::read_yaml(strConfig)
  }
  if(length(args)>3) batchKey <- args[4]
  
  source(paste0(dirname(gsub("--file=","",grep("file=",commandArgs(),value=T))),"/readH5ad.R"))
  print(system.time(processH5ad(strH5ad,batchKey,strOut,config$expScaler,config$PrepSCTFindMarkers)))
}

main()