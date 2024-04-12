PKGloading <- function(){
  require(Seurat)
  require(SeuratObject)
  require(SeuratDisk)
  require(sctransform)
  require(rhdf5)
  require(Matrix)
  #require(future) #no effect on SCT
  #plan("multiprocess", workers = 8)
  require(peakRAM)
  options(future.globals.maxSize=8000*1024^2)
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
    message("\t\t--> SCT normalization is selected <--")
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
        #expScale <- D@assays$SCT@SCTModel.list[[1]]@median_umi
      }
    }
  }else{
    message("\t\t--> LogNormal normalization is selected <--")
    message("\tScale: ",expScale)
    D <- NormalizeData(D,assay="RNA",normalization.method = "LogNormalize", scale.factor = expScale,verbose = FALSE)
  }
  saveX(D,strOut,gID)
}
saveX <- function(D,strH5,ggID){
  message("\tsaving expression: ",strH5)
  saveRDS(D,paste0(strH5,".rds"))
  if("SCT"%in%names(D)){
    X <- Matrix::t(Matrix::Matrix(D@assays$SCT@layers$data,sparse=T))
    cID <- D@assays$SCT@cells[['data']]
    gID <- D@assays$SCT@features[['data']]
  }else{
    X <- Matrix::t(Matrix::Matrix(D@assays$RNA@layers$data,sparse=T))
    cID <- D@assays$RNA@cells[['data']]
    gID <- D@assays$RNA@features[['data']]
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
  cat(paste(cID,collapse="\n"),sep="",file=gsub("h5$","cID",strH5))
  cat(paste(ggID[gID],collapse="\n"),sep="",file=gsub("h5$","gID",strH5))
  #cat(expScale,file=gsub("h5$","scaleF",strH5))
}
mergeAllbatches <- function(strRDS,strOut,batchKey){
  D <- NULL
  message("\tmerge all normalized seruat objects")
  for(one in strRDS){
    oneD <- readRDS(one)
    message("\t\t",oneD@meta.data[1,batchKey],"\t",which(strRDS==one),"/",length(strRDS))
    if(is.null(D)) D <- oneD
    else D <- merge(D,oneD)
  }
  if("SCT" %in% names(D)){
    message("\tPrepSCTFindMarkers might take a while")
    D <- PrepSCTFindMarkers(D)
    cat(D@assays$SCT@SCTModel.list[[1]]@median_umi,file=gsub("rds$","scaleF",strOut))
  }
  saveRDS(D,strOut)
}
main <- function(){
  suppressMessages(suppressWarnings(PKGloading()))
  batchKey="library_id" #"batch"
  args = commandArgs(trailingOnly=TRUE)
  if(length(args)<2) stop("Path to h5ad file, the output file and config file are required!")
  if(length(args)==2){
    strRDS <- paste0(unlist(strsplit(args[1],",")),".rds")
    print(peakRAM(mergeAllbatches(strRDS,args[2],batchKey)))
  }else{
    strH5ad <- args[1]
    if(!file.exists(strH5ad)) stop(paste0("H5ad file (",strH5ad,") does not exist!"))
    strOut <- args[2]
    strConfig <- args[3]
    if(!file.exists(strConfig)){
      stop(paste0("Config file (",strConfig,") does not exist!"))
    }else{
      config <- yaml::read_yaml(strConfig)
    }
    scaleF=as.numeric(args[4])
    if(length(args)>4) batchKey <- args[5]
    
    source(paste0(dirname(gsub("--file=","",grep("file=",commandArgs(),value=T))),"/readH5ad.R"))
    print(peakRAM(processH5ad(strH5ad,batchKey,strOut,scaleF,config$PrepSCTFindMarkers)))
  }
}

main()