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
  require(BiocParallel)
  options(future.globals.maxSize=8000*1024^2)
}

processH5ad <- function(strH5ad,batch,strOut,expScale,bPrepSCT,core=10){
  if(is.null(bPrepSCT)) bPrepSCT <- F
  X <- getX(strH5ad,batchID=batch,core=core)
  gID <- setNames(rownames(X[[1]]),gsub("_","-",rownames(X[[1]])))
  X <- lapply(X,function(one){
    rownames(one) <- names(gID)
    return(one)
  })
  #rownames(X) <- names(gID) #gsub("_","-",rownames(X))
  D <- CreateSeuratObject(counts=X,
                          project="SCT",
                          meta.data=getobs(strH5ad))
  rm(X)
  gc()
  if(is.null(expScale) || expScale==0){
    message("\t\t--> SCT normalization is selected <--")
    Dmedian <- NA
    if(!bPrepSCT){
      message("\t***median UMI is used to normalize across samples by scale_factor from vst***")
      expScale <- Dmedian <- median(colSums(D@assays$RNA@counts))#min(colSums(D@assays$RNA@counts))#
    }
    D <- SCTransform(D,vst.flavor="v2",
                     return.only.var.genes = FALSE,
                     scale_factor=medianUMI)
    if(!is.null(bPrepSCT) && bPrepSCT){
      message("\t***PrepSCTFindMarkers***\n\t\tMight take a while ...")
      D <- PrepSCTFindMarkers(D)
    }
    D[["SCT"]] <- split(D[["SCT"]],f=unlist(DDnorm[[batch]],use.names=F))
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
    D[["SCT"]] <- JoinLayers(D[["SCT"]])
    X <- Matrix::t(Matrix::Matrix(D@assays$SCT@layers$data,sparse=T))
    cID <- D@assays$SCT@cells[['data']]
    gID <- D@assays$SCT@features[['data']]
  }else{
    D[["RNA"]] <- JoinLayers(D[["RNA"]])
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
mergeAllbatches <- function(strRDS,strOut,core){
  D <- NULL
  message("\tmerge all normalized seruat objects")
  Dlist <- bplapply(strRDS,function(one){
    message("\t\tReading ",basename(one),"\t",which(strRDS==one),"/",length(strRDS)," @",Sys.time())
    return(readRDS(one))
  },
  BPPARAM = MulticoreParam(workers=min(core,length(strRDS),max(1,parallelly::availableCores()-2))))
  # tasks in MulticoreParam seems exectute as batches, each batch every tasks needs to be all completed until the next batch
  gc()
  message("\t\tMerging ...@",Sys.time())
  D <- merge(Dlist[[1]],Dlist[-1])
  rm(Dlist)
  gc()
  message("\t\tFinished @",Sys.time())
  if("SCT" %in% names(D)){
    message("\tPrepSCTFindMarkers might take a while")
    D <- PrepSCTFindMarkers(D)
    cat(D@assays$SCT@SCTModel.list[[1]]@median_umi,file=gsub("rds$","scaleF",strOut))
  }
  message("\tSaving ...")
  saveRDS(D,strOut)
}
main <- function(){
  suppressMessages(suppressWarnings(PKGloading()))
  batchKey="library_id" #"batch"
  args = commandArgs(trailingOnly=TRUE)
  #if(length(args)<2) stop("Path to h5ad file, the output file and config file are required!")
  task <- args[1]
  if(task=="MERGE"){
    strRDS <- paste0(unlist(strsplit(args[2],",")),".rds")
    print(peakRAM(mergeAllbatches(strRDS,args[3],as.numeric(args[4]))))
  }else if(task=="NORM"){
    strH5ad <- args[2]
    if(!file.exists(strH5ad)) stop(paste0("H5ad file (",strH5ad,") does not exist!"))
    strOut <- args[3]
    strConfig <- args[4]
    if(!file.exists(strConfig)){
      stop(paste0("Config file (",strConfig,") does not exist!"))
    }else{
      config <- yaml::read_yaml(strConfig)
    }
    scaleF=as.numeric(args[5])
    if(length(args)>5) batchKey <- args[6]
    
    source(paste0(dirname(gsub("--file=","",grep("file=",commandArgs(),value=T))),"/readH5ad.R"))
    print(peakRAM(processH5ad(strH5ad,batchKey,strOut,scaleF,config$PrepSCTFindMarkers,
                              ifelse(is.null(config$subprocess),5,config$subprocess))))
  }else{
    stop("Unknown task in SCT.R: ",task)
  }
}

main()