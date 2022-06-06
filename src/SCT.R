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

processH5ad <- function(strH5ad,batch,strOut){
  X <- getX(strH5ad)
  gID <- setNames(rownames(X),gsub("_","-",rownames(X)))
  rownames(X) <- gsub(rownames(X),"_","-")
  D <- CreateSeuratObject(counts=X,
                          project="SCT",
                          meta.data=getobs(strH5ad))
  Dmedian <- median(colSums(D@assays$RNA@counts))
  Dlist <- SplitObject(D,split.by=batch)
  Dlist <- sapply(Dlist,function(one,medianUMI){
    message("\t\tSCT one ...")
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
    SCT <- Dlist[[1]]
  }else{
    SCT <- merge(Dlist[[1]], y=Dlist[-1])
  }
  saveX(SCT,strOut,gID)
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
saveX <- function(SCT,strH5,gID){
  message("\tsaving SCT ...")
  saveRDS(SCT,paste0(strH5,".rds"))
  X <- Matrix::t(Matrix::Matrix(SCT@assays$SCT@data,sparse=T))
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
  cat(paste(rownames(X),collapse="\n"),sep="",file=gsub("h5$","cID",strH5))
  cat(paste(gID[colnames(X)],collapse="\n"),sep="",file=gsub("h5$","gID",strH5))
}
main <- function(){
  suppressMessages(suppressWarnings(PKGloading()))
  batchKey="library_id" #"batch"
  args = commandArgs(trailingOnly=TRUE)
  if(length(args)<2) stop("Path to h5ad file and the output file are required!")
  strH5ad <- args[1]
  if(!file.exists(strH5ad)) stop(paste0("H5ad file (",strH5ad,") does not exist!"))
  strOut <- args[2]
  if(length(args)>2) batchKey <- args[3]
  print(system.time(processH5ad(strH5ad,batchKey,strOut)))
}

main()