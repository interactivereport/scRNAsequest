#!/usr/bin/env Rscript

strPipePath <- ""
##
loadingPKG <- function(){
  require(Seurat)
  require(cowplot)
  require(patchwork)
  require(sctransform)
  require(rhdf5)
  require(Matrix)
  require(Azimuth)
  require(scales)
  options(future.globals.maxsize=3145728000,stringsAsFactors=F)
  #source(paste0(strPipePath,"/src/Azimuth_create.R"))
}
run_cmd <- function(cmd){
  cmdR <- tryCatch(system(cmd,intern=T),
                   error=function(e){
                     message(e)
                     return("")
                   })
  return(cmdR)
}

## msg
MsgExit <- function(msg=""){
  if(nchar(msg)>3) message("ERROR: ",msg)
  MsgPower()
  q()
}
MsgPower <- function(){
  message("\nPowered by the Computational Biology Group [zhengyu.ouyang@biogen.com;kejie.li@biogen.com]")
  message("------------")
}
MsgHelp <- function(){
  message("\nscRef /path/to/a/output/folder === or === scRef /path/to/a/Ref/config/file\n")
  message("The folder will be created if it does not exist.")
  message("The Ref config file will be generated automatically when a path is provided")
  message("===== CAUTION =====")
  message("\t1. This process will add a seurat reference data into the scRNAsequest pipeline PERMANENTLY!")
  message("\t2. Make sure the data provided for reference building is SCT transformed!")
  MsgPower()
  q()
}
MsgInit <- function(){
  cmdURL=paste0("cd ",strPipePath,";git config --get remote.origin.url")
  message("###########\n## scAnalyzer: ",run_cmd(cmdURL))
  message("## Pipeline Path: ",strPipePath)
  cmdHEAD=paste0("cd ",strPipePath,";git rev-parse HEAD")
  message("## git HEAD: ",run_cmd(cmdHEAD),"\n###########")
  message("\nLoading resources")
}

main <- function(){
  args = commandArgs()
  strPipePath <<- normalizePath(dirname(dirname(sapply(strsplit(grep("file=",args,value=T),"="),tail,1))))

  strInput <- commandArgs(trailingOnly=T)
  if(length(strInput)==0){
    MsgHelp()
  }
  if(file.exists(strInput[1]) && !dir.exists(strInput[1])){
    createRef(strInput[1])
  }else{
    initRef(strInput[1])
  }

}

initRef <- function(strDir){
  strDir <- normalizePath(strDir)
  dir.create(strDir,showWarnings=F)
  config <- readLines(paste0(strPipePath,"/src/ref.yml"))
  config <- gsub("OUTPUTdir",strDir,config)
  writeLines(config,paste0(strDir,"/refConfig.yml"))
  message(paste0("\n---> Please update the config file @",paste0(strDir,"/refConfig.yml")))
  MsgExit()
}

checkExist <- function(strF){
  if(is.null(strF) || !file.exists(strF))
    return(F)
  return(T)
}
checkConfig <- function(strConfig,sysConfig){
  message("Checking config file ...")
  config <- yaml::read_yaml(strConfig)
  sapply(c("ref_name","ref_link","ref_src","ref_platform","ref_PCA","ref_label"),
         function(x,sInfo){
           a <- sInfo[[x]]
           if(length(a)==0 || nchar(a[1])==0) MsgExit(paste0("Please provided required infomraiton @",x))
         },config)
  
  if(config$ref_name%in%sysConfig$ref && !config$overwrite){
    MsgExit(paste0("ref_name (",config$ref_name,") exists, please either change it or set 'overwrite: True'"))
  }
  if(!checkExist(config$ref_rds))
    if(!checkExist(config$ref_h5ad_raw))
      MsgExit(paste0("Either one seurat object rds file or raw h5ad files is required to be existed"))

  return(config)
}
getobs <- function(strH5ad){
  message("\tobtainning obs ...")
  obs <- h5read(strH5ad,"obs")
  meta <- do.call(cbind.data.frame, obs[grep("^_",names(obs),invert=T)])
  if(!"_index"%in%names(obs)) dimnames(meta) <- list(obs[[1]],grep("^_",names(obs),invert=T,value=T))
  else dimnames(meta) <- list(obs[["_index"]],grep("^_",names(obs),invert=T,value=T))
  
  if("__categories"%in%names(obs)){
    for(one in names(obs[["__categories"]])){
      meta[,one] <- obs[["__categories"]][[one]][1+meta[,one]]
    }
  }
  return(meta)
}
getX <- function(strH5ad,cID){
  message("\tobtainning X ...")
  X <- h5read(strH5ad,"X")
  gID <- h5read(strH5ad,"var/_index")
  cID <- h5read(strH5ad,"obs/_index")
  if((max(X$indices)+1)==length(gID)){ # CSR sparse matrix
    M <- sparseMatrix(i=X$indices+1,p=X$indptr,x=as.numeric(X$data),
                      dims=c(length(gID),length(cID)),
                      dimnames=list(gID,cID))
  }else if((max(X$indices)+1)==length(cID)){#CSC sparse matrix
    M <- sparseMatrix(j=X$indices+1,p=X$indptr,x=as.numeric(X$data),
                      dims=c(length(gID),length(cID)),
                      dimnames=list(gID,cID))
  }
  return(M)
}
getobsm <- function(strH5ad,key){
  k <- h5ls(strH5ad,recursive=2)
  k <- k[grepl("obsm",k[,1]),2]
  if(!key%in%k) return(NULL)
  X <- h5read(strH5ad,paste0("obsm/",key))
  return(t(X))
}
getSCT <- function(strH5ad,batch){
  D <- CreateSeuratObject(counts=getX(strH5ad),
                          project="SCT",
                          meta.data=getobs(strH5ad))
  Dmedian <- median(colSums(D@assays$RNA@counts))
  if(!batch%in%colnames(D@meta.data)) MsgExit("ref_batch (",batch,") is not one of the annotations")
  Dlist <- SplitObject(D,split.by=batch)
  rm("D")
  gc()
  message("\tSCTransform ...")
  Dlist <- sapply(Dlist,function(one,medianUMI){
    message("\t\tSCT ",unique(unlist(one[[batch]],use.names=F))," ...")
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
  # obtain one unified SCT model for Azimuth, and the SCT model is not used for mapping which cause problem of normalizing
  Dtmp <- SCTransform(SCT,method = 'glmGamPoi',new.assay.name="SCT",return.only.var.genes = FALSE,verbose = FALSE)
  SCT[["SCT"]]@SCTModel.list <- Dtmp[["SCT"]]@SCTModel.list
  VariableFeatures(SCT) <- SelectIntegrationFeatures(Dlist,nfeatures=3000)
  return(SCT)
}
saveRef <- function(D,config,sysConfig){
  strRef <- file.path(sysConfig$refDir,
                       gsub("[ [:punct:]]","_",paste("azimuth",config$ref_name,config$ref_src,config$ref_platform,"ref")),
                       "ref.Rds")
  message("ref saved @",strRef)
  dir.create(dirname(strRef),showWarnings=F)
  saveRDS(D,strRef)
  
  strReadme <- file.path(dirname(strRef),"readme")
  cat("This is created by scRef based on a project @",config$output,"\n\n",
      sep="",file=strReadme)
  conn <- file(strReadme,"a")
  sink(conn, type="message")
  MsgInit()
  sink(type="message")
  close(conn)
  # update system config
  addOne <- setNames(list(list(ref_file=strRef,
                               ref_link=config$ref_link,
                               ref_src=config$ref_src,
                               ref_platform=config$ref_platform,
                               ref_assay="refAssay",
                               ref_neighbors="refdr.annoy.neighbors",
                               ref_reduction="refDR",
                               ref_reduction.model="refUMAP",
                               ref_label=config$ref_label)),
                     config$ref_name)
  sysConfig$ref <- c(sysConfig$ref,config$ref_name)
  sysConfig <- c(sysConfig[names(sysConfig)!="methods"],
                 addOne,
                 sysConfig["methods"])
  strSys <- paste0(strPipePath,"/src/sys.yml")
  comments <- grep("^#",readLines(strSys),value=T)
  yaml::write_yaml(sysConfig,strSys)
  cat(paste(comments,collapse="\n"),"\n",
      sep="",file=strSys,append=T)
  message("\nA new reference (",config$ref_name,") is added into the scAnalyzer!")
  MsgPower()
}
createRef <- function(strConfig){
  MsgInit()
  sysConfig <- yaml::read_yaml(paste0(strPipePath,"/src/sys.yml"))
  config <- checkConfig(strConfig,sysConfig)
  suppressMessages(suppressWarnings(loadingPKG()))
  #"/camhpc/ngs/projects/TST11837/dnanexus/20220311155416_maria.zavodszky/sc20220403_0/TST11837_SCT.h5ad"
  if(!is.null(config$ref_rds) && file.exists(config$ref_rds)){
    D <- readRDS(config$ref_rds)
  }else if(!is.null(config$ref_h5ad_raw) && file.exists(config$ref_h5ad_raw)){
    D <- getSCT(config$ref_h5ad_raw,config$ref_batch)
    DefaultAssay(D) <- "SCT"
    xy <- getobsm(config$ref_h5ad,paste0("X_",config$ref_PCA))
    if(is.null(xy)) MsgExit(paste0("The 'ref_PCA' (",config$ref_PCA,") is not in the h5ad file"))
    rownames(xy) <- colnames(D)
    D[[config$ref_PCA]] <- CreateDimReducObject(embeddings=xy[colnames(D),],
                                                key="PC_",
                                                assay="SCT")
  }else{
    MsgExit(paste0("Either one seurat object rds file or two h5ad files are required to be existed"))
  }
  if(sum(!config$ref_label%in%colnames(D@meta.data)))
    MsgExit(paste0("The following annotation labels are not in the h5ad file:\n",
                   paste(config$ref_label[!config$ref_label%in%colnames(D@meta.data)],collapse=", ")))
  message("Processing ...")
  D <- FindNeighbors(D, dims = 1:30, reduction=config$ref_PCA,verbose = FALSE)
  D <- RunSPCA(D, npcs=ncol(D[[config$ref_PCA]]), graph = 'SCT_snn')
  D <- RunUMAP(D, dims = 1:30, reduction="spca",return.model=TRUE)
  saveRDS(D,file.path(config$output,"ref_full.rds"))
  D_Azimuth <- AzimuthReference(
    D,
    refUMAP = "umap",
    refDR = "spca",
    refAssay = "SCT",
    dims = 1:ncol(D[[config$ref_PCA]]),
    plotref = "umap",
    metadata = config$ref_label #c('genotype','age','sex','cluster',"celltype")
  )
  
  saveRef(D_Azimuth,config,sysConfig)
  
}
main()
