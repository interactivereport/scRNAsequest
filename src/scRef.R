#!/usr/bin/env Rscript

strPipePath <- ""
## loading packages -----
loadingPKG <- function(){
  require(Seurat)
  require(cowplot)
  require(patchwork)
  require(sctransform)
  require(rhdf5)
  require(Matrix)
  #require(Azimuth)
  require(scales)
  options(future.globals.maxsize=3145728000,stringsAsFactors=F)
  source(paste0(strPipePath,"/src/azimuth.R"))
}
## msg -----
MsgExit <- function(msg=""){
  if(nchar(msg)>3) message("ERROR: ",msg)
  MsgPower()
  q()
}
MsgPower <- function(){
  message("\nPowered by the Research Data Sciences group [zhengyu.ouyang@biogen.com;kejie.li@biogen.com]")
  message("------------")
}
MsgHelp <- function(){
  message("\nscRef /path/to/a/output/folder === or === scRef /path/to/a/Ref/config/file\n")
  message("The folder has to be existed.")
  message("The Ref config file will be generated automatically when a path is provided")
  message("===== CAUTION =====")
  message("\t1. This process will add a seurat reference data into the scRNAsequest pipeline PERMANENTLY!")
  message("\t2. Make sure the data provided for reference building is SCT transformed!")
  MsgPower()
  q()
}
MsgInit <- function(){
  message("\n\n***** ",Sys.time()," *****")
  if(dir.exists(file.path(strPipePath,".git"))){
    gitConfig <- readLines(file.path(strPipePath,".git","config"))
    pos <- grep("^\\[.*\\]$",gitConfig)
    sel <- (1:length(pos))[grepl("remote",gitConfig[pos])&grepl("origin",gitConfig[pos])]
    if(length(sel)==0) return()
    pos <- c(pos,length(gitConfig)+1)
    url <- sapply(strsplit(grep("^url",trimws(gitConfig[pos[sel]:(pos[sel+1]-1)]),value=T),
                           " "),tail,1)
    gitLog <- unlist(strsplit(unlist(tail(data.table::fread(file.path(strPipePath,".git","logs","HEAD"),header=F),1))[1]," "))
    message("###########\n## scRNAsequest: ",url)
    message("## Pipeline Path: ",strPipePath)
    message("## Pipeline Date: ",
            format(as.POSIXct(as.numeric(tail(gitLog,2)[1]),
                              origin="1970-01-01"),
                   format="%Y-%m-%d %H:%M:%S"),
            " ",tail(gitLog,1))
    message("## git HEAD: ",gitLog[2],"\n###########\n")
  }
  message("\nLoading resources")
}
## ----
run_cmd <- function(cmd){
  cmdR <- tryCatch(system(cmd,intern=T),
                   error=function(e){
                     message(e)
                     return("")
                   })
  return(cmdR)
}

## main functions -----
main <- function(){
  args = commandArgs()
  strPipePath <<- normalizePath(dirname(dirname(sapply(strsplit(grep("file=",args,value=T),"="),tail,1))))
  MsgInit()
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
  dir.create(strDir,showWarnings=F)
  strDir <- normalizePath(strDir)
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
checkH5adRefSetting <- function(config){
  xy <- getobsm(config$ref_h5ad,paste0("X_",config$ref_PCA))
  if(is.null(xy)) MsgExit(paste0("The 'ref_PCA' (",config$ref_PCA,") is not in the h5ad file"))
  
  meta <- getobs(config$ref_h5ad)
  if(sum(!config$ref_label%in%colnames(meta))>0)
    MsgExit(paste0("The following annotation labels (case sensitive) are not in the h5ad file:\n",
                   paste(config$ref_label[!config$ref_label%in%colnames(meta)],collapse=", ")))
}
checkRDSRefSeting <- function(config,D){
  if(!"SCT" %in% names(D@assays)) MsgExit(paste0("The 'SCT' assay is required but not in the rds file"))
  if(!config$ref_PCA %in% names(D@reductions)) MsgExit(paste0("The 'ref_PCA' (",config$ref_PCA,") is not in the rds file"))
  if(sum(!config$ref_label%in%colnames(D@meta.data))>0)
    MsgExit(paste0("The following annotation labels (case sensitive) are not in the rds file:\n",
                   paste(config$ref_label[!config$ref_label%in%colnames(D@meta.data)],collapse=", ")))
  #https://github.com/satijalab/azimuth/wiki/Azimuth-Reference-Format
  if(ncol(D[[config$ref_PCA]])<50)
    MsgExit(paste0(config$ref_PCA," should contain at least 50 dimensions for use with Azimuth"))
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
  VariableFeatures(SCT) <- SelectIntegrationFeatures(Dlist,nfeatures=3000)
  return(SCT)
}
unifySCTmodel <- function(SCT){
  if(length(levels(SCT[["SCT"]])) == 1) return(SCT)
  # obtain one unified SCT model for Azimuth, and the SCT model is not used for mapping which cause problem of normalizing
  Dtmp <- SCTransform(SCT,method = 'glmGamPoi',new.assay.name="SCT",return.only.var.genes = FALSE,verbose = FALSE)
  if(is.null(levels(SCT[["SCT"]]))){
    Dtmp[["SCT"]]$data <- SCT[["SCT"]]$data
    Dtmp[["SCT"]]$scale.data <- SCT[["SCT"]]$scale.data
    Dtmp[["SCT"]]$var.features <- SCT[["SCT"]]$var.features
    SCT <- Dtmp
  }else{
    SCT[["SCT"]]@SCTModel.list <- Dtmp[["SCT"]]@SCTModel.list
  }
  return(SCT)
}
saveRef <- function(D,config,sysConfig,strAzimuth="azimuth"){
  message("saving to scAnalyzer ...")
  strRef <- file.path(sysConfig$refDir,
                       gsub("[ [:punct:]]","_",paste(strAzimuth,config$ref_name,config$ref_src,config$ref_platform,"ref")),
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
  system(paste("chmod a+r -R",dirname(strRef)))
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
  sysConfig <- yaml::read_yaml(paste0(strPipePath,"/src/sys.yml"))
  config <- checkConfig(strConfig,sysConfig)
  suppressMessages(suppressWarnings(loadingPKG()))
  #"/camhpc/ngs/projects/TST11837/dnanexus/20220311155416_maria.zavodszky/sc20220403_0/TST11837_SCT.h5ad"
  strRef <- file.path(config$output,paste0(config$ref_name,"_for_scAnalyzer.rds"))
  if(!file.exists(strRef)){
    strTemp <- file.path(config$output,"ref_notFor_scAnalyzer.rds")
    if(!file.exists(strTemp)){
      if(!is.null(config$ref_rds) && file.exists(config$ref_rds)){
        message("Reading rds file @",config$ref_rds," ...")
        D <- readRDS(config$ref_rds)
        DefaultAssay(D) <- "SCT"
        checkRDSRefSeting(config,D)
        if(length(D@assays$SCT@var.features)==0){
          VariableFeatures(D) <- D@assays[[D@reductions[[config$ref_PCA]]@assay.used]]@var.features
        }
        D@reductions[[config$ref_PCA]]@assay.used <- "SCT"
      }else if(!is.null(config$ref_h5ad_raw) && file.exists(config$ref_h5ad_raw)){
        message("Reading h5ad file ...")
        checkH5adRefSetting(config)
        D <- getSCT(config$ref_h5ad_raw,config$ref_batch)
        DefaultAssay(D) <- "SCT"
        xy <- getobsm(config$ref_h5ad,paste0("X_",config$ref_PCA))
        rownames(xy) <- colnames(D)
        D[[config$ref_PCA]] <- CreateDimReducObject(embeddings=xy[colnames(D),],
                                                    key="PC_",
                                                    assay="SCT")
      }else{
        MsgExit(paste0("Either one seurat object rds file or two h5ad files are required to be existed"))
      }
      
      message("Processing ...")
      D <- unifySCTmodel(D)
      D <- FindNeighbors(D, dims = 1:30, reduction=config$ref_PCA,verbose = FALSE)
      D <- RunSPCA(D, npcs=ncol(D[[config$ref_PCA]]), graph = 'SCT_snn')
      D <- RunUMAP(D, dims = 1:30, reduction="spca",return.model=TRUE)
      saveRDS(D,strTemp)
    }else{
      message("Previous tmp file found @",strTemp,"\nPlease remove it if all new reference is needed.\nLoading ...")
      D <- readRDS(strTemp)
    }
    
    message("Creating Azimuth reference ...")
    D_ref <- AzimuthReference(
      D,
      refUMAP = "umap",
      refDR = "spca",
      refAssay = "SCT",
      dims = 1:ncol(D[[config$ref_PCA]]),
      plotref = "umap",
      metadata = config$ref_label #c('genotype','age','sex','cluster',"celltype")
    )
    saveRDS(D_ref,strRef)
  }else{
    message("Found the existed ref @",strRef,"\n\tPlease remove/rename it rerun is prefered!")
  }
  
  if(config$publish){
    D_Azimuth <- readRDS(strRef)
    saveRef(D_Azimuth,config,sysConfig)
  }else{
    message("The private reference could be used by provide the following full path to 'ref_name' in scAnalyzer config file:")
    message("\t",strRef)
  }
    
  
}
main()
