#!/usr/bin/env Rscript

strPipePath <- ""
global_feature_n <- 5000
## loading packages -----
loadingPKG <- function(){
  require(Seurat)
  require(Azimuth)
  require(future)
  require(BiocParallel)
  require(dplyr)
  options(future.globals.maxSize=8000*1024^2,stringsAsFactors=F)
  #require(cowplot)
  #require(patchwork)
  #require(sctransform)
  #require(rhdf5)
  #require(Matrix)
  #require(scales)
  #source(paste0(strPipePath,"/src/azimuth.R"))
}
## msg -----
MsgExit <- function(...){
  msg <- paste0(...)
  if(length(msg)>0 && nchar(msg)>3) message("ERROR: ",msg)
  MsgPower()
  q()
}
MsgPower <- function(){
  message("\nPowered by ", yaml::read_yaml(paste0(strPipePath,"/src/sys.yml"))$powerby)
  message("------------")
}
MsgHelp <- function(){
  message("\nscRef /path/to/a/output/folder === or === scRef /path/to/a/Ref/config/file\n")
  message("The folder has to be existed.")
  message("The Ref config file will be generated automatically when a path is provided")
  message("===== CAUTION =====")
  message("\t1. If 'publish: True' is set in config. This process will add a seurat reference data into the scRNAsequest pipeline PERMANENTLY!")
  message("\t2. If seurat object is provided, make sure the data provided for reference building is SCT transformed!")
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
createRef <- function(strConfig){
  customRef <- file.path(strPipePath,"src","sys_ref.csv")
  sysRefDir <- yaml::read_yaml(paste0(strPipePath,"/src/sys.yml"))$refDir
  source(paste0(dirname(gsub("--file=","",grep("file=",commandArgs(),value=T))),"/readH5ad.R"))
  config <- checkConfig(strConfig,sysRefDir)
  suppressMessages(suppressWarnings(loadingPKG()))
  source(paste0(dirname(gsub("--file=","",grep("file=",commandArgs(),value=T))),"/azimuth.R"))
  config$core = ifelse(is.null(config$core),4,config$core)
  
  strRef <- file.path(config$output,config$ref_name)
  if(!dir.exists(strRef) || !file.exists(file.path(strRef,'ref.Rds')) || !file.exists(file.path(strRef,'idx.annoy'))){
    succ <- tryCatch({
      suppressWarnings(
      switch(config$ref_integration,
             Layer=sctIntegrationLayer(config$ref_raw,
                                       config$ref_batch,
                                       strRef,
                                       config$ref_label,
                                       config$core,
                                       config$ref_reduction),
             RPCA=sctIntegrationRPCA(config$ref_raw,
                                     config$ref_batch,
                                     strRef,
                                     config$ref_label,
                                     config$core,
                                     config$ref_reduction),
             CCA=sctIntegrationCCA(config$ref_raw,
                                   config$ref_batch,
                                   strRef,
                                   config$ref_label,
                                   config$core,
                                   config$ref_reduction)
             ))
    },error=function(eMsg){
      return(F)
    })
    if(is.null(succ))MsgExit("unknown SCT integration method: ",config$ref_integration)
    if(!succ) MsgExit("Failed ",config$ref_integration," integration on SCT! Please reduce the number of batch samples or try to use other integration method!")
    saveInfo(strConfig,strRef)
  }else{
    message("Found the existed ref @",strRef,"\n\tPlease remove/rename it if a rerun is prefered!")
  }
  if(config$publish){
    D_ref <- readRDS(file.path(strRef,"ref.Rds"))
    saveRef(strRef,config,sysRefDir,dim(D_ref)[2])
  }else{
    message("The private reference can be used by provide the following full path to 'ref_name' in scAnalyzer config file:")
    message("\t",strRef)
  }
  MsgPower()
}

checkExist <- function(strF){
  if(is.null(strF) || !file.exists(strF))
    return(F)
  return(T)
}
checkConfig <- function(strConfig,refPath){
  message("Checking config file ...")
  config <- yaml::read_yaml(strConfig)
  stopifnot(file.exists(config$ref_raw))
  stopifnot(!is.null(config$ref_batch))
  stopifnot(length(config$ref_label)>0)
  stopifnot(length(config$ref_name)==1)
  stopifnot(nchar(config$ref_name)>3)
  
  if(config$publish){
    sapply(c("ref_version","ref_summary","ref_species","ref_system","ref_tech"),
           function(x,sInfo){
             a <- sInfo[[x]]
             if(length(a)==0 || nchar(a[1])==0) MsgExit("Please provided required infomraiton @",x)
           },config)
    seuratRef <- SeuratData::AvailableData()$Dataset
    if(config$ref_name%in%SeuratData::AvailableData()$Dataset)
      MsgExit("Seurat ref (",config$ref_name,") exists and cannot be overwrriten")
    strSysRef <- file.path(refPath,"scRNAsequest_ref.csv")
    if(!file.exists(strSysRef))
      cat("Dataset,Version,Summary,species,system,ncells,tech\n",sep="",file=strSysRef)
    allRef <- data.table::fread(strSysRef)
    if(config$ref_name%in%allRef$Dataset && !config$overwrite){
      MsgExit("Public ref (",config$ref_name,") exists!")
    }
  }
  checkH5adRefSetting(config$ref_raw,config$ref_reduction,config$ref_batch,config$ref_label)
  return(config)
}
checkH5adRefSetting <- function(ref_h5ad,ref_reduction,ref_batch,ref_label){
  if(!grepl("h5ad$",ref_h5ad)) return()
  if(!is.null(ref_reduction)){
    xy <- getobsm(ref_h5ad,paste0("X_",ref_reduction))
    if(is.null(xy)) MsgExit("The 'ref_reduction' (",ref_reduction,") is not in the h5ad file")
    if(ncol(xy)<50) MsgExit("At least 50 dimentions are required in reduction where",ref_reduction,"only contains",ncol(xy))
  }
  
  meta <- getobs(ref_h5ad)
  if(!is.null(ref_batch) && !ref_batch%in%colnames(meta))
    MsgExit("The following annotation labels (case sensitive) are not in the h5ad file:\n",
                   paste(ref_label[!ref_label%in%colnames(meta)],collapse=", "))
  if(sum(!ref_label%in%colnames(meta))>0)
    MsgExit("The following annotation labels (case sensitive) are not in the h5ad file:\n",
                   paste(ref_label[!ref_label%in%colnames(meta)],collapse=", "))
  
}
saveRef <- function(strRef,config,refDir,nCell){
  message("saving to the scRNAsequest ...")
  strSysRef <- file.path(refDir,"scRNAsequest_ref.csv")
  allRef <- data.table::fread(strSysRef)
  
  if(config$ref_name%in%list.dirs(refDir,full.names=F)){
    file.rename(file.path(refDir,config$ref_name),
                file.path(refDir,paste0(config$ref_name,"_archived",format(Sys.time(),"%Y%m%d"),"_",Sys.getenv("USER")))) 
  }
  file.copy(strRef,refDir,recursive = TRUE)
  
  allRef <- data.frame(data.table::fread(strSysRef))
  if(config$ref_name%in%allRef$Dataset) allRef <- allRef[allRef$Dataset!=config$ref_name,]
  allRef <- rbind(allRef,
                  data.frame(Dataset=config$ref_name,
                             Version=config$ref_version,
                             Summary=config$ref_summary,
                             species=config$ref_species,
                             system=config$ref_system,
                             ncells=nCell,
                             tech=config$ref_tech))
  data.table::fwrite(allRef,strSysRef)
  message("\nA new reference (",config$ref_name,") is added into the scRNAsequest!")
}
saveInfo <- function(strConfig,strOut){
  file.copy(strConfig,strOut)
  conn <- file(file.path(strOut,"readme"),"w")
  sink(file=conn,type='message')
  MsgInit()
  sink(type='message')
  close(conn)
}

# put the whole process in one function avoiding memory copy in R
sctIntegrationCCA <- function(strRaw,batch,strRef,ref_label,core=4,ref_reduction=NULL){
  message("\n\n=== CCA integration on SCT assays")
  #plan("default")
  refReduct <- NULL
  if(grepl("rds$",strRaw)){
    D <- readRDS(strRaw)
    if(!is.null(ref_reduction))
      refReduct <- D@reductions[[ref_reduction]]@cell.embeddings
  }else if(grepl("h5seurat$",strRaw)){
    D <- LoadH5Seurat(strRaw)
    if(!is.null(ref_reduction))
      refReduct <- D@reductions[[ref_reduction]]@cell.embeddings
  }else if(grepl("h5ad$",strRaw)){
    D <- CreateSeuratObject(counts=getX(strRaw),
                            project="scRef",
                            meta.data=getobs(strRaw))
    if(!is.null(ref_reduction)){
      refReduct <- getobsm(strRaw,paste0("X_",ref_reduction))
      rownames(refReduct) <- colnames(D)
    }
  }else{
    stop("unknown format: ",strRaw)
  }
  if("SCTModel.list" %in% slotNames(D[[DefaultAssay(D)]]) && length(D[[DefaultAssay(D)]]@SCTModel.list)==1){
    message("*** Directly use ",DefaultAssay(D)," assay with only 1 SCT model!")
  }else{
    Dlist <- SplitObject(D,split.by=batch)
    rm(D)
    message("\tSCT (",length(Dlist),")...")
    Dlist <- bplapply(1:length(Dlist),function(i){
      message("\t\t",Dlist[[i]][[batch]][1,1]," ...")
      return(suppressMessages(suppressWarnings(
        SCTransform(Dlist[[i]],vst.flavor="v2",
                    return.only.var.genes = FALSE,
                    verbose = FALSE)
      )))
    },BPPARAM = MulticoreParam(workers=min(core,length(Dlist)),tasks=length(Dlist)))
    message("\tPrepare Integration ...")
    features <- SelectIntegrationFeatures(Dlist,nfeatures=global_feature_n,verbose = FALSE)
    Dlist <- PrepSCTIntegration(Dlist,assay="SCT",anchor.features=features,verbose = FALSE)
    #plan("multisession", workers = core)
    message("\tFind Anchors ...")
    anchors <- FindIntegrationAnchors(
      object.list = Dlist,
      anchor.features = features,
      normalization.method = "SCT",
      dims = 1:30,
      verbose = FALSE
    )
    rm(Dlist)
    message("\tCCA Integration ...")
    D <- IntegrateData(
      anchorset = anchors,
      normalization.method = "SCT",
      verbose = FALSE
    )
  }
  reductName <- "scRNASequest"
  if(!is.null(refReduct)){
    message("\t***Using the reduction: ",ref_reduction)
    D[[reductName]] <- CreateDimReducObject(embeddings=refReduct[colnames(D),],
                                                key="PC_",
                                                assay=DefaultAssay(D))
  }else{
    message("\tPCA ...")
    D <- RunPCA(D,npcs=50,verbose = FALSE,reduction.name=reductName)
  }
  dim_n <- ncol(D[[reductName]])
  graph <- paste0(DefaultAssay(D),"_snn")
  message("\tFind neighbors ...")
  if(!graph%in%names(D)) D <- FindNeighbors(D,dims=1:dim_n,
                                            reduction=reductName,
                                            verbose=FALSE)
  message("\tSPCA ...")
  D <- RunSPCA(D, npcs=dim_n,graph=graph,verbose=F)
  message("\tUMAP ...")
  suppressMessages(
    D <- RunUMAP(D,dims=1:dim_n,reduction="spca",
                 umap.method="umap-learn",metric="correlation",
                 return.model=TRUE,verbose=F))
  message("\tSave tmp")
  saveRDS(D,paste0(strRef,".rds"))

  message("Creating Azimuth reference ...")
  D_ref <- AzimuthReference(D,refAssay=DefaultAssay(D),
                            metadata=ref_label)
  dir.create(strRef,showWarnings=F,recursive=T)
  SaveAzimuthReference(D_ref,paste0(strRef,.Platform$file.sep))
  return(T)
}
sctIntegrationRPCA <- function(strRaw,batch,strRef,ref_label,core=4,ref_reduction=NULL){
  message("\n\n=== RPCA integration on SCT assays")
  #plan("default")
  refReduct <- NULL
  if(grepl("rds$",strRaw)){
    D <- readRDS(strRaw)
    if(!is.null(ref_reduction))
      refReduct <- D@reductions[[ref_reduction]]@cell.embeddings
  }else if(grepl("h5seurat$",strRaw)){
    D <- LoadH5Seurat(strRaw)
    if(!is.null(ref_reduction))
      refReduct <- D@reductions[[ref_reduction]]@cell.embeddings
  }else if(grepl("h5ad$",strRaw)){
    D <- CreateSeuratObject(counts=getX(strRaw),
                            project="scRef",
                            meta.data=getobs(strRaw))
    if(!is.null(ref_reduction)){
      refReduct <- getobsm(strRaw,paste0("X_",ref_reduction))
      rownames(refReduct) <- colnames(D)
    }
  }else{
    stop("unknown format: ",strRaw)
  }
  if("SCTModel.list" %in% slotNames(D[[DefaultAssay(D)]]) && length(D[[DefaultAssay(D)]]@SCTModel.list)==1){
    message("*** Directly use ",DefaultAssay(D)," assay with only 1 SCT model!")
  }else{
    Dlist <- SplitObject(D,split.by=batch)
    rm(D)
    message("\tSCT (",length(Dlist),")...")
    Dlist <- bplapply(1:length(Dlist),function(i){
      message("\t\t",Dlist[[i]][[batch]][1,1]," ...")
      return(suppressMessages(suppressWarnings(
        SCTransform(Dlist[[i]],vst.flavor="v2",
                    return.only.var.genes = FALSE,
                    verbose = FALSE)
      )))
    },BPPARAM = MulticoreParam(workers=min(core,length(Dlist)),tasks=length(Dlist)))
    message("\tPrepare Integration ...")
    features <- SelectIntegrationFeatures(Dlist,nfeatures=global_feature_n,verbose = FALSE)
    Dlist <- PrepSCTIntegration(Dlist,assay="SCT",anchor.features=features,verbose = FALSE)
    Dlist <- bplapply(1:length(Dlist),function(i){
      return(suppressMessages(suppressWarnings(
        RunPCA(Dlist[[i]],features=features,verbose = FALSE)
      )))
    },BPPARAM = MulticoreParam(workers=min(core,length(Dlist)),tasks=length(Dlist)))
    #plan("multisession", workers = core)
    message("\tFind RPCA Anchors ...")
    anchors <- FindIntegrationAnchors(
      object.list = Dlist,
      anchor.features = features,
      normalization.method = "SCT",
      dims = 1:30,
      reduction = "rpca",
      verbose = FALSE
    )
    rm(Dlist)
    message("\tRPCA Integration ...")
    D <- IntegrateData(
      anchorset = anchors,
      normalization.method = "SCT",
      features=features,
      verbose = FALSE
    )
  }
  reductName <- "scRNASequest"
  if(!is.null(refReduct)){
    message("\t*** Using the reduction: ",ref_reduction)
    D[[reductName]] <- CreateDimReducObject(embeddings=refReduct[colnames(D),],
                                            key="PC_",
                                            assay=DefaultAssay(D))
  }else{
    message("\tPCA ...")
    D <- RunPCA(D,npcs=50,verbose = FALSE,reduction.name=reductName)
    
  }
  dim_n <- ncol(D[[reductName]])
  graph <- paste0(DefaultAssay(D),"_snn")
  message("\tFind neighbors ...")
  if(!graph%in%names(D)) D <- FindNeighbors(D,dims=1:dim_n,
                                            reduction=reductName,
                                            verbose=FALSE)
  message("\tSPCA ...")
  D <- RunSPCA(D, npcs=dim_n,graph=graph,verbose=F)
  message("\tUMAP ...")
  suppressMessages(
    D <- RunUMAP(D,dims=1:dim_n,reduction="spca",
                 umap.method="umap-learn",metric="correlation",
                 return.model=TRUE,verbose=F))
  
  message("\tSave tmp")
  saveRDS(D,paste0(strRef,".rds"))
  
  message("Creating Azimuth reference ...")
  browser()
  D_ref <- AzimuthReference(D,refAssay=DefaultAssay(D),
                            metadata=ref_label)
  dir.create(strRef,showWarnings=F,recursive=T)
  SaveAzimuthReference(D_ref,paste0(strRef,.Platform$file.sep))
  return(T)
}
sctIntegrationLayer <- function(strRaw,batch,strRef,ref_label,core=4,ref_reduction=NULL){
  message("\n\n=== Harmony integration on SCT assays")
  #plan("multisession", workers = core)
  refReduct <- NULL
  if(grepl("rds$",strRaw)){
    D <- readRDS(strRaw)
    if(!is.null(ref_reduction))
      refReduct <- D@reductions[[ref_reduction]]@cell.embeddings
  }else if(grepl("h5seurat$",strRaw)){
    D <- LoadH5Seurat(strRaw)
    if(!is.null(ref_reduction))
      refReduct <- D@reductions[[ref_reduction]]@cell.embeddings
  }else if(grepl("h5ad$",strRaw)){
    D <- CreateSeuratObject(counts=getX(strRaw),
                            project="scRef",
                            meta.data=getobs(strRaw))
    if(!is.null(ref_reduction)){
      refReduct <- getobsm(strRaw,paste0("X_",ref_reduction))
      #colnames(refReduct) <- paste
    }
  }else{
    stop("unknown format: ",strRaw)
  }
  refModel <- NULL
  if("SCTModel.list" %in% slotNames(D[[DefaultAssay(D)]]) && length(D[[DefaultAssay(D)]]@SCTModel.list)==1){
    message("*** Directly use ",DefaultAssay(D)," assay with only 1 SCT model!")
  }else{
    message("\tSCT ...")
    nBatches <- length(unique(unlist(D[[batch]],use.names=F)))
    if(nBatches==1){
    	D <- SCTransform(D,return.only.var.genes=F)#,vst.flavor="v2",verbose=F
    }else{
    	message("\tFind unify SCT model")
    	A <- D@meta.data %>% dplyr::group_by(across(all_of(batch))) %>% dplyr::count()
    	ncells <- ceiling(max(A$n)/1000)*1000
    	selID <- A[order(A$n,decreasing=T),][[batch]][1:min(3,nrow(A))]
    	# select samples with most annotations
    	A <- D@meta.data %>% 
    		dplyr::filter(!!sym(batch)%in%selID) %>% 
    		dplyr::group_by(across(all_of(c(batch,ref_label)))) %>% 
    		dplyr::count()
    	selID <- names(table(A[batch]))[table(A[batch])==max(table(A[batch]))]
    	# select samples with max cells in the smallest annotation
    	if(length(selID)>1){
    		A <- A %>% dplyr::filter(!!sym(batch)%in%selID)%>% dplyr::group_by(across(all_of(batch))) %>% dplyr::summarise(MIN=min(n))
    		selID <- A[[batch]][order(A$MIN,decreasing=T)][1]
    	}
    	subD <- subset(D,batch==selID)
    	subD <- SCTransform(subD,return.only.var.genes=F,variable.features.n=global_feature_n,ncells=ncells,verbose=F)
    	message("\tSplit batches ...")
    	D[["RNA"]] <- split(D[["RNA"]],f=unlist(D[[batch]],use.names=F))
    	refModel <- subD@assays$SCT@SCTModel.list[[1]]
    	D <- SCTransform(D,vst.flavor="v2",
    									 variable.features.n=global_feature_n,
    									 ncells=ncells,
    									 reference.SCT.model=refModel,
    									 return.only.var.genes=F,
    									 verbose=F)
    	#refModel <- D@assays$SCT@SCTModel.list[[selID]]
    	cA <- lapply(D[["SCT"]]@SCTModel.list,function(x)return(x@cell.attributes))
    	names(cA) <- NULL
    	refModel@cell.attributes <- do.call(rbind,cA)
    	#D[["SCT"]]@SCTModel.list <- list(refModel=refModel)
    }
  }
  reductName <- "scRNASequest"
  if(!is.null(refReduct)){
    message("\t*** Using the reduction: ",ref_reduction)
    D[[reductName]] <- CreateDimReducObject(embeddings=refReduct[colnames(D),],
                                                key="PC_",
                                                assay=DefaultAssay(D))
  }else{
    message("\tPCA ...")
    D <- RunPCA(D,npcs=50,verbose = FALSE)
    message("\tIntegration by HarmonyIntegration...")
    D <- IntegrateLayers(D,method=HarmonyIntegration,
                         orig.reduction="pca",new.reduction=reductName,
                         normalization.method = "SCT",
                         assay="SCT",verbose=FALSE)
  }
  dim_n <- ncol(D[[reductName]])
  graph <- paste0(DefaultAssay(D),"_snn")
  message("\tFind neighbors")
  if(!graph%in%names(D)) D <- FindNeighbors(D,dims=1:dim_n,
                                            reduction=reductName,
                                            verbose=FALSE)
  message("\tSPCA")
  D <- RunSPCA(D, npcs=dim_n,graph=graph,verbose=F)
  message("\tUMAP")
  suppressMessages(
    D <- RunUMAP(D,dims=1:dim_n,reduction="spca",
                 umap.method="umap-learn",metric="correlation",
                 return.model=TRUE,verbose=F))
  if(!is.null(refModel))
  	D[["SCT"]]@SCTModel.list <- list(refModel=refModel)
  message("\tSave tmp")
  saveRDS(D,paste0(strRef,"_ref.rds"))
  
  message("Creating Azimuth reference ...")
  D_ref <- suppressMessages(suppressWarnings(
    AzimuthReference(D,refAssay=DefaultAssay(D),metadata=ref_label)))
  dir.create(strRef,showWarnings=F,recursive=T)
  suppressMessages(suppressWarnings(
    SaveAzimuthReference(D_ref,paste0(strRef,.Platform$file.sep))))
  return(T)
}

main()
