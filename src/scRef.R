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
  require(Azimuth)
  require(scales)
  options(future.globals.maxSize=8000*1024^2,stringsAsFactors=F)
  #source(paste0(strPipePath,"/src/azimuth.R"))
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
checkConfig <- function(strConfig,refPath){
  message("Checking config file ...")
  config <- yaml::read_yaml(strConfig)
  stopifnot(file.exists(config$ref_raw))
  stopifnot(length(config$ref_reduction)==1)
  stopifnot(length(config$ref_label)>0)
  stopifnot(length(config$ref_name)==1)
  stopifnot(nchar(config$ref_name)>3)
  
  if(config$publish){
    sapply(c("ref_version","ref_summary","ref_species","ref_system","ref_tech"),
           function(x,sInfo){
             a <- sInfo[[x]]
             if(length(a)==0 || nchar(a[1])==0) MsgExit(paste0("Please provided required infomraiton @",x))
           },config)
    seuratRef <- SeuratData::AvailableData()$Dataset
    if(config$ref_name%in%SeuratData::AvailableData()$Dataset)
      MsgExit(paste0("Seurat ref (",config$ref_name,") exists and cannot be overwrriten"))
    strSysRef <- file.path(refPath,"scRNAsequest_ref.csv")
    if(!file.exists(strSysRef))
      cat("Dataset,Version,Summary,species,system,ncells,tech\n",sep="",file=strSysRef)
    allRef <- data.table::fread(strSysRef)
    if(config$ref_name%in%allRef$Dataset && !config$overwrite){
      MsgExit(paste0("Public ref (",config$ref_name,") exists!"))
    }
  }
  return(config)
}
checkH5adRefSetting <- function(config){
  xy <- getobsm(config$ref_h5ad,paste0("X_",config$ref_reduction))
  if(is.null(xy)) MsgExit(paste0("The 'ref_reduction' (",config$ref_reduction,") is not in the h5ad file"))
  
  if(ncol(xy)<50) MsgExit(paste("At least 50 dimentions are required in reduction where",config$ref_reduction,"only contains",ncol(xy)))
  
  meta <- getobs(config$ref_h5ad)
  if(sum(!config$ref_label%in%colnames(meta))>0)
    MsgExit(paste0("The following annotation labels (case sensitive) are not in the h5ad file:\n",
                   paste(config$ref_label[!config$ref_label%in%colnames(meta)],collapse=", ")))
}
checkRDSRefSeting <- function(config,D){
  refAssay <- DefaultAssay(D)
  message("Default (active) seurat assay is: ",refAssay)
  #if(!"SCT" %in% names(D)) MsgExit(paste0("The 'SCT' assay is required but not in the rds file"))
  if(!"SCTModel.list" %in% slotNames(D[[refAssay]]) || length(D[[refAssay]]@SCTModel.list)==0) 
    MsgExit("The SCTModel.list is reqired in the activated assay ",refAssay)
  if(!config$ref_reduction %in% names(D)) MsgExit(paste0("The 'ref_reduction' (",config$ref_reduction,") is not in the rds file"))
  if(ncol(D[[config$ref_reduction]])<50) MsgExit(paste0("At least 50 dimentions are required in reduction where ",config$ref_reduction," only contains ",ncol(D[[config$ref_reduction]])))
  #if(D@reductions[[config$ref_reduction]]@assay.used!="SCT") MsgExit(paste0("The specified PCA (",config$ref_reduction,") did NOT derived from SCT but ",D@reductions$pca@assay.used))
  if(sum(dim(D@assays$SCT@data)==dim(D@assays$SCT@scale.data))!=2)
    MsgExit(paste0("Both of 'data' and 'scale.data' are required in SCT assay with the same dimensions"))
  if(sum(!config$ref_label%in%colnames(D@meta.data))>0)
    MsgExit(paste0("The following annotation labels (case sensitive) are not in the rds file:\n",
                   paste(config$ref_label[!config$ref_label%in%colnames(D@meta.data)],collapse=", ")))
  #https://github.com/satijalab/azimuth/wiki/Azimuth-Reference-Format
  #if(!"SCT_snn"%in%names(D)) warning(paste("SCT_snn neighber is not available, will be calculated based on SCT assay and provided reduction",config$ref_reduction))
}
prepareSeurat <- function(D,batch){
  refAssay <- DefaultAssay(D)
  if(length(D[[refAssay]]@SCTModel.list)==1){
    if(length(VariableFeatures(D))<1)
      D <- FindVariableFeatures(D)
  }else{
    message("Multiple SCT models detected!\nUsing ",batch," to integrate SCT models")
    Dlist <- SplitObject(D,split.by=batch)
    D <- integrateSCT(Dlist)
  }
  return(D)
}
integrateSCT <- function(Dlist){
  if(length(Dlist)==1) return(Dlist[[1]])
  features <- SelectIntegrationFeatures(Dlist,nfeatures=3000)
  Dlist <- PrepSCTIntegration(Dlist,anchor.features=features)
  anchors <- FindIntegrationAnchors(
    object.list = Dlist,
    anchor.features = features,
    normalization.method = "SCT",
    dims = 1:50
  )
  SCT <- IntegrateData(
    anchorset = anchors,
    normalization.method = "SCT",
    dims = 1:50
  )
  return(SCT)
}
getSCT <- function(strH5ad,batch){
  D <- CreateSeuratObject(counts=getX(strH5ad),
                          project="SCT",
                          meta.data=getobs(strH5ad))
  if(!batch%in%colnames(D@meta.data)) MsgExit("ref_batch (",batch,") is not one of the annotations")
  Dlist <- SplitObject(D,split.by=batch)
  rm("D")
  gc()
  message("\tSCTransform ...")
  Dlist <- sapply(Dlist,function(one){
    message("\t\tSCT ",unique(unlist(one[[batch]],use.names=F))," ...")
    return(suppressMessages(suppressWarnings(
      SCTransform(one,method = 'glmGamPoi',
                  new.assay.name="SCT",
                  return.only.var.genes = FALSE,
                  verbose = FALSE)
    )))
  })
  return(integrateSCT(Dlist))
}
unifySCTmodel <- function(SCT){
  if(length(levels(SCT[["SCT"]])) == 1) return(SCT)
  # obtain one unified SCT model for Azimuth, and the SCT model is not used for mapping which cause problem of normalizing
  
  
  
  
  
  
  Dtmp <- SCTransform(SCT,method = 'glmGamPoi',new.assay.name="SCT",return.only.var.genes = FALSE,verbose = FALSE)
  if(is.null(levels(SCT[["SCT"]]))){
    if(min(dim(SCT[["SCT"]]@data))>100) Dtmp[["SCT"]]@data <- SCT[["SCT"]]@data
    if(sum(dim(SCT[["SCT"]]@scale.data)==dim(SCT[["SCT"]]@data))==2) Dtmp[["SCT"]]@scale.data <- SCT[["SCT"]]@scale.data
    if(length(SCT[["SCT"]]@var.features)>100) Dtmp[["SCT"]]@var.features <- SCT[["SCT"]]@var.features
    SCT <- Dtmp
  }else{
    SCT[["SCT"]]@SCTModel.list <- Dtmp[["SCT"]]@SCTModel.list
  }
  return(SCT)
}
saveRef <- function(strRef,config,refDir,nCell){
  message("saving to the scRNAsequest ...")
  strSysRef <- file.path(refPath,"scRNAsequest_ref.csv")
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
  data.table::fwrite(strSysRef)
  message("\nA new reference (",config$ref_name,") is added into the scRNAsequest!")
  MsgPower()
}
createRef <- function(strConfig){
  customRef <- file.path(strPipePath,"src","sys_ref.csv")
  sysRefDir <- yaml::read_yaml(paste0(strPipePath,"/src/sys.yml"))$refDir
  config <- checkConfig(strConfig,sysRefDir)
  suppressMessages(suppressWarnings(loadingPKG()))
  source(paste0(dirname(gsub("--file=","",grep("file=",commandArgs(),value=T))),"/readH5ad.R"))
  #"/camhpc/ngs/projects/TST11837/dnanexus/20220311155416_maria.zavodszky/sc20220403_0/TST11837_SCT.h5ad"
  strRef <- file.path(config$output,config$ref_name)

  if(!dir.exists(strRef) || !file.exists(file.path(strRef,'ref.Rds')) || !file.exists(file.path(strRef,'idx.annoy'))){
    strTemp <- file.path(config$output,"ref_notFor_scAnalyzer.rds")
    if(!file.exists(strTemp)){
      if(!is.null(config$ref_rds) && file.exists(config$ref_rds)){
        message("Reading rds file @",config$ref_rds)
        D <- readRDS(config$ref_rds)
        checkRDSRefSeting(config,D)
        D <- prepareSeurat(D)
        D@reductions[[config$ref_reduction]]@assay.used <- DefaultAssay(D)
      }else if(!is.null(config$ref_h5ad_raw) && file.exists(config$ref_h5ad_raw)){
        message("Reading h5ad file ...")
        checkH5adRefSetting(config)
        D <- getSCT(config$ref_h5ad_raw,config$ref_batch)
        xy <- getobsm(config$ref_h5ad,paste0("X_",config$ref_reduction))
        rownames(xy) <- colnames(D)
        D[[config$ref_reduction]] <- CreateDimReducObject(embeddings=xy[colnames(D),],
                                                    key="PC_",
                                                    assay=DefaultAssay(D))
      }else{
        MsgExit(paste0("Either one seurat object rds file or two h5ad files are required to be existed"))
      }
      
      message("Azimuth Processing ...")
      #https://github.com/satijalab/azimuth/wiki/Azimuth-Reference-Format
      #D <- unifySCTmodel(D)
      dim_n <- ncol(D[[config$ref_reduction]])
      graph <- paste0(DefaultAssay(D),"_snn")
      if(!graph%in%names(D)) D <- FindNeighbors(D, dims = 1:dim_n, reduction=config$ref_reduction,verbose = FALSE)
      D <- RunSPCA(D, npcs=dim_n, graph = graph)
      D <- RunUMAP(D, dims = 1:dim_n, reduction="spca",umap.method="umap-learn",metric = "correlation",return.model=TRUE)
      saveRDS(D,strTemp)
    }else{
      message("Previous tmp file found @",strTemp,"\nPlease remove it if all new reference is needed.\nLoading ...")
      D <- readRDS(strTemp)
    }
    message("Creating Azimuth reference ...")
    D_ref <- AzimuthReference(D,refAssay=DefaultAssay(D),
                              metadata=config$ref_label)
    dir.create(strRef,showWarnings=F,recursive=T)
    SaveAzimuthReference(D_ref,strRef)
  }else{
    message("Found the existed ref @",strRef,"\n\tPlease remove/rename it rerun is prefered!")
    D_ref <- readRDS(file.path(strRef,"ref.Rds"))
  }
  
  if(config$publish){
    saveRef(strRef,config,sysRefDir,dim(D_ref)[2])
  }else{
    message("The private reference can be used by provide the following full path to 'ref_name' in scAnalyzer config file:")
    message("\t",strRef)
  }
}
main()
