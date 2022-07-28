strPipePath <- ""
loadingPKG <- function(){
  require(Seurat)
  require(cowplot)
  require(patchwork)
  require(sctransform)
  require(rhdf5)
  require(Matrix)
}

# message -----
MsgExit <- function(msg=""){
  if(nchar(msg)>3) message("ERROR: ",msg)
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
    message("###########\n## ExpressionAnalysis: ",url)
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
MsgHelp <- function(){
  message("\nsc2celldepot /path/to/a/output/folder === or === sc2celldepot /path/to/a/config/file\n")
  message("Please create the folder before running sc2celldepot.")
  message("The data config file will be generated automatically when a path is provided")
  MsgPower()
  q()
}
MsgComplete <- function(strh5ad){
  message(strh5ad," is created!")
  message("To add it into the CellDepot:\n1. Please copy the above h5ad file into the celldepot folder")
  message("\tThe celldepot folder can be found in celldepot website under 'Create Projects'")
  message("2. Fill the 'Create Projects' page, and click 'Save'")
}
MsgPower <- function(){
  sysConfig <- yaml::read_yaml(paste0(strPipePath,"/src/sys.yml"))
  message("\nPowered by ",sysConfig$powerby)
  message("------------")
}


# create initial config file -----
initConfig <- function(strD){
  dir.create(strD,showWarnings=F)
  config <- readLines(file.path(strPipePath,"src/external2h5ad.yml"))
  config <- gsub("OUTPUT_DIR",normalizePath(strD),config)
  writeLines(config,file.path(strD,"sc2celldepot.yml"))
  message("Please complete the config file @",file.path(strD,"sc2celldepot.yml"))
}

# process --------
process2H5ad <- function(strConfig){
  suppressMessages(suppressWarnings(loadingPKG()))
  config <- yaml::read_yaml(strConfig)
  if(is.null(config$prefix)||nchar(config$prefix)<2) MsgExit(paste("Please specify the prefix of result h5ad file"))
  if(!is.null(config$seuratObj)) return(seurat2h5ad(config))
  return(createH5ad(config))
  
}
## seurat to h5ad ------
seurat2h5ad <- function(config){
  message("Creating h5ad from processed Seurat Object")
  if(!file.exists(config$seuratObj)) MsgExit(paste("Missing seurat RDS @",config$seuratObj))
  message("reading seurat RDS")
  D <- readRDS(config$seuratObj)
  if(!config$seuratSCT %in% names(D@assays) || !config$seuratUMI %in% names(D@assays))
    MsgExit(paste0("assay (",config$seuratUMI," and/or ",config$seuratSCT,") is missing from ",config$seuratObj,"!"))
  if(length(config$seuratMeta)>0 && sum(!config$seuratMeta%in%colnames(D@meta.data)))
    MsgExit(paste("unknown meta column selected:",
                  paste(config$seuratMeta[!config$seuratMeta%in%colnames(D@meta.data)],collapse="; ")))
  return(transformH5ad(D,config$output,config$prefix,config$seuratSCT,
                       config$seuratUMI,config$seuratMeta))
}
## To h5ad -----
transformH5ad <- function(D,strOut,prefix,expAssay,umiAssay=NULL,selMeta=NULL){
  strTemp <- file.path(strOut,paste0("tmp",sample(100000,1)))
  if(!dir.exists(strTemp)) dir.create(strTemp)
  strH5ad <- file.path(strOut,paste0(prefix,".h5ad"))
  message("\tH5ad trasforming ...")
  a <- tryCatch(saveH5AD(D,strTemp,strH5ad,
                         selAssay=expAssay,
                         umiAssay=umiAssay,
                         selMeta=selMeta),
                error=function(cond){
                  message(cond$message)
                })
  unlink(strTemp,recursive=T)
  return(strH5ad)
}
saveH5AD <- function(D,strTemp,strH5ad,selAssay="SCT",umiAssay=NULL,selMeta=NULL,maxAnnoCat=50){
  message("\texpression matrix")
  X <- Matrix::t(Matrix::Matrix(D@assays[[selAssay]]@data,sparse=T))
  strH5 <- file.path(strTemp,"X.hdf5")
  saveHdf5(X,strH5)
  
  ## meta ----
  message("\tmeta information")
  tmp <- D@meta.data
  if(length(selMeta)>0) tmp <-tmp[,selMeta,drop=F]
  meta <- data.frame(row.names=rownames(tmp))
  for(i in 1:ncol(tmp)){
    if(is.numeric(tmp[,i])){
      tmp[is.na(tmp[,i]),i] <- 0
    }else if(length(unique(tmp[,i]))<maxAnnoCat){
      tmp[is.na(tmp[,i])|is.nan(tmp[,i])|tmp[,i]=="nan",i] <- "unknown"
    }else{
      next
    }
    meta <- cbind(meta,tmp[,i,drop=F])  
  }
  colnames(meta) <- gsub("\\.","_",colnames(meta))
  write.csv(meta,file=file.path(strTemp,"meta.csv"))
  
  ## reduction -----
  message("\treduction coordinates")
  for(i in names(D@reductions)){
    write.csv(D@reductions[[i]]@cell.embeddings,
              file=file.path(strTemp,paste0("layout.",i,".csv")))
  }
  
  ## to H5ad ----
  system(paste0(strPipePath,"/src/toH5ad.py ",
                strTemp," ",strH5ad))
  
  ## save the raw counts ----
  if(!is.null(umiAssay)){
    message("\traw count matrix")
    X <- Matrix::t(Matrix::Matrix(D@assays[[umiAssay]]@counts,sparse=T))
    saveHdf5(X,strH5)
    system(paste0(strPipePath,"/src/toH5ad.py ",
                  strTemp," ",gsub("h5ad$","raw_added.h5ad",strH5ad)))
  }
  
  
}
saveHdf5 <- function(X,strH5){
  if(file.exists(strH5)) file.remove(strH5)
  suppressMessages(suppressWarnings({
    h5write(X@x,file=strH5,name="data")
    h5write(X@i,file=strH5,name="indices")
  }))
  h5write(X@p,file=strH5,name="indptr")
  h5write(dim(X),file=strH5,name="shape")
  h5write(rownames(X),file=strH5,name="row_names")
  h5write(colnames(X),file=strH5,name="col_names")
  h5closeAll()
}
## seperated files (express/meta/reduction) to h5ad ------
createH5ad <- function(config){
  strRDS <- paste0(config$output,"/",config$prefix,".rds")
  if(file.exists(strRDS)){
    message("\n*** Using existing processed data found: ",strRDS,"! ***")
    D <- readRDS(strRDS)
  }else{
    message("Creating h5ad from expression and cell meta inforamtion")
    checkConfig(config)
    message("\tObtain meta information")
    meta <- getData(config$annotation)
    selMeta <- checkUse(meta,config$annotationUse)
    X <- getX(config$expression)
    reduction <- getData(config$reduction$files)
    D <- createSeuratObj(X,meta,selMeta,config)
    
    rawAssay <- NULL
    expAssay <- "RNA"
    if(config$dataUMI){
      if(is.null(config$sample_column))
        config$sample_column <- "orig.ident"
      D <- processSCT(D,config$sample_column)
      rawAssay <- "RNA"
      expAssay <- "SCT"
    }
    D <- processLayout(D,reduction,config)
    saveRDS(D,file=strRDS)
  }
  
  return(transformH5ad(D,config$output,config$prefix,
                       expAssay,rawAssay))
}
checkConfig <- function(config){
  message("")
  if(is.null(config$expression) || length(config$expression)==0) MsgExit(paste("Missing Expression file"))
  if(is.null(config$annotation) || length(config$annotation)==0) MsgExit(paste("Missing cell annotation file"))
  for(one in config$expression){
    #message(one)
    if(!file.exists(one)) MsgExit(paste("Missing Expression file:\n",one))
  }
    
  for(one in config$annotation){
    #message(one)
    if(!file.exists(one)) MsgExit(paste("Missing Expression file:\n",one))
  }
  if(!is.null(config$reduction$files) && length(config$reduction$files)>0)
    for(one in config$reduction$files)
      if(!file.exists(one)) MsgExit(paste("Missing reduction file:\n",one))
}
getX <- function(strFs){
  message("\tObtain expression")
  X <- list()
  for(one in strFs){
    if(dir.exists(one)){
      X[[length(X)+1]] <- Read10X(one)
    }else if(grepl("h5$",one)){
      X[[length(X)+1]] <- Read10X_h5(one)
    }else{
      tmp <- data.table::fread(one)
      X[[length(X)+1]] <- Matrix(as.matrix(data.table::fread(one),rownames=1),sparse=T)
    }
  }
  return(X)
}
getData <- function(strFs){
  meta <- list()
  for(one in strFs){
    tmp <- data.table::fread(one)
    meta[[length(meta)+1]] <- data.frame(tmp[,-1],row.names=unlist(tmp[,1]))
  }
  return(meta)
}
checkUse <- function(meta,selMeta){
  comm <- Reduce(intersect,lapply(meta,colnames))
  if(length(comm)==0) MsgExit(paste("No common columns among all annotation files"))
  if(is.null(selMeta) || length(selMeta)==0) return(comm)
  selMeta <- interaction(selMeta,comm)
  if(length(selMeta)==0) MsgExit(paste("No sepesified annotations are common columns among all annotation files"))
  return(selMeta)
}
createSeuratObj <- function(X,meta,selMeta,config){
  #save(X,meta,selMeta,config,file="data.RData")
  message("\tCreating seurat object from expression and annotation files...")
  meta <- lapply(meta, function(x)return(x[,selMeta]))
  if(length(X)==length(meta)){
    Dlist <- list()
    for(i in 1:length(X)){
      cID <- intersect(colnames(X[[i]]),rownames(meta[[i]]))
      message("\t\t",length(cID)," cells overlapped between expression and annotation for ",basename(config$expression[i]))
      if(length(cID)<10){
        message("\t\t\tLess than 10 cells! SKIP!")
        next
      }
      Dlist[[i]] <- CreateSeuratObject(counts=X[[i]][,cID],
                                       project=basename(config$expression[i]),
                                       meta.data=meta[[i]][cID,,drop=F])
      #UMI <- X[[i]][,colnames(X[[i]])%in%cID]
      #Dlist[[i]] <- CreateSeuratObject(counts=UMI,
      #                                 project=basename(config$expression[i]),
      #                                 meta.data=meta[[i]][colnames(UMI),,drop=F])
    }
    if(is.null(Dlist)) MsgExit(paste("Too few matching cell IDs"))
    if(length(Dlist)==1){
      D <- Dlist[[1]]
    }else{
      D <- merge(Dlist[[1]], y=Dlist[-1],add.cell.ids=paste0("S",1:length(Dlist)))
    }
  }else{
    gID <- Reduce(intersect,lapply(X,rownames))
    X <- lapply(X,function(x)return(x[gID,]))
    X <- Reduce(cbind,X)
    
    meta <- Reduce(rbind,meta)
    cID <- intersect(rownames(meta),colnames(X))
    message("\t\t",length(cID)," cells overlapped between expression and annotation")
    if(length(cID)<10){
      MsgExit(paste("Too few matching cell IDs"))
    }
    D <- CreateSeuratObject(counts=X[,cID],
                            project=basename(config$expression[i]),
                            meta.data=meta[cID,])
  }
  return(D)
}
processSCT <- function(D,batch=NULL){
  message("\tUsing SCT to normalize and scale the data ...")
  Dmedian <- median(colSums(D@assays$RNA@counts))
  if(!is.null(batch)) Dlist <- SplitObject(D,split.by=batch)
  else Dlist <- list(D)
  Dlist <- lapply(Dlist,function(one,medianUMI){
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
    D <- Dlist[[1]]
  }else{
    D <- merge(Dlist[[1]], y=Dlist[-1])
    VariableFeatures(D) <- SelectIntegrationFeatures(Dlist,nfeatures=3000)
  }
  return(D)
}
processLayout <- function(D,reduction,config){
  message("\tobtain the layout reduction")
  if(length(reduction)==0){
    stop("Empty reduction layout")
    if(length(VariableFeatures(D))<50) D <- FindVariableFeatures(D, nfeatures = 3000)
    D <- RunPCA(D, features = VariableFeatures(D))
    D <- RunUMAP(D,dims=1:50)
  }else{
    commKey <- Reduce(intersect,lapply(reduction,colnames))
    if(length(config$expression)==length(reduction)){
      message("\t\tmatching cell ID ...")
      cID <- rownames(D@meta.data)
      for(i in 1:length(reduction)){
        cID1 <- rownames(reduction[[i]])
        cID2 <- paste0("S",i,"_",rownames(reduction[[i]]))
        if(sum(cID1%in%cID)<sum(cID2%in%cID)){
          rownames(reduction[[i]]) <- cID2
          message("\t\t\tupdating cell name to:": paste(head(cID2),collapse=", "))
        }
      }
    }
    reduction <- lapply(reduction,function(x)return(x[,commKey]))
    reduction <- Reduce(rbind,reduction)
    
    for(one in names(config$reduction)){
      if(one=="files") next
      if(sum(config$reduction[[one]]%in%commKey)<2){
        message("\t\t",one,": Less than 2 dimension --- skip ---")
        next
      }
      oneReduction <- as.matrix(reduction[rownames(reduction)%in%rownames(D@meta.data),
                                    colnames(reduction)%in%config$reduction[[one]]])
      message("\t\t",one,": ",nrow(oneReduction)," cells overlap with expression and annotation")
      if(nrow(oneReduction)<10){
        message("\t\t\t--- skip ---")
        next
      }
      if(sum(!rownames(D@meta.data)%in%rownames(oneReduction))>0){
        message("adding 0 for missing cells in reduction ",one)
        cID <- rownames(D@meta.data)[!rownames(D@meta.data)%in%rownames(oneReduction)]
        oneReduction <- rbind(oneReduction,matrix(0,nrow=length(cID),ncol=ncol(oneReduction),
                                            dimnames=list(cID,colnames(oneReduction))))
      }
      colnames(oneReduction) <- paste(one,1:ncol(oneReduction),sep="_")
      D[[one]] <- CreateDimReducObject(embeddings=oneReduction,
                                       key=paste0(one,"_"),
                                       assay=D@active.assay)
    }
  }
  if(length(D@reductions)==0)
    stop("No reduction layout is extracted successfully")
  return(D)
}
## main ---------
main <- function(){
  args = commandArgs()
  strPipePath <<- normalizePath(dirname(dirname(sapply(strsplit(grep("file=",args,value=T),"="),tail,1))))
  MsgInit()
  strInput <- commandArgs(trailingOnly=T)
  if(length(strInput)==0){
    MsgHelp()
  }
  if(file.exists(strInput[1]) && !dir.exists(strInput[1])){
    strH5ad = process2H5ad(strInput[1])
    MsgComplete(strH5ad)
  }else{
    initConfig(strInput[1])
  }
  MsgPower()
}



main()
