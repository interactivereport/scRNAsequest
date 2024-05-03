PKGloading <- function(){
  require(Seurat)
  require(SeuratObject)
  require(SeuratDisk)
  require(sctransform)
  require(ggplot2)
  require(reshape2)
  require(peakRAM)
  require(BiocParallel)
  require(Azimuth)
  options(future.globals.maxSize=8000*1024^2,check.names=F) 
  #SeuratData::AvailableData()
}
refTmpName <- paste(c("A",sample(c(LETTERS[1:20],letters[1:20],0:9),15,replace=T)),collapse="")
updateRef <- function(ref,config){
  checkList <- list(ref_assay="refAssay",
                    ref_neighbors="refdr.annoy.neighbors",
                    ref_reduction="refDR",
                    ref_reduction.model="refUMAP",
                    ref_label=NA)  
  for(one in names(checkList)){
    if(!one %in% names(config)){
      if(!is.na(checkList[[one]])){
        tmp <- tryCatch({
          ref[[checkList[[one]]]]
        },error=function(e){
          stop(paste(checkList[[one]],"cannot be found!\n\tPlease use scRef to generate reference!"))
        })
        rm("tmp")
        config[[one]] <- checkList[[one]]
      }else{
        config[[one]] <- names(ref@tools$AzimuthReference@colormap)
      }
    }
  }
  
  for(one in config[["ref_label"]])
    if(sum(nchar(levels(ref@meta.data[[one]]))==0)>0)
      levels(ref@meta.data[[one]])[nchar(levels(ref@meta.data[[one]]))==0] <- "unknown"

  return(list(config=config,ref=ref))
}
extractMap <- function(mapD,rName){
  X <- mapD@meta.data[,c("mapping.score",grep("^predicted",colnames(mapD@meta.data),value=T))]
  firstOne <- F
  for(oneReduc in names(mapD@reductions)){
    oneX <- mapD[[oneReduc]]@cell.embeddings
    if(!firstOne && ncol(oneX)>=50){
      colnames(oneX) <- paste0("pca_",1:ncol(oneX))
      firstOne <- T
    }else{
      colnames(oneX) <- paste0(ifelse(grepl("umap",oneReduc,ignore.case=T),"umap",oneReduc),"_",1:ncol(oneX))
    }
    X <- cbind(X,oneX)
  }
  if(rName!=refTmpName)
    colnames(X) <- paste0(rName,"_",colnames(X))
  X$cID <- rownames(X)
  return(X)
}
mergeMap <- function(mapL){
  D <- do.call(plyr::rbind.fill,mapL)
  meta <- data.frame(D[,grep("cID",colnames(D),invert=T,value=T)],
                     check.names=F)
  meta <- dplyr::bind_cols(lapply(meta,function(x){
    if(is.numeric(x)){
      x[is.na(x)] <- min(x,na.rm=T)
    }else{
      x <- as.character(x)
      x[is.na(x)] <- "missing"
    }
    return(x)
  }))
  #meta[colnames(meta)] <- do.call(cbind,)
  return(data.frame(row.names=D[,"cID"],
                    meta,check.names=F))
}
processAzimuth5 <- function(strH5ad,batch,refList,subCores){
  message("\tCreating seurat object ...")
  D <- CreateSeuratObject(counts=getX(strH5ad,batchID=batch,core=subCores),
                          project="seuratRef",
                          meta.data=getobs(strH5ad))
  message("\tMapping each reference ...")
  maplist <- bplapply(names(refList),function(oneRef){
    oneD <- suppressMessages(suppressWarnings(
      RunAzimuth(D,reference=refList[oneRef],verbose=F)
    ))
    #print(oneD)
    return(extractMap(oneD,oneRef))
  },
  BPPARAM = MulticoreParam(workers=min(subCores,length(refList),max(1,parallelly::availableCores()-2)),
                             tasks=length(refList)))
  meta <- do.call(merge,c(maplist,by="cID",all=T,sort=F))
  meta <- dplyr::bind_cols(lapply(meta,function(x){
    if(is.numeric(x)){
      x[is.na(x)] <- min(x,na.rm=T)
    }else{
      x <- as.character(x)
      x[is.na(x)] <- "missing"
    }
    return(x)
  }))
  return(data.frame(row.names=meta$cID,
                    meta[,-1],check.names=F))
}
processSCTref <- function(strH5ad,batch,refList){
  message("\tCreating seurat object ...")
  D <- CreateSeuratObject(counts=getX(strH5ad,batchID=batch),
                          project="seuratRef",
                          meta.data=getobs(strH5ad))
  Dlist <- SplitObject(D,split.by=batch)
  #Dmedian <- median(colSums(D@assays$RNA@counts))
  #cellN <- dim(D)[2]
  #message("\tmemory usage before mapping: ",sum(sapply(ls(),function(x){object.size(get(x))})),"B for ",cellN," cells")
  rm(D)
  gc()
  maplist <- bplapply(1:length(Dlist),function(i){#,medianUMI
    bID <- Dlist[[i]]@meta.data[1,batch]
    message("\tmapping ",bID)
    #after checking/testing, the below would return proper SCT nomralized data
    #/edgehpc/dept/compbio/users/zouyang/process/PRJNA544731/src/SCT_scale_batch.ipynb
    #https://github.com/satijalab/sctransform/issues/128
    #Azimuth v5 mapping
    meta <- data.frame(cID=colnames(Dlist[[i]]))
    for(oneRef in names(refList)){
      message("\t\t",oneRef)
      oneD <- tryCatch({
        suppressMessages(suppressWarnings(
          RunAzimuth(Dlist[[i]],reference=refList[oneRef],verbose=F)
        ))
      },error=function(err){
        return(NULL)
      })
      if(is.null(oneD)) next
      #saveRDS(oneD,gsub(".h5ad",paste0("_",oneRef,".rds"),strH5ad))
      oneM <- extractMap(oneD,oneRef)
      meta <- merge(meta,oneM,by="cID",all=T)
    }
    return(meta)
  },BPPARAM = MulticoreParam(workers=min(5,length(Dlist),max(1,parallelly::availableCores()-2)),
                             tasks=length(Dlist)))
  return(mergeMap(maplist))
}
plotCrossAnno <- function(D,refID,strOut){
  if(length(refID)==1) return()
  pdf(paste0(strOut,".pdf"),width=8,height=8)
  for(i in 1:(length(refID)-1)){
    for(j in (i+1):length(refID)){
      sel1 <- colnames(D)[!grepl("score$",colnames(D))&grepl(paste0("^",refID[i],"_predicted"),colnames(D))]
      sel2 <- colnames(D)[!grepl("score$",colnames(D))&grepl(paste0("^",refID[j],"_predicted"),colnames(D))]
      for(selA in sel1){
        for(selB in sel2){
          X <- melt(table(D[,c(selA,selB)]))
          print(ggplot(X,aes(.data[[selA]],.data[[selB]]))+
                  geom_tile(aes(fill = value),show.legend=F)+
                  geom_text(aes(label = value))+
                  scale_fill_gradient(low = "white", high = "red")+
                  theme_minimal()+
                  theme(axis.text.x = element_text(angle=90)))
        }
      }
    }
  }
  a <- dev.off()
}
checkRef <- function(refNames,sysConfig){
  if(!is.list(refNames)){
    refNames <- setNames(list(refNames),refTmpName)
  }
  
  seuratRef <- SeuratData::AvailableData()
  seuratRef <- seuratRef$Dataset[grepl("Azimuth Reference",seuratRef$Summary)]
  sysRef <- NULL
  strSysRef = file.path(sysConfig$refDir,"scRNAsequest_ref.csv")
  if(file.exists(strSysRef)){
    sysRef <- data.table::fread(strSysRef)$Dataset
  }
  strRef <- c()
  for(one in names(refNames)){
    if(is.null(one) || nchar(one)==0)
      stop("missing ref name!")
    if(grepl(":",one) && is.null(refNames[[one]]))
      stop("A space is required after ':' in defined reference: ",one)
    
    if(!is.null(refNames[[one]]) && refNames[[one]]%in%sysRef)
       refNames[[one]] <- file.path(sysConfig$refDir,refNames[[one]])
    if(!is.null(refNames[[one]]) && (refNames[[one]]%in%seuratRef || 
                                     (dir.exists(refNames[[one]]) && 
                                      file.exists(file.path(refNames[[one]],"ref.Rds")) &&
                                      file.exists(file.path(refNames[[one]],"idx.annoy"))))){
      strRef <- append(strRef,setNames(refNames[[one]],one))
    }else{
      stop("Unknown reference: ",refNames[[one]],"\n\tIf system reference was used, please contact Admin!")
    }

  }
  if(length(strRef)==0){
    stop("unknown reference: ",paste(refNames,collapse="; "))
  }
  return(strRef)
}
main <- function(){
  selfPath <- gsub("^--file=","",grep("^--file=",commandArgs(),value=T)[1])
  suppressMessages(suppressWarnings(PKGloading()))
  batchKey="library_id" #"batch"
  args = commandArgs(trailingOnly=TRUE)
  if(length(args)<2) stop("Path to h5ad file and the output file are required!")
  strH5ad <- args[1]
  if(!file.exists(strH5ad)) stop(paste0("H5ad file (",strH5ad,") does not exist!"))
  strConfig <- args[2]
  if(!file.exists(strConfig)) stop(paste0("H5ad file (",strConfig,") does not exist!"))
  
  config <- yaml::read_yaml(strConfig)
  sysConfig <- yaml::read_yaml(paste0(dirname(selfPath),"/sys.yml"))
  refList <- checkRef(config$ref_name,sysConfig)

  strOut <- args[3]
  subCores <- 5
  if(length(args)>3) subCores <- as.numeric(args[4])
  if(length(args)>4) batchKey <- args[5]
  
  source(paste0(dirname(selfPath),"/readH5ad.R"))
  #print(peakRAM(D <- processSCTref(strH5ad,batchKey,refList)))
  print(peakRAM(D <- processAzimuth5(strH5ad,batchKey,refList,
                                     ifelse(is.null(config$subprocess),5,config$subprocess))))
  saveRDS(D,strOut)
  plotCrossAnno(D,names(refList),strOut)
}

main()