PKGloading <- function(){
  require(Seurat)
  require(SeuratObject)
  require(SeuratDisk)
  require(sctransform)
  require(ggplot2)
  require(reshape2)
  require(peakRAM)
  options(future.globals.maxSize=8000*1024^2) 
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
processSCTrefOne <- function(strH5ad,batch,config){
  if(grepl("^http",config$ref_file))
    reference <- readRDS(url(config$ref_file))
  else
    reference <- readRDS(config$ref_file)
  res <- updateRef(reference,config)
  config <- res$config
  reference <- res$ref
  rm(res)
  print(peakRAM({
    message("\t\tCreating seurat object ...")
    D <- CreateSeuratObject(counts=getX(strH5ad),
                            project="SCT",
                            meta.data=getobs(strH5ad))
    cellN <- dim(D)[2]
    Dmedian <- median(colSums(D@assays$RNA@counts))
    Dlist <- SplitObject(D,split.by=batch)
    #message("\tmemory usage before mapping: ",sum(sapply(ls(),function(x){object.size(get(x))})),"B for ",cellN," cells")
    rm(D)
    gc()
  }))
  print(peakRAM({
    Dlist <- sapply(Dlist,function(one,medianUMI){
      bID <- one@meta.data[1,batch]
      message("\t\tmapping ",bID)
      #after checking/testing, the below would return proper SCT nomralized data
      #/edgehpc/dept/compbio/users/zouyang/process/PRJNA544731/src/SCT_scale_batch.ipynb
      #https://github.com/satijalab/sctransform/issues/128
      oneD <- suppressMessages(suppressWarnings(
        SCTransform(one,method = 'glmGamPoi',
                    new.assay.name="SCT",
                    return.only.var.genes = FALSE,
                    scale_factor=medianUMI,
                    verbose = FALSE)
      ))
      anchors <- suppressMessages(suppressWarnings(
        FindTransferAnchors(
          reference = reference,
          query = oneD,
          k.filter = NA,
          reference.neighbors = config$ref_neighbors,
          reference.assay = config$ref_assay,
          query.assay = "SCT",
          reference.reduction = config$ref_reduction,
          normalization.method = "SCT",
          features = intersect(rownames(x = reference), VariableFeatures(object = oneD)),
          dims = 1:50,
          mapping.score.k = 100,
          verbose=F
        )
      ))
      oneD <- suppressMessages(suppressWarnings(
        MapQuery(
          anchorset = anchors,
          query = oneD,
          reference = reference,
          refdata = setNames(as.list(config$ref_label),config$ref_label),
          reference.reduction = config$ref_reduction, 
          reduction.model = config$ref_reduction.model,
          verbose=F
        )
      ))
      return(oneD)
    },Dmedian)
  }))
  #message("\tmemory usage after mapping: ",sum(sapply(ls(),function(x){object.size(get(x))})),"B for ",cellN," cells")
  print(peakRAM({
    meta = extractX(Dlist,batch)
  }))
  return(meta)
}
extractX <- function(Dlist,batch){
  meta <- NULL
  for(i in 1:length(Dlist)){
    message("\t\tmerging ",Dlist[[i]]@meta.data[1,batch])
    X <- Dlist[[i]]@meta.data[,grep("^predicted",colnames(Dlist[[i]]@meta.data)),drop=F]
    layout <- c()
    for(j in names(Dlist[[i]]@reductions)){
      layout <- cbind(layout,Dlist[[i]]@reductions[[j]]@cell.embeddings[,1:2])
    }
    colnames(layout) <- sapply(strsplit(colnames(layout),"_"),function(x){
      ix <- grep("umap",x,ignore.case = T)
      if(length(ix)>0) x[ix] <- "umap"
      return(paste(x,collapse="_"))
    })
    X <- cbind(X,layout)
    X <- cbind(cID=rownames(X),X)
    meta <- rbind(meta,X)
  }
  return(meta)
}
processSCTref <- function(strH5ad,batch,refList,strOut){
  D <- NULL
  for(one in names(refList)){
    message("*** Mapping ",one)
    X <- processSCTrefOne(strH5ad,batch,refList[[one]])
    if(one!=refTmpName){
      cNames <- gsub("^predicted",paste0("predicted.",one),colnames(X))
      cNames[!grepl("^cID$|^predicted",cNames)] <- 
        sapply(strsplit(cNames[!grepl("^cID$|^predicted",cNames)],"_"),
               function(x){return(paste(c(head(x,-1),one,tail(x,1)),collapse="_"))})
      colnames(X) <- cNames
    }
    if(is.null(D)) D <- X
    else{
      D <- merge(D,X,by="cID",all=T)
    }
  }
  data.table::fwrite(D,strOut)
  plotCrossAnno(D,names(refList),strOut)
}
plotCrossAnno <- function(D,refID,strOut){
  if(length(refID)==1) return()
  pdf(gsub("csv$","pdf",strOut),width=8,height=8)
  for(i in 1:(length(refID)-1)){
    for(j in (i+1):length(refID)){
      sel1 <- colnames(D)[!grepl("score$",colnames(D))&grepl(paste0("predicted.",refID[i]),colnames(D))]
      sel2 <- colnames(D)[!grepl("score$",colnames(D))&grepl(paste0("predicted.",refID[j]),colnames(D))]
      for(selA in sel1){
        for(selB in sel2){
          X <- melt(table(D[,c(selA,selB)]))
          print(ggplot(X,aes_string(selA,selB))+
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
  strRef <- list()
  for(one in names(refNames)){
    if(is.null(one) || nchar(one)==0)
      stop("missing ref name!")
    if(grepl(":",one) && is.null(refNames[[one]]))
      stop(paste("A space is required after ':' in define reference:",one))
    
    if(grepl("rds$",refNames[[one]]) && file.exists(refNames[[one]])){
      strRef <- c(strRef,setNames(list(list(ref_file=refNames[[one]])),one))
    }else if(!is.null(refNames[[one]]) && refNames[[one]]%in%names(sysConfig)){
      strRef <- c(strRef,setNames(list(sysConfig[[refNames[[one]]]]),one))
    }
  }
  if(length(strRef)==0){
    stop(paste("unknown reference:",paste(refNames,collapse=", ")))
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
  strRef <- checkRef(config$ref_name,sysConfig)

  strOut <- args[3]
  if(length(args)>3) batchKey <- args[4]
  
  source(paste0(dirname(selfPath),"/readH5ad.R"))
  processSCTref(strH5ad,batchKey,strRef,strOut)
}

main()