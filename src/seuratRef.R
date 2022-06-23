PKGloading <- function(){
  require(Seurat)
  require(SeuratObject)
  require(SeuratDisk)
  require(sctransform)
  require(rhdf5)
  require(Matrix)
  options(future.globals.maxsize=3145728000)
}
mapOne <- function(oneD){
  
}
checkRef <- function(ref,config){
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
  return(config)
}
processSCTref <- function(strH5ad,batch,config,strOut){
  if(grepl("^http",config$ref_file))
    reference <- readRDS(url(config$ref_file))
  else
    reference <- readRDS(config$ref_file)
  config <- checkRef(reference,config)
  D <- CreateSeuratObject(counts=getX(strH5ad),
                          project="SCT",
                          meta.data=getobs(strH5ad))
  Dmedian <- median(colSums(D@assays$RNA@counts))
  Dlist <- SplitObject(D,split.by=batch)
  Dlist <- sapply(Dlist,function(one,medianUMI){
    message("\t\tmapping one ...")
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
        mapping.score.k = 100
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
  
  if(length(Dlist)==1){
    SCT <- Dlist[[1]]
  }else{
    SCT <- merge(Dlist[[1]], y=Dlist[-1],merge.dr = names(Dlist[[1]]@reductions))
  }
  
  ## save layout and annotation (not iterative mapping yet!)
  meta <- SCT@meta.data[,grep("^predicted",colnames(SCT@meta.data)),drop=F]
  layout <- c()
  for(i in names(SCT@reductions)){
    layout <- cbind(layout,SCT@reductions[[i]]@cell.embeddings[,1:2])
  }
  #colnames(layout) <- paste0("seuratRef_",colnames(layout))
  colnames(layout) <- sapply(strsplit(colnames(layout),"_"),function(x){
    ix <- grep("umap",x,ignore.case = T)
    if(length(ix)>0) x[ix] <- "umap"
    return(paste(x,collapse="_"))
    })#paste0("seuratRef_",colnames(layout))
  X <- cbind(meta,layout)
  X <- cbind(cID=rownames(X),X)
  data.table::fwrite(X,strOut)
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
  if(file.exists(config$ref_name) && grepl("rds$",config$ref_name)){
    oneRef <- list(ref_file=config$ref_name)
  }
  else if(!is.null(config$ref_name) && config$ref_name%in%names(sysConfig)){
    oneRef <- sysConfig[[config$ref_name]]
  }else{
    stop(paste("unknown reference:",config$ref_name))
  }

  strOut <- args[3]
  if(length(args)>3) batchKey <- args[4]
  processSCTref(strH5ad,batchKey,oneRef,strOut)
  
}

main()