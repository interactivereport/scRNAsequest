PKGloading <- function(){
  require(rhdf5)
  require(Matrix)
  require(dplyr)
  require(BiocParallel)
}
suppressMessages(suppressWarnings(PKGloading()))

getobs <- function(strH5ad){
  message("\tobtainning obs ...")
  obs <- h5read(strH5ad,"obs")
  sel <- names(obs)[sapply(obs,function(x)return(is.null(names(x))))&!grepl("^_|index$|^barcode$",names(obs))]
  meta <- do.call(cbind.data.frame, obs[sel]) %>% mutate_all(~ifelse(is.nan(.),NA,.))
  #meta <- do.call(cbind.data.frame, obs[grep("^_",names(obs),invert=T)])
  #dimnames(meta) <- list(obs[["_index"]],grep("^_",names(obs),invert=T,value=T))
  if(ncol(meta)==0) meta <- data.frame(row.names=obs[[grep("index$|^barcode$",names(obs))]])
  else rownames(meta) <- obs[[grep("index$|^barcode$",names(obs))]]
  for(one in names(obs[["__categories"]])){
    if(min(meta[,one])<0){
      ann <- meta[,one]+1
      ann[ann<1] <- max(ann)+1
      annLable <- obs[["__categories"]][[one]]
      annLable <- c(annLable,"NAN")
      meta[,one] <- annLable[ann]
    }else{
      meta[,one] <- as.vector(obs[["__categories"]][[one]])[1+meta[,one]]
    }
  }
  # for anndata v 0.8
  for(one in names(obs)[sapply(obs,function(x)return(!is.null(names(x))))&!grepl("^_|index$",names(obs))]){
    if(sum(c("categories","codes")%in%names(obs[[one]]))!=2) next
    meta[[one]]<- rep(NA,length(obs[[one]]$codes))
    selNA <- obs[[one]]$codes<0
    meta[[one]][!selNA] <- obs[[one]]$categories[1+obs[[one]]$codes[!selNA]]
    #if(sum(obs[[one]]$codes<0)>0)
    #  stop(paste0("NaN or NA found in ",one,", please use a string such as 'unknown' or 'undetermined' instead!"))
    #meta[[one]] <- obs[[one]]$categories[1+obs[[one]]$codes]
  }
  return(meta)
}
getID <- function(strH5ad,keys,grp){
  if("_index" %in% keys$name[grepl(grp,keys$group)]){
    return(h5read(strH5ad,paste0(grp,"/_index")))
  }else if("index" %in% keys$name[grepl(grp,keys$group)]){
    return(h5read(strH5ad,paste0(grp,"/index")))
  }else if("barcode" %in% keys$name[grepl(grp,keys$group)]){
    return(h5read(strH5ad,paste0(grp,"/barcode")))
  }else if("feature_name" %in% keys$name[grepl(grp,keys$group)]){
    gName <- h5read(strH5ad,paste0(grp,"/feature_name"))
    if(sum(c("categories","codes")%in%names(gName))==2){
      gName <- gName$categories[1+gName$codes]
    }
    return(gName)
  }else{
    stop(paste("unknown adata format: Neither index or _index exists in group",grp))
  }
}
getBatchX <- function(strH5ad,batchID,useRaw=T,core=Inf){
  meta <- getobs(strH5ad)
  message("\tobtainning X ...")
  keys <- h5ls(strH5ad)
  message("\t\textracting counts by batch ID: ",batchID)
  if(min(diff(as.numeric(factor(meta[[batchID]],levels=unique(meta[[batchID]])))))<0)
    stop("Samples are NOT ordered by batch ID in the give h5ad!")
  message("\t\textracting gene name")
  if(useRaw && sum(grepl("/raw/var",keys$group))>0){
    message("\t\t\tFound .raw.var")
    gID <- getID(strH5ad,keys,"/raw/var") #h5read(strH5ad,"/raw/var/_index")
  }else{
    gID <- getID(strH5ad,keys,"/var") #h5read(strH5ad,"/var/_index")
  }
  message("\t\textracting cell name")
  if(useRaw && sum(grepl("/raw/obs",keys$group))>0){
    message("\t\tFound .raw.obs")
    cID <- getID(strH5ad,keys,"/raw/obs") #h5read(strH5ad,"/raw/obs/_index")
  }else{
    cID <- getID(strH5ad,keys,"/obs") #h5read(strH5ad,"/obs/_index")
  }
  gID <- as.vector(gID)
  cID <- as.vector(cID)
  if(useRaw && sum(grepl("/raw/X",keys$group))>0){
    message("\t\t\tFound .raw.X")
    selGroup <- "/raw/X"
  }else{
    selGroup <- "/X"
  }
  indptr <- h5read(strH5ad,paste0(selGroup,"/indptr"),bit64conversion="double")
  if((length(indptr)-1)!=length(cID))
    stop("Given h5ad is NOT CSC, cannot get batch X!")
  bNames <- unique(meta[[batchID]])#[1:5]
  X <- bplapply(setNames(bNames,bNames),
                function(one){
                  message("\t\tsample: ",one,"\t",which(bNames==one),"/",length(bNames),appendLF=F)
                  cIx <- seq_along(meta[[batchID]])[meta[[batchID]]==one]
                  iStart <- indptr[min(cIx)]+1
                  iEnd <- indptr[max(cIx)+1]
                  if((iEnd-iStart+1)>(2^31-1)) stop("Max elements in sparse matrix is more than 2^31-1")
                  x <- h5read(strH5ad,paste0(selGroup,"/data"),index=list(iStart:iEnd))
                  i <- h5read(strH5ad,paste0(selGroup,"/indices"),index=list(iStart:iEnd))+1
                  p <- indptr[min(cIx):(max(cIx)+1)]-indptr[min(cIx)]
                  M <- sparseMatrix(i=i,p=p,x=as.numeric(x),
                                    dims=c(length(gID),length(cIx)),
                                    dimnames=list(gID,cID[cIx]))
                  rm(x,i,p)
                  gc()
                  return(M)
                },BPPARAM = MulticoreParam(workers=min(core,length(bNames),max(1,parallelly::availableCores()-2)),
                                           tasks=length(bNames)))
  if(F){
    X <- lapply(setNames(unique(meta[[batchID]]),unique(meta[[batchID]])),
                function(one){
                  message("\t\tsample: ",one)
                  cIx <- seq_along(meta[[batchID]])[meta[[batchID]]==one]
                  iStart <- indptr[min(cIx)]+1
                  iEnd <- indptr[max(cIx)+1]
                  if((iEnd-iStart+1)>(2^31-1)) stop("Max elements in sparse matrix is 2^31-1")
                  x <- h5read(strH5ad,paste0(selGroup,"/data"),index=list(iStart:iEnd))
                  i <- h5read(strH5ad,paste0(selGroup,"/indices"),index=list(iStart:iEnd))+1
                  p <- indptr[min(cIx):(max(cIx)+1)]-indptr[min(cIx)]
                  M <- sparseMatrix(i=i,p=p,x=as.numeric(x),
                                    dims=c(length(gID),length(cIx)),
                                    dimnames=list(gID,cID[cIx]))
                  return(M)
                })
  }
  return(X)
}
getX <- function(strH5ad,batchID=NULL,useRaw=T,core=5){
  if(!is.null(batchID))
    return(getBatchX(strH5ad,batchID,useRaw,core))
  message("\tobtainning X ...")
  keys <- h5ls(strH5ad)
  message("\t\textracting counts")
  if(useRaw && sum(grepl("/raw/X",keys$group))>0){
    message("\t\t\tFound .raw.X")
  	Xdim <- as.numeric(keys[grepl("/raw/X",keys$group) & grepl("data",keys$name),'dim'])
    if(is.na(Xdim) || Xdim>(2^31-1))
      stop(paste("Max allowed element length in sparse matrix is 2^31-1, input:",Xdim))
    X <- h5read(strH5ad,"/raw/X")
  }else if(sum(grepl("/X",keys$group))>0){
  	Xdim <- as.numeric(keys[grepl("/X",keys$group) & grepl("data",keys$name),'dim'])
    if(is.na(Xdim) || Xdim>(2^31-1))
      stop(paste("Max allowed element length in sparse matrix is 2^31-1, input:",Xdim))
    X <- h5read(strH5ad,"X")
  }else{
  	stop("Sparse matrix is required for either .X or .raw.X!")
  }
  message("\t\textracting gene name")
  if(useRaw && sum(grepl("/raw/var",keys$group))>0){
    message("\t\t\tFound .raw.var")
    gID <- getID(strH5ad,keys,"/raw/var") #h5read(strH5ad,"/raw/var/_index")
  }else{
    gID <- getID(strH5ad,keys,"/var") #h5read(strH5ad,"/var/_index")
  }
  message("\t\textracting cell name")
  if(useRaw && sum(grepl("/raw/obs",keys$group))>0){
    message("\t\tFound .raw.obs")
    cID <- getID(strH5ad,keys,"/raw/obs") #h5read(strH5ad,"/raw/obs/_index")
  }else{
    cID <- getID(strH5ad,keys,"/obs") #h5read(strH5ad,"/obs/_index")
  }
  gID <- as.vector(gID)
  cID <- as.vector(cID)
  message("\t\tcreating sparse matrix")
  if((max(X$indices)+1)==length(gID) || (length(X$indptr)-1)==length(cID)){ # CSC sparse matrix
    M <- sparseMatrix(i=X$indices+1,p=X$indptr,x=as.numeric(X$data),
                      dims=c(length(gID),length(cID)),
                      dimnames=list(gID,cID))
  }else if((max(X$indices)+1)==length(cID) || (length(X$indptr)-1)==length(gID)){#CSR sparse matrix
    M <- sparseMatrix(j=X$indices+1,p=X$indptr,x=as.numeric(X$data),
                      dims=c(length(gID),length(cID)),
                      dimnames=list(gID,cID))
  }else{
    stop(paste("Error in reading X from h5ad:",strH5ad))
  }
  return(M)
}
getobsm <- function(strH5ad,key){
  AllK <- h5ls(strH5ad,recursive=2)
  k <- AllK[grepl("obsm",AllK[,1]),2]
  if(!key%in%k) return(NULL)
  X <- h5read(strH5ad,paste0("obsm/",key))
  colnames(X) <- getID(strH5ad,AllK,"/obs")    #h5read(strH5ad,"/obs/_index")
  return(t(X))
}
getobsmKey <- function(strH5ad){
  k <- h5ls(strH5ad,recursive=2)
  k <- k[grepl("obsm",k[,1]),2]
  return(k)
}