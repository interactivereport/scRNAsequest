PKGloading <- function(){
  require(rhdf5)
  require(Matrix)
}
suppressMessages(suppressWarnings(PKGloading()))

getobs <- function(strH5ad){
  message("\tobtainning obs ...")
  obs <- h5read(strH5ad,"obs")
  sel <- names(obs)[sapply(obs,function(x)return(is.null(names(x))))&!grepl("^_|index$|^barcode$",names(obs))]
  meta <- do.call(cbind.data.frame, obs[sel])
  #meta <- do.call(cbind.data.frame, obs[grep("^_",names(obs),invert=T)])
  #dimnames(meta) <- list(obs[["_index"]],grep("^_",names(obs),invert=T,value=T))
  rownames(meta) <- obs[[grep("index$|^barcode$",names(obs))]]
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
    if(sum(obs[[one]]$codes<0)>0)
      stop(paste0("NaN or NA found in ",one,", please use a string such as 'unknown' or 'undetermined' instead!"))
    meta[[one]] <- obs[[one]]$categories[1+obs[[one]]$codes]
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
getX <- function(strH5ad,useRaw=T){
  message("\tobtainning X ...")
  keys <- h5ls(strH5ad)
  message("\t\textracting counts")
  if(useRaw && sum(grepl("/raw/X",keys$group))>0){
    message("\t\t\tFound .raw.X")
    X <- h5read(strH5ad,"/raw/X")
  }else{
    X <- h5read(strH5ad,"X")
  }
  if(length(X$data)>(2^31-1))
    stop(paste("Max elements in sparse matrix is 2^31-1, input:",length(X$data)))
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
  if((max(X$indices)+1)==length(gID) || (length(X$indptr)-1)==length(cID)){ # CSR sparse matrix
    M <- sparseMatrix(i=X$indices+1,p=X$indptr,x=as.numeric(X$data),
                      dims=c(length(gID),length(cID)),
                      dimnames=list(gID,cID))
  }else if((max(X$indices)+1)==length(cID) || (length(X$indptr)-1)==length(gID)){#CSC sparse matrix
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