PKGloading <- function(){
  require(rhdf5)
  require(Matrix)
}
PKGloading()
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
