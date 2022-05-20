
PKGloading <- function(){
  require(Seurat)
  require(RcppCNPy)
  require(readr)
  require(rhdf5)
  require(Matrix)
  options(stringsAsFactors = FALSE)
  
  rlang::env_unlock(env = asNamespace('base'))
  rlang::env_binding_unlock(env = asNamespace('base'))
  message <<- function(...,domain = NULL, appendLF = TRUE){
    cat(...,"\n",sep="")
  }
  rlang::env_binding_lock(env = asNamespace('base'))
  rlang::env_lock(asNamespace('base'))
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
  suppressMessages(suppressWarnings(PKGloading()))
  batchKey="library_id" #"batch"#
  args = commandArgs(trailingOnly=TRUE)
  if(length(args)<3) stop("Path to 3 files are required!")
  strH5ad <- args[1]
  ## expr
  #adata_expr = npyLoad(args[1])
  #adata_expr = t(adata_expr)
  ## batch
  #adata_batch = read_csv(args[2])
  #adata_batch <- as.data.frame(adata_batch)
  ## gene names
  #adata_var = read_csv(args[3])
  
  #head(adata_batch)
  #head(adata_var)
  #colnames(adata_expr) <- adata_batch$index
  #row.names(adata_expr) <- adata_var$index
  #row.names(adata_batch) <- adata_batch$index
  
  print(Sys.time())
  print(system.time({
    #concat <- CreateSeuratObject(counts = adata_expr, meta.data = adata_batch)
    concat <- CreateSeuratObject(counts=getX(strH5ad),
                                 project="SCT",
                                 meta.data=getobs(strH5ad))
    # split cells based on project id
    concatlst <- SplitObject(concat, split.by = batchKey)
    concat.anchors <- FindIntegrationAnchors(object.list = concatlst)
    concat.integrated <- IntegrateData(anchorset = concat.anchors, preserve.order = TRUE)
    #saveRDS(concat.integrated, file = "working_data/concat_seurat3_integrated.RDS")
    
    ## return integrated expr matrix
    adata_expr_integrated = data.frame(concat.integrated@assays$integrated@data)
    # adata_expr_integrated_genes = row.names(adata_expr_integrated)
    write_csv(adata_expr_integrated, file = args[2])
    
    ## also export gene names, cell names are in the matrix
    write_csv(data.frame(concat.integrated@assays$integrated@var.features), file = args[3])
  }))
  print(Sys.time())
}

main()