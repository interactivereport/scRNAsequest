
PKGloading <- function(){
  require(Seurat)
  require(RcppCNPy)
  require(readr)
  options(stringsAsFactors = FALSE)
  if(F){
    rlang::env_unlock(env = asNamespace('base'))
    rlang::env_binding_unlock(env = asNamespace('base'))
    message <<- function(...,domain = NULL, appendLF = TRUE){
      cat(paste(...,collapse=", "),"\n",sep="")
    }
    rlang::env_binding_lock(env = asNamespace('base'))
    rlang::env_lock(asNamespace('base'))
  }
  
}

main <- function(){
  suppressMessages(suppressWarnings(PKGloading()))
  batchKey="library_id" #"batch"#
  args = commandArgs(trailingOnly=TRUE)
  if(length(args)<2) stop("Path to 2 files are required!")
  strH5ad <- args[1]
  clusterResolution <- as.numeric(args[3])
  cluster_method <- args[4]
  source(paste0(dirname(gsub("--file=","",grep("file=",commandArgs(),value=T))),"/readH5ad.R"))
  print(system.time({
    D <- CreateSeuratObject(counts=getX(strH5ad),
                                 project="SCT",
                                 meta.data=getobs(strH5ad))
    # split cells based on project id
    Dlist <- SplitObject(D, split.by = batchKey)
    ## following https://satijalab.org/seurat/articles/integration_rpca.html (with SCTransform)
    #Dmedian <- median(colSums(D@assays$RNA@counts))
    rm(D)
    gc()
    Dlist <- sapply(Dlist,function(one,medianUMI){
      message("\t\tSCT one ...")
      #after checking/testing, the above would return proper SCT nomralized data
      #/edgehpc/dept/compbio/users/zouyang/process/PRJNA544731/src/SCT_scale_batch.ipynb
      #https://github.com/satijalab/sctransform/issues/128
      return(suppressMessages(suppressWarnings(
        SCTransform(one,method = 'glmGamPoi',
                    new.assay.name="SCT",
                    return.only.var.genes = FALSE,
                    #scale_factor=medianUMI,
                    verbose = FALSE)
      )))
    },Dmedian)
    if(length(Dlist)>1){
      features <- SelectIntegrationFeatures(object.list = Dlist, nfeatures = 3000)
      Dlist <- PrepSCTIntegration(object.list = Dlist, anchor.features = features)
      Dlist <- lapply(Dlist, FUN = RunPCA, features = features)
      
      anchors <- FindIntegrationAnchors(Dlist, normalization.method = "SCT",
                                        anchor.features = features, dims = 1:50,
                                        reduction = "rpca", k.anchor = 20)
      rm(Dlist)
      gc()
      Dsct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:50)
      rm(anchors)
      gc()
    }else{
      Dsct <- Dlist[[1]]
    }
    
    Dsct <- RunPCA(Dsct, verbose = FALSE)
    Dsct <- RunUMAP(Dsct, reduction = "pca", dims = 1:50)
    Dsct <- FindNeighbors(Dsct,dims = 1:50)
    
    message("Clustering ",cluster_method," (",clusterResolution,") ...")
    if(grepl("Leiden",cluster_method,ignore.case=T)){
      Dsct <- FindClusters(Dsct,verbose = FALSE,resolution=clusterResolution,algprithem=4)
    }else{
      Dsct <- FindClusters(Dsct,verbose = FALSE,resolution=clusterResolution)#default Algorithm: 1 = original Louvain algorithm
    }
    # prepare saving
    D <- cbind(Dsct[["seurat_clusters"]],
               Dsct[["pca"]]@cell.embeddings,
               Dsct[["umap"]]@cell.embeddings)
    colnames(D) <- gsub("^PC","pca",colnames(D))
    colnames(D) <- gsub("^UMAP","umap",colnames(D))
    colnames(D) <- gsub("seurat_clusters","seuratRPCA_cluster",colnames(D))
    D <- cbind(cID=rownames(D),D)
    write_csv(D,args[2])
  }))
}

main()