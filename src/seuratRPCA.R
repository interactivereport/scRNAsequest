
PKGloading <- function(){
  require(Seurat)
  require(peakRAM)
  #require(RcppCNPy)
  #require(readr)
  options(stringsAsFactors = FALSE)
}

runRPCA <- function(strH5ad,batch,strOut,cluster_method,clusterResolution){
  D <- CreateSeuratObject(counts=getX(strH5ad),
                          project="RPCA",
                          meta.data=getobs(strH5ad))
  
  D[["RNA"]] <- split(D[["RNA"]],f=unlist(D[[batch]],use.names=F))
  message("\tSCT ...")
  D <- SCTransform(D,vst.flavor="v2",
                   return.only.var.genes=F,
                   verbose=F)
  if(length(D[["RNA"]]@layers)>1){
    D <- RunPCA(D,npcs=50,verbose = FALSE)
    D <- IntegrateLayers(D,method=RPCAIntegration,
                         orig.reduction = "pca",
                         normalization.method="SCT",
                         verbose = F)
  }
  D <- RunPCA(D,npcs=50,verbose = FALSE)
  D <- RunUMAP(D,dims=1:50,
               reduction="pca",
               umap.method="umap-learn",
               metric = "correlation",
               verbose=F)
  D <- FindNeighbors(D,dims = 1:50,verbose=F)
  
  message("Clustering ",cluster_method," (",clusterResolution,") ...")
  if(grepl("Leiden",cluster_method,ignore.case=T)){
    D <- FindClusters(D,verbose = FALSE,resolution=clusterResolution,method="igraph",algorithm=4)
  }else{
    D <- FindClusters(D,verbose = FALSE,resolution=clusterResolution)#default Algorithm: 1 = original Louvain algorithm
  }
  # prepare saving
  X <- cbind(D[["seurat_clusters"]],
             D[["pca"]]@cell.embeddings,
             D[["umap"]]@cell.embeddings)
  colnames(X) <- gsub("^PC","pca",colnames(X))
  colnames(X) <- gsub("^UMAP","umap",colnames(X))
  colnames(X) <- gsub("seurat_clusters","seuratRPCA_cluster",colnames(X))
  X <- cbind(cID=rownames(X),X)
  write_csv(X,strOut)
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
  print(peakRam(
    runRPCA(strH5ad,batchKey,args[2],cluster_method,clusterResolution)))
}

main()