PKGloading <- function(){
  require(Seurat)
  require(peakRAM)
  #require(RcppCNPy)
  #require(readr)
  options(stringsAsFactors = FALSE,
          future.globals.maxSize=50*1024^3,
          future.seed=TRUE)
}
runRPCA <- function(strH5ad,batch,strOut,cluster_method,clusterResolution,subcore=5,geneN=3000){
  D <- getSCTransform(strH5ad,batch,subcore,geneN=geneN,sctAssay=T,only.var.genes=T)
  plan("multisession",workers=subcore)
  message("\tPCA ...")
  D <- RunPCA(D,npcs=50,verbose = FALSE)
  dr <- "pca"
  if(length(D[["RNA"]]@layers)>1){
    message("\tRPCA integration ...")
    dr <- "rpca"
    D <- IntegrateLayers(D,method=RPCAIntegration,
                         orig.reduction = "pca",
                         normalization.method="SCT",
                         new.reduction=dr,
                         verbose = F)
  }
  message("\tFind neighbors ...")
  D <- FindNeighbors(D,reduction=dr,dims = 1:50,verbose=F)
  message("\tClustering ",cluster_method," (",clusterResolution,") ...")
  if(grepl("Leiden",cluster_method,ignore.case=T)){
    D <- FindClusters(D,verbose = FALSE,resolution=clusterResolution,method="igraph",algorithm=4)
  }else{
    D <- FindClusters(D,verbose = FALSE,resolution=clusterResolution)#default Algorithm: 1 = original Louvain algorithm
  }
  message("\tUMAP ...")
  D <- suppressWarnings(suppressMessages(
    RunUMAP(D,dims=1:50,
               reduction=dr,
               umap.method="umap-learn",
               metric = "correlation",
               verbose=F)))
  # prepare saving
  X <- cbind(D[["seurat_clusters"]],
             D[["pca"]]@cell.embeddings,
             D[["umap"]]@cell.embeddings)
  colnames(X) <- gsub("^PC","pca",colnames(X))
  colnames(X) <- gsub("^UMAP","umap",colnames(X))
  colnames(X) <- gsub("seurat_clusters","seuratRPCA_cluster",colnames(X))
  #X <- cbind(cID=rownames(X),X)
  data.table::fwrite(cbind(cID=rownames(X),X),strOut)
}
main <- function(){
  suppressMessages(suppressWarnings(PKGloading()))
  batchKey="library_id" #"batch"#
  args = commandArgs(trailingOnly=TRUE)
  if(length(args)<2) stop("Path to 2 files are required!")
  strH5ad <- args[1]
  strConfig <- args[2]
  if(!file.exists(strConfig)){
  	stop(paste0("Config file (",strConfig,") does not exist!"))
  }else{
  	config <- yaml::read_yaml(strConfig)
  }
  clusterResolution <- ifelse(is.null(config$clusterResolution),0.8,config$clusterResolution)
  cluster_method <- ifelse(is.null(config$clusterMethod),"Louvain",config$clusterMethod)
  subcore <- ifelse(is.null(config$subprocess),5,config$subprocess)
  geneN <- ifelse(is.null(config$harmonyBatchGene),3000,config$harmonyBatchGene)
  source(paste0(dirname(gsub("--file=","",grep("file=",commandArgs(),value=T))),"/SCTransform.R"),chdir=T)
  print(peakRAM(
    runRPCA(strH5ad,batchKey,args[2],cluster_method,clusterResolution,subcore,geneN)))
}

main()