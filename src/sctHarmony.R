PKGloading <- function(){
  require(Seurat)
  require(harmony)
  options(stringsAsFactors = FALSE)
}

# following the discussion @https://github.com/immunogenomics/harmony/issues/41
processH5ad <- function(strH5ad,batch,strOut,bPrepSCT){
  assayName <- "sctHarmony"
  dimN <- 50
  X <- getX(strH5ad)
  gID <- setNames(rownames(X),gsub("_","-",rownames(X)))
  rownames(X) <- names(gID) #gsub("_","-",rownames(X))
  D <- CreateSeuratObject(counts=X,
                          project="sctHarmony",
                          meta.data=getobs(strH5ad))
  rm(X)
  gc()
  Dlist <- SplitObject(D,split.by=batch)
  Dlist <- sapply(Dlist,function(one){
    bID <- one@meta.data[1,batch]
    message("\t\t",bID)
    one <- suppressMessages(suppressWarnings(
      SCTransform(one,method = 'glmGamPoi',
                  new.assay.name=assayName,
                  return.only.var.genes = FALSE,
                  verbose = FALSE)
    ))
    one <- FindVariableFeatures(one, selection.method = "vst", nfeatures = 3000)
    return(one)
  })
  
  selFN <- 3000
  features <- SelectIntegrationFeatures(Dlist,nfeatures=selFN)
  D <- merge(Dlist[[1]], y=Dlist[-1],project="sctHarmony")
  rm(Dlist)
  gc()
  DefaultAssay(D) <- assayName
  VariableFeatures(D) <- features
  #D <- PrepSCTFindMarkers(D)
  D <- RunPCA(D, npcs = dimN)
  D <- RunHarmony(D,
                      assay.use=assayName,
                      reduction = 'pca',
                      group.by.vars = batch)
  
  D <- FindNeighbors(D, dims = 1:dimN,reduction="harmony")
  D <- RunSPCA(D,reduction.key="harmonySPC",
               graph = paste0(assayName,'_snn'))
  D <- RunUMAP(D,reduction="harmony",
               dims = 1:dimN)
  D <- RunTSNE(D,reduction="harmony",
               dims = 1:dimN)
  D <- FindClusters(object = D)

  X <- cbind(cID=row.names(D@meta.data),
             D@reductions$harmony@cell.embeddings,
             D@reductions$umap@cell.embeddings,
             D@reductions$tsne@cell.embeddings,
             D@meta.data[,c("seurat_clusters",batch)])
  colnames(X) <- gsub("UMAP","umap",gsub("harmony","pca",gsub("seurat_clusters","sctHarmony_cluster",colnames(X))))
  #save(X,D,file=gsub("csv$","RData",strOut))
  data.table::fwrite(X,strOut)
}

main <- function(){
  suppressMessages(suppressWarnings(PKGloading()))
  batchKey="library_id" #"batch"
  args = commandArgs(trailingOnly=TRUE)
  if(length(args)<3) stop("Path to h5ad file, the output file and config file are required!")
  strH5ad <- args[1]
  if(!file.exists(strH5ad)) stop(paste0("H5ad file (",strH5ad,") does not exist!"))
  strOut <- args[2]
  strConfig <- args[3]
  if(!file.exists(strConfig)){
    stop(paste0("Config file (",strConfig,") does not exist!"))
  }else{
    config <- yaml::read_yaml(strConfig)
  }
  if(length(args)>3) batchKey <- args[4]
  
  
  source(paste0(dirname(gsub("--file=","",grep("file=",commandArgs(),value=T))),"/readH5ad.R"))
  print(system.time(processH5ad(strH5ad,batchKey,strOut,config$PrepSCTFindMarkers)))
}

main()