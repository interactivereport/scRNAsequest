#https://github.com/immunogenomics/harmony/issues/41
PKGloading <- function(){
  require(Seurat)
  require(harmony)
  require(peakRAM)
  require(BiocParallel)#parallelly::availableCores()-2
  #register(MulticoreParam())#parallelly::availableCores()-2
  require(future)
  options(stringsAsFactors = FALSE,
          future.globals.maxSize=12*1024^3,
          future.seed=TRUE)
  plan("multicore")
}
# using seurat5
processV5 <- function(strH5ad,batch,strOut,geneN,core=10){
  assayName <- "sctHarmony"
  dimN <- 50
  X <- getX(strH5ad,batchID=batch,core=core)
  gID <- setNames(rownames(X[[1]]),gsub("_","-",rownames(X[[1]])))
  X <- lapply(X,function(one){
    rownames(one) <- names(gID)
    return(one)
  })
  meta <- getobs(strH5ad)
  Dlist <- bplapply(names(X),
                function(id){
                  message("\tSCT ",id,"\t",which(names(X)==id),"/",length(X),appendLF=F)
                  one <- tryCatch({
                    DD <- CreateSeuratObject(counts=X[id],
                                            project=assayName,
                                            meta.data=meta[colnames(X[[id]]),])
                    SCTransform(DD,vst.flavor="v2",
                                new.assay.name=assayName,
                                return.only.var.genes=FALSE,
                                verbose=FALSE)
                  },error=function(ee){
                    message("\t\tError: ",conditionMessage(ee))
                    DD <- SCTransform(D,vst.flavor="v1",
                                     return.only.var.genes = FALSE,
                                     new.assay.name=assayName,
                                     verbose=FALSE)
                    return(DD)
                    })
                  rm(DD)
                  gc()
                  return(one)
                },
                BPPARAM = MulticoreParam(workers=min(core,length(X),max(1,parallelly::availableCores()-2)),
                                         tasks=length(X)))
  rm(X)
  gc()
  strHVG <- paste0(dirname(strH5ad),"/hvg.csv")
  if(file.exists(strHVG)){
    message("\tUsing Highly Variable Features from scanpy.highly_variable_genes")
    hvg <- data.table::fread(strHVG,header=T)
    selGene <- intersect(unlist(hvg[[1]]),Reduce(intersect,lapply(Dlist,function(one)return(rownames(one[[assayName]])))))
  }else{
    message("\tFinding Highly Variable Features by Seurat.SelectIntegrationFeatures ...")
    selGene <- SelectIntegrationFeatures(Dlist,nfeatures=geneN)
  }
  message("\t\t",length(selGene)," features")
  if(length(selGene)<100) stop(paste0("Two few features (",length(selGene),")! Please increase the number of harmonyBatchGene (remove tmp folder)!"))
  #saveRDS(Dlist,file="sctHarmony.rds")
  message("\tsaving ...")
  D <- NULL
  for(i in 1:length(Dlist)){
    message("\t\t",Dlist[[1]]@meta.data[1,batch])
    oneD <- data.frame(t(Dlist[[1]]@assays[[assayName]]@scale.data[selGene,]),check.names=F)
    if(is.null(D)) D <- oneD
    else D <- rbind(D,oneD)
    Dlist[[1]] <- NULL
    gc()
  }
  saveRDS(D,strOut)
}

# following the discussion @https://github.com/immunogenomics/harmony/issues/41
processSCT <- function(strH5ad,batch,strOut,geneN){
  assayName <- "sctHarmony"
  dimN <- 50
  X <- getX(strH5ad)
  rownames(X) <- gsub("_","-",rownames(X))
  D <- CreateSeuratObject(counts=X,
                          project=assayName,
                          meta.data=getobs(strH5ad))
  cellN <- dim(D)[2]
  rm(X)
  gc()
  Dlist <- SplitObject(D,split.by=batch)
  #message("memory usage before SCT: ",sum(sapply(ls(),function(x){object.size(get(x))})),"B for ",cellN," cells")
  rm(D)
  gc()
  message("\tSCT ",length(Dlist)," samples ...")
  if(is.null(geneN)) geneN <- 3000
  if(F){
    Dlist <- lapply(Dlist,function(one){
      bID <- one@meta.data[1,batch]
      message("\t\t",bID)
      one <- suppressMessages(suppressWarnings(
        SCTransform(one,method = 'glmGamPoi',
                    new.assay.name=assayName,
                    return.only.var.genes = FALSE,
                    verbose = FALSE)
      ))
      #one <- FindVariableFeatures(one, selection.method = "vst", nfeatures = 5000,verbose=F)
      return(one)
    })#,BPPARAM = MulticoreParam(workers=2)
  }else{
    Dlist <- bplapply(1:length(Dlist),function(i){
      bID <- Dlist[[i]]@meta.data[1,batch]
      message("\t\t",bID)
      oneD <- tryCatch({
        suppressMessages(suppressWarnings(
          SCTransform(Dlist[[i]],method = 'glmGamPoi',
                      new.assay.name=assayName,
                      return.only.var.genes = FALSE,
                      verbose = FALSE)
        ))
      },error=function(cond){
        message("\t\twithout glmGamPoi for ",bID)
        return(suppressMessages(suppressWarnings(
          SCTransform(Dlist[[i]],
                      new.assay.name=assayName,
                      return.only.var.genes = FALSE,
                      verbose = FALSE)
        )))
      })
      oneD <- FindVariableFeatures(oneD, selection.method = "vst", nfeatures = geneN,verbose=F)
      return(oneD)
    },BPPARAM = MulticoreParam(workers=min(5,length(Dlist),max(1,parallelly::availableCores()-2)),
                               tasks=length(Dlist)))#min(length(Dlist),parallelly::availableCores()-2)
  }
  strHVG <- paste0(dirname(strH5ad),"/hvg.csv")
  if(file.exists(strHVG)){
    message("\tUsing Highly Variable Features from scanpy.highly_variable_genes")
    hvg <- data.table::fread(strHVG,header=T)
    selGene <- intersect(unlist(hvg[[1]]),Reduce(intersect,lapply(Dlist,function(one)return(rownames(one[[assayName]])))))
  }else{
    message("\tFinding Highly Variable Features by Seurat.SelectIntegrationFeatures ...")
    selGene <- SelectIntegrationFeatures(Dlist,nfeatures=geneN)
  }
  message("\t\t",length(selGene)," features")
  if(length(selGene)<100) stop(paste0("Two few features (",length(selGene),")! Please increase the number of harmonyBatchGene (remove tmp folder)!"))
  #saveRDS(Dlist,file="sctHarmony.rds")
  message("\tsaving ...")
  D <- NULL
  for(i in 1:length(Dlist)){
    message("\t\t",Dlist[[1]]@meta.data[1,batch])
    oneD <- data.frame(t(Dlist[[1]]@assays[[assayName]]@scale.data[selGene,]),check.names=F)
    if(is.null(D)) D <- oneD
    else D <- rbind(D,oneD)
    Dlist[[1]] <- NULL
  }
  rm(Dlist)
  saveRDS(D,strOut)
}
processPCA <- function(strPCA,strOut,batch,clusterResolution,cluster_method){
  message("starting Harmony ...")
  PCA <- getobsm(strPCA,"X_pca")
  meta <- getobs(strPCA)
  if(length(unique(meta[[batch]]))>1){
    PCA <- harmony::RunHarmony(PCA,meta,batch,do_pca=FALSE,verbose=FALSE)
  }
  message("Finishing Harmony ...")
  dimN <- ncol(PCA)
  colnames(PCA) <- paste0("pca_",1:dimN)
  D <- CreateSeuratObject(counts=data.frame(matrix(0,nrow=2,ncol=nrow(PCA),
                                                   dimnames=list(paste0("G",1:2),rownames(PCA))),
                                            check.names=F),
                          meta.data=meta)
  D[['harmony']] <- CreateDimReducObject(embeddings=PCA,key="pca_",assay="RNA")
  #saveRDS(D,"sctHarmony.rds")
  message("Find Neighbor ...")
  D <- FindNeighbors(D, dims = 1:dimN,reduction="harmony",verbose = FALSE)
  message("UMAP ...")
  D <- RunUMAP(D,reduction="harmony",dims = 1:dimN,
               umap.method ="umap-learn",metric = "correlation",
               verbose = FALSE)
  #message("tSNE ...")
  #tryCatch({D <- RunTSNE(D,reduction="harmony",dims = 1:dimN,verbose = FALSE)},
  #         error=function(cond){
  #           message("\terror:")
  #           message(cond)
  #           message()
  #         })
  #D <- RunTSNE(D,reduction="harmony",dims = 1:dimN,verbose = FALSE)
  message("Clustering ",cluster_method," (",clusterResolution,") ...")
  if(grepl("Leiden",cluster_method,ignore.case=T)){
    D <- FindClusters(D,verbose = FALSE,resolution=clusterResolution,algprithem=4)
  }else{
    D <- FindClusters(D,verbose = FALSE,resolution=clusterResolution)
  }
  message("saving ...")
  X <- cbind(#cID=row.names(D@meta.data),
    D@reductions$harmony@cell.embeddings,
    D@reductions$umap@cell.embeddings,
    #D@reductions$tsne@cell.embeddings,
    D@meta.data[,c("seurat_clusters",batch)])
  colnames(X) <- gsub("UMAP","umap",gsub("seurat_clusters","sctHarmony_cluster",colnames(X)))
  saveRDS(X,strOut)
}

main <- function(){
  suppressMessages(suppressWarnings(PKGloading()))
  batchKey="library_id" #"batch"
  args = commandArgs(trailingOnly=TRUE)
  
  
  source(paste0(dirname(gsub("--file=","",grep("file=",commandArgs(),value=T))),"/readH5ad.R"))
  if(length(args)<2) stop("Path to h5ad file, the output file and config file are required!")
  
  task <- args[1]
  if(task=='SCT'){
    strH5ad <- args[2]
    if(!file.exists(strH5ad)) stop(paste0("H5ad file (",strH5ad,") does not exist!"))
    strOut <- args[3]
    strConfig <- args[4]
    
    if(!file.exists(strConfig)){
      stop(paste0("Config file (",strConfig,") does not exist!"))
    }else{
      config <- yaml::read_yaml(strConfig)
    }
    if(length(args)>4) batchKey <- args[5]
    print(peakRAM(processV5(strH5ad,batchKey,strOut,config$harmonyBatchGene,
                            ifelse(is.null(config$subprocess),5,config$subprocess))))
    #print(peakRAM(processSCT(strH5ad,batchKey,strOut,config$harmonyBatchGene)))
  }else if(task=="Harmony"){
    strPCA <- args[2]
    if(!file.exists(strPCA)) stop(paste0("PCA file (",strPCA,") does not exist!"))
    strOut <- args[3]
    clu_reso <- as.numeric(args[4])
    clu_method <- args[5]
    if(length(args)>5) batchKey <- args[6]
    print(peakRAM(processPCA(strPCA,strOut,batchKey,clu_reso,clu_method)))
  }else{
    stop(paste("Unknown sctHarmony task:",task))
  }
}

main()