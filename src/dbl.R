PKGloading <- function(){
  require(Seurat)
  require(scDblFinder)
  require(ggplot2)
  require(Matrix)
  require(BiocParallel)
  source(paste0(dirname(gsub("--file=","",grep("file=",commandArgs(),value=T))),"/readH5ad.R"))
}

dbl <- function(strH5ad,batch,strOut){
  X <- getX(strH5ad)
  gID <- setNames(rownames(X),gsub("_","-",rownames(X)))
  rownames(X) <- names(gID)
  meta <- getobs(strH5ad)
  D <- CreateSeuratObject(counts=X,
                          project="SCT",
                          meta.data=meta)
  rm(X,meta)
  gc()
  Dlist <- SplitObject(D,split.by=batch)
  rm(D)
  gc()
  DBL <- lapply(Dlist,function(one,medianUMI){
    cat("\t",one@meta.data[1,batch],"\n")
    message("-----",one@meta.data[1,batch],"\n")
    Xdbl <- scDblFinder(one@assays$RNA@counts)
    message("\n")
    return(cbind(data.frame(Xdbl@colData),one@meta.data[colnames(Xdbl),]))
  })
  rm(Dlist)
  gc()
  DBL <- dplyr::bind_rows(DBL)
  dbl_plot(DBL,strOut,batch)
  saveRDS(DBL,paste0(strOut,".rds"))
  write.csv(DBL[,c("scDblFinder.class","scDblFinder.score")],file=strOut)
}
dbl_single <- function(strUMI,strOut,strBarcode="library_id"){
  if(dir.exists(strUMI)){
    X <- Read10X(strUMI)
  }else if(grepl("h5$",strUMI)){
    suppressWarnings(X <- Read10X_h5(strUMI))
  }else if(grepl("csv$|tsv$",strUMI)){
    X <- data.table::fread(strUMI)
    X <- as(as.matrix(data.frame(row.names=unlist(X[,1]),X[,-1])),"sparseMatrix")
  }else{
    stop(paste("Unsupported UMI format:",strUMI))
  }
  if(file.exists(strBarcode)){
    barcodes <- unlist(data.table::fread(strBarcode,header=F)[,1],use.names=F)
    X <- X[,barcodes]
  }
  message(ncol(X),"\tcells for scDblFinder")
  DBL <- tryCatch(
    {
      Xdbl <- scDblFinder(X,BPPARAM=MulticoreParam(max(1,parallelly::availableCores()-2)))
      data.frame(Xdbl@colData)
    },
    error=function(cond){
      message("*** Error in scDblFinder for ",strUMI," ***")
      return(data.frame(row.names=barcodes,
                        scDblFinder.class=rep('Failed',length(barcodes)),
                        scDblFinder.score=rep(0,length(barcodes))))
    }
  )
  DBL <- cbind(DBL,nCount_RNA=colSums(X),nFeature_RNA=diff(X@p))
  saveRDS(DBL,paste0(strOut,".rds"))
  dbl_plot(DBL,strOut)
  write.csv(DBL[,c("scDblFinder.class","scDblFinder.score")],file=strOut)
}
dbl_plot <- function(D,strF,batch=NULL){
  pdf(paste0(strF,".pdf"))
  D <- D[order(D$scDblFinder.score),]
  if(is.null(batch)){
    dbl_oneplot(D)
  }else{
    for(one in unique(D[,batch])){
      dbl_oneplot(D[D[,batch]==one,],one)
    }
  }
  a<- dev.off()
}
dbl_oneplot <- function(DD,gTitle=""){
  dbNum <- sum(DD[,"scDblFinder.class"]=="doublet")
  print(ggplot(DD,aes(nCount_RNA,nFeature_RNA,color=scDblFinder.class))+
          geom_point(size=1)+
          ggtitle(gTitle)+
          annotate(geom="text",label=paste0(dbNum," (",round(dbNum/nrow(DD)*100,2),"%) doublets"),
                   x=Inf,y=-Inf,hjust=1,vjust=0)+
          theme_light())
  print(ggplot(DD,aes(nCount_RNA,nFeature_RNA,color=scDblFinder.score))+
          geom_point(size=1)+
          scale_color_continuous(type="viridis")+
          ggtitle(gTitle)+
          theme_light())
}
#source("../../src/dbl.R");dbl("raw/TST11837_oyoung_raw_prefilter.h5ad","library_id","dbl/TST11837_oyoung.csv")
main <- function(){
  suppressMessages(suppressWarnings(PKGloading()))
  batchKey="library_id" #"batch"
  args = commandArgs(trailingOnly=TRUE)
  if(length(args)<2) stop("Path to h5ad file and the output file are required!")
  strH5ad <- args[1]
  if(!file.exists(strH5ad)) stop(paste0("H5ad file (",strH5ad,") does not exist!"))
  strOut <- args[2]
  if(length(args)>2) batchKey <- args[3]
  
  if(grepl("h5ad$",strH5ad)){
    print(system.time(dbl(strH5ad,batchKey,strOut)))
  }else{
    print(system.time(dbl_single(strH5ad,strOut,batchKey)))
  }
  
}

main()
  