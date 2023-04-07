PKGloading <- function(){
  require(Seurat)
  require(scDblFinder)
  require(ggplot2)
  source(paste0(dirname(gsub("--file=","",grep("file=",commandArgs(),value=T))),"/readH5ad.R"))
}

dbl <- function(strH5ad,batch,strOut){
  conn <- file(paste0(strOut,".log"),"w")
  sink(conn,type="message")
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
  sink(type="message")
  close(conn)
  rm(Dlist)
  gc()
  DBL <- dplyr::bind_rows(DBL)
  dbl_plot(DBL,batch,strOut)
  saveRDS(DBL,paste0(strOut,".rds"))
  write.csv(DBL[,c("scDblFinder.class","scDblFinder.score")],file=strOut)
}
dbl_plot <- function(D,batch,strF){
  pdf(paste0(strF,".pdf"))
  D <- D[order(D$scDblFinder.score),]
  for(one in unique(D[,batch])){
    DD <- D[D[,batch]==one,]
    dbNum <- sum(DD[,"scDblFinder.class"]=="doublet")
    print(ggplot(DD,aes(nCount_RNA,nFeature_RNA,color=scDblFinder.class))+
            geom_point(size=1)+
            ggtitle(one)+
            annotate(geom="text",label=paste0(dbNum," (",round(dbNum/nrow(DD)*100,2),"%) doublets"),
                     x=Inf,y=-Inf,hjust=1,vjust=0)+
            theme_light())
    print(ggplot(DD,aes(nCount_RNA,nFeature_RNA,color=scDblFinder.score))+
            geom_point(size=1)+
            scale_color_continuous(type="viridis")+
            ggtitle(one)+
            theme_light())
  }
  a<- dev.off()
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
  
  print(system.time(dbl(strH5ad,batchKey,strOut)))
}

main()
  