PKGloading <- function(){
  require(kBET)
  require(ggplot2)
  require(dplyr)
  require(BiocParallel)
  options(stringsAsFactors = FALSE)
}

estimateBatch <- function(i,umapList,batchKey){
  tryCatch({
    oneM <- names(umapList)[i]
    message("\t working on ",oneM)
    umap <- umapList[[i]]
    batch <- unlist(umap[,batchKey])
    umap <- umap[,colnames(umap)!=batchKey]
    batch.estimate <- kBET::kBET(umap, batch, do.pca=FALSE, plot = FALSE, k0 = 100) ## look at 100 neighbors
    kbet <- data.frame(kBET = batch.estimate$stats$kBET.observed,
                       method = oneM, stringsAsFactors = FALSE)
    return(kbet)
  },error=function(e){
    print(e)
    return(NULL)
  })
}
plotBatch <- function(kBET_all,strKBET){
  kBET_all_rank <- kBET_all  %>% 
    group_by(method)  %>% 
    summarise(kBET_mean = mean(kBET))  %>% 
    arrange(kBET_mean)  %>% 
    filter(method != "nocorrection")
  
  kBET_all <- kBET_all  %>% 
    mutate(method = factor(method, levels = c(kBET_all_rank$method, "nocorrection")))
  
  p <- ggplot(kBET_all, aes(x = method, y = kBET)) + geom_boxplot() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = -45, hjust = 0), axis.title.x = element_blank(), plot.margin = margin(0.2, 0.8, 0.2 , 0.2, unit = 'in')) +
    ylim(0, 1) 

  dir.create(dirname(strKBET),showWarnings=F)
  ggsave(filename=strKBET,plot=p, width = 5, height = 3)
}
getUMAP <- function(strH5ad,batchKey){
  meta <- getobs(strH5ad)
  umapKey <- getobsmKey(strH5ad)
  uList <- list()
  for(one in grep("umap",umapKey,value=T)){
    mID <- gsub("^X_|_umap","",one)
    if(mID=="umap") mID <- "raw"
    uList[[mID]] <- cbind(data.frame(getobsm(strH5ad,one)[rownames(meta),]),meta[,batchKey])
    colnames(uList[[mID]]) <- c(paste0("umap_",1:(ncol(uList[[mID]])-1)),batchKey)
  }
  return(uList)
}
main <- function(){
  selfPath <- gsub("^--file=","",grep("^--file=",commandArgs(),value=T)[1])
  suppressMessages(suppressWarnings(PKGloading()))
  source(paste0(dirname(selfPath),"/readH5ad.R"))
  batchKey <- "library_id"
  
  args = commandArgs(trailingOnly=TRUE)
  if(length(args)<2) stop("the final h5ad file and the output file are required!")
  strH5ad <- args[1]
  strOut <- args[2]
  umapList <- getUMAP(strH5ad,batchKey)
  kBETall <- bind_rows(bplapply(1:length(umapList),estimateBatch,
                                umapList,batchKey,
                                BPPARAM = MulticoreParam(workers=min(length(umapList),max(1,parallelly::availableCores()-2)),
                                                         tasks=length(umapList))))
  plotBatch(kBETall,strOut)
}
main()