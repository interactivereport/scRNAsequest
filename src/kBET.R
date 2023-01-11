PKGloading <- function(){
  require(kBET)
  require(hdf5r)
  require(ggplot2)
  require(dplyr)
  require(rhdf5)
  #require(foreach)
  #require(doParallel)
  require(BiocParallel)
  register(MulticoreParam(workers=8))
  #options(MulticoreParam=quote(MulticoreParam(workers=8)))
}
getobs <- function(strH5ad){
  #message("\tobtainning obs ...")
  obs <- h5read(strH5ad,"obs")
  meta <- do.call(cbind.data.frame, obs[grep("^_",names(obs),invert=T)])
  if(!"_index"%in%names(obs)) dimnames(meta) <- list(obs[[1]],grep("^_",names(obs),invert=T,value=T))
  else dimnames(meta) <- list(obs[["_index"]],grep("^_",names(obs),invert=T,value=T))
  
  if("__categories"%in%names(obs)){
    for(one in names(obs[["__categories"]])){
      meta[,one] <- obs[["__categories"]][[one]][1+meta[,one]]
    }
  }
  return(meta)
}
getobsm <- function(strH5ad,key){
  k <- h5ls(strH5ad,recursive=2)
  k <- k[grepl("obsm",k[,1]),2]
  if(!key%in%k) return(NULL)
  X <- h5read(strH5ad,paste0("obsm/",key))
  return(t(X))
}

estimateBatch <- function(oneM,prefix,batchKey){
  #batchRaw <- unlist(getobs(paste0(prefix,"_raw.h5ad"))[,batchKey,drop=F])
  tryCatch({
    message("\t working on ",oneM)
    
    #strH5ad <- paste0(prefix,"_",oneM,".h5ad")
    strH5ad <- paste0(file.path(dirname(prefix),oneM,basename(prefix)),".h5ad")
    
    batch <- unlist(getobs(strH5ad)[,batchKey],use.names=F) #as.vector(batchRaw[rownames(getobs(strH5ad))]) #
    umap <- getobsm(strH5ad,"X_umap")
    if(is.null(umap)) stop(paste(oneM,"doesn't contain X_umap"))
    
    #adata = hdf5r::h5file(paste("working_data/", fh, sep = ""), mode = 'r')
    ## get batch labels
    ## accessing the fields in the h5 file is learned from Seurat codes
    #ncells = adata[['obs']]$dims
    
    ## batch is stored as factor, so when you get it out it'll be integer
    #batch = sapply(1:ncells, function(x) return(adata[['obs']][x][['library_id']]), simplify = FALSE)
    #batch = unlist(batch,use.names=F)
    
    ## X_umap
    #umap = sapply(1:ncells, function(x) return(adata[['obsm']][['X_umap']][,x]), simplify = FALSE)
    #umap = as.data.frame(do.call(rbind, umap))
    batch.estimate <- kBET::kBET(umap, batch, do.pca=FALSE, plot = FALSE, k0 = 100) ## look at 100 neighbors
    
    #saveRDS(batch.estimate, file = paste("working_data/", fh, ".batch.estimate.k0_100.RDS", sep = ""))
    kbet <- data.frame(kBET = batch.estimate$stats$kBET.observed, method = oneM, stringsAsFactors = FALSE)
    return(kbet)
  },error=function(e){
    print(e)
    return(NULL)
  })
}
plotBatch <- function(kBET_all,prefix){
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
  
  #ggplot(kBET_all, aes(x = method, y = kBET)) + geom_boxplot() +
  #    theme_minimal() +
  #    theme(axis.text.x = element_text(angle = -45, hjust = 0), axis.title.x = element_blank(), plot.margin = margin(0.2, 0.8, 0.2 , 0.2, unit = 'in')) +
  #    ylim(0, 1) #+ scale_x_discrete(labels = c("DESC", "Harmony", "Seurat3", "combat", "Scanorama(svd)", "linear regression", "mnn", "bbknn", "bbknn(trim)", "Scanorama",  "no correction"))
  #ggsave(filename = "figures/kBET_boxplot_across_methods_umap_k0_100.pdf", width = 5, height = 3)
  strKBET <- paste0(file.path(dirname(prefix),"evaluation",basename(prefix)),"_kBET_umap_k0_100.pdf")
  dir.create(dirname(strKBET),showWarnings=F)
  ggsave(filename=strKBET,plot=p, width = 5, height = 3)
}

main <- function(){
  suppressMessages(suppressWarnings(PKGloading()))
  batchKey <- "library_id"
  args = commandArgs(trailingOnly=TRUE)
  if(length(args)<2) stop("h5ad prefix and the method list are required!")
  prefix <- args[1]
  methods <- unlist(strsplit(args[2],","))
  
  kBETall <- bind_rows(bplapply(methods,estimateBatch,prefix,batchKey))
  plotBatch(kBETall,prefix)
}
main()