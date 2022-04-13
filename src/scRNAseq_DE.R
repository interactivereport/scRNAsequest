## wrapper function for pipeline_class.R
## to have a better user experience, easier usage
PKGloading <- function(){
  require(rhdf5)
  require(Matrix)
  options(future.globals.maxsize=3145728000,stringsAsFactors=F)
}
#
checkFileExist <- function(strF,msg="file"){
  if(!file.exists(strF)){
    message("ERROR: ",msg," does not exist:\n",strF,"\n")
    q()
  }
  return(strF)
}
main <- function(strConfig){
  suppressMessages(suppressWarnings(PKGloading()))
  message("\tscDEG starting ...")
  config <- yaml::read_yaml(strConfig)
  prefix <- paste0(config$output,"/",config$prj_name)
  strDEG <- checkFileExist(config$DEG_desp,"DEG description file")
  strH5adraw <- checkFileExist(paste0(prefix,"_raw.h5ad"),"raw count h5ad file")
  strH5ad <- checkFileExist(paste0(prefix,".h5ad"),"pipeline output h5ad file")
  
  #sample,cluster,group,alt,ref,covars[comma separated],method[default NEBULA]
  compInfo <- as.data.frame(data.table::fread(strDEG))
  if(nrow(compInfo)==0){
    message("Empty comparison description file!\nDONE")
    q()
  }
  colnames(compInfo) <- sapply(strsplit(colnames(compInfo),"[[:punct:]]"),head,1)
  meta <- getobs(strH5ad)
  ## check compInfo against meta
  message("\tComparison description file checking ...")
  selCol <- apply(compInfo,1,function(x){
    hh <- unlist(x[c("sample","cluster","group")])
    if(!is.na(x["covars"]) && nchar(x["covars"])>2){
      coV <- unlist(strsplit(x["covars"],"\\+"),use.name=F)
      if(sum(hh%in%coV)>0) stop("Overlap covars with other setups!")
      hh <- c(hh,coV)
    }
    if(sum(!hh%in%colnames(meta))>0){
      message(paste(hh[!hh%in%colnames(meta)],collapse=", "))
      stop("The above header defined in DEG description file are NOT in the obs of h5ad file")
    }
    v <- x[c("alt","ref")]
    if(sum(!v%in%unique(meta[,x["group"]]))>0){
      message(paste(v[!v%in%unique(meta[,x["group"]])],collapse=", "))
      stop("The above entris defined in DEG description file are NOT in the obs column")
    }
    return(hh)
  })
  meta <- meta[,unique(unlist(selCol,use.names=F))]
  # obtain the raw counts
  X <- getX(strH5adraw)
  X <- X[,rownames(meta)]
  
  ## create scDEG folder
  message("\tsaving counts and meta ...")
  scDEGpath <- paste0(prefix,"_scDEG/")
  dir.create(scDEGpath,showWarnings=F)
  strCounts <- paste0(scDEGpath,"counts.rds")
  strMeta <- paste0(scDEGpath,"meta.rds")
  saveRDS(X,strCounts)
  saveRDS(meta,strMeta)
  
  ## enumerate all comparisons
  message("\tcreating scDEG tasks ...")
  cmds <- apply(compInfo,1,function(x){
    if(is.na(x["method"])){
      x["method"] <- "NEBULA"
      if(is.na(x["model"])) x["model"] <- "HL"
    }
    coV <- NULL
    if(!is.na(x["covars"]) && nchar(x["covars"])>2) coV <- unlist(strsplit(x["covars"],"\\+"),use.name=F)
    strOut <- file.path(scDEGpath,paste(gsub("_",".",x[c("cluster","group","method","model")]),collapse="_"))
    return(scRNAseq_DE(strCounts,strMeta,strOut,x["method"],
                x["sample"],x["cluster"],
                x["group"],x["ref"],x["alt"],
                x["model"],coV,
                
                min.cells.per.gene = config$min.cells.per.gene,
                min.genes.per.cell = config$min.genes.per.cell,
                min.perc.cells.per.gene = config$min.perc.cells.per.gene,
                perc_filter = config$perc_filter,
                
                R6_min.cells.per.gene = config$R6_min.cells.per.gene,
                R6_min.perc.cells.per.gene = config$R6_min.perc.cells.per.gene,
                R6_min.cells.per.gene.type = config$R6_min.cells.per.gene.type,
                R6_cells.per.gene.filter = config$R6_cells.per.gene.filter,
                R6_perc.cells.filter = config$R6_perc.cells.filter,
                R6_perc.filter = config$R6_perc.filter,
                R6_perc.filter.type = config$R6_perc.filter.type,
                R6_perc_threshold = config$R6_perc_threshold,
                R6_min.ave.pseudo.bulk.cpm = config$R6_min.ave.pseudo.bulk.cpm,
                R6_pseudo.bulk.cpm.filter = config$R6_pseudo.bulk.cpm.filter,
                R6_min.cells.per.subj = config$R6_min.cells.per.subj))
  })
  cmds <- unlist(cmds)#,use.names=F
  write(rjson::toJSON(cmds),paste0(prefix,"_scDEG.cmd.json"))
  #print(head(cmds))
  #writeLines(paste(cmds,collapse="\n"),paste0(prefix,"_scDEG.cmd"))
  cat("scDEG task creation completed")
  
}
getobs <- function(strH5ad){
  message("\tobtainning obs ...")
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
getX <- function(strH5ad){
  message("\tobtainning X ...")
  X <- h5read(strH5ad,"X")
  gID <- h5read(strH5ad,"var/_index")
  cID <- h5read(strH5ad,"obs/_index")
  if((max(X$indices)+1)==length(gID)){ # CSR sparse matrix
    M <- sparseMatrix(i=X$indices+1,p=X$indptr,x=as.numeric(X$data),
                      dims=c(length(gID),length(cID)),
                      dimnames=list(gID,cID))
  }else if((max(X$indices)+1)==length(cID)){#CSC sparse matrix
    M <- sparseMatrix(j=X$indices+1,p=X$indptr,x=as.numeric(X$data),
                      dims=c(length(gID),length(cID)),
                      dimnames=list(gID,cID))
  }
  return(M)
}


scRNAseq_DE <- function(
    countRDS,
    metaRDS,
    output,
    method,
    column_sample,
    column_cluster,
    column_group=NULL,
    grp_ref=NULL,
    grp_alt=NULL,
    method_model=NULL,
    column_covars=NULL,
    
    min.cells.per.gene = 3,
    min.genes.per.cell = 250,
    min.perc.cells.per.gene = 0.00,
    perc_filter = TRUE,

    R6_min.cells.per.gene = 3,
    R6_min.perc.cells.per.gene = 0.1,
    R6_min.cells.per.gene.type = "or",
    R6_cells.per.gene.filter = TRUE,
    R6_perc.cells.filter = TRUE,
    R6_perc.filter = FALSE,
    R6_perc.filter.type = "and",
    R6_perc_threshold = 0.90,
    R6_min.ave.pseudo.bulk.cpm = 1,
    R6_pseudo.bulk.cpm.filter = FALSE,
    R6_min.cells.per.subj = 3,
    
    core=1,
    parallel=FALSE,
    addSRC=NULL
){
    env <- as.list(environment())
    checkInput(env)
    saveRDS(env,file=file.path(output,"env.rds"))

    meta <- readRDS(metaRDS)
    cmd <- c()
    for(one in unique(meta[,column_cluster])){
      #"Rscript cmdPath/scRNAseq_DE.R cmdPath grpInterest" paste(jID,one,addSRC)
      cmd <- c(cmd,paste("cd",output,";Rscript",paste0(scRNAseq_DE_path,"/scRNAseq_DE.R"),
                         scRNAseq_DE_path,one))
    }
    names(cmd) <- paste(basename(output),unique(meta[,column_cluster]),grp_alt,grp_ref,sep="_")
    return(cmd)
}

checkInput <- function(env){
    if(!file.exists(env$countRDS) || !file.exists(env$metaRDS)){
        stop("Either count RDS file or meta RDS file is missing!")
    }
    meta <- readRDS(env$metaRDS)
    for(one in c(env$column_sample,env$column_cluster,env$column_group,env$column_covars)){
        if(!one%in%colnames(meta))
            stop(paste(one,"is not in the sample meta table!"))
    }
    if(!is.null(env$column_group)){
        if(is.null(env$grp_ref) || is.null(env$grp_alt))stop(paste("grp_ref and grp_alt are required for",env$column_group))
        for(one in c(env$grp_ref,env$grp_alt)){
            if(!one%in%unique(meta[,env$column_group]))
                stop(paste(one,"is not in the column",env$column_group,"from sample meta table!"))
        }
    }
    allMethods <- c("t_test", "u_test","edgeR","limma","DESeq2","MAST","limma_cell_level","glmmTMB","NEBULA","ancova")
    if(!env$method%in%allMethods)
        stop(paste0("method (",env$method,")is not supported!\nSelect from ",
                    paste(allMethods,collapse=", ")))
    if(env$method=="NEBULA"){
        if(!env$method_model%in%c("LN", "HL")) stop("method_model has to be LN or HL for NEBULA method!")
        if(env$method_model!="HL") warning("method_model is recommended to be HL for NEBULA method!")
    }
    else if(env$method=="glmmTMB"){
        if(!env$method_model%in% c("nbinom2", "nbinom1", "poisson", "nbinom2zi", "nbinom1zi"))
            stop("method_model has to be nbinom2, nbinom1, poisson, nbinom2zi or nbinom1zi for glmmTMB method!")
        if(env$method_model!="nbinom2") warning("method_model is recommended to be nbinom2 for NEBULA method!")
    }
    system(paste("mkdir -p",env$output))
}

scRNAseq_DE_one <- function(
    env,
    cluster_interest,
    strSrc=NULL
){
    if(!is.null(strSrc)) source(paste0(strSrc,"/pipeline_class.R"))
    message("===== read counts and meta information =====")
    counts <- readRDS(env$countRDS)
    allMeta <- readRDS(env$metaRDS)[colnames(counts),]
    allMeta$cell <- rownames(allMeta)
    message("===== ",env$method,":",cluster_interest," =====")
    #strOut <- paste0(env$output,"/",env$method,"_",env$column_cluster,"/")
    #system(paste("mkdir -p",strOut))
    strOut <- env$output
    if(!is.null(env$column_group)){
        strF <- file.path(strOut,paste0(env$grp_alt,".vs.",env$grp_ref,"_",gsub("_",".",cluster_interest),".QC.pdf"))
        sce <- BiostatsSingleCell$new(count_data = counts,
                                      meta_data = allMeta,
                                      sampleId_col = env$column_sample,
                                      cluster_col = env$column_cluster,
                                      treatment_col = env$column_group)
        sce$set_group_mode(cluster_of_interest = cluster_interest, ref_group = env$grp_ref, alt_group =env$grp_alt)

    }else{
        strF <- file.path(strOut,paste0(cluster_interest,".vs.Rest","_",gsub("_",".",env$column_cluster),".QC.pdf"))
        intrestGrp <- as.character(allMeta[,env$column_cluster])
        intrestGrp[intrestGrp!=cluster_interest] <- "others"
        meta <- cbind(allMeta,all="all",intrestGrp=intrestGrp)
        
        sce <- BiostatsSingleCell$new(count_data = counts,
                                      meta_data = meta,
                                      sampleId_col = env$column_sample,
                                      cluster_col = "all",
                                      treatment_col = "intrestGrp")
        sce$set_group_mode(cluster_of_interest = "all", ref_group = "others", alt_group = cluster_interest)
        
    }
    sce$make_QCplots(strF)
    sce$apply_filter(min.cells.per.gene = env$min.cells.per.gene, min.genes.per.cell = env$min.genes.per.cell,
                     min.perc.cells.per.gene = env$min.perc.cells.per.gene,perc_filter = env$perc_filter) # 0% expression requirement
    strF <- gsub("QC.pdf","csv",strF)
    
    if(tryCatch({
        sce_qc <- sce$apply_filter_contrasts_R6(min.cells.per.gene = env$R6_min.cells.per.gene,
                                                min.perc.cells.per.gene = env$R6_min.perc.cells.per.gene,
                                                perc.cells.filter = env$R6_perc.cells.filter,
                                                min.cells.per.gene.type = env$R6_min.cells.per.gene.type,
                                                cells.per.gene.filter = env$R6_cells.per.gene.filter,
                                                perc.filter = env$R6_perc.filter,
                                                perc.filter.type = env$R6_perc.filter.type,
                                                perc_threshold = env$R6_perc_threshold,
                                                min.ave.pseudo.bulk.cpm = env$R6_min.ave.pseudo.bulk.cpm,
                                                pseudo.bulk.cpm.filter = env$R6_pseudo.bulk.cpm.filter,
                                                min.cells.per.subj = env$R6_min.cells.per.subj)
        T
    },error=function(err){
        print(err)
        F
    })){
        covars <- env$column_covars
        system.time(de <- switch(env$method,
                                  't_test' = sce_qc$t_test_pipeline(),
                                  'u_test' = sce_qc$u_test_pipeline(),
                                  'edgeR'= sce_qc$edgeR_pipeline(covs = covars),
                                  'limma'= sce_qc$limma_pipeline(covs = covars),
                                  'DESeq2'= sce_qc$DESeq2_pipeline(covs = covars),
                                  'glmmTMB'= sce_qc$glmmTMB_pipeline(covs = covars, family = env$method_model, cores=env$core,detection_rate = FALSE),
                                  'MAST'= sce_qc$MAST_pipeline(covs = covars, detection_rate = TRUE),
                                  'limma_cell_level'= sce_qc$limma_cell_level_pipeline(covs = covars),
                                  'ancova'= sce_qc$ancova_pipeline(covs = covars),
                                  'NEBULA'= sce_qc$nebula_pipeline(covs = covars,method=env$method_model)))
        write.csv(de, file=strF, row.names = FALSE)
        p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = env$method)
        ggsave(gsub("csv","png",strF))
    }
    message("Successful!")
}

args <- commandArgs(trailingOnly=TRUE)
if(length(args)==1){
  a <- commandArgs()
  strPath <- gsub("^--file=","",grep("^--file",a,value=T)[1])
  scRNAseq_DE_path <<- dirname(normalizePath(strPath))
  main(checkFileExist(args[1],"config file"))
}
if(length(args)>1){
  print(system.time(
    scRNAseq_DE_one(readRDS("env.rds"),
                    paste(args[2]),
                    args[1])
  ))
}

