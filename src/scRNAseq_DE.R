## wrapper function for pipeline_class.R
## to have a better user experience, easier usage
PKGloading <- function(){
  require(rhdf5)
  require(Matrix)
  require(peakRAM)
  options(future.globals.maxSize=16*1024^3,stringsAsFactors=F)
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
  message("scDEG starting ...")
  config <- yaml::read_yaml(strConfig)
  strDEG <- checkFileExist(config$DEG_desp,"DEG description file")
  if(!is.null(config$prj_name)){
    prefix <- paste0(config$output,"/",config$prj_name)
    strH5adraw <- checkFileExist(paste0(file.path(config$output,"raw",config$prj_name),".h5ad"),"raw count h5ad file")
    strH5ad <- checkFileExist(paste0(file.path(config$output,config$prj_name),".h5ad"),"pipeline output h5ad contains cell annotation file")
  }else{
    prefix <- paste0(config$output,"/",config$DBname)
    strH5adraw <- checkFileExist(config$UMI,"raw count h5ad file")
    strH5ad <- checkFileExist(config$meta,"cell annotation file")
  }
  
  
  #sample,cluster,group,alt,ref,covars[comma separated],method[default NEBULA]
  compInfo <- as.data.frame(data.table::fread(strDEG))
  if(nrow(compInfo)==0){
    message("Empty comparison description file!\nDONE")
    q()
  }
  colnames(compInfo) <- sapply(strsplit(colnames(compInfo),"[[:punct:]]"),head,1)
  meta <- getMeta(strH5ad)
  ## check compInfo against meta
  message("\tComparison description file checking ...")
  if("comparisonName"%in%colnames(compInfo) && length(unique(compInfo$comparisonName))!=nrow(compInfo))
    stop("Unique comparison names are required!")
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

  ## create scDEG folder
  message("\tcreating scDEG folder ...")
  scDEGpath <- paste0(prefix,"_scDEG/")
  dir.create(scDEGpath,showWarnings=F)

  ## enumerate all comparisons
  message("\tcreating scDEG tasks ...")
  cmds <- apply(compInfo,1,function(x){
    if(is.na(x["method"]) || nchar(x["method"])<1){
      x["method"] <- "NEBULA"
      if(is.na(x["model"])) x["model"] <- "HL"
    }
    coV <- NULL
    if(!is.na(x["covars"]) && nchar(x["covars"])>2) coV <- unlist(strsplit(x["covars"],"\\+"),use.name=F)
    if("comparisonName"%in%names(x)){
      message(x["comparisonName"])
      strOut <- file.path(scDEGpath,x["comparisonName"])
    }else{
      strD <- paste(c(gsub("_",".",c(paste(x[c("alt","ref")],collapse=".vs."),
                                     x[c("cluster","group","method","model")]))),collapse="_")
      strOut <- list.files(scDEGpath,strD)
      if(length(strOut)>0) strOut <- file.path(scDEGpath,strOut[1])
      else strOut <- file.path(scDEGpath,strD)
    }
    if(dir.exists(strOut)) message("\tUsing existing: ",strOut)
    #system(paste("rm -rf",strOut))
    return(scRNAseq_DE(strH5adraw,strH5ad,strOut,x["method"],
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
  #message(paste(paste(names(cmds),cmds,sep=": "),collapse="\n"))
  write(rjson::toJSON(cmds),paste0(prefix,"_scDEG.cmd.json"))
  #print(head(cmds))
  #writeLines(paste(cmds,collapse="\n"),paste0(prefix,"_scDEG.cmd"))
  cat("scDEG task creation completed")
  
}

scRNAseq_DE <- function(
    strCount,
    strMeta,
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

    meta <- getMeta(strMeta)
    cmd <- c()
    for(one in unique(meta[,column_cluster])){
      #"Rscript cmdPath/scRNAseq_DE.R cmdPath grpInterest" paste(jID,one,addSRC)
      cmd <- c(cmd,paste0("cd ",output,"; Rscript ",scRNAseq_DE_path,"/scRNAseq_DE.R ",
                         scRNAseq_DE_path,' "',one,'"'))
    }
    names(cmd) <- gsub(" ",".",paste(basename(output),unique(meta[,column_cluster]),grp_alt,grp_ref,sep="_"))
    
    return(list(cmd))
}

getMeta <- function(strMeta){
  if(grepl("rds$",strMeta)){
    meta <- readRDS(strMeta)
  }else if(grepl("h5ad$",strMeta)){
    meta <- getobs(strMeta)
  }else{
    stop(paste("unknown meta format:",strMeta))
  }
  return(meta)
}
getUMI <- function(strUMI){
  if(grepl("rds$",strUMI)){
    UMI <- readRDS(strUMI)
  }else if(grepl("h5ad$",strUMI)){
    UMI <- getX(strUMI)
  }else{
    stop(paste("unknown UMI format:",strUMI))
  }
  return(UMI)
}

checkInput <- function(env){
    if(!file.exists(env$strCount) || !file.exists(env$strMeta)){
        stop("Either count RDS file or meta RDS file is missing!")
    }
    meta <- getMeta(env$strMeta)
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
checkCellMeta <- function(meta,env,cluster_interest){
  sel <- meta[,env$column_cluster]==cluster_interest
  # check total cells
  if(sum(sel)<10){
    message("\tSkip: Too few cells (",sum(sel),")")
    return()
  }
  # check if any subject in both groups
  subjGroup <- ftable(meta[sel,c(env$column_sample,env$column_group)])
  print(subjGroup)
  if(ncol(subjGroup)<2){
    message("Skip: missing one group!")
    return()
  }
  # check if at least 2 subjects in each group
  message("Checking minimal ",env$R6_min.cells.per.subj," cells per subject")
  subjGroupCell <- apply(subjGroup,2,function(x)return(sum(x>=env$R6_min.cells.per.subj)))
  if(sum(subjGroupCell[c(env$grp_ref,env$grp_ref)]>1)<2){
    message("Skip: One group contains less than 2 subjects!")
    return()
  }
  # check if pseudo bulk test, the cell level cov is not supported
  if(env$method %in% c("t_test","u_test","edgeR","limma","DESeq2")){
    if(length(unique(meta[,env$column_sample])) != nrow(unique(meta[,c(env$column_sample,env$column_group,env$column_covars)]))){
      message("Skip: cell level meta (e.g. pct_MT, library_size) is NOT supported in pseudo bulk comparison!")
      return()
    }
  }
  subjCell <- table(meta[sel,env$column_sample])
  sel <- sel & (meta[,env$column_sample]%in%names(subjCell)[subjCell>=env$R6_min.cells.per.subj])
  subjMeta <- meta[sel,c(env$column_sample,env$column_group,env$column_covars)]
  for(one in colnames(subjMeta)){
    if(sum(is.na(subjMeta[,one]))>0){
      message("Skip: NA found in ",one)
      return()
    }
  }
  subjMeta <- subjMeta[!duplicated(subjMeta[,env$column_sample]),]
  rownames(subjMeta) <- NULL
  print(subjMeta)
  return(as.vector(sel))
}
scRNAseq_DE_one <- function(
    env,
    cluster_interest,
    strSrc=NULL
){
    if(!is.null(strSrc)) source(paste0(strSrc,"/pipeline_class.R"))
    message("===== ",env$method,":",cluster_interest," =====")
    strOut <- env$output
    if(!is.null(env$column_group)){
      strF <- file.path(strOut,paste0(gsub("__","_",paste0(env$grp_alt,".vs.",env$grp_ref)),"__",
                                      gsub("__","_",paste(env$column_cluster,cluster_interest,sep=":")),"__",
                                      gsub("__","_",env$method),
                                      ".csv"))
    }else{
      strF <- file.path(strOut,paste0(cluster_interest,".vs.Rest","_",gsub("_",".",env$column_cluster),".csv"))
      intrestGrp <- as.character(allMeta[,env$column_cluster])
      intrestGrp[intrestGrp!=cluster_interest] <- "Rest"
      allMeta <- cbind(allMeta,all="all",intrestGrp=intrestGrp)
      env$column_cluster <- "all"
      env$column_group <- "intrestGrp"
      env$grp_ref <- "Rest"
      env$grp_alt <- cluster_interest
      cluster_interest <- "all"
    }
    if(file.exists(strF)){
      message("\tSkip: ",strF," exists!")
      return()
    }

    message("===== read counts and meta information =====")
    counts <- getUMI(env$strCount)
    allMeta <- getMeta(env$strMeta)[colnames(counts),]
    allMeta$cell <- rownames(allMeta)
    ## subset only the interested cluster cells would be included
    sel <- checkCellMeta(allMeta,env,cluster_interest)
    if(is.null(sel)) return()
    allMeta <- allMeta[sel,]
    counts <- counts[,sel]
    message("meta: ",nrow(allMeta),"; UMI:",ncol(counts))
    #strOut <- paste0(env$output,"/",env$method,"_",env$column_cluster,"/")
    #system(paste("mkdir -p",strOut))
    sce <- BiostatsSingleCell$new(count_data = counts,
                                  meta_data = allMeta,
                                  sampleId_col = env$column_sample,
                                  cluster_col = env$column_cluster,
                                  treatment_col = env$column_group)
    sce$set_group_mode(cluster_of_interest = cluster_interest, ref_group = env$grp_ref, alt_group =env$grp_alt)
    sce$make_QCplots(gsub("csv$","QC.pdf",strF))
    sce$apply_filter(min.cells.per.gene = env$min.cells.per.gene, min.genes.per.cell = env$min.genes.per.cell,
                     min.perc.cells.per.gene = env$min.perc.cells.per.gene,perc_filter = env$perc_filter) # 0% expression requirement

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
        saveRDS(de,file=paste0(strF,".rds"))
        selCol <- c("ID","log2FC","Pvalue","FDR")
        if("res.tab" %in% names(de))
          de <- de$res.tab
        if(sum(selCol %in% names(de))==4)
          write.csv(de[selCol], file=strF, row.names = FALSE)
        else
          message("\n\n*** MISSING required columns (",
                  paste(selCol,collapse=","),
                  ") in DE (",
                  paste0(strF,".rds"),
                  ")results ***\n")
          
        p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = env$method)
        ggsave(gsub("csv","png",strF))
        return(c(file.info(env$strCount)$size/(1024*1024),sum(allMeta[,env$column_cluster]==cluster_interest)))
    }
    return()
}

args <- commandArgs(trailingOnly=TRUE)
if(length(args)==1){
  a <- commandArgs()
  strPath <- gsub("^--file=","",grep("^--file",a,value=T)[1])
  scRNAseq_DE_path <<- dirname(normalizePath(strPath))
  source(paste0(scRNAseq_DE_path,"/readH5ad.R"))
  main(checkFileExist(args[1],"config file"))
}
if(length(args)>1){
  selGrp <- paste(args[-1],collapse=" ")
  #message("\n\n\n",args[-1],": ",selGrp,"\n\n\n")
  suppressMessages(suppressWarnings(PKGloading()))
  strAppPath <- dirname(gsub("--file=","",grep("file=",commandArgs(),value=T)))
  source(file.path(strAppPath,"readH5ad.R"))
  memUse <- peakRAM(Finfo <- scRNAseq_DE_one(readRDS("env.rds"),
                                            selGrp,
                                            args[1]))
  print(memUse)
  strLog <- file.path(dirname(strAppPath),"log","scDEG_memory.log")
  if(!file.exists(strLog)){
    dir.create(dirname(strLog),showWarnings=F)
    cat("File_size_MB\tCell_number\tTime_s\tRam_MiB\tPeak_Ram_MiB\n",sep="",file=strLog)
    Sys.chmod(strLog,"666",use_umask = F)
  }
  if(!is.null(Finfo)){
    message("Recording memory usage")
    cat(paste(round(c(Finfo,unlist(memUse[2:4])),4),collapse="\t"),"\n",sep="",file=strLog,append=T)
  }
  message("Successful!")
}

