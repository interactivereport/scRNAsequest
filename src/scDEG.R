require(R6)
require(cli)
require(BiocParallel)

sparse_batch_apply <- function(X,MARGIN,FUN,shortLab="",...){
  step <- 500
  while(step>10){
    Xres <- c()
    ii <- 0
    Xlength <- ifelse(MARGIN==1,nrow(X),ncol(X))
    tryCatch({
      steps <- unique(c(seq(0,Xlength,step),Xlength))
      stepPair <- sapply(2:length(steps),function(i)return(list(c(steps[i-1]+1,steps[i]))))
      message(shortLab,": step ",step)
      Xres <- bplapply(stepPair,function(step,marg,oneFun){
        if(marg==1){
          A <- apply(X[step[1]:step[2],],MARGIN,FUN)#,...
        }else if(marg==2){
          A = apply(X[,step[1]:step[2]],MARGIN,FUN)#,...
        }else{
          stop("unknown MARGIN: ",MARGIN)
        }
        return(A)
      },marg=MARGIN,oneFun=FUN,#...,
      BPPARAM = MulticoreParam(workers=min(30,length(stepPair),max(1,parallelly::availableCores()-1)),
                               tasks=length(stepPair)))#,progressbar=T
      if(sum(!is.null(unlist(sapply(Xres,ncol))))==0){
        Xres <- unlist(Xres)
        #if(MARGIN==1) Xres <- Xres[rownames(X)]
        #else Xres <- Xres[colnames(X)]
      }else{
        Xres <- dplyr::bind_cols(Xres)
        #if(MARGIN==1) Xres <- Xres[,rownames(X)]
        #else Xres <- Xres[,colnames(X)]
      }
      step <- -1
    },
    error=function(cond){
      step <- floor(step/2)
      message("reduce step: ",step)
    })
  }
  #print(head(Xres))
  if(step==-1) return(Xres)
  stop("sparse_batch_apply")
}
# R6 class, save memory by keep one copy of the data matrix
# print out message for user to remove data outside R6 once the initialization process is finished
scDEG <- R6Class("scDEG",
                 public=list(
                   initialize=function(X=NULL,meta=NULL,
                                       id_col=NULL,cluster_col=NULL,grp_col=NULL,ctrl_value=NULL,alt_value=NULL,
                                       strX=NULL,strMeta=NULL,
                                       pipelinePath=""){
                     message("Loading ...")
                     stopifnot(!is.null(id_col))
                     stopifnot(!is.null(cluster_col))
                     stopifnot(!is.null(grp_col))
                     stopifnot(!is.null(ctrl_value))
                     stopifnot(!is.null(alt_value))
                     stopifnot(!is.null(X) || (!is.null(strX) && file.exists(strX)))
                     stopifnot(!is.null(meta) || (!is.null(strMeta) && file.exists(strMeta)))
                     suppressMessages(suppressWarnings(private$loadPKG(pipelinePath)))
                     
                     message("Reading ...")
                     private$id_col <- id_col
                     private$cluster_col <- cluster_col
                     private$grp_col <- grp_col
                     private$ctrl_value <- ifelse(is.na(ctrl_value)||nchar(ctrl_value)<2,NA,ctrl_value)
                     private$alt_value <- ifelse(is.na(alt_value)||nchar(alt_value)<2,NA,alt_value)
                     if(is.null(meta)) private$getMeta(strMeta)
                     else private$meta <- meta
                     private$check_meta()
                     
                     if(is.null(X)){
                       private$getUMI(strX)
                     }else{
                       private$X <- X
                       message("***** Please considering remove the original data matrix before 'run'! *****")
                     }
                     
                     if(sum(!rownames(private$meta)%in%colnames(private$X))>0){
                       stop("There are ",sum(!rownames(private$meta)%in%colnames(private$X))," cells defined in the meta but not in UMI!")
                     }
                     private$X <- private$X[,rownames(private$meta)]
                     
                     message("The scDEG initialization is completed successfully!")
                   },
                   run=function(alg,fliter_list,clusters=NULL,covar=NULL,prefix=NULL,...){
                     list2env(fliter_list,envir=environment())
                     #The following should be in the fitler_list
                     stopifnot(exists('rmGene'))
                     stopifnot(exists('min.cells.per.gene'))
                     stopifnot(exists('min.perc.cells.per.gene'))
                     stopifnot(exists('min.cells.per.gene.type'))
                     stopifnot(exists('perc_threshold'))
                     stopifnot(exists('perc.filter.type'))
                     stopifnot(exists('lib_size_low'))
                     stopifnot(exists('lib_size_high'))
                     stopifnot(exists('min.genes.per.cell'))
                     stopifnot(exists('min.cells.per.subj'))
                     stopifnot(exists('min.ave.pseudo.bulk.cpm'))
                     
                     if(is.null(clusters)) clusters <- unique(private$meta[,private$cluster_col])
                     
                     allDEG <- NULL
                     for(one in clusters){
                       private$g_index <- rep(T,nrow(private$X))
                       private$c_index <- private$meta[,private$cluster_col]==one
                       private$scDEG_apply_filter_contrast(rmGene,
                                                   min.cells.per.gene,min.perc.cells.per.gene,min.cells.per.gene.type,
                                                   #set both number to be 0 to avoid any gene filtering
                                                   perc_threshold, perc.filter.type,
                                                   lib_size_low, lib_size_high,
                                                   min.genes.per.cell, min.cells.per.subj)
                       sInfo <- private$scDEG_check_model(covar)
                       if(!is.null(prefix)) data.table::fwrite(sInfo,paste0(prefix,"_",one,"_sampleInfo.csv"))
                       
                       if(alg %in%c('DESeq2')){
                         private$createPseudo()
                         private$scDEG_apply_filter_pseudoBulk(min.ave.pseudo.bulk.cpm)
                         message("pseudo bulk model: ",paste(colnames(private$sInfo),collapse="+"))
                       }
                       
                       oneDEG <- function(){
                         de <- switch (alg,
                           'NEBULA' = private$scDEG_Nebula(covar,...),
                           'DESeq2' = private$scDEG_DESeq2(...)
                         )
                         return(de)
                       }
                       de <- cbind(oneDEG(),cluster=one)
                       if(is.null(allDEG)) allDEG <- de
                       else allDEG <- rbind(allDEG,de)
                     }
                     return(allDEG)

                   }
                 ),
                 private=list(
                   X=NULL,meta=NULL,id_col=NULL,cluster_col=NULL,grp_col=NULL,ctrl_value=NULL,alt_value=NULL,
                   c_index=NULL,g_index=NULL,pseudoX=NULL,sInfo=NULL,
                   loadPKG=function(strPath){
                     stopifnot(require(data.table))
                     stopifnot(require(dplyr))
                     stopifnot(require(glue))
                     stopifnot(require(SingleCellExperiment))
                     stopifnot(require(edgeR))
                     stopifnot(require(DESeq2))
                     stopifnot(require(apeglm))
                     stopifnot(require(glmmTMB))
                     stopifnot(require(BiocParallel))
                     stopifnot(require(Seurat))
                     stopifnot(require(MAST))
                     stopifnot(require(Matrix))
                     stopifnot(require(slam))
                     stopifnot(require(foreach))
                     stopifnot(require(doMC))
                     stopifnot(require(biomaRt))
                     stopifnot(require(emmeans))
                     stopifnot(require(nebula))
                     stopifnot(require(scater))
                     #stopifnot(require(peakRAM))
                     source(file.path(strPath,"readH5ad.R"))
                   },
                   getUMI=function(strUMI){
                     if(grepl("rds$",strUMI)){
                       private$X <- readRDS(strUMI)
                     }else if(grepl("h5ad$",strUMI)){
                       private$X <- getX(strUMI)
                     }else{
                       stop(paste("unknown UMI format:",strUMI))
                     }
                     gc()
                   },
                   getMeta=function(strMeta){
                     if(grepl("rds$",strMeta)){
                       private$meta <- readRDS(strMeta)
                     }else if(grepl("h5ad$",strMeta)){
                       private$meta <- getobs(strMeta)
                     }else{
                       stop(paste("unknown meta format:",strMeta))
                     }
                   },
                   check_meta=function(){
                     if(!private$id_col %in% colnames(private$meta)) stop(paste(private$id_col,"NOT in meta"))
                     if(!private$cluster_col %in% colnames(private$meta)) stop(paste(private$cluster_col,"NOT in meta"))
                     if(!private$grp_col %in% colnames(private$meta)) stop(paste(private$grp_col,"NOT in meta"))
                     if(!is.na(private$ctrl_value) && !private$ctrl_value %in% private$meta[,private$grp_col]) stop(paste(private$ctrl_value,"NOT in",private$grp_col))
                     if(!is.na(private$alt_value) && !private$alt_value %in% private$meta[,private$grp_col]) stop(paste(private$alt_value,"NOT in",private$grp_col))
                     private$meta <- private$meta[order(private$meta[,private$id_col]),]
                   },
                   scDEG_apply_filter=function(rmGene=c("Mt-","MT-","mt-"), lib_size_low = 200, lib_size_high = 20*10^6,
                                               min.cells.per.gene = 50,min.perc.cells.per.gene = 0.01,
                                               min.genes.per.cell = 500){
                     # remove low expressed genes
                     if(length(rmGene)){
                       message("Filtering genes which start with ",paste(rmGene,collapse=", "))
                       private$g_index <- private$g_index & !grepl(paste(paste0("^",rmGene),collapse="|"),rownames(private$X))
                     }
                     if(min.perc.cells.per.gene>0 && min.perc.cells.per.gene<1){
                       message("\tCalculating minimal cell number of a gene based on  min.perc.cells.per.gene: ",min.perc.cells.per.gene)
                       min.cells.per.gene <- max(min.cells.per.gene,ceiling(min.perc.cells.per.gene*ncol(private$X)))
                     }
                     message("Filtering genes by minimal cell number: ",min.cells.per.gene)
                     private$g_index <- private$g_index & sparse_batch_apply(private$X,1,function(x)return(sum(x>0)),shortLab="CellN per gene")>min.cells.per.gene
                     message("--- Total of ",sum(private$g_index)," genes")
                     # remove low sequence depth cells
                     libSize <- colSums(private$X[as.vector(private$g_index),])
                     nGene <- sparse_batch_apply(private$X[as.vector(private$g_index),],2,function(x)return(sum(x>0,na.rm=T)),shortLab="GeneN per cell")
                     private$cellFiltering(libSize>lib_size_low & libSize<lib_size_high,
                                           paste0("Filtering cells by sequence depth: ",lib_size_low," ~ ",lib_size_high))
                     private$cellFiltering(nGene>min.genes.per.cell,
                                           paste0("Filtering cells by min genes: ",min.genes.per.cell))
                     message("--- Total of ",sum(private$c_index)," cells")
                   },
                   scDEG_apply_filter_contrast=function(rmGene=c("Mt-","MT-","mt-"),
                                                        min.cells.per.gene = 50,min.perc.cells.per.gene = 0.10,min.cells.per.gene.type = "and",
                                                        #set both number to be 0 to avoid any gene filtering
                                                        perc_threshold = 0.9, perc.filter.type = "and",
                                                        lib_size_low = 200, lib_size_high = 20*10^6,
                                                        min.genes.per.cell = 500, min.cells.per.subj = 5){
                     if(length(rmGene)>0){
                       message("Filtering genes which start with ",paste(rmGene,collapse=", "))
                       private$g_index <- private$g_index & !grepl(paste(paste0("^",rmGene),collapse="|"),rownames(private$X))
                     }
                     if(!is.na(private$ctrl_value) && !is.na(private$alt_value)){
                       ctrl <- private$c_index & private$meta[,private$grp_col]==private$ctrl_value
                       alter <- private$c_index & private$meta[,private$grp_col]==private$alt_value
                       allIndex <- ctrl|alter
                       ctrlG <- sparse_batch_apply(private$X[,as.vector(ctrl)],1,function(x)return(sum(x>0,na.rm=T)),shortLab="Ctrl cellN per gene")
                       alterG <- sparse_batch_apply(private$X[,as.vector(alter)],1,function(x)return(sum(x>0)),shortLab="Alt cellN per gene")
                       if(min.cells.per.gene>0){
                         message("Filtering genes by minimal cell number: ",min.cells.per.gene)
                         if(min.perc.cells.per.gene>0 && min.perc.cells.per.gene<1){
                           message("\tCalculating minimal cell number of a gene based on min.perc.cells.per.gene: ",min.perc.cells.per.gene)
                           min.cells.per.gene <- max(min.cells.per.gene,ceiling(min.perc.cells.per.gene*min(sum(ctrl),sum(alter))))
                         }
                         if(grepl('and',min.cells.per.gene.type,ignore.case=T)){
                           private$g_index <- private$g_index & (ctrlG>min.cells.per.gene & alterG>min.cells.per.gene)
                         }else{
                           private$g_index <- private$g_index & (ctrlG>min.cells.per.gene | alterG>min.cells.per.gene)
                         }
                         message("\t",sum(private$g_index)," genes")
                       }
                       
                       if(perc_threshold<1){
                         message("Filtering genes by maximal zero-expression percentile within a group: ",perc_threshold)
                         if(grepl('and',perc.filter.type,ignore.case=T)){
                           private$g_index <- private$g_index & (ctrlG/sum(ctrl) > (1-perc_threshold) & 
                                                                   alterG/sum(alter) > (1-perc_threshold))
                         }else{
                           private$g_index <- private$g_index & (ctrlG/sum(ctrl) > (1-perc_threshold) | 
                                                                   alterG/sum(alter) > (1-perc_threshold))
                         }
                       }
                     }else{
                       allIndex <- private$c_index
                       allG <- sparse_batch_apply(private$X[,as.vector(allIndex)],1,function(x)return(sum(x>0,na.rm=T)),shortLab="All cellN per gene")
                       if(min.cells.per.gene>0){
                         message("Filtering genes by minimal cell number: ",min.cells.per.gene)
                         if(min.perc.cells.per.gene>0 && min.perc.cells.per.gene<1){
                           message("\tCalculating minimal cell number of a gene based on min.perc.cells.per.gene: ",min.perc.cells.per.gene)
                           min.cells.per.gene <- max(min.cells.per.gene,ceiling(min.perc.cells.per.gene*sum(allIndex)))
                         }
                         private$g_index <- private$g_index & allG>min.cells.per.gene
                         message("\t",sum(private$g_index)," genes")
                       }
                       if(perc_threshold<1){
                         message("Filtering genes by maximal zero-expression percentile within a group: ",perc_threshold)
                         private$g_index <- private$g_index & (allG/sum(allIndex) > (1-perc_threshold))
                       }
                     }
                     message("--- Total of ",sum(private$g_index)," genes")
                     stopifnot(sum(private$g_index)>100)
                     
                     # remove low sequence depth cells
                     message("remove low sequence depth cells")
                     libSize <- nGene <- rep(0,length(private$c_index))
                     libSize[allIndex] <- colSums(private$X[as.vector(private$g_index),as.vector(allIndex)])
                     nGene[allIndex] <- sparse_batch_apply(private$X[as.vector(private$g_index),as.vector(allIndex)],2,function(x)return(sum(x>0,na.rm=T)),shortLab="GeneN per cells")
                     private$cellFiltering(libSize>lib_size_low & libSize<lib_size_high,
                                   paste0("Filtering cells by sequence depth: ",lib_size_low," ~ ",lib_size_high))
                     private$cellFiltering(nGene>min.genes.per.cell,
                                   paste0("Filtering cells by min genes: ",min.genes.per.cell))
                     
                     # remove the samples with low cell number
                     cNum <- table(private$meta[private$c_index,private$id_col])
                     #if(sum(cNum<min.cells.per.subj)>0) message("removing: ",paste(names(cNum)[cNum<min.cells.per.subj],collapse=","))
                     private$cellFiltering(private$meta[,private$id_col]%in%names(cNum[cNum>=min.cells.per.subj]),
                                   paste0("Filtering samples with low cell number: ",min.cells.per.subj))
                     
                     message("--- Total of ",sum(private$c_index)," cells")
                   },
                   scDEG_check_model=function(covar){
                     if(sum(!covar%in%colnames(private$meta))>0) stop("Undefinied covar in meta!")
                     message("Please check the following sample meta information:")
                     meta <- private$meta[private$c_index,c(private$id_col,private$cluster_col,private$grp_col,covar)] %>% 
                       group_by_at(private$id_col) %>% 
                       group_modify(function(x,k){
                         oneInfo <- data.frame(row.names=k,cellN=nrow(x))
                         endCol <- c()
                         for(one in colnames(x)){
                           if(is.numeric(x[[one]])){
                             oneInfo <- cbind(oneInfo,setNames(data.frame(mean(x[[one]])),paste0(one,"_mean")))
                             oneInfo <- cbind(oneInfo,setNames(data.frame(sd(x[[one]])),paste0(one,"_sd")))
                           }else{
                             nC <- table(x[[one]])
                             if(length(nC)==1){
                               oneInfo <- cbind(oneInfo,setNames(data.frame(names(nC)),one))
                             }else if(length(nC)>(nrow(x)/4)){
                               warning("\tSkip! There are more than a quater of unique values of ",one," in ",k)
                               oneInfo <- cbind(oneInfo,setNames(data.frame("Skip"),one))
                             }else{
                               oneInfo <- cbind(oneInfo,setNames(data.frame(paste(paste(names(nC),nC,sep=": "),collapse="; ")),one))
                               endCol <- c(endCol,one)
                             }
                           }
                         }
                         if(length(endCol)>0) oneInfo <- cbind(oneInfo[,!colnames(oneInfo)%in%endCol],oneInfo[,endCol,drop=F])
                         return(oneInfo)
                       })
                     print(data.frame(meta))
                     
                     private$sInfo <- as.data.frame(meta[,!(grepl('_sd$',colnames(meta)))])
                     colnames(private$sInfo) <- gsub("_mean$","",colnames(private$sInfo))
                     rownames(private$sInfo) <- private$sInfo[,private$id_col]
                     private$sInfo <- private$sInfo[,c(private$grp_col,covar),drop=F]
                     private$sInfo <- private$sInfo[,sapply(private$sInfo,function(x)return(sum(grepl("^Skip$",x)|(grepl(":",x)&grepl(";",x)))==0)),drop=F]
                     message("return scDEG_check_model")
                     return(meta)
                   },
                   scDEG_apply_filter_pseudoBulk=function(min.ave.pseudo.bulk.cpm = 1){
                     sTotal <- colSums(private$pseudoX)
                     g_avg_cpm <- apply(private$pseudoX,1,function(x)return(mean(x/sTotal*1e6)))
                     message("Filtering genes by minimal average pseudo CPM: ",min.ave.pseudo.bulk.cpm)
                     private$pseudoX <- private$pseudoX[g_avg_cpm>min.ave.pseudo.bulk.cpm,]
                     message("--- Total of ",nrow(private$pseudoX)," genes")
                   },
                   # DEG by NEBULA pipeline -------
                   scDEG_Nebula=function(covar,saveRaw=NULL,...){
                     meta <-  private$meta[private$c_index,c(private$grp_col,covar),drop=F]
                     if(is.na(private$ctrl_value) && is.na(private$alt_value)){
                       message("\t\tas.numeric for ",private$grp_col)
                       meta[,private$grp_col] <- as.numeric(meta[,private$grp_col])
                     }else{
                       meta[,private$grp_col] <- factor(meta[,private$grp_col],levels=c(private$ctrl_value,private$alt_value))
                     }
                     df = model.matrix(as.formula(paste0("~",paste(colnames(meta),collapse="+"))),data=meta)
                     #print(dim(meta))
                     #print(dim(df))
                     #message('sample size:',sum(private$c_index))
                     re = nebula(private$X[as.vector(private$g_index),as.vector(private$c_index)],
                                 private$meta[private$c_index,private$id_col],
                                 pred=df,offset=colSums(private$X[as.vector(private$g_index),as.vector(private$c_index)]),...)
                     #message("***** ",saveRaw," *****")
                     if(!is.null(saveRaw)) saveRDS(re,saveRaw)
                     final_table <- data.frame(ID=re$summary[,"gene"],re$summary[,grep(private$grp_col,colnames(re$summary))])
                     colnames(final_table) <- c("ID","logFC","se","Pvalue")
                     final_table$log2FC <- log2(exp(final_table$logFC))
                     final_table$FDR <- p.adjust(final_table[,"Pvalue"],method="fdr")
                     final_table$method <- "NEBULA"
                     final_table$algorithm <- re$algorithm
                     first4 <- c("ID","log2FC","Pvalue","FDR")
                     final_table <- final_table[,c(first4,colnames(final_table)[!colnames(final_table)%in%first4])]
                     return(final_table)
                   },
                   # DEG by DESeq2 ---------
                   scDEG_DESeq2=function(...){
                     private$X <- NULL
                     gc()
                     register(MulticoreParam(max(1,parallelly::availableCores()-1)))
                     meta <- private$sInfo
                     meta[,private$grp_col] <- factor(meta[,private$grp_col],levels=c(private$ctrl_value,private$alt_value))
                     #X <- private$pseudoX
                     #save(X,meta,file="test.RData")
                     dds <- DESeqDataSetFromMatrix(countData=private$pseudoX,
                                                   colData=meta,
                                                   design=as.formula(paste0("~",paste(colnames(meta),collapse="+"))))
                     dds <- DESeq(dds,parallel=T)
                     res = lfcShrink(dds, coef = paste(private$grp_col, private$alt_value, "vs", private$ctrl_value, sep = "_"),
                                     format="DataFrame",type = "apeglm", parallel = T)
                     final_table <- data.frame(row.names=NULL,
                                               ID=rownames(res),
                                               log2FC=res$log2FoldChange,
                                               Pvalue=res$pvalue,
                                               FDR=res$padj,
                                               se=res$lfcSE,
                                               method="DESeq2",
                                               algorithm='lfcShrink')
                     return(final_table)
                   },
                   #create pseudo -------
                   createPseudo=function(){
                     #private$pseudoX <- NULL
                     message("Creating pseudo bulk")
                     private$pseudoX <- bplapply(rownames(private$sInfo),function(one){
                       return(setNames(data.frame(rowSums(private$X[as.vector(private$g_index),as.vector(private$c_index&private$meta[,private$id_col]==one)])),
                                       one))
                     },BPPARAM = MulticoreParam(workers=min(30,nrow(private$sInfo),max(1,parallelly::availableCores()-1)),
                                              tasks=nrow(private$sInfo),progressbar=T))
                     private$pseudoX <- as.matrix(dplyr::bind_cols(private$pseudoX))
                     #for(one in rownames(private$sInfo)) private$pseudoX <- cbind(private$pseudoX,rowSums(private$X[private$g_index,private$c_index&private$meta[,private$id_col]==one]))
                     #colnames(private$pseudoX) <- rownames(private$sInfo)
                     mode(private$pseudoX) <- 'integer'
                   },
                   #cell filtering -------
                   cellFiltering=function(sel,msg){
                     message(msg)
                     private$c_index <- private$c_index & sel
                     message("\t",sum(private$c_index)," cells")
                   }
                 ))
