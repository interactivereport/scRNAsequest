suppressMessages(suppressWarnings({
	require(Seurat)
	require(SeuratObject)
	require(SeuratDisk)
	require(sctransform)
	source("readH5ad.R")
}))

getSCTransform <- function(strH5ad,batch,core,assayName="SCT",geneN=3000,hvg=NULL,sctAssay=F,only.var.genes=F){
	#message("coreN: ",core)
	X <- getX(strH5ad,batchID=batch,core=core)
	gID <- setNames(rownames(X[[1]]),gsub("_","-",rownames(X[[1]])))
	X <- lapply(X,function(one){
		rownames(one) <- names(gID)
		return(one)
	})
	meta <- getobs(strH5ad)
	message("SCTransform ...")#
	Dlist <- bplapply(names(X),function(id){
											message("\tSCT ",id,"\t",which(names(X)==id),"/",length(X),appendLF=F)
											DD <- CreateSeuratObject(counts=X[[id]],#meta.data=meta[unlist(sapply(X[sID],colnames),use.names=F),]
																							 project=assayName,
																							 meta.data=meta[colnames(X[[id]]),])
											one <- tryCatch({
												SCTransform(DD,vst.flavor="v2",
																		new.assay.name=assayName,
																		n_genes=geneN,
																		variable.features.n=geneN,
																		return.only.var.genes=only.var.genes,
																		verbose=FALSE)
											},error=function(ee){
												message("\t\tError: ",conditionMessage(ee))
												DD <- SCTransform(DD,vst.flavor="v1",
																					return.only.var.genes = only.var.genes,
																					new.assay.name=assayName,
																					verbose=FALSE)
												return(DD)
											})
											rm(DD)
											gc()
											suppressWarnings({
												one[["RNA"]] <- split(one[["RNA"]],f=unlist(one[[batch]],use.names=F))
												if(!sctAssay)
													one[[assayName]] <- split(one[[assayName]],f=unlist(one[[batch]],use.names=F))
											})
											return(one)
										},
										BPPARAM = MulticoreParam(workers=min(core,length(X),max(1,parallelly::availableCores()-2)),
																						 tasks=length(X)))
	intG <- SelectIntegrationFeatures(Dlist,nfeatures=geneN)
	if(!is.null(hvg) && !sctAssay){
		if(length(hvg)>10) intG <- intersect(hvg,Reduce(intersect,lapply(Dlist,function(one)return(gID[rownames(one[[assayName]])]))))
		if(length(intG)<100) stop(paste0("Two few features (",length(intG),")! Please increase the number of harmonyBatchGene (remove tmp folder)!"))
		message("\t\t",length(intG)," features")
		if(!sctAssay)
		Dlist <- lapply(Dlist,function(one){
			return(one[intG])
			})
	}
	message("\tmerging ...")
	if(length(Dlist)==1){
		D <- Dlist[[1]]
	}else{
		D <- merge(Dlist[[1]],Dlist[-1],project="SCT")
	}
	rm(Dlist)
	gc()
	D@misc <- list(gID=gID)
	VariableFeatures(D) <- intG
	return(D)
}