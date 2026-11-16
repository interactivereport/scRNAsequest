# Parse Biosciences output folder: count matrix (DGE.mtx in older, count_matrix.mtx in newer versions),
# all_genes.csv and cell_metadata.csv. Kept in line with findParse_mtx/getParse_mtx in obtainData.py.
findParse <- function(strDir){
  strMTX <- file.path(strDir,c("DGE.mtx","count_matrix.mtx"))
  strMTX <- strMTX[file.exists(strMTX)]
  if(length(strMTX)==0) return(NULL)
  return(strMTX[1])
}
# Returns a genes x cells sparse matrix. On disk:
# the matrix is cells x genes (used without transpose check beyond dimensions), cell names are sample_bc_wells, gene names are gene_name (made unique)
readParse <- function(strMTX){
  strDir <- dirname(strMTX)
  msg <- c()
  for(one in c("all_genes.csv","cell_metadata.csv")) if(!file.exists(file.path(strDir,one))) msg <- c(msg,paste("missing file:",one))
  if(length(msg)>0) stop(paste0("Incompatible Parse Biosciences input folder ",strDir,":\n\t",paste(msg,collapse="\n\t")))
  gInfo <- data.table::fread(file.path(strDir,"all_genes.csv"),data.table=F)
  cInfo <- data.table::fread(file.path(strDir,"cell_metadata.csv"),data.table=F)
  if(!"gene_name"%in%colnames(gInfo)) msg <- c(msg,paste0("all_genes.csv: missing column 'gene_name' (found: ",paste(colnames(gInfo),collapse=", "),")"))
  missCol <- setdiff(c("sample","bc_wells"),colnames(cInfo))
  if(length(missCol)>0) msg <- c(msg,paste0("cell_metadata.csv: missing column(s) '",paste(missCol,collapse="', '"),"' (found: ",paste(colnames(cInfo),collapse=", "),")"))
  if(length(msg)>0) stop(paste0("Incompatible Parse Biosciences input folder ",strDir,":\n\t",paste(msg,collapse="\n\t")))
  X <- Matrix::readMM(strMTX)
  if(nrow(X)!=nrow(cInfo) || ncol(X)!=nrow(gInfo))
    stop(paste0(basename(strMTX)," is ",nrow(X)," x ",ncol(X),", but ",nrow(cInfo)," cells (cell_metadata.csv) x ",
                nrow(gInfo)," genes (all_genes.csv) is expected"))
  X <- as(Matrix::t(X),"CsparseMatrix")
  dimnames(X) <- list(make.unique(as.character(gInfo$gene_name),sep="-"),
                      paste0(cInfo$sample,"_",cInfo$bc_wells))
  return(X)
}
