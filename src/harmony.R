library(harmony)

main <- function(){
  args = commandArgs(trailingOnly=TRUE)
  if(length(args)<3) stop("ERROR in harmony.R: at least 3 input arguments/files needed: PCA, batch, output")
  strPCA <- args[1]
  strBatch <- args[2]
  strOut <- args[3]
  
  load(strPCA)
  pca <- my_df
  
  load(strBatch)
  batch <- as.vector(my_df)
  
  print(Sys.time())
  # pca_harmony <- HarmonyMatrix(pca, batch, do_pca = FALSE, theta=4)
  print(system.time(pca_harmony <- HarmonyMatrix(pca, unlist(batch), do_pca = FALSE)))
  print(Sys.time())
  save(pca_harmony, file=strOut)
}


main()