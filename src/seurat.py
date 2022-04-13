#!/usr/bin/env python

import subprocess, os, h5py, sys, warnings
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import time
import rpy2
from rpy2 import robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()

sc.settings.verbosity = 3
#sc.logging.print_versions()
sc.settings.n_jobs = 20

def msgError(msg):
  print(msg)
  exit()
# read .RData file as a pandas dataframe
def load_rdata_file(filename):
    r_data = robjects.r['get'](robjects.r['load'](filename))
    #df = pandas2ri.rpy2py(r_data)
    df = r_data
    return df

# write pandas dataframe to an .RData file
def save_rdata_file(df, filename):
    r_data = pandas2ri.py2rpy(df)
    robjects.r.assign("my_df", r_data)
    robjects.r("save(my_df, file='{}')".format(filename))

def main():
  print("starting seurat process ...")
  if len(sys.argv)<1:
    msgError("ERROR: raw h5ad file required!")
  strH5ad = sys.argv[1]
  if not os.path.isfile(strH5ad):
    msgError("ERROR: %s does not exist!"%strH5ad)
  if not strH5ad.endswith("raw.h5ad"):
    msgError("ERROR: %s is not raw h5ad file required!"%strH5ad)

  #strExp = strH5ad.replace("raw.h5ad","seurat_py2r_adata_expr.npy")
  #strBatch = strH5ad.replace("raw.h5ad","seurat_py2r_adata_batch.csv.gz")
  #strVar = strH5ad.replace("raw.h5ad","seurat_py2r_adata_var.csv.gz")
  
  strOut = strH5ad.replace("raw.h5ad","seurat_r2py_integrated.csv.gz")
  strGene = strH5ad.replace("raw.h5ad","seurat_r2py_integrated.genes.csv.gz")

  adata = sc.read_h5ad(strH5ad)
  #adata_expr = pd.DataFrame.sparse.from_spmatrix(data=adata.X.transpose(), index=adata.var_names, columns=adata.obs_names)
  #adata_batch = pd.DataFrame(data=adata.obs.library_id,index=adata.obs_names)
  #adata_batch = adata_batch.astype('str') ## need to convert to string or rpy2 will get error "TypeError: Parameter 'categories' must be list-like", when running at rpy2=2.9.4, pandas=0.24.2, rpy2 3.0 will solve the problem according to here: https://bitbucket.org/rpy2/rpy2/issues/509/latest-pandas-breaking-rpy2-conversion, however wasn't able to update because of dependencies
  ## write to disk for Seurat in R
  #np.save(strExp, adata.X.astype(np.float64))
  # adata_expr.to_csv("working_data/concat_4seurat_adata_expr.csv.gz") ## this is very slow
  #adata_batch.to_csv(strBatch)
  ## also write gene names
  #adata.var.to_csv(strVar)
  # the following is NOT used O'Young
  #save_rdata_file(pd.DataFrame(adata_expr), "seurat_py2r_adata_expr.Rdata")
  #save_rdata_file(pd.DataFrame(adata_batch), "seurat_py2r_adata_batch.Rdata")
  #cmd = "Rscript %s %s %s %s %s %s"%(os.path.join(os.path.dirname(os.path.realpath(__file__)),"seurat.R"),
  #                          strExp,strBatch,strVar,strOut,strGene)
  print("\tintegration ...")
  cmd = "Rscript %s %s %s %s"%(os.path.join(os.path.dirname(os.path.realpath(__file__)),"seurat.R"),
                            strH5ad,strOut,strGene)
  subprocess.run(cmd,shell=True,check=True)
  
  adata_expr_integrated = pd.read_csv(strOut, header = 0, index_col=False)
  print("\tintegration completed")
  adata_expr_integrated.columns = [i.replace(".", "-") for i in adata_expr_integrated.columns]
  set(adata_expr_integrated.columns) == set(adata.obs_names)
  ## by default suerat only integrate / correct the expr of 2000 genes, we have set to 20000 in the above block to let it integrate all genes
  ## read gene names
  adata_expr_integrated_genes = pd.read_csv(strGene, header = 0, index_col = False)
  adata_expr_integrated['gene'] = adata_expr_integrated_genes.iloc[:, 0]
  adata_expr_integrated = adata_expr_integrated.set_index('gene')
  adata_expr_integrated = adata_expr_integrated.transpose()
  #adata_expr_integrated
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ## make a dummy obj and run scale and PCA, then stuck the PCA embeddings to adata
  
    bdata = sc.AnnData(X=adata_expr_integrated, obs=adata.obs)
    ## scale each gene to unit variance, clip values exceeding sd 10
    sc.pp.scale(bdata, max_value=10)
    ## PCA
    ## this is multi-processed
    print("\tPCA ...")
    sc.tl.pca(bdata, svd_solver='arpack', n_comps = 100)
  
    ## inspect the variance explained by each PC
    #sc.pl.pca_variance_ratio(bdata, log=False, n_pcs=100)
  
    ## copy these two fields to adata
    adata.uns = bdata.uns
    adata.obsm = bdata.obsm

    ## do the clustering
    ## clustering
    ## this is using the regular neighbors function from scanpy
    npc = 50
    print("\tclustering ...")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=npc)
    sc.tl.louvain(adata,key_added="Seurat.louvain")
    print("\tembedding ...")
    sc.tl.umap(adata)
    sc.tl.rank_genes_groups(adata, 'Seurat.louvain')
    #sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
    sc.tl.tsne(adata, n_pcs=npc)
    print("\tsaving ...")
    adata.write_h5ad(strH5ad.replace("raw.h5ad","Seurat.h5ad"))
  print("=== seurat integration process completed! ===")

if __name__ == "__main__":
  main()
