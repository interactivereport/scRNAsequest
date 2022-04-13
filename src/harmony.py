#!/usr/bin/env python

import sys, os, warnings, subprocess
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import rpy2
from rpy2 import robjects
from rpy2.robjects import pandas2ri

pandas2ri.activate()
sc.settings.n_jobs = 20

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
  print("starting Harmony ...")
  if len(sys.argv)<1:
    msgError("ERROR: raw h5ad file required!")
  strH5ad = sys.argv[1]
  if not os.path.isfile(strH5ad):
    msgError("ERROR: %s does not exist!"%strH5ad)
  if not strH5ad.endswith("raw.h5ad"):
    msgError("ERROR: %s is not raw h5ad file required!"%strH5ad)
  
  adata = sc.read_h5ad(strH5ad)
  ## scale each gene to unit variance, clip values exceeding sd 10
  print("\tscaling ...")
  sc.pp.scale(adata, max_value=10)
  ## PCA
  ## this is multi-processed (not multi-processed by O'Young)
  print("\tPCA ...")
  sc.tl.pca(adata, svd_solver='arpack', n_comps = 100)
  npcs = 50
  ## correction
  pca = adata.obsm['X_pca']
  batch = adata.obs['library_id']
  
  strPCA = strH5ad.replace("raw.h5ad","harmony_py2r_pca.Rdata")
  strLib = strH5ad.replace("raw.h5ad","harmony_py2r_batch.Rdata")
  strOut = strH5ad.replace("raw.h5ad","harmony_r2py.Rdata")
  save_rdata_file(pd.DataFrame(pca), strPCA)
  save_rdata_file(pd.DataFrame(batch), strLib)
  
  print("\tharmonization ...")
  cmd = "Rscript %s %s %s %s"%(os.path.join(os.path.dirname(os.path.realpath(__file__)),"harmony.R"),
                            strPCA,strLib,strOut)
  subprocess.run(cmd,shell=True,check=True)
  if not os.path.isfile(strOut):
    msgError("".join(cmdR.stdout.decode("utf-8"))+"\nERROR: harmony failed!")
  
  pca_harmony = load_rdata_file(strOut)
  adata.obsm['X_pca'] = np.array(pca_harmony)
  ## clustering
  print("\tclustering ...")
  sc.pp.neighbors(adata, n_neighbors=10, n_pcs=npcs)
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    sc.tl.louvain(adata,key_added="Harmony.louvain")
    ## umap embedding 
    print("\tembedding ...")
    sc.tl.umap(adata, init_pos='spectral')
    sc.tl.rank_genes_groups(adata, 'Harmony.louvain')
    sc.tl.tsne(adata, n_pcs=npcs)
    print("\tsaving ...")
    adata.write_h5ad(strH5ad.replace("raw.h5ad","Harmony.h5ad"))
  print("=== harmony process completed! ===")

if __name__ == "__main__":
  main()
