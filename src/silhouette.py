#!/usr/bin/env python
# coding: utf-8
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.metrics import silhouette_samples
import os, re, sys
import matplotlib.pyplot as plt

def msgError(msg):
  print(msg)
  exit()

def calc_sil_pca50(prefix,oneM):
  """
  this will be using euclidean distances in the PCA space
  """
  strH5ad = "%s_%s.h5ad"%(prefix,oneM)
  print(strH5ad)
  adata = sc.read(strH5ad,backed=True)
  sil_coeff = None
  if not 'X_pca' in adata.obsm.keys():
    print("%s is missing obsm 'X_pca'"%strH5ad)
    return oneM,sil_coeff
  cKey = [one for one in adata.obs.keys() if 'louvain' in one]
  if len(cKey)==0:
    cKey = [one for one in adata.obs.keys() if 'cluster' in one]
  if len(cKey)>1:
    tmp = [one for one in cKey if oneM in one]
    if len(tmp)==1:
      cKey = tmp
  print("\t%s: %s"%(oneM,cKey[0]))
  sil_coeff = silhouette_samples(X=adata.obsm['X_pca'][:, :50], labels=np.array(adata.obs[cKey[0]].values))
  return oneM,sil_coeff

def make_df(i):
  ## make individual coeff_pca50 to a pandas df
  df0 = pd.DataFrame(i[1])
  df0.columns = ['Silhouette_coefficients']
  #method = i[0].replace("concat12_", "").replace("clustered.h5ad", "")
  df0['method'] = i[0]
  return(df0)

def main():
  print("starting silhouette evaluation ...")
  if len(sys.argv)<3:
    msgError("ERROR: project prefix and methods required!")
  prefix = sys.argv[1]
  methods = sys.argv[2].split(",")
  
  coeff_pca50 = [one for one in [calc_sil_pca50(prefix,m) for m in methods] if one[1] is not None]
  coeff_pca50_dfs = list(map(make_df, coeff_pca50))
  ## bind
  coeff_pca50_df = pd.concat(coeff_pca50_dfs)
  # boxplot
  ax = coeff_pca50_df.boxplot(by='method',rot=90)
  ax.set_title("Silhouette_coefficients")
  plt.grid()
  plt.savefig("%s_Silhouette_boxplot_pc50.pdf"%prefix,bbox_inches="tight")
  plt.close()

if __name__ == "__main__":
  main()
