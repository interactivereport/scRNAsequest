#!/usr/bin/env python
# coding: utf-8
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.metrics import silhouette_samples
import os, re, sys, multiprocessing
import matplotlib.pyplot as plt

def msgError(msg):
  print(msg)
  exit()

def calc_sil_pca50(prefix,oneM):
  """
  this will be using euclidean distances in the PCA space
  """
  #strH5ad = "%s_%s.h5ad"%(prefix,oneM)
  strH5ad="%s.h5ad"%os.path.join(os.path.dirname(prefix),oneM,os.path.basename(prefix))
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

def cal_sil_pca50_one(pcaKey,strH5ad):
  k = re.sub("^X_|_pca","",pcaKey)
  if k=='pca':
    k="raw"
  sil_coeff = None
  adata = sc.read(strH5ad,backed=True)
  obs = adata.obs.copy()
  X = adata.obsm[pcaKey].copy()
  if np.count_nonzero(np.abs(X).sum(axis=1)<0.001) > X.shape[0]/4:
    print("--> More than a quarter cells with zero in %s PCA: SKIP! <--"%k)
  else:
    cKey=[one for one in obs.columns if re.search('^%s.*cluster$|^%s.*louvain$'%(k,k),one,re.IGNORECASE)]
    if len(cKey)>0:
      sil_coeff = silhouette_samples(X=X[:, :50], labels=np.array(obs[cKey[0]].values))
    print(print("\tFinishing %s: %s"%(pcaKey,cKey[0])))
  return k,sil_coeff

def make_df(i):
  ## make individual coeff_pca50 to a pandas df
  df0 = pd.DataFrame(i[1])
  df0.columns = ['Silhouette_coefficients']
  #method = i[0].replace("concat12_", "").replace("clustered.h5ad", "")
  df0['method'] = i[0]
  return(df0)
def setupDir(strOut):
  try:
    os.makedirs(strOut)
  except FileExistsError:
    pass

def main():
  print("starting silhouette evaluation ...")
  if len(sys.argv)<3:
    msgError("ERROR: project final h5ad and output are required!")
  strH5ad = sys.argv[1]
  strPDF = sys.argv[2]
  
  PCAkey = [(one,strH5ad) for one in sc.read(strH5ad,backed=True).obsm.keys() if '_pca' in one]
  print("\tworking on",", ".join([a[0] for a in PCAkey]))
  with multiprocessing.Pool(processes=len(PCAkey)) as pool:
    coeff_pca50 = pool.starmap(cal_sil_pca50_one,PCAkey)
  coeff_pca50 = [one for one in coeff_pca50 if one[1] is not None]

  coeff_pca50_dfs = list(map(make_df, coeff_pca50))
  ## bind
  coeff_pca50_df = pd.concat(coeff_pca50_dfs)
  # boxplot
  ax = coeff_pca50_df.boxplot(by='method',rot=90)
  ax.set_title("Silhouette_coefficients")
  plt.grid()
  #strPDF = "%s_Silhouette_boxplot_pc50.pdf"%os.path.join(os.path.dirname(prefix),'evaluation',os.path.basename(prefix))
  setupDir(os.path.dirname(strPDF))
  plt.savefig(strPDF,bbox_inches="tight")
  plt.close()

if __name__ == "__main__":
  main()
