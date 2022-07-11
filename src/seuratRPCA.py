#!/usr/bin/env python

import subprocess, os, h5py, sys, warnings
import anndata as ad
import pandas as pd
import numpy as np
import seaborn as sns
import time

def msgError(msg):
  print(msg)
  exit()

def main():
  strPipePath=os.path.dirname(os.path.realpath(__file__))
  print("starting seurat fast integration (RPCA) ...")
  if len(sys.argv)<1:
    msgError("ERROR: raw h5ad file required!")
  strH5ad = sys.argv[1]
  if not os.path.isfile(strH5ad):
    msgError("ERROR: %s does not exist!"%strH5ad)
  if not strH5ad.endswith("raw.h5ad"):
    msgError("ERROR: %s is not raw h5ad file required!"%strH5ad)

  strOut = strH5ad.replace("raw.h5ad","seurat_rpca.csv.gz")

  print("\tRPCA integration ...")
  cmd = "Rscript %s %s %s"%(os.path.join(strPipePath,"seuratRPCA.R"),
                            strH5ad,strOut)
  subprocess.run(cmd,shell=True,check=True)
  
  D = ad.read_h5ad(strH5ad,backed="r")
  Dbatch = D.obs["library_id"].copy()
  meta = pd.read_csv(strOut,index_col=0,header=0)
  meta.index = list(meta.index)
  for one in meta.columns:
    if meta[one].nunique()<100:
      meta[one]=meta[one].astype("category")
  print("\tRPCA integration completed")
  FakeD = pd.DataFrame(0,columns=["FakeG1","FakeG2"],index=meta.index)
  D = ad.AnnData(FakeD)
  D.obs = pd.concat([meta[[one for one in meta.columns if one.startswith("seuratRPCA")]],Dbatch],axis=1,join="inner")
  for one in set([one.rsplit("_",1)[0] for one in [one for one in meta.columns if not one.startswith("seuratRPCA")]]):
    D.obsm['X_%s'%one] = meta[[a for a in meta.columns if a.startswith(one)]].values
  print("\tsaving ...")
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    D.write(strH5ad.replace("raw.h5ad","SeuratRPCA.h5ad"))
  print("=== seurat RPCA integration process completed! ===")

if __name__ == "__main__":
  main()
