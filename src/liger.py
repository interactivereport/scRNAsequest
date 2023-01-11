#!/usr/bin/env python

import subprocess, os, h5py, sys, warnings, math, re, yaml
import anndata as ad
import scanpy as sc
import numpy as np
from scipy import sparse
from scipy.sparse import csc_matrix
import pandas as pd

def msgError(msg):
  print(msg)
  exit()

def main():
  print("starting liger ...")
  if len(sys.argv)<1:
    msgError("ERROR: raw h5ad file required!")
  strH5ad = sys.argv[1]
  if not os.path.isfile(strH5ad):
    msgError("ERROR: %s does not exist!"%strH5ad)
  #if not strH5ad.endswith("raw.h5ad"):
  #  msgError("ERROR: %s is not raw h5ad file required!"%strH5ad)
  
  strConfig = sys.argv[2]
  if not os.path.isfile(strConfig):
    msgError("ERROR: %s does not exist!"%strConfig)
  with open(strConfig,"r") as f:
    config = yaml.safe_load(f)

  print("\tcalculate highly variable genes")
  strHVG= "%s_hvg.csv"%os.path.join(config["output"],"Liger",config["prj_name"])#strH5ad.replace("raw.h5ad","liger_hvg.csv")
  D = sc.read_h5ad(strH5ad)
  # 95 percentile to normalize
  sc.pp.normalize_total(D,target_sum=math.ceil(np.percentile(D.X.sum(axis=1).transpose().tolist()[0],95)))
  sc.pp.log1p(D)
  sc.pp.highly_variable_genes(D, min_mean=0.01, max_mean=3, min_disp=0.5)
  D.var.highly_variable.to_csv(strHVG)
  Dbatch = D.obs["library_id"].copy()
  
  strCSV = re.sub("_hvg.csv$",".csv",strHVG)#strH5ad.replace("raw.h5ad","liger.csv")
  cmd = "Rscript %s %s %s %s"%(os.path.join(os.path.dirname(os.path.realpath(__file__)),"liger.R"),
                            strH5ad,strHVG,strCSV)
  subprocess.run(cmd,shell=True,check=True)

  meta = pd.read_csv(strCSV,index_col=0,header=0)
  meta.index = list(meta.index)
  for one in meta.columns:
    if meta[one].nunique()<100:
      meta[one]=meta[one].astype("category")
    
  print("\tliger R completed")
  FakeD = pd.DataFrame({"FakeG%d"%i:[0 for j in range(meta.shape[0])] for i in range(2)},
                      index=meta.index)
  D = ad.AnnData(FakeD)
  D.obs = pd.concat([meta[[one for one in meta.columns if one.startswith("liger")]],Dbatch],axis=1,join="inner")
  for one in set([one.rsplit("_",1)[0] for one in [one for one in meta.columns if not one.startswith("liger")]]):
    D.obsm['X_%s'%one] = meta[[a for a in meta.columns if a.startswith(one)]].values
  print("\tsaving ...")
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    D.write(re.sub("_hvg.csv$",".h5ad",strHVG))
  print("liger process completed!")

if __name__ == "__main__":
  main()
