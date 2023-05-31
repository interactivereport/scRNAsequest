#!/usr/bin/env python

import subprocess, os, h5py, sys, warnings, yaml, re
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
  #if not strH5ad.endswith("raw.h5ad"):
  #  msgError("ERROR: %s is not raw h5ad file required!"%strH5ad)

  strConfig = sys.argv[2]
  if not os.path.isfile(strConfig):
    msgError("ERROR: %s does not exist!"%strConfig)
  with open(strConfig,"r") as f:
    config = yaml.safe_load(f)
  strOut = "%s.csv.gz"%os.path.join(config["output"],"SeuratRPCA",config["prj_name"])#strH5ad.replace("raw.h5ad","seurat_rpca.csv.gz")

  cmd = "Rscript %s %s %s |& tee %s/SeuratRPCA.log"%(os.path.join(strPipePath,"seuratRPCA.R"),
                            strH5ad,strOut,os.path.dirname(strOut))
  if os.path.isfile(strOut):
    print("Using previous SeuratRPCA results: %s\n***=== Important: If a new run is desired, please remove/rename the above file "%strOut)
  else:  
    subprocess.run(cmd,shell=True,check=True)#,stdout=subprocess.PIPE
  if not os.path.isfile(strOut):
    msgError("\tERROR: SeuratRPCA failed!")

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
    D.obsm['X_seuratRPCA_%s'%one] = meta[[a for a in meta.columns if a.startswith(one)]].values
  print("\tsaving ...")
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    D.write(re.sub("csv.gz","h5ad",strOut))
  print("=== seurat RPCA integration process completed! ===")

if __name__ == "__main__":
  main()
