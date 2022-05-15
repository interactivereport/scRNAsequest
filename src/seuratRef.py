#!/usr/bin/env python

import subprocess, os, h5py, sys, warnings, logging, yaml
import anndata as ad
from scipy import sparse
from scipy.sparse import csc_matrix
import pandas as pd

logging.disable()

def msgError(msg):
  print(msg)
  exit()

def main():
  print("starting seurat reference mapping ...")
  if len(sys.argv)<2:
    msgError("ERROR: raw h5ad file and the config file are required!")
  strH5ad = sys.argv[1]
  if not os.path.isfile(strH5ad):
    msgError("ERROR: %s does not exist!"%strH5ad)
  if not strH5ad.endswith("raw.h5ad"):
    msgError("ERROR: %s is not raw h5ad file required!"%strH5ad)
  strConfig = sys.argv[2]
  if not os.path.isfile(strConfig):
    msgError("ERROR: %s does not exist!"%strConfig)
  with open(strConfig,"r") as f:
    config = yaml.safe_load(f)
  if config['ref_name'] is None:
    print("No reference specified, END")
    return
  print(strH5ad)
  D = ad.read_h5ad(strH5ad) #,backed=True
  Dbatch = D.obs["library_id"].copy()

  strCSV = strH5ad.replace("raw.h5ad","seuratRef.csv")
  cmd = "Rscript %s %s %s %s"%(os.path.join(os.path.dirname(os.path.realpath(__file__)),"seuratRef.R"),
                            strH5ad,strConfig,strCSV)
  subprocess.run(cmd,shell=True,check=True)
  
  meta = pd.read_csv(strCSV,index_col=0,header=0)
  meta.index = list(meta.index)
  print("\tmapping completed")
  FakeD = pd.DataFrame({"FakeG%d"%i:[0 for j in range(meta.shape[0])] for i in range(2)},
                      index=meta.index)
  D = ad.AnnData(FakeD)
  D.obs = pd.concat([meta[[one for one in meta.columns if one.startswith("predicted")]],Dbatch],axis=1,join="inner")
  for one in set([one.rsplit("_",1)[0] for one in [one for one in meta.columns if not one.startswith("predicted")]]):
    D.obsm['X_%s'%one] = meta[[a for a in meta.columns if a.startswith(one)]].values
  print("\tsaving ...")
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    D.write(strH5ad.replace("raw.h5ad","SeuratRef.h5ad"))
  print("Mapping to reference completed!")

if __name__ == "__main__":
  main()
