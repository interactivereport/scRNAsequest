#!/usr/bin/env python

import subprocess, os, h5py, sys, warnings, re, yaml
from scipy import sparse
from scipy.sparse import csc_matrix
import pandas as pd
import anndata as ad

def msgError(msg):
  print(msg)
  exit()

def main():
  print("starting SCT+Harmony ...")
  batchKey="library_id"
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

  strCSV = "%s.csv"%os.path.join(config["output"],"sctHarmony",config["prj_name"])#strH5ad.replace("raw.h5ad","sctHarmony.csv")
  cmd = "Rscript %s %s %s %s"%(os.path.join(os.path.dirname(os.path.realpath(__file__)),"sctHarmony.R"),
                            strH5ad,strCSV,strConfig)
  subprocess.run(cmd,shell=True,check=True)
  if not os.path.isfile(strCSV):
    msgError("".join(cmdR.stdout.decode("utf-8"))+"\nERROR: sctHarmony failed!")

  meta = pd.read_csv(strCSV,index_col=0,header=0)
  meta.index = list(meta.index)
  for one in meta.columns:
    if meta[one].nunique()<100:
      meta[one]=meta[one].astype("category")
    
  print("\tsctHarmony R completed")
  FakeD = pd.DataFrame({"FakeG%d"%i:[0 for j in range(meta.shape[0])] for i in range(2)},
                      index=meta.index)
  D = ad.AnnData(FakeD)
  selMeta = [one for one in meta.columns if 'cluster' in one.lower()]+[batchKey]
  D.obs = meta[selMeta]
  for one in set([one.rsplit("_",1)[0] for one in [one for one in meta.columns if not one in selMeta]]):
    D.obsm['X_%s'%one] = meta[[a for a in meta.columns if a.startswith(one)]].values
  print("\tsaving ...")
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    D.write(re.sub("csv$","h5ad",strCSV))
  print("sctHarmony process completed!")

if __name__ == "__main__":
  main()
