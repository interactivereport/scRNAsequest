#!/usr/bin/env python

import subprocess, os, h5py, sys, warnings, logging, yaml, re, datetime, glob, functools
import cmdUtility as cU
import anndata as ad
from scipy import sparse
from scipy.sparse import csc_matrix
import batchUtility as bU
import pyarrow.feather as feather
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
readRDS = robjects.r['readRDS']

print=functools.partial(print, flush=True)
logging.disable(level=logging.INFO)
warnings.simplefilter(action='ignore', category=FutureWarning)
strPipePath=os.path.dirname(os.path.realpath(__file__))
batchKey="library_id"

def msgError(msg):
  print(msg)
  exit()

def runOneBatch(oneH5ad,strConfig,oneMeta):
  print("***** processing: %s *****"%os.path.basename(oneH5ad))
  cmd = "Rscript %s %s %s %s |& tee %s"%(os.path.join(strPipePath,"seuratRef.R"),
                            oneH5ad,strConfig,oneMeta,re.sub("rds$","log",oneMeta))
  #print(cmd)
  subprocess.run(cmd,shell=True,check=True)
  return(oneMeta)
  
def batchRef(strH5ad,strConfig,strMeta,batchCell,subCore=5):
  if os.path.isfile(strMeta):
    print("Using previous SeuratRef results: %s\n***=== Important: If a new run is desired, please remove/rename the above file "%strMeta)
    meta = feather.read_feather(strMeta)
    return(meta)
  h5adList = sorted(bU.splitBatch(strH5ad,os.path.join(os.path.dirname(strMeta),"tmp"),batchCell,batchKey))
  if len(h5adList)==0:
    msgError("No h5ad!")
  print("There are total of %d batches"%len(h5adList))
  parallelCMD = []
  for oneH5ad in h5adList:
    oneMeta = re.sub("h5ad$","rds",oneH5ad)
    if not os.path.isfile(oneMeta):
      parallelCMD.append(functools.partial(runOneBatch,oneH5ad,strConfig,oneMeta))
  if len(parallelCMD)>0:
    metaList=cU.parallel_cmd(parallelCMD,min(subCore,len(parallelCMD)))
  mapInfo=[]
  print("Reading batch results ...")
  for oneH5ad in h5adList:
    print("\t%s"%os.path.basename(oneH5ad))
    oneMeta = re.sub("h5ad$","rds",oneH5ad)
    if not os.path.isfile(oneMeta):
      msgError("\tERROR: %s SeuratRef failed!"%os.path.basename(oneH5ad))
    oneMap = pandas2ri.rpy2py_dataframe(readRDS(oneMeta))
    mapInfo.append(oneMap)
  meta = pd.concat(mapInfo)
  meta.index = list(meta.index)
  for col in meta.columns:
    if meta[col].dtype=='float64':
      meta[col]=meta[col].astype('float32')
  feather.write_feather(meta,strMeta)
  print("mapping completed")
  return(meta)

def main():
  print("starting seurat reference mapping ...")
  if len(sys.argv)<2:
    msgError("ERROR: raw h5ad file and the config file are required!")
  config = bU.inputCheck(sys.argv)
  if config == False :
    return()
  strH5ad = sys.argv[1]
  strConfig=sys.argv[2]

  D = ad.read_h5ad(strH5ad,backed="r") #,backed=True
  Dbatch = D.obs[batchKey].copy()
  strMeta = "%s.feather"%os.path.join(config["output"],"SeuratRef",config["prj_name"])#strH5ad.replace("raw.h5ad","seuratRef.csv")
  meta = batchRef(strH5ad,strConfig,strMeta,config.get('batchCell'),subCore=5 if config.get('subprocess') is None else config.get('subprocess'))
  
  FakeD = pd.DataFrame({"FakeG%d"%i:[0 for j in range(meta.shape[0])] for i in range(2)},
                      index=meta.index)
  D = ad.AnnData(FakeD)
  annoCol = [one for one in meta.columns if "_predicted." in one or "_mapping." in one]
  D.obs = pd.concat([meta[annoCol],Dbatch],axis=1,join="inner")
  for one in set([one.rsplit("_",1)[0] for one in list(set(meta.columns)-set(annoCol))]):
    D.obsm['X_SeuratRef_%s'%one] = meta[[a for a in meta.columns if a.startswith(one)]].values
  print("\tsaving ...")
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    D.write(re.sub("feather","h5ad",strMeta))
  print("Mapping to reference completed!")

if __name__ == "__main__":
  main()
