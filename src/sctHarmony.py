#!/usr/bin/env python

import subprocess, os, h5py, sys, warnings, re, yaml, logging, glob, functools, random, time,resource,gc
import cmdUtility as cU
from scipy import sparse
from scipy.sparse import csc_matrix
import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np
import batchUtility as bU
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
readRDS = robjects.r['readRDS']

print=functools.partial(print, flush=True)
warnings.simplefilter(action='ignore', category=FutureWarning)
logging.disable(level=logging.INFO)
strPipePath=os.path.dirname(os.path.realpath(__file__))
batchKey="library_id"
def msgError(msg):
  print(msg)
  exit()

def runOneSCT(oneH5ad,strConfig,strSCT):
  print("***** process: %s *****"%os.path.basename(oneH5ad))
  cmd = "Rscript %s SCT %s %s %s |& tee %s/sctHarmony.log"%(os.path.join(strPipePath,"sctHarmony.R"),
                              oneH5ad,strSCT,strConfig,os.path.dirname(strSCT))
  subprocess.run(cmd,shell=True,check=True)

def sct(strH5ad,strConfig,strPCA,batchCell,hvgN,subCore=5):
  if os.path.isfile(strPCA):
    print("\tUsing previous sct PCA results: %s\n***=== Important: If a new run is desired, please remove/rename the above file "%strPCA)
    return()
  h5adList = sorted(bU.splitBatch(strH5ad,os.path.join(os.path.dirname(strPCA),"tmp"),batchCell,batchKey,hvgN))
  if len(h5adList)==0:
    msgError("No h5ad!")
  print("There are total of %d batches"%len(h5adList))
  Dlist=[]
  cN=0
  sctD=None
  parallelCMD = []
  strSCTsuffix = ".h5" # ".rds"
  for oneH5ad in h5adList:
    strSCT = re.sub(".h5ad$",strSCTsuffix,oneH5ad)
    if not os.path.isfile(strSCT):
      parallelCMD.append(functools.partial(runOneSCT,oneH5ad,strConfig,strSCT))
  if len(parallelCMD)>0:
    metaList=cU.parallel_cmd(parallelCMD,min(subCore,len(parallelCMD)))
  print("Reading batches ...")
  for oneH5ad in h5adList:
    strSCT = re.sub(".h5ad$",strSCTsuffix,oneH5ad)
    print("\t",os.path.basename(strSCT))
    if not os.path.isfile(strSCT):
      msgError("\tERROR: %s sctHarmony failed in SCT step!"%os.path.basename(oneH5ad))
    if strSCT.endswith('rds'):
    	oneD=ad.AnnData(pandas2ri.rpy2py_dataframe(readRDS(strSCT)))
    elif strSCT.endswith('h5'):
    	with h5py.File(strSCT,'r') as f:
    		oneD=ad.AnnData(pd.DataFrame(np.array(f['X']),columns=np.array(f['var_name'],dtype='str'),index=np.array(f['obs_name'],dtype='str')))
    print("***** finishing  %d cells and %d genes *****"%(oneD.shape[0],oneD.shape[1]))
    Dlist.append(oneD)
    cN += oneD.shape[0]
    del oneD
    gc.collect()
    print("\tTotal %d cells\tPeak memory %.2fG"%(cN,resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2))
  print("Merging SCT batches ...")
  sctD=ad.concat(Dlist,join="outer")
  del Dlist
  gc.collect()
  print("\tPeak memory %.2fG"%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2))

  sctD.X[np.isnan(sctD.X)] = 0
  print("Total: %d cells and %d genes"%(sctD.shape[0],sctD.shape[1]))
  batchV=sc.read_h5ad(strH5ad,backed="r").obs[batchKey].copy()
  sctD.obs[batchKey]=batchV[sctD.obs.index]
  print("PCA ...")
  sc.tl.pca(sctD,n_comps=50,svd_solver='arpack')
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    sctD.write(strPCA)
  print("\tThe sct PCA step for all samples are completed with peak memory %.2fG!"%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2))
  return None

def runRharmony(strPCA,strMeta,strConfig):
  clusterResolution=0.8 if clusterResolution is None else clusterResolution
  clusterMethod='Louvain' if clusterMethod is None else clusterMethod
  cmd = "Rscript %s Harmony %s %s %s |& tee %s/sctHarmony.log"%(os.path.join(strPipePath,"sctHarmony.R"),
                              strPCA,strMeta,strConfig,os.path.dirname(strMeta))
  subprocess.run(cmd,shell=True,check=True)

def sctHarmony(strH5ad,strConfig,strMeta,batchCell,hvgN,clusterResolution,clusterMethod,subCore):
  if os.path.isfile(strMeta):
    print("Using previous sctHarmony results: %s\n***=== Important: If a new run is desired, please remove/rename the above file "%strMeta)
    meta = pandas2ri.rpy2py_dataframe(readRDS(strMeta))
    return(meta)
  strPCA = re.sub("rds$","pca.h5ad",strMeta)
  sct(strH5ad,strConfig,strPCA,batchCell,hvgN,subCore)
  runRharmony(strPCA,strMeta,strConfig)
  if not os.path.isfile(strMeta):
    msgError("\tERROR: sctHarmony failed in final harmony step!")
  meta = pandas2ri.rpy2py_dataframe(readRDS(strMeta))
  print("sctHarmony completed")
  return(meta)

def main():
  print("starting SCT+Harmony ...")
  if len(sys.argv)<2:
    msgError("ERROR: raw h5ad file and the config file are required!")
  config = bU.inputCheck(sys.argv)
  strH5ad = sys.argv[1]
  strConfig=sys.argv[2]

  strMeta = "%s.rds"%os.path.join(config["output"],"sctHarmony",config["prj_name"])#strH5ad.replace("raw.h5ad","sctHarmony.csv")
  meta = sctHarmony(strH5ad,strConfig,strMeta,config.get('batchCell'),config.get('harmonyBatchGene'),
    config.get('clusterResolution'),config.get('clusterMethod'),
    subCore=5 if config.get('subprocess') is None else config.get('subprocess'))

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
    D.obsm['X_sctHarmony_%s'%one] = meta[[a for a in meta.columns if a.startswith(one)]].values
  print("\tsaving ...")
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    D.write(re.sub("rds$","h5ad",strMeta))
  print("sctHarmony process completed!")

if __name__ == "__main__":
  start_time = time.time()
  main()
  print("---sctHarmony: total time %s seconds ---" % (time.time() - start_time))
