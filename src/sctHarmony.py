#!/usr/bin/env python

import subprocess, os, h5py, sys, warnings, re, yaml, logging, glob, functools, random, time,resource,gc,datetime
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
  cmd = "Rscript %s SCT %s %s %s |& tee %s.log"%(os.path.join(strPipePath,"sctHarmony.R"),
                              oneH5ad,strSCT,strConfig,strSCT)
  subprocess.run(cmd,shell=True,check=True)

def sct(strH5ad,strConfig,strPCA,batchCell,hvgN,subCore=5):
  if os.path.isfile(strPCA):
    print("\tUsing previous sct PCA results: %s\n***=== Important: If a new run is desired, please remove/rename the above file "%strPCA)
    return()
  h5adList = sorted(bU.splitBatch(strH5ad,os.path.join(os.path.dirname(strPCA),"tmp"),batchCell,batchKey,hvgN))
  if len(h5adList)==0:
    msgError("No h5ad!")
  print(datetime.datetime.now(),": There are total of %d batches"%len(h5adList))
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
  print(datetime.datetime.now(),": Reading batches ...")
  g=[]
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
    		g.append(np.array(f['var_name'],dtype='str'))
    print(datetime.datetime.now(),": ***** finishing  %d cells and %d genes *****"%(oneD.shape[0],oneD.shape[1]))
    Dlist.append(oneD)
    cN += oneD.shape[0]
    del oneD
    gc.collect()
    print(datetime.datetime.now(),": \tTotal %d cells\tPeak memory %.2fG"%(cN,resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2))
  print(datetime.datetime.now(),": Merging SCT batches ...")
  sctD=ad.concat(Dlist,join="outer")
  del Dlist
  gc.collect()
  g = list(functools.reduce(set.intersection, map(set, g)))
  g = random.sample(g,k=min(len(g),int((2**31-1)/sctD.shape[0])))
  tmpD = ad.AnnData(sctD[:,g].to_df())
  tmpD.X = csc_matrix(tmpD.X)
  sctD.raw = tmpD
  del tmpD
  gc.collect()
  print(datetime.datetime.now(),": \tPeak memory %.2fG"%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2))

  sctD.X[np.isnan(sctD.X)] = 0
  print(datetime.datetime.now(),": Total: %d cells and %d genes"%(sctD.shape[0],sctD.shape[1]))
  batchV=sc.read_h5ad(strH5ad,backed="r").obs[batchKey].copy()
  sctD.obs[batchKey]=batchV[sctD.obs.index]
  print(datetime.datetime.now(),": PCA ...")
  sc.tl.pca(sctD,n_comps=50,svd_solver='arpack')
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    sctD.write(strPCA)
  print(datetime.datetime.now(),": \tThe sct PCA step for all samples are completed with peak memory %.2fG!"%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2))
  return None

def runRharmony(strPCA,strH5,strConfig):
  cmd = "Rscript %s Harmony %s %s %s |& tee %s/sctHarmony.log"%(os.path.join(strPipePath,"sctHarmony.R"),
                              strPCA,strH5,strConfig,os.path.dirname(strH5))
  subprocess.run(cmd,shell=True,check=True)

def sctHarmony(strH5ad,strConfig,config):
	batchCell = config.get('batchCell')
	hvgN=config.get('harmonyBatchGene')
	subCore=5 if config.get('subprocess') is None else config.get('subprocess')
	strHarmony = "%s.h5ad"%os.path.join(config["output"],"sctHarmony",config["prj_name"])#strH5ad.replace("raw.h5ad","sctHarmony.csv")
	if os.path.isfile(strHarmony):
		print(datetime.datetime.now(),": Using previous harmony results: %s\n\tIf a new run is needed, please rename/remove the above file"%strHarmony)
		return
	strPCA = re.sub("h5ad$","pca.h5ad",strHarmony)
	sct(strH5ad,strConfig,strPCA,batchCell,hvgN,subCore)
	strHarmonyPCA = re.sub("h5ad","pca.h5",strHarmony)
	runRharmony(strPCA,strHarmonyPCA,strConfig)
	if not os.path.isfile(strHarmonyPCA):
		msgError("Harmony failed in R")
	with h5py.File(strHarmonyPCA,'r') as f:
		PCA = np.array(f['PCA']).transpose()
		cID = np.array(f['obs_name'])
	print(datetime.datetime.now(),": \tsctHarmony R completed")
	D = ad.AnnData(pd.DataFrame({"FakeG%d"%i:[0 for j in range(PCA.shape[0])] for i in range(2)},
		index=cID))
	D.obsm['X_sctHarmony_pca'] = PCA
	print(datetime.datetime.now(),": \tFinding Neighbors")
	sc.pp.neighbors(D,n_neighbors=20,use_rep='X_sctHarmony_pca')
	print(datetime.datetime.now(),": \tUMAP")
	sc.tl.umap(D)
	cMethod = "louvain" if config.get('clusterMethod') is None else config.get('clusterMethod')
	cResolution = 0.8 if config.get('clusterResolution') is None else config.get('clusterResolution')
	print(datetime.datetime.now(),": Clustering %s (%f)"%(cMethod,cResolution))
	if bool(re.search('leiden',cMethod,re.IGNORECASE)):
		cMethod = 'leiden'
		sc.tl.leiden(D,resolution=cResolution)
	elif bool(re.search('louvain',cMethod,re.IGNORECASE)):
		cMethod = 'louvain'
		sc.tl.louvain(D,resolution=cResolution)
	else:
		msgError("unknow clustering method (leiden or louvain): %s"%cMethod)
	D.obs.columns = [re.sub(cMethod,'sctHarmony_cluster',_,flags=re.IGNORECASE) for _ in D.obs.columns]
	D.obsm['X_sctHarmony_umap'] = D.obsm["X_umap"]
	del D.obsm["X_umap"]
	D1 = ad.read_h5ad(strPCA,backed='r')
	D.obs[batchKey] = D1.obs[batchKey].copy()
	D.write(strHarmony)
	print(datetime.datetime.now(),": sctHarmony process completed!")

def main():
  print(datetime.datetime.now(),": starting SCT+Harmony ...")
  if len(sys.argv)<2:
    msgError("ERROR: raw h5ad file and the config file are required!")
  config = bU.inputCheck(sys.argv)
  strH5ad = sys.argv[1]
  strConfig=sys.argv[2]
  sctHarmony(strH5ad,strConfig,config)

if __name__ == "__main__":
  start_time = time.time()
  main()
  print("---sctHarmony: total time %s seconds ---" % (time.time() - start_time))
