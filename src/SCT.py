#!/usr/bin/env python

import subprocess, os, h5py, sys, warnings, re, yaml,functools,logging,math
import scanpy as sc
from scipy import sparse
from scipy.sparse import csc_matrix
import pandas as pd
import batchUtility as bU

warnings.simplefilter(action='ignore', category=FutureWarning)
logging.disable(level=logging.INFO)
batchKey="library_id"
strPipePath=os.path.dirname(os.path.realpath(__file__))
print=functools.partial(print, flush=True)

def msgError(msg):
  print(msg)
  exit()
def runOneBatch(oneH5ad,strConfig,strExp,scaleF):
  cmd = "Rscript %s %s %s %s %d |& tee %s"%(os.path.join(strPipePath,"SCT.R"),
                            oneH5ad,strExp,strConfig,scaleF,re.sub("h5$","log",strExp))
  subprocess.run(cmd,shell=True,check=True)
def runMergeSeurat(h5List,seuratF):
  cmd = "Rscript %s %s %s |& tee %s"%(os.path.join(strPipePath,"SCT.R"),
                            ",".join(h5List),seuratF,re.sub("rds$","log",seuratF))
  subprocess.run(cmd,shell=True,check=True)
def saveScale(strH5ad,strConfig,newH5ad):
  with open(strConfig,"r") as f:
    config = yaml.safe_load(f)
  scaleF = 0 if config.get("expScaler") is None else config.get("expScaler")
  D = sc.read_h5ad(strH5ad,backed="r")
  if scaleF<=100 and scaleF>0:
    scaleF = math.ceil(np.percentile(D.obs.n_counts,scaleF)/1000)*1000
  with open(re.sub("h5ad$","scaleF",newH5ad),"w") as f:
    f.write(str(scaleF))
  return scaleF
def batchNorm(strH5ad,strConfig,newH5ad,batchCell):
  if os.path.isfile(newH5ad):
    print("Using previous Normalized Expression: %s\n***=== Important: If a new run is desired, please remove/rename the above file "%newH5ad)
    return
  scaleF = saveScale(strH5ad,strConfig,newH5ad)
  h5adList = sorted(bU.splitBatch(strH5ad,os.path.join(os.path.dirname(newH5ad),"tmp"),batchCell,batchKey))
  if len(h5adList)==0:
    msgError("No h5ad!")
  print("There are total of %d batches"%len(h5adList))
  Exp=None
  h5List = []
  for oneH5ad in h5adList:
    print("***** batch: %s *****"%os.path.basename(oneH5ad))
    oneH5 = re.sub("h5ad$","h5",oneH5ad)
    if not os.path.isfile(oneH5):
      runOneBatch(oneH5ad,strConfig,oneH5,scaleF)
    if not os.path.isfile(oneH5):
      msgError("\tERROR: %s Expression normalization failed!"%os.path.basename(oneH5ad))
    print("\tcreating h5ad from normlized h5 ...")
    h5List.append(oneH5)
    f = h5py.File(oneH5,'r')
    X = sparse.csc_matrix((f['data'],f['indices'],f['indptr']),f['shape'])
    f.close()
    cID = pd.read_csv(re.sub("h5$","cID",oneH5),header=None,index_col=0)
    gID = pd.read_csv(re.sub("h5$","gID",oneH5),header=None,index_col=0)
    oneD = sc.AnnData(X)
    oneD.obs_names = list(cID.index)
    oneD.var_names = list(gID.index)
    print("***** finishing  %d cells and %d genes *****"%(oneD.shape[0],oneD.shape[1]))
    if Exp is None:
      Exp = oneD
    else:
      Exp = Exp.concatenate(oneD,batch_key=None,index_unique=None)#,join='outer'
    print("After merge: %d cells %d genes\n\n"%(Exp.shape[0],Exp.shape[1]))
  runMergeSeurat(h5List,re.sub("h5ad$","rds",newH5ad))
  Exp.uns['scaleF'] = pd.read_csv(re.sub("h5ad$","scaleF",newH5ad),header=None,index_col=0).index[0]
  rawD=sc.read_h5ad(strH5ad,backed="r")
  Exp.obs=rawD.obs.copy().loc[Exp.obs.index,:]
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    print("\thighly_variable_genes ...")
    sc.pp.highly_variable_genes(Exp,n_top_genes=3000,batch_key=batchKey)
    Exp1 = Exp[:, Exp.var.highly_variable].copy()
    sc.pp.scale(Exp1, max_value=10)
    print("\tPCA ...")
    npc = 50
    clusterKey = "normalized_cluster"
    sc.tl.pca(Exp1, svd_solver='arpack', n_comps = npc)
    sc.pp.neighbors(Exp1, n_neighbors=10, n_pcs=npc)
    sc.tl.leiden(Exp1,key_added=clusterKey)
    print("\tUMAP ...")
    sc.tl.umap(Exp1)
    print("\ttSNE ...")
    sc.tl.rank_genes_groups(Exp1,clusterKey)
    sc.tl.tsne(Exp1, n_pcs=npc)
    print("\tsaving ...")
    Exp1 = Exp1[Exp.obs_names]
    Exp.obs[clusterKey]=Exp1.obs[clusterKey]
    Exp.obsm['X_normalized_umap'] = Exp1.obsm["X_umap"]
    Exp.obsm['X_normalized_tsne'] = Exp1.obsm["X_tsne"]
    Exp.obsm['X_normalized_pca'] = Exp1.obsm["X_pca"]
    Exp.write(newH5ad)
  return
  
def main():
  print("starting expression normalization ...")
  if len(sys.argv)<1:
    msgError("ERROR: raw h5ad file required!")
  config = bU.inputCheck(sys.argv)
  strH5ad = sys.argv[1]
  strConfig = sys.argv[2]
  newH5ad = "%s.h5ad"%os.path.join(config["output"],"SCT",config["prj_name"])
  
  batchNorm(strH5ad,strConfig,newH5ad,config.get("batchCell"))
  if not os.path.isfile(newH5ad):
    msgError("Error in normalization step!")
  print("Express normalization completed")

if __name__ == "__main__":
  main()
