#!/usr/bin/env python

import subprocess, os, h5py, sys, warnings, re
import scanpy as sc
from scipy import sparse
from scipy.sparse import csc_matrix
import pandas as pd

def msgError(msg):
  print(msg)
  exit()
  
def main():
  print("starting SCTransform ...")
  if len(sys.argv)<1:
    msgError("ERROR: raw h5ad file required!")
  strH5ad = sys.argv[1]
  if not os.path.isfile(strH5ad):
    msgError("ERROR: %s does not exist!"%strH5ad)
  if not strH5ad.endswith("raw.h5ad"):
    msgError("ERROR: %s is not raw h5ad file required!"%strH5ad)
  strH5 = strH5ad.replace("raw.h5ad","SCT.h5")
  
  cmd = "Rscript %s %s %s"%(os.path.join(os.path.dirname(os.path.realpath(__file__)),"SCT.R"),
                            strH5ad,strH5)
  subprocess.run(cmd,shell=True,check=True)
  if not os.path.isfile(strH5):
    msgError("".join(cmdR.stdout.decode("utf-8"))+"\nERROR: SCTransform failed!")
  print("\tcreating h5ad from SCT h5 ...")
  f = h5py.File(strH5,'r')
  X = sparse.csc_matrix((f['data'],f['indices'],f['indptr']),f['shape'])
  #cID = [one.decode('utf-8') for one in f['row_names']]
  #gID = [one.decode('utf-8') for one in f['col_names']]
  f.close()
  cID = pd.read_csv(re.sub("h5$","cID",strH5),header=None,index_col=0)
  gID = pd.read_csv(re.sub("h5$","gID",strH5),header=None,index_col=0)

  D = sc.AnnData(X)
  #D.obs_names = pd.Index(cID)
  #D.var_names = pd.Index(gID)
  D.obs_names = list(cID.index)
  D.var_names = list(gID.index)
  
  Draw = sc.read_h5ad(strH5ad)
  D.obs=pd.concat([D.obs,Draw.obs],axis=1,join="inner")
  D.var=pd.concat([D.var,Draw.var],axis=1,join='inner')
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    print("\tPCA ...")
    sc.tl.pca(D, svd_solver='arpack', n_comps = 100)
    npc = 50
    sc.pp.neighbors(D, n_neighbors=10, n_pcs=npc)
    sc.tl.louvain(D,key_added="SCT.louvain")
    print("\tembedding ...")
    sc.tl.umap(D)
    sc.tl.rank_genes_groups(D, 'SCT.louvain')
    sc.tl.tsne(D, n_pcs=npc)
    print("\tsaving ...")
    D.write(strH5ad.replace("raw.h5ad","SCT.h5ad"))
  
  print("Completed SCTransform ...")

if __name__ == "__main__":
  main()
