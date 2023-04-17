#!/usr/bin/env python

import subprocess, os, h5py, sys, warnings, re, yaml
import anndata as ad
from scipy import sparse
from scipy.sparse import csc_matrix
import pandas as pd

def msgError(msg):
  print("***** "+msg)
  exit()

def dbl(config,strH5ad,adata,filterRes):
  print("starting doublet finding ...")
  if not os.path.isfile(strH5ad):
    msgError("ERROR: %s does not exist!"%strH5ad)
  # run dbl process
  os.makedirs(os.path.join(config["output"],"dbl"),exist_ok=True)
  strDBL = "%s.csv"%os.path.join(config["output"],"dbl",config["prj_name"])
  cmd = "Rscript %s %s %s |& tee %s/scDblFinder.log"%(os.path.join(os.path.dirname(os.path.realpath(__file__)),"dbl.R"),
                            strH5ad,strDBL,os.path.dirname(strDBL))
  if os.path.isfile(strDBL):
    print("Using previous scDBL results: %s\n***=== Important: If a new run is desired, please remove/rename the above file "%strDBL)
  else:  
    subprocess.run(cmd,shell=True,check=True)#,stdout=subprocess.PIPE
  if not os.path.isfile(strDBL):
    msgError("\tERROR: doublet finding failed!")

  print("=== doublet finding process is completed! ===")
  # merge dbl results
  meta = pd.read_csv(strDBL,index_col=0,header=0)
  adata.obs = adata.obs.merge(meta,how='left',left_index=True,right_index=True)
  # filtering
  if config.get("dbl_filter") is True:
    adata = adata[adata.obs['scDblFinder.class'].isin(['singlet'])]
    print("Filtering by dbl class: ",adata.shape[0]," cells")
    filterRes.append("doublet,byClass,%d,%d\n"%(adata.shape[0],adata.shape[1]))
  if type(config.get("dbl_filter")) is float:
    adata = adata[adata.obs['scDblFinder.score']<=config.get("dbl_filter")]
    print("Filtering by dbl score (",config.get("dbl_filter"),"): ",adata.shape[0]," cells")
    filterRes.append("doublet,byScore(%f),%d,%d\n"%(config.get("dbl_filter"),adata.shape[0],adata.shape[1]))
  # return
  return adata,filterRes

