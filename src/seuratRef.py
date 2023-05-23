#!/usr/bin/env python

import subprocess, os, h5py, sys, warnings, logging, yaml, re, datetime, glob, functools
import anndata as ad
from scipy import sparse
from scipy.sparse import csc_matrix
import pandas as pd
import pyarrow.feather as feather
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
readRDS = robjects.r['readRDS']

print=functools.partial(print, flush=True)
logging.disable()
strPipePath=os.path.dirname(os.path.realpath(__file__))

def msgError(msg):
  print(msg)
  exit()

def inputCheck(args):
  strH5ad = args[1]
  if not os.path.isfile(strH5ad):
    msgError("ERROR: %s does not exist!"%strH5ad)
  strConfig = sys.argv[2]
  if not os.path.isfile(strConfig):
    msgError("ERROR: %s does not exist!"%strConfig)
  with open(strConfig,"r") as f:
    config = yaml.safe_load(f)
  if config['ref_name'] is None:
    print("No reference specified, END")
    return False
  return config

def splitBatch(strH5ad,strMeta,batchCell):
  if batchCell is None:
    print("Batch process (batchCell) is not set in sys.yml, large amount of memory might be required")
    h5adList=[strH5ad]
  else:
    print("Batch process")
    strOut=os.path.dirname(strMeta)+"/tmp/"
    os.makedirs(strOut,exist_ok=True)
    h5adList= glob.glob(strOut+"*.h5ad")
    if len(h5adList)==0:
      print("Reading ...")
      D=ad.read_h5ad(strH5ad)
      sampleCellN = D.obs.library_id.value_counts()
      #print(sampleCellN)
      sID=[]
      cellN=0
      batchN=0
      for one in list(sampleCellN.index):
        sID.append(one)
        cellN+=sampleCellN[one]
        #print(cellN)
        if cellN>batchCell or one==list(sampleCellN.index)[-1]:
          print("batch %d: %d samples"%(batchN,len(sID)))
          strH5ad=strOut+"tmp_%d.h5ad"%batchN
          batchN +=1
          with open(re.sub("h5ad$","txt",strH5ad),"w") as f:
            f.write("\n".join(sID))
          with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            D1=D[D.obs.library_id.isin(sID),:]
            D1.write(strH5ad)
          sID=[]
          cellN=0
          h5adList.append(strH5ad)
      del D
  return(h5adList)

def runOneBatch(oneH5ad,strConfig,oneMeta):
  cmd = "Rscript %s %s %s %s |& tee %s"%(os.path.join(strPipePath,"seuratRef.R"),
                            oneH5ad,strConfig,oneMeta,re.sub("rds$","log",oneMeta))
  subprocess.run(cmd,shell=True,check=True)
  
def batchRef(strH5ad,strConfig,strMeta,batchCell):
  if os.path.isfile(strMeta):
    print("Using previous SeuratRef results: %s\n***=== Important: If a new run is desired, please remove/rename the above file "%strMeta)
    meta = feather.read_feather(strMeta)
    return(meta)
  h5adList = splitBatch(strH5ad,strMeta,batchCell)
  if len(h5adList)==0:
    msgError("No h5ad!")
  print("There are total of %d batches"%len(h5adList))
  mapInfo=[]
  for oneH5ad in h5adList:
    print("***** batch: %s *****"%os.path.basename(oneH5ad))
    oneMeta = re.sub("h5ad$","rds",oneH5ad)
    if not os.path.isfile(oneMeta):
      runOneBatch(oneH5ad,strConfig,oneMeta)
    if not os.path.isfile(oneMeta):
      msgError("\tERROR: %s SeuratRef failed!"%os.path.basename(oneH5ad))
    oneMap = pandas2ri.rpy2py_dataframe(readRDS(oneMeta))
    mapInfo.append(oneMap)
    print("\n\n")
  meta = pd.concat(mapInfo)
  meta.index = list(meta.index)
  feather.write_feather(meta,strMeta)
  print("mapping completed")
  return(meta)

def main():
  print("starting seurat reference mapping ...")
  if len(sys.argv)<2:
    msgError("ERROR: raw h5ad file and the config file are required!")
  config = inputCheck(sys.argv)
  if config == False :
    return()
  strH5ad = sys.argv[1]
  strConfig=sys.argv[2]

  D = ad.read_h5ad(strH5ad,backed="r") #,backed=True
  Dbatch = D.obs["library_id"].copy()
  strMeta = "%s.feather"%os.path.join(config["output"],"SeuratRef",config["prj_name"])#strH5ad.replace("raw.h5ad","seuratRef.csv")
  meta = batchRef(strH5ad,strConfig,strMeta,config.get('batchCell'))
  
  FakeD = pd.DataFrame({"FakeG%d"%i:[0 for j in range(meta.shape[0])] for i in range(2)},
                      index=meta.index)
  D = ad.AnnData(FakeD)
  D.obs = pd.concat([meta[[one for one in meta.columns if one.startswith("predicted")]],Dbatch],axis=1,join="inner")
  for one in set([one.rsplit("_",1)[0] for one in [one for one in meta.columns if not one.startswith("predicted")]]):
    D.obsm['X_%s'%one] = meta[[a for a in meta.columns if a.startswith(one)]].values
  print("\tsaving ...")
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    D.write(re.sub("feather","h5ad",strMeta))
  print("Mapping to reference completed!")

if __name__ == "__main__":
  main()
