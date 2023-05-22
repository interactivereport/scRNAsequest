#!/usr/bin/env python

import subprocess, os, h5py, sys, warnings, logging, yaml, re, datetime, glob, functools
import anndata as ad
from scipy import sparse
from scipy.sparse import csc_matrix
import pandas as pd

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

def splitBatch(strH5ad,strCSV):
  with open(strPipePath+"/sys.yml","r") as f:
    sysCon = yaml.safe_load(f)
  batchCell= sysCon.get("batchCell")
  if batchCell is None:
    print("Batch process (batchCell) is not set in sys.yml, large amount of memory might be required")
    h5adList=[strH5ad]
  else:
    print("Batch process")
    strOut=os.path.dirname(strCSV)+"/tmp/"
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

def runOneBatch(oneH5ad,strConfig,oneCSV):
  cmd = "Rscript %s %s %s %s |& tee %s"%(os.path.join(strPipePath,"seuratRef.R"),
                            oneH5ad,strConfig,oneCSV,re.sub("csv$","log",oneCSV))
  subprocess.run(cmd,shell=True,check=True)
  
def batchRef(strH5ad,strConfig,strCSV):
  if os.path.isfile(strCSV):
    print("Using previous SeuratRef results: %s\n***=== Important: If a new run is desired, please remove/rename the above file "%strCSV)
    meta = pd.read_csv(strCSV,index_col=0,header=0)
    meta.index = list(meta.index)
    return(meta)
  h5adList = splitBatch(strH5ad,strCSV)
  if len(h5adList)==0:
    msgError("No h5ad!")
  print("There are total of %d batches"%len(h5adList))
  mapInfo=[]
  for oneH5ad in h5adList:
    print("***** batch: %s *****"%os.path.basename(oneH5ad))
    oneCSV = re.sub("h5ad$","csv",oneH5ad)
    if not os.path.isfile(oneCSV):
      runOneBatch(oneH5ad,strConfig,oneCSV)
    if not os.path.isfile(oneCSV):
      msgError("\tERROR: %s SeuratRef failed!"%os.path.basename(oneH5ad))
    oneMap = pd.read_csv(oneCSV,index_col=0,header=0)
    mapInfo.append(oneMap)
    print("\n\n")
  meta = pd.concat(mapInfo)
  meta.index = list(meta.index)
  meta.to_csv(strCSV)
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
  strCSV = "%s.csv"%os.path.join(config["output"],"SeuratRef",config["prj_name"])#strH5ad.replace("raw.h5ad","seuratRef.csv")
  meta = batchRef(strH5ad,strConfig,strCSV)
  
  FakeD = pd.DataFrame({"FakeG%d"%i:[0 for j in range(meta.shape[0])] for i in range(2)},
                      index=meta.index)
  D = ad.AnnData(FakeD)
  D.obs = pd.concat([meta[[one for one in meta.columns if one.startswith("predicted")]],Dbatch],axis=1,join="inner")
  for one in set([one.rsplit("_",1)[0] for one in [one for one in meta.columns if not one.startswith("predicted")]]):
    D.obsm['X_%s'%one] = meta[[a for a in meta.columns if a.startswith(one)]].values
  print("\tsaving ...")
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    D.write(re.sub("csv$","h5ad",strCSV))
  print("Mapping to reference completed!")

if __name__ == "__main__":
  main()
