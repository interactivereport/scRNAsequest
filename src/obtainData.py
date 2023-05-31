import os,multiprocessing,warnings,sys,logging,glob,re,random, yaml
from datetime import datetime
import scanpy as sc
import pandas as pd
from scipy.sparse import csc_matrix
import cmdUtility as cU
import dbl

warnings.simplefilter(action='ignore', category=FutureWarning)
logging.disable(level=logging.INFO)

strPipePath=os.path.dirname(os.path.realpath(__file__))
UMIcol="h5path"
ANNcol="metapath"
IntronExon="intron_exon_count_path"
batchKey="library_id"
#meta = pd.read_csv("/mnt/depts/dept04/compbio/users/zouyang/projects/AD_sc_atlas/data/cellbender_final.csv")
#meta = meta[0:3]
#sID="Sample_Name"
def msgError(msg):
  print("***** "+msg)
  exit()
def getSysConfig():
  with open(os.path.join(strPipePath,"sys.yml"),"r") as f:
    config = yaml.safe_load(f)
  return(config)
def getIntronExon(strF,cID):
  IEcount = pd.read_csv(strF,sep=None,index_col=0)
  if len(list(set(IEcount.index) & set(cID)))<len(cID):
    IEcount.index = list(IEcount.index+"-1")
  IEcount = IEcount.loc[cID,:]
  IErate = IEcount.apply(lambda x:x/sum(x),axis=1)
  IErate.columns = [a.replace("count","rate") for a in IErate.columns]
  IE = IEcount.merge(IErate,left_index=True,right_index=True)
  return(IE)

def getData_one(oneMeta,sID,strOut):
  print("\t%s"%oneMeta[sID])
  if os.path.isdir(oneMeta[UMIcol]):
    for one in glob.glob("%s/*mtx"%oneMeta[UMIcol])+glob.glob("%s/*tsv"%oneMeta[UMIcol]):
      if not os.path.isfile("%s.gz"%one):
        print("\t\tsave %s as gz"%one)
        with open(one,"rb") as Fin:
          with gzip.open("%s.gz"%one,"wb") as Fout:
            shutil.copyfileobj(Fin, Fout)
    adata = sc.read_10x_mtx(oneMeta[UMIcol])
  elif oneMeta[UMIcol].endswith('.h5'):
    adata = sc.read_10x_h5(oneMeta[UMIcol])
  elif oneMeta[UMIcol].endswith('.csv'):
    adata = sc.read_csv(oneMeta[UMIcol]).transpose()
  elif oneMeta[UMIcol].endswith('.tsv'):
    adata = sc.read_csv(oneMeta[UMIcol],'\t').transpose()
  else:
    Exit("Unsupported UMI format: %s"%oneMeta[UMIcol])
  adata.var_names_make_unique()
  sc.pp.filter_cells(adata,min_counts=1)
  adata.X = csc_matrix(adata.X)
  ## add intro/exon counts/ratio if exists
  if IntronExon in oneMeta.keys() and os.path.isfile(oneMeta[IntronExon]):
    IE = getIntronExon(oneMeta[IntronExon],adata.obs_names)
    adata.obs = adata.obs.merge(IE,left_index=True,right_index=True)
  if ANNcol in oneMeta.keys() and not pd.isna(oneMeta[ANNcol]):
    if not os.path.exists(oneMeta[ANNcol]):
      print("\t\tWarning: cell level meta file does not exist!\n\t\t\t%s"%oneMeta[ANNcol])
    else:
      annMeta = pd.read_csv(oneMeta[ANNcol],index_col=0)
      adata = adata[adata.obs.index.isin(list(annMeta.index))]
      adata.obs = pd.merge(adata.obs,annMeta,left_index=True,right_index=True)
      print("\t\tCell level meta available, cell number: %d"%adata.shape[0])
  for one in oneMeta.keys():
    if not 'path' in one and not one==sID:
      adata.obs[one]=oneMeta[one]
  # add dbl detection
  adata=dbl.singleDBL(strOut,oneMeta[UMIcol],oneMeta[sID],adata)
  return adata

def getData_batch(strMeta,sID,strOut):
  print("processing one sample batch UMI ...")
  strH5ad = re.sub("csv$","h5ad",strMeta)
  if os.path.isfile(strH5ad):
    print("\tUsing previous batch reading result: %s\nPlease remove/rename above file if new batch reading is required!"%strH5ad)
    return
  meta = pd.read_csv(strMeta)
  adatals = []
  with multiprocessing.Pool(processes=min(5,meta.shape[0])) as pool:
    adatals = pool.starmap(getData_one,[[row,sID,strOut] for idx, row in meta.iterrows()])
  print("\tmerging samples ...")
  if len(adatals)==1:
    adata = adatals[0]
    adata.obs[batchKey]=meta[sID][0]
  else:   
    adata = sc.AnnData.concatenate(*adatals,
      join="outer",
      batch_categories=meta[sID],
      batch_key=batchKey)
  print("\tmerging samples completed!")
  ## remove duplicated columns in var
  varCol = [one.split("-")[0] for one in adata.var.columns]
  varInx = [i for i,v in enumerate(varCol) if not v in varCol[:i]]
  adata.var = adata.var.iloc[:,varInx]
  adata.var.columns=[varCol[i] for i in varInx]
  print("\tsaving adata")
  adata.write(strH5ad)

def splitMeta(meta,strOut):
  os.makedirs(strOut,exist_ok=True)
  ix = list(meta.index)
  random.shuffle(ix)
  sysConfig=getSysConfig()
  split_size=50 if sysConfig.get("readSampleChunk") is None else sysConfig.get("readSampleChunk")
  strMeta = glob.glob(os.path.join(strOut,"tmp_*.csv"))
  if len(strMeta)==0:
    nbatch = 0
    for one in [meta.iloc[ix[i:i+split_size]] for i in range(0,len(ix),split_size)]:
      strMeta.append(os.path.join(strOut,"tmp_%d.csv"%nbatch))
      one.to_csv(strMeta[-1],index=False)
      nbatch += 1
  return strMeta

def mergeBatch(strMetas,strH5ad):
  D = None
  for one in strMetas.split(","):
    oneH5ad = re.sub("csv$","h5ad",one)
    print("\tMerging ",os.path.basename(oneH5ad))
    if not os.path.isfile(oneH5ad):
      msgError("Reading batch failed: missing %s"%oneH5ad)
    D1=sc.read_h5ad(oneH5ad)
    if D is None:
      D = D1
    else:
      D = D.concatenate(D1,join="outer",batch_key=None,index_unique=None,)
    print("\t+++ %d cells +++"%D.shape[0])
    del D1
  print("\tFiltering genes by minimal cell 1")
  sc.pp.filter_genes(D,min_cells=1)
  print("\tsaving",os.path.basename(strH5ad))
  D.write(strH5ad)

def getData(meta,config,strH5ad):
  print("starting Reading ...")
  strOut = os.path.join(os.path.dirname(strH5ad),"tmp")
  if os.path.isdir(strOut) and config["newProcess"]:
    outTmp = strOut+"_"+datetime.today().strftime('%Y%m%d')
    os.rename(strOut,outTmp)
  if os.path.isfile(strH5ad):
    if config["newProcess"]:
      h5adTmp = re.sub("h5ad$","%s.h5ad"%datetime.today().strftime('%Y%m%d'),strH5ad)
      print("--> New process is set to 'True' <--\n\tRename %s to %s"%(strH5ad,h5adTmp))
      os.rename(strH5ad,h5adTmp)
    else:
      print("\tPrevious read data is used: %s\n\tPlease remove/rename the above file if new reading is desired!"%strH5ad)
      return
    #D = sc.read_h5ad(strH5ad,backed="r")
    #return D
  metaBatch = sorted(splitMeta(meta,strOut))
  sID =config["sample_name"]
  cmd=[]
  for one in metaBatch:
    oneH5ad = re.sub("csv$","h5ad",one)
    if os.path.isfile(oneH5ad):
      print("\tUsing previous batch reading result: %s\nPlease remove/rename above file if new batch reading is required!"%oneH5ad)
      continue
    cmd.append("python -u %s/obtainData.py %s %s %s"%(strPipePath,one,sID,config["output"]))
  if len(cmd)>0:
    cU.submit_cmd({"ReadB%d"%i:cmd[i] for i in range(len(cmd))},config)
  print("merging reading batches")
  cU.submit_cmd({"ReadBm":"python -u %s/obtainData.py %s %s"%(strPipePath,",".join(metaBatch),strH5ad)},config)
  if not os.path.isfile(strH5ad):
    msgError("Reading batch failed at merging batches")
  #D = sc.read_h5ad(strH5ad,backed="r")
  #return D

def main():
  print("starting Reading ...")
  if len(sys.argv)<3:
    msgError("ERROR: sample meta files, project output path and sample name column header are required!")
  if len(sys.argv)==3:
    return mergeBatch(sys.argv[1],sys.argv[2])
  getData_batch(sys.argv[1],sys.argv[2],sys.argv[3])

if __name__ == "__main__":
  main()
