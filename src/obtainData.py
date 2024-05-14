import os,multiprocessing,warnings,sys,logging,glob,re,random, yaml,functools,gc,resource
from datetime import datetime
import scanpy as sc
import anndata as ad
import pandas as pd
from scipy.sparse import csc_matrix
import cmdUtility as cU
import dbl

print=functools.partial(print, flush=True)
warnings.simplefilter(action='ignore', category=FutureWarning)
logging.disable(level=logging.INFO)

strPipePath=os.path.dirname(os.path.realpath(__file__))
UMIcol="h5path"
ANNcol="metapath" #first column is the cell name/id 
GENEcol='genepath' # first column with 'name' in the header is gene name otherwise the first column will be used
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
  IEcount = pd.read_csv(strF,sep=" ",index_col=0)
  if len(list(set(IEcount.index) & set(cID)))<len(cID):
    IEcount.index = list(IEcount.index+"-1")
  IEcount = IEcount.loc[cID,:]
  IErate = IEcount.apply(lambda x:x/sum(x),axis=1)
  IErate.columns = [a.replace("count","rate") for a in IErate.columns]
  IE = IEcount.merge(IErate,left_index=True,right_index=True)
  return(IE)

def getParse_mtx(strDir):
  D = sc.read_mtx(os.path.join(strDir,"DGE.mtx"))
  gInfo = pd.read_csv(os.path.join(strDir,"all_genes.csv"))
  gInfo.index = ad.utils.make_index_unique(pd.Index(gInfo.gene_name.values))
  D.var = gInfo
  cInfo = pd.read_csv(os.path.join(strDir,"cell_metadata.csv"))
  cInfo.index = cInfo['sample'] + "_" + cInfo['bc_wells']
  D.obs=cInfo
  return D
def Check10Xmtx(strDir):
  for one in glob.glob("%s/*mtx"%strDir)+glob.glob("%s/*tsv"%strDir):
      if not os.path.isfile("%s.gz"%one):
        print("\t\tsave %s as gz"%one)
        with open(one,"rb") as Fin:
          with gzip.open("%s.gz"%one,"wb") as Fout:
            shutil.copyfileobj(Fin, Fout)
def getData_one(oneMeta,sID,strOut,dblScore=True):
  print("\t%s"%oneMeta[sID])
  if os.path.isdir(oneMeta[UMIcol]):
    if os.path.isfile(os.path.join(oneMeta[UMIcol],"DGE.mtx")):#Parse biosciences
      adata = getParse_mtx(oneMeta[UMIcol])
    else:
      Check10Xmtx(oneMeta[UMIcol])
      adata = sc.read_10x_mtx(oneMeta[UMIcol])
  elif oneMeta[UMIcol].endswith('.h5'):
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      adata = sc.read_10x_h5(oneMeta[UMIcol])
  elif oneMeta[UMIcol].endswith('.h5ad'):
    adata = sc.read_h5ad(oneMeta[UMIcol])
  elif oneMeta[UMIcol].endswith('.csv'):
    adata = sc.read_csv(oneMeta[UMIcol]).transpose()
  elif oneMeta[UMIcol].endswith('.tsv'):
    adata = sc.read_csv(oneMeta[UMIcol],'\t').transpose()
  elif oneMeta[UMIcol].endswith('.mtx') or oneMeta[UMIcol].endswith('.mtx.gz'):
    if GENEcol in oneMeta.keys() and ANNcol in oneMeta.keys():
      adata = sc.read_mtx(oneMeta[UMIcol])
      gInfo = pd.read_csv(oneMeta[GENEcol])
      gName = [s for s in gInfo.columns if re.search('name',s,re.IGNORECASE)]
      gName = list(gInfo.columns)[0] if len(gName)==0 else gName[0]
      gInfo.index = list(gInfo[gName])
      if gInfo.shape[0] == adata.shape[0]:
        adata = adata.T
      adata.var = gInfo
      cInfo = pd.read_csv(oneMeta[ANNcol])
      adata.obs_names= cInfo.iloc[:,0]
    else:
      msgError("%s and %s columns are required for mtx file provided in % column in sample meta information!"%(GENEcol,ANNcol,UMIcol))
  else:
    msgError("Unsupported UMI format: %s"%oneMeta[UMIcol])
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
  if dblScore:
    adata=dbl.singleDBL(strOut,oneMeta[UMIcol],oneMeta[sID],adata)
  return adata

def getData_batch(strMeta,sID,strOut,dblScore=True):
  print("processing one sample batch UMI ...")
  strH5ad = re.sub("csv$","h5ad",strMeta)
  if os.path.isfile(strH5ad):
    print("\tUsing previous batch reading result: %s\nPlease remove/rename above file if new batch reading is required!"%strH5ad)
    return
  meta = pd.read_csv(strMeta)
  adatals = []
  with multiprocessing.Pool(processes=min(5,meta.shape[0])) as pool:
    adatals = pool.starmap(getData_one,[[row,sID,strOut,dblScore] for idx, row in meta.iterrows()])
  print("\tmerging samples ...")
  if len(adatals)==1:
    adata = adatals[0]
    adata.obs[batchKey]=meta[sID][0]
    adata.obs.index += ("-" + meta[sID][0])
  else:   
    adata = ad.concat(adatals,
      join="outer",
      keys=meta[sID],
      label=batchKey,
      index_unique="-")
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
      strMeta.append(os.path.join(strOut,"tmp_{0:03}.csv".format(nbatch)))
      one.to_csv(strMeta[-1],index=False)
      nbatch += 1
  return strMeta

def mergeBatch(strMetas,strH5ad):
  Dlist = []
  cellN = 0
  for one in strMetas.split(","):
    oneH5ad = re.sub("csv$","h5ad",one)
    print("\tReading ",os.path.basename(oneH5ad))
    if not os.path.isfile(oneH5ad):
      msgError("Reading batch failed: missing %s"%oneH5ad)
    Dlist.append(ad.read_h5ad(oneH5ad))
    # the following does NOT save more memory
    #D1=ad.read_h5ad(oneH5ad)
    #if D is None:
    #  D = D1
    #else:
    #  #D = ad.concat([D,D1],join="outer")
    #  #D = D.concatenate(D1,join="outer",batch_key=None,index_unique=None)
    #del D1
    #gc.collect()
    cellN += Dlist[-1].shape[0]
    print("\t+++ %d cells +++ peak memory %.2fG"%(cellN,resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2))
  print("\tMerging ...")
  D = ad.concat(Dlist,join="outer")
  del Dlist
  gc.collect()
  print("\t+++ %d cells +++ peak memory %.2fG"%(D.shape[0],resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2))
  print("\tFiltering genes by minimal cell 1")
  sc.pp.filter_genes(D,min_cells=1)
  print("\tsaving",os.path.basename(strH5ad))
  D.write(strH5ad)
  print("\tFinal peak memory %.2fG"%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2))

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
  dbl = True if config.get('dblscore') is None else config.get('dblscore')
  metaBatch = sorted(splitMeta(meta,strOut))
  sID =config["sample_name"]
  cmd=[]
  for one in metaBatch:
    oneH5ad = re.sub("csv$","h5ad",one)
    if os.path.isfile(oneH5ad):
      print("\tUsing previous batch reading result: %s\nPlease remove/rename above file if new batch reading is required!"%oneH5ad)
      continue
    cmd.append("python -u %s/obtainData.py %s %s %s %s"%(strPipePath,one,sID,config["output"],dbl))
  if len(cmd)>0:
    cU.submit_cmd({"ReadB{0:03}".format(i):cmd[i] for i in range(len(cmd))},config)#"ReadB{0:03}".format(i)
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
  getData_batch(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]=='True')

if __name__ == "__main__":
  main()
