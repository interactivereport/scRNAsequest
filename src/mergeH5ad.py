import yaml,io,os,sys,subprocess,re,logging,warnings,glob
import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
import cmdUtility as cU

warnings.simplefilter(action='ignore', category=FutureWarning)
logging.disable(level=logging.INFO)
tempDir="raw"
strPipePath=os.path.dirname(os.path.realpath(__file__))

def msgError(msg):
  print("***** "+msg)
  exit()
def combine(mIDs,prefix,config):
  CKmethods = [one for one in mIDs if os.path.isfile("%s.h5ad"%os.path.join(config['output'],one,config["prj_name"]))]+['raw']
  if not "SCT" in CKmethods:
    msgError("SCT is missing! and it is required expression for visualization!")
  
  clusterR = 2 if config.get("major_cluster_rate") is None else config.get("major_cluster_rate")
  refName = str(config.get("ref_name"))
  print("combining all methods results ...")
  strH5ad = prefix+".h5ad"
  strRaw = prefix+"_raw_obsAdd.h5ad"
  cmds={}
  if not os.path.isfile(strH5ad):
    cmds["inNorm"] = "python -u %s/mergeH5ad.py %s %s %s %s %.2f %s"%(strPipePath,strH5ad,
        "%s.h5ad"%os.path.join(config['output'],"SCT",config["prj_name"]),",".join(CKmethods),prefix,clusterR,refName)
  else:
    print("\tUsing previous merged h5ad: %s\n\t\tPlease rename/remove the above file if a new merge is needed!"%strH5ad)
  if not os.path.isfile(strRaw):
    cmds['inRaw']="python -u %s/mergeH5ad.py %s %s %s %s %.2f %s"%(strPipePath,strRaw,
        "%s.h5ad"%os.path.join(config['output'],tempDir,config["prj_name"]),",".join(CKmethods),prefix,clusterR,refName)
  else:
    print("\tUsing previous merged h5ad: %s\n\t\tPlease rename/remove the above file if a new merge is needed!"%strRaw)
  if len(cmds)>0:
    cU.submit_cmd(cmds,config)
  if not os.path.isfile(strH5ad) or not os.path.isfile(strRaw):
    msgError("Error in integration all methods!")
    
  print("kBET & silhouette evaluation ...")
  strKbet = "%s/evaluation/kBet.pdf"%os.path.dirname(prefix)
  strSil = "%s/evaluation/silhouette.pdf"%os.path.dirname(prefix)
  strSeurat = prefix+".h5seurat"
  cmds={}
  if os.path.isfile(strKbet):
    print("\tkBET result exists: %s\n\t\tPlease rename/remove the above file if a new kBET evaluation is needed!"%strKbet)
  else:
    cmds["kBET"]="Rscript %s/kBET.R %s.h5ad %s"%(strPipePath,prefix,strKbet)
  if os.path.isfile(strSil):
    print("\tSilhouette result exists: %s\n\t\tPlease rename/remove the above file if a new silhouette evaluation is needed!"%strSil)
  else:
    cmds["silhouette"]="python -u %s/silhouette.py %s.h5ad %s"%(strPipePath,prefix,strSil)
  if os.path.isfile(strSeurat):
    print("\th5seurat exists: %s\n\t\tPlease rename/remove the above file if a new h5seurat merge is needed!"%strSeurat)
  else:
    cmds["inSeurat"]="python -u %s/mergeH5ad.py %s"%(strPipePath,prefix)
  if len(cmds)>0:
    cU.submit_cmd(cmds,config)
  if not os.path.isfile(strKbet) or not os.path.isfile(strSil) or not os.path.isfile(strSeurat):
    print("WARNING: evaluation encountered issues!")

  return ad.read_h5ad(strH5ad,backed="r").uns.get('scaleF')

def integrateH5ad(strH5ad,methods,prefix,outH5ad,majorR=None,ref_name=None):
  print("reading h5ad for integration")
  D = ad.read_h5ad(strH5ad)
  obsm = pd.DataFrame(index=D.obs.index)
  seuratRefLab = None
  seuratRefCluster= None
  for one in methods:
    if one in ["SCT"]:#,"raw"
      continue
    oneH5ad = "%s.h5ad"%os.path.join(os.path.dirname(prefix),one,os.path.basename(prefix))#"%s_%s.h5ad"%(prefix,one)
    if not os.path.isfile(oneH5ad):
      print("Warning: ignore missing h5ad (%s) for method %s"%(oneH5ad,one))
      continue
    print("\tmerging ",one)
    D1 = ad.read_h5ad(oneH5ad,backed=True)
    obs = D1.obs.copy()
    addObs = [one for one in obs.columns if not one in D.obs.columns]
    if one == "SeuratRef":
      seuratRefCluster=obs.columns
      prefix_ref=""
      if isinstance(ref_name,dict):
        prefix_ref=list(ref_name.keys())[0]
      seuratRefLab = [i for i in addObs if not i.endswith("score") and i.startswith("predicted."+prefix_ref)]
    if len(addObs)>0:
      D.obs=D.obs.merge(obs[addObs],how="left",left_index=True,right_index=True)
    # check the order of cells
    for k in D1.obsm.keys():
      #kname = k.replace("X_","X_%s_"%one)
      obsm1 = obsm.merge(pd.DataFrame(D1.obsm[k],index=D1.obs.index),how="left",left_index=True,right_index=True)
      D.obsm[k] = obsm1.fillna(0).to_numpy()
  ## assign seurat labels to other integration clusters
  if majorR is not None and majorR<=1 and seuratRefLab is not None and ref_name is not None:
    print("----- Reassign cluster/leiden to seurat label transfer")
    for aCluster in [aCluster for aCluster in D.obs.columns if aCluster.endswith('louvain') or aCluster.endswith('cluster')]:
      if aCluster in seuratRefCluster:
        continue
      for aLabel in seuratRefLab:
        D.obs["%s_%s"%(aCluster,aLabel)] = findMajor(D.obs[[aCluster,aLabel]],majorR)
  ## remove NA/nan
  for one in D.obs.columns:
    if D.obs[one].isna().sum()>0:
      print("Fix NA: %s"%one)
      if D.obs[one].dtype == 'category':
        D.obs[one] = list(D.obs[one])
        D.obs[one].fillna("NA")
  D.write(outH5ad)
def findMajor(X,majorR):
  # first column is the cluster number, second column is the label
  colName = list(X.columns)
  Xsize = X.groupby(colName).size().reset_index().rename(columns={0:'count'})
  Xmax = (Xsize[Xsize.groupby(colName[0])['count'].transform(max)==Xsize['count']]
    .merge(Xsize.groupby(colName[0])['count'].sum().reset_index().rename(columns={'count':'sum'}),on=colName[0]))
  Xmax["name"] = Xmax.apply(lambda x:x[colName[1]] if x['count']/x['sum']>majorR else x[colName[0]],axis=1)
  Xsel = Xmax.set_index(colName[0])['name'].to_dict()
  return(X[colName[0]].map(Xsel))
def saveSeuratObj(prefix):
  strH5ad = "%s.h5ad"%prefix
  strRDS = glob.glob('%s*.rds'%os.path.join(os.path.dirname(prefix),"SCT",os.path.basename(prefix)))
  if len(strRDS)>0 and os.path.isfile(strRDS[0]):
    cU.run_cmd("Rscript %s/seuratObj.R %s %s"%(strPipePath,strRDS[0],strH5ad))

def main():
  print(sys.argv)
  if len(sys.argv)==2:
    saveSeuratObj(sys.argv[1])
    return()
  if len(sys.argv)<7:
    msgError("ERROR: 6 arguments are required for mergeH5ad!")
  outH5ad=sys.argv[1]
  inH5ad=sys.argv[2]
  methods=sys.argv[3].split(",")
  prefix=sys.argv[4]
  clusterR=float(sys.argv[5])
  refName=None if sys.argv[6]=="None" else sys.argv[6]
  integrateH5ad(inH5ad,methods,prefix,outH5ad,clusterR,refName)

if __name__ == "__main__":
  main()
