import yaml,io,os,sys,subprocess,re,logging,warnings,glob,functools,gc
import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
import cmdUtility as cU
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import scipy.stats as ss


print=functools.partial(print, flush=True)
warnings.simplefilter(action='ignore', category=FutureWarning)
logging.disable(level=logging.INFO)
tempDir="raw"
strPipePath=os.path.dirname(os.path.realpath(__file__))
rawSuffix = "_raw_obsAdd.h5ad"

def msgError(msg):
  print("***** "+msg)
  exit()
def getRefName(config):
  # can be dict, name of a reference, or a path to a rds file
  refName = config.get("ref_name") # can be dict, name of a reference, or a rds
  if refName is None:
    refName='None'
  elif isinstance(refName,dict):
    refName=",".join(refName.keys())
  else:
    refName=','
  return refName
def combine(mIDs,prefix,config):
  CKmethods = [one for one in mIDs if os.path.isfile("%s.h5ad"%os.path.join(config['output'],one,config["prj_name"]))]+['raw']
  if not "SCT" in CKmethods:
    msgError("SCT is missing! and it is required expression for visualization!")
  
  clusterR = 2 if config.get("major_cluster_rate") is None else config.get("major_cluster_rate")
  refName = getRefName(config)
  print("combining all methods results ...")
  strH5ad = prefix+".h5ad"
  strRaw = prefix+rawSuffix
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
def obsCOR(obs,strH5ad):
  if not strH5ad.endswith(rawSuffix):
    return
  print("Correlation among obs for scDEG consideration")
  obs = obs.loc[:,obs.apply(lambda x: x.nunique()>1)]
  dt = {}
  for cName in obs.columns:
    if not pd.api.types.is_categorical_dtype(obs[cName]) and obs[cName].nunique()<100:
      dt[cName]='category'
  obs = obs.astype(dt).copy()
  pVal = pd.DataFrame(1,index=obs.columns,columns=obs.columns)
  strPDF = os.path.join(os.path.dirname(strH5ad),"obsCor.pdf")
  with PdfPages(strPDF) as pdf:
    for a in obs.columns:
      print("\t%d/%d: %s"%(list(obs.columns).index(a),obs.shape[1],a))
      if obs[a].nunique()<2:
        continue
      obsT = obs.copy()
      for b in obs.columns:
        bbox_inches="tight"
        if obs[b].nunique()<2:
          continue
        if a==b:
          pVal.loc[a,b]=0
          continue
        if pd.api.types.is_numeric_dtype(obsT[a]) and pd.api.types.is_numeric_dtype(obsT[b]):
          s,pVal.loc[a,b]=ss.pearsonr(obsT[a].values,obsT[b].values)
          ax = obsT.plot.scatter(x=a,y=b,s=1,zorder=0)
        elif pd.api.types.is_numeric_dtype(obsT[a]) and pd.api.types.is_categorical_dtype(obsT[b]):
          sList = [v[1].values for v in obsT[[a,b]].groupby([b])[a]]
          s,pVal.loc[a,b] = ss.f_oneway(*sList)
          ax = obsT.boxplot(column=a,by=b,flierprops={'marker': '.','markersize':2},zorder=0,rot=90)
          ax.set_xlabel(b)
          ax.set_ylabel(a)
          ax.grid(linestyle='dotted')
          if obsT[a].max()>1000:
            ax.set_yscale('log')
        elif pd.api.types.is_categorical_dtype(obsT[a]) and pd.api.types.is_numeric_dtype(obsT[b]):
          sList = [v[1].values for v in obsT[[a,b]].groupby([a])[b]]
          s,pVal.loc[a,b] = ss.f_oneway(*sList)
          ax = obsT.boxplot(column=b,by=a,vert=False,flierprops={'marker': '.','markersize':2},zorder=0)  
          ax.set_xlabel(b)
          ax.set_ylabel(a)
          ax.grid(linestyle='dotted')
          if obsT[b].max()>1000:
            ax.set_xscale('log')
        elif pd.api.types.is_categorical_dtype(obsT[a]) and pd.api.types.is_categorical_dtype(obsT[b]):
          A = obsT[[a,b]]#.drop_duplicates()
          ct = pd.crosstab(A[a],A[b])
          if ct.nunique().nunique()==1:
            s,pVal.loc[a,b]=(0,1)
          else:
            s,pVal.loc[a,b],d,expF = ss.chi2_contingency(ct)
          plt.figure(figsize=(max(6,ct.shape[1]/2),max(3,ct.shape[0]/2)))#.set_figheight(max(2,ct.shape[0]/2))
          the_table = plt.table(ct.values,rowLabels=ct.index,colLabels=ct.columns,rowLoc='right',cellLoc='center',loc='center',zorder=2)#
          the_table.auto_set_font_size(False)
          the_table.set_fontsize(10)
          the_table.scale(1, 1.5)
          for k,cell in the_table.get_celld().items():#_cells.keys():
              if k[0]==0:
                  cell.set_text_props(rotation=45,horizontalalignment='left',verticalalignment='baseline')#
                  cell.set_linewidth(0)
          plt.axis('off')
          plt.figtext(0.5,0.93,b, horizontalalignment='center', size=10, weight='bold')
          plt.figtext(0.04,0.1,a,size=10, weight='bold',rotation='vertical',verticalalignment='bottom')#,horizontalalignment='right', verticalalignment='bottom'#,rotation_mode='anchor'
          plt.tight_layout(rect=(0.08,0,0.98,0.90),pad=0.1)
          ax = plt.gca()
          bbox_inches=None
        else:
          print("unknown type for %s or %s"%(a,b))
          continue
        ax.set_rasterization_zorder(1)
        #ax.set_title(label="stat=%.3f pvalue: %.3f"%(s,pVal.loc[a,b]),loc='left')
        plt.figtext(0.02,0.93,"stat=%.3f pvalue: %.3f"%(s,pVal.loc[a,b]), horizontalalignment='left',verticalalignment='bottom',size=15, weight='bold')
        res=plt.title('')
        res=plt.suptitle('')
        pdf.savefig(bbox_inches=bbox_inches)
        res=plt.close()
      del obsT
      gc.collect()
  pVal.to_csv(re.sub('pdf$','csv',strPDF))

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
    D1 = ad.read_h5ad(oneH5ad,backed="r")
    obs = D1.obs.copy()
    addObs = [one for one in obs.columns if not one in D.obs.columns]
    if one == "SeuratRef":
      seuratRefCluster=obs.columns
      prefix_ref="" if ref_name is None else ref_name.split(',')[0]
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
  obsCOR(D.obs,outH5ad)
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
  #print(sys.argv)
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
