import yaml, io, os, sys, re, logging, warnings, glob, math, functools, time, resource
from datetime import datetime
import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from scipy.sparse import csc_matrix
from natsort import natsorted
from timeit import default_timer as timer
import cmdUtility as cU
import obtainData as oD
import dbl
import timeoutFun as tF

print=functools.partial(print, flush=True)
warnings.simplefilter(action='ignore', category=FutureWarning)
logging.disable(level=logging.INFO)
batchKey="library_id"
UMIcol="h5path"
seqMetrics="sequence_metrics"
Rmarkdown="Rmarkdown"
tempDir="raw"
qcDir="QC"
beRaster=True
strPipePath=os.path.dirname(os.path.realpath(__file__))

def getSystemConfig():
  with open(os.path.join(strPipePath,"sys.yml"),"r") as f:
    config = yaml.safe_load(f)
  return config
def msgError(msg):
  print("***** "+msg)
  exit()
def updateRmarkdown(config):
  global Rmarkdown
  Rmarkdown = os.path.join(config["output"],Rmarkdown)
  os.makedirs(Rmarkdown,exist_ok=True)
def runQCmsg(config):
  print("Please check the following QC files @%s:\n\tsequencingQC.csv\n\tsequencingQC.pdf\n\tprefilter.QC.pdf\n\tpostfilter.QC.pdf"%os.path.join(config["output"],qcDir))
  print("And then:")
  print("\t1. Remove any outlier sample from the sample meta table %s"%config['sample_meta'])
  print("\t2. Update config file on cutoff values for cell filtering")
  print("\t3. After making sure the cell filtering setting (might several iteration), set 'runAnalysis: True' in the config file.")
  #print("\t4. (Optional) consider to enable parallel by setting: 'parallel: sge' for CAMHPC or 'parallel: slurm' for EdgeHPC.")
def formatFileName(strF):
  return re.sub('_$','',re.sub('[^A-Za-z0-9]+','_',strF))
def plotSeqQC(meta,sID,strOut,grp=None,redo=None):
  print("plotting sequence QC ...")
  strOut = os.path.join(strOut,qcDir)
  os.makedirs(strOut,exist_ok=True)
  seqQC = []
  newFormat = False
  for i in range(meta.shape[0]):
    strF = meta[seqMetrics][i]#glob.glob(os.path.join(os.path.dirname(meta[UMIcol][i]),"%s*metrics_summary.csv"%meta[sID][i]))
    #strF = os.path.join(os.path.dirname(meta[UMIcol][i]),"%s.metrics_summary.csv"%meta[sID][i])
    #if os.path.isfile(strF):
    if os.path.isfile(strF):
      print("\tQC: %s"%strF)
      one = pd.read_csv(strF,thousands=",")
      if one.shape[0]>3:
      	newFormat=True
      	mergeKeys = list(one.columns)
      	a=mergeKeys.pop()
      	one.columns = [re.sub('Metric Value',meta[sID][i],_) for _ in one.columns]
      else:
      	one.index=[meta[sID][i]]
      seqQC.append(one)
    else:
      print("\tMissing QC: ",meta[sID][i])
      #return
  if len(seqQC)<1:
    print("***NO sequence QC***")
    return
  if newFormat:
  	QC = functools.reduce(lambda l,r: pd.merge(l,r,on=mergeKeys,how='left'),seqQC)
  	QC.to_csv("%s/sequencingQC.csv"%strOut)
  	return
  else:
  	QC = pd.concat(seqQC)
  k=list(QC.columns)
  for i,one in enumerate(k):
    if pd.api.types.is_string_dtype(QC[one]) and QC[one][0].endswith("%"):
      QC[one] = QC[one].str.rstrip('%').astype('float')
      k[i]=one+"%"
  QC.columns=k
  QC.to_csv("%s/sequencingQC.csv"%strOut)
  QC['sample'] = QC.index
  h = 4.8
  w = max(6.4,QC.shape[0]*2/10)

  with PdfPages("%s/sequencingQC.pdf"%strOut) as pdf:
    for one in k:
      ax = QC.plot.bar(x='sample',y=one,rot=90,legend=False,figsize=(w,h))
      ax.set_title(one)
      plt.grid()
      pdf.savefig(bbox_inches="tight")
      plt.savefig(os.path.join(Rmarkdown,"sequencingQC_%s.png"%formatFileName(one)),bbox_inches="tight")
      plt.close()
    if not grp==None:
      for oneG in grp:
        if oneG in meta.columns:
          QC[oneG] = list(meta[oneG])
          w = max(4,meta[oneG].nunique()*2/10)
          for one in k:
            ax = QC[[one,oneG]].boxplot(by=oneG,rot=90,figsize=(w,h))
            ax.set_title(one)
            plt.grid()
            pdf.savefig(bbox_inches="tight")
            plt.savefig(os.path.join(Rmarkdown,"sequencingQC_%s_%s.png"%(oneG,formatFileName(one))),bbox_inches="tight")
            plt.close()
def preprocess(adata,config):
  print("preprocessing ... ")
  varKey=[]
  rmGene=np.full(adata.shape[1],False)
  if "MTstring" in config.keys():#older version
    MTstring=config["MTstring"]
    #get MT genes
    if MTstring is None or len(MTstring)<2:
      for one in ["MT-","Mt-","mt-"]:
        mito_genes = adata.var_names.str.startswith(one)
        if (mito_genes).sum()>0:
          break
    else:
      mito_genes = adata.var_names.str.startswith(MTstring)
    adata.obs['pct_counts_mt'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1 * 100
    if config["rmMT"]:
      print("\tall mitochondrial genes removed")
      adata._inplace_subset_var(np.invert(mito_genes))
  elif "gene_group" in config.keys():
    for k in config['gene_group']:
      if type(config['gene_group'][k]['startwith']) is not list:
        msgError("config error with %s, 'startwith' has to be a list"%k)
      gList = np.full(adata.shape[1],False)
      for one in config['gene_group'][k]['startwith']:
        if len(one)>1:
          gList |= adata.var_names.str.startswith(one)
      if gList.sum()==0:
        print("\tNo genes found for %s"%k)
        continue
      print("\tGene group %s contains %d genes"%(k,gList.sum()))
      adata.var[k]=gList
      varKey += [k]
      if config['gene_group'][k]['rm']:
        print("\t%s genes will be removed"%k)
        rmGene |= gList
  else:
    msgError("Unknown config format! Either 'MTstring' or 'gene_group' is required")
  if len(varKey)>0:
    sc.pp.calculate_qc_metrics(adata,qc_vars=varKey,inplace=True)
    for k in varKey:
      adata.obs.drop("total_counts_%s"%k,axis=1,inplace=True)
      adata.obs.drop("log1p_total_counts_%s"%k,axis=1,inplace=True)
  else:
    sc.pp.calculate_qc_metrics(adata, inplace=True)
  if rmGene.sum()>0:
    adata._inplace_subset_var(np.invert(rmGene))
    print("\tTotal of %d genes are removed",rmGene.sum())
  print("\tcompleted!")
  #return adata
def filtering(adata,config,filterRes):
  print("filtering ...")
  min_cells=config["min.cells"]
  min_features=config["min.features"]
  highCount_cutoff=config["highCount.cutoff"]
  highGene_cutoff=config["highGene.cutoff"]

  print("\tfiltering cells and genes")
  selObs = np.full(adata.shape[0],True)
  selVar = np.full(adata.shape[1],True)
  if "mt.cutoff" in config.keys():
    mt_cutoff=config["mt.cutoff"]
    sTime = timer()
    selObs = np.logical_and(selObs,adata.obs.pct_counts_mt<mt_cutoff)
    filterRes.append("MT cutoff,%f,%d,%d\n"%(mt_cutoff,adata.shape[0],adata.shape[1]))
    print("\t\tfiltered cells with mt.cutoff %d left %d cells"%(mt_cutoff,np.sum(selObs)))
  elif "gene_group" in config.keys():
    for k in config['gene_group']:
      if not "pct_counts_%s"%k in adata.obs.columns:
        print("\t\tskip %s"%k)
        continue
      selObs = np.logical_and(selObs,adata.obs["pct_counts_%s"%k]<config['gene_group'][k]["cutoff"])
      filterRes.append("%s,%f,%d,%d\n"%(k,config['gene_group'][k]["cutoff"],np.sum(selObs),np.sum(selVar)))
      print("\t\tfiltered cells with %s<%d%% left %d cells"%(k,config['gene_group'][k]["cutoff"],np.sum(selObs)))
  else:
    msgError("Unknown config format! Either 'mt.cutoff' or 'gene_group' is required")
    
  ## filtering low content cells and low genes
  selVar=np.logical_and(selVar,adata.var['n_cells_by_counts']>=min_cells)
  filterRes.append("min cell,%d,%d,%d\n"%(min_cells,np.sum(selObs),np.sum(selVar)))
  print("\t\tfiltered genes with min.cells %d left %d genes"%(min_cells,np.sum(selVar)))
  selObs = np.logical_and(selObs,adata.obs.n_genes_by_counts>=min_features)
  filterRes.append("min gene,%d,%d,%d\n"%(min_features,np.sum(selObs),np.sum(selVar)))
  print("\t\tfiltered cells with min.features %d left %d cells"%(min_features,np.sum(selObs)))
  selObs = np.logical_and(selObs,adata.obs.n_genes_by_counts<highGene_cutoff)
  filterRes.append("max gene,%d,%d,%d\n"%(highGene_cutoff,np.sum(selObs),np.sum(selVar)))
  print("\t\tfiltered cells with highGene.cutoff %d left %d cells"%(highGene_cutoff,np.sum(selObs)))
  selObs = np.logical_and(selObs,adata.obs.total_counts<highCount_cutoff)
  filterRes.append("max UMI,%d,%d,%d\n"%(highCount_cutoff,np.sum(selObs),np.sum(selVar)))
  print("\t\tfiltered cells with highCount.cutoff %d left %d cells"%(highCount_cutoff,np.sum(selObs)))
  introSel = [item for item in adata.obs.columns if re.search('intron.*rate$', item,re.IGNORECASE)]
  if len(introSel)==1:
    sel = introSel[0]
    if config.get('intron.cutoff.min') is not None and config['intron.cutoff.min']>0:
      selObs = np.logical_and(selObs,adata.obs[sel]>config["intron.cutoff.min"])
      filterRes.append("intron rate min,%f,%d,%d\n"%(config["intron.cutoff.min"],np.sum(selObs),np.sum(selVar)))
      print("\t\tfiltered cells with intron.cutoff.min %f left %d cells"%(config["intron.cutoff.min"],np.sum(selObs)))
    if config.get('intron.cutoff.max') is not None and config['intron.cutoff.max']<1:
      selObs = np.logical_and(selObs,adata.obs[sel]<config["intron.cutoff.max"])
      filterRes.append("intron rate max,%f,%d,%d\n"%(config["intron.cutoff.max"],np.sum(selObs),np.sum(selVar)))
      print("\t\tfiltered cells with intron.cutoff.max %f left %d cells"%(config["intron.cutoff.max"],np.sum(selObs)))

  with open(os.path.join(Rmarkdown,"filter.csv"),"w") as f:
    f.writelines(filterRes)

  if np.sum(selObs)<10:
    msgError("Few cells (%d<10) left after filtering, please check the filtering setting in config to contitue!"%np.sum(selObs))
  print("\tSubsetting")
  sTime=timer()
  adata._inplace_subset_var(selVar)
  adata._inplace_subset_obs(selObs)
  print("\tCompleted in %.f seconds"%(timer()-sTime))
def plotQC(adata,strPDF,grp=None):
  print("plotting UMI QC ...")
  strRmark = os.path.join(Rmarkdown,os.path.splitext(os.path.basename(strPDF))[0])
  with PdfPages(strPDF) as pdf:
    savePDF_PNG(sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',color=batchKey,alpha=0.5,show=False),
      pdf,"%s_couts_genes.png"%strRmark)
    savePDF_PNG(sc.pl.highest_expr_genes(adata, n_top=20,show=False),
      pdf,"%s_topGene.png"%strRmark)
    w = max(6.4,adata.obs[batchKey].nunique()*2/10)
    plt.rcParams["figure.figsize"] = (w,4.8)
    
    pltUMI = sc.pl.violin(adata, keys = 'total_counts',groupby=batchKey,rotation=90,show=False)
    adataM = adata.obs['total_counts'].median()
    pltUMI.hlines(adataM,xmin=-1,xmax=adata.obs[batchKey].nunique(),color='b')#"Median: %d"%round(adataM)
    pltUMI.set_title("Median: %d"%round(adataM))
    savePDF_PNG(pltUMI,
      pdf,"%s_counts.png"%strRmark)
    savePDF_PNG(sc.pl.violin(adata, keys = 'n_genes_by_counts',groupby=batchKey,rotation=90,show=False),
      pdf,"%s_genes.png"%strRmark)

    for k in [one for one in adata.obs.columns if one.startswith('pct_counts_')]:
      savePDF_PNG(sc.pl.violin(adata, keys =k,groupby=batchKey,rotation=90,show=False),
        pdf,"%s_%s.png"%(strRmark,k))

    plt.rcParams['figure.figsize'] = plt.rcParamsDefault['figure.figsize']

    if not grp==None:
      for oneG in grp:
        if oneG in adata.obs.columns:
          savePDF_PNG(sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',color=oneG,alpha=0.5,show=False),
                      pdf,"%s_%s_couts_genes.png"%(strRmark,oneG))

          w = max(6.4,adata.obs[oneG].nunique()*2/10)
          plt.rcParams["figure.figsize"] = (w,4.8)

          savePDF_PNG(sc.pl.violin(adata, keys = 'n_genes_by_counts',groupby=oneG,rotation=90,show=False,order=list(adata.obs[oneG].unique())),
                      pdf,"%s_%s_genes.png"%(strRmark,oneG))
          savePDF_PNG(sc.pl.violin(adata, keys = 'total_counts',groupby=oneG,rotation=90,show=False,order=list(adata.obs[oneG].unique())),
                      pdf,"%s_%s_counts.png"%(strRmark,oneG))

          for k in [one for one in adata.obs.columns if one.startswith('pct_counts_')]:
            savePDF_PNG(sc.pl.violin(adata, keys =k,groupby=oneG,rotation=90,show=False,order=list(adata.obs[oneG].unique())),
                        pdf,"%s_%s_%s.png"%(strRmark,oneG,k))
          plt.rcParams['figure.figsize'] = plt.rcParamsDefault['figure.figsize']
def savePDF_PNG(ax,pdf,strPNG):
  if beRaster:
    for one in ax.get_children():
      if 'PathCollection' in str(one):
        one.set_rasterized(True)
  pdf.savefig(bbox_inches="tight")
  plt.savefig(strPNG,bbox_inches="tight")
  plt.close()
def checkCells(adata):
  sysConfig = getSystemConfig()
  cellN = adata.obs[batchKey].value_counts()
  cellN = cellN[cellN<sysConfig["minCell"]]
  if cellN.shape[0]==0:
    return
  msgError("The following samples contains less than %d cells, please either relax the filtering or remove them:\n\t%s"%(sysConfig["minCell"],", ".join(list(cellN.index))))
def obtainRAWobsm(D,cluster_reso,cluster_method,reg=None):
  # 95 percentile to normalize
  print("Finding the obs and obsm for the raw")
  #print("\tReading post-filtered h5ad ...")
  #D = ad.read_h5ad(strH5ad)
  print("\thighly_variable_genes seurat_v3 ...")
  # span to removing low expressed genes?: https://github.com/scverse/scanpy/issues/1504
  span=0.3
  print("\t\thttps://github.com/scverse/scanpy/issues/1504")
  while span<=1:
    print("\t\ttry with span: %.2f"%span)
    try:
      sc.pp.highly_variable_genes(D,n_top_genes=5000,flavor="seurat_v3",batch_key=batchKey,span=span)
      break
    except:
      span += 0.1
  D = D[:, D.var.highly_variable].copy()
  print("\tnormalization ...")
  sc.pp.normalize_total(D,target_sum=math.ceil(np.percentile(D.obs.n_counts,80)))
  sc.pp.log1p(D)
  if not reg is None:
    sc.pp.regress_out(D, reg)#['total_counts', 'pct_counts_mt']
  print("\tscale ...")
  sc.pp.scale(D, max_value=10)
  print("\tpca ...")
  sc.tl.pca(D, svd_solver='arpack',n_comps = 100)
  npcs = 50
  print("\tfinding neighbors ...")
  sc.pp.neighbors(D, n_neighbors=10, n_pcs=npcs)
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    clusterKey="raw_cluster"
    print("\tclustering %s (%.2f) ..."%(cluster_method,cluster_reso))
    if bool(re.search('Leiden',cluster_method)):
      sc.tl.leiden(D,resolution=cluster_reso,key_added=clusterKey)
    else:
      sc.tl.louvain(D,resolution=cluster_reso,key_added=clusterKey)
    ## umap embedding
    print("\tumap ...")
    sTime = time.time()
    sc.tl.umap(D, init_pos='spectral')
    sc.tl.rank_genes_groups(D,clusterKey)
    #print("\ttSNE ...")
    #try:
    #  with tF.time_limit(max(3600*10,5*int(time.time()-sTime))):
    #    sc.tl.tsne(D, n_pcs=npcs)
    #except tF.TimeoutException as e:
    #  print("\t\tTime out! NO tSNE!")
    for one in D.obsm_keys():
      D.obsm[re.sub("^X_","X_raw_",one)] = D.obsm.pop(one)
  return D.obs,D.obsm #, D.var.highly_variable
def main():
  task = sys.argv[1]
  if task == 'QC':
    preH5ad = sys.argv[2]
    postH5ad = sys.argv[3]
    strConfig = sys.argv[4]
    with open(strConfig,"r") as f:
      config=yaml.safe_load(f)
    updateRmarkdown(config)
    global beRaster
    beRaster = beRaster if config.get("rasterizeFig") is None else config.get("rasterizeFig")
    print("Reading pre-filtered h5ad ...")
    sTime = timer()
    adata = sc.read_h5ad(preH5ad)
    print("\t%.2f seconds"%(timer()-sTime))
    print("Prefilter Report: %d cells with %d genes"%(adata.shape[0],adata.shape[1]))
    filterRes=["Filter,cutoff,cell_number,gene_number\n"]
    filterRes.append("Prefilter,,%d,%d\n"%(adata.shape[0],adata.shape[1]))
    preprocess(adata,config)
    strPrefilterQC = os.path.join(config["output"],qcDir,"prefilter.QC.pdf")
    if not os.path.isfile(strPrefilterQC):
      plotQC(adata,strPrefilterQC,config["group"])
    if config["filter_step"]:
      #adata = dbl.filterDBL(adata,config,filterRes)
      dbl.filterDBL(adata,config,filterRes)
      #adata = filtering(adata,config,filterRes)
      filtering(adata,config,filterRes)
      plotQC(adata,os.path.join(config["output"],qcDir,"postfilter.QC.pdf"),config["group"])
    checkCells(adata)
    adata.write(postH5ad)
    print("\t%.2f seconds"%(timer()-sTime))
    print("\tQC final peak memory %.2fG"%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2))
  elif task == 'rawOBSM':
    strH5ad = sys.argv[2]
    strOut = sys.argv[3]
    cluster_reso = float(sys.argv[4])
    cluster_method = sys.argv[5]
    if os.path.isfile(strOut):
      print("Using previous raw obsm & obs: %s\n\tPlease rename/remove above file to rerun!"%strOut)
      return
    print("Reading post-filtered h5ad ...")
    adata = sc.read_h5ad(strH5ad)
    adata.obs,adata.obsm=obtainRAWobsm(adata,cluster_reso,cluster_method)
    adata.write(strOut)
    print("\trawOBSM final peak memory %.2fG"%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2))
  else:
    msgError("Unknown task for processQC: ",task)

def QC(preH5ad,postH5ad,config,strConfig):
  reRunQC = True if config.get('reRunQC') is None else config.get('reRunQC')
  if os.path.isfile(postH5ad):
    if config["newProcess"] or reRunQC:
      os.remove(postH5ad)
    else:
      print("\tPrevious QC result is used: %s\n\tPlease remove/rename the above file if a new QC step is desired!"%postH5ad)
      return
  print("process QC/Filtering ...")
  cU.submit_cmd({"QC":"python -u %s/processQC.py QC %s %s %s"%(strPipePath,preH5ad,postH5ad,strConfig)},config)
  if not os.path.isfile(postH5ad):
    msgError("Failed in QC/filtering step!")
  runQCmsg(config)
def runQC(config,meta,strConfig):
  updateRmarkdown(config)
  prefix = os.path.join(config["output"],tempDir,config["prj_name"])
  plotSeqQC(meta,config["sample_name"],config["output"],config["group"])
  os.makedirs(os.path.dirname(prefix),exist_ok=True)
  strPrefilter = "%s_raw_prefilter.h5ad"%prefix
  oD.getData(meta,config,strPrefilter)
  strPostfilter = "%s_raw_postfilter.h5ad"%prefix
  QC(strPrefilter,strPostfilter,config,strConfig)
  #create cmd for rawOBSM
  return rawOBSM(config,strPostfilter,"%s.h5ad"%prefix)

def rawOBSM(config,postH5ad,strH5ad):
  if os.path.isfile(strH5ad):
    if config["newProcess"]:
      h5adTmp = re.sub("h5ad$","%s.h5ad"%datetime.today().strftime('%Y%m%d'),strH5ad)
      print("--> New process is set to 'True' <--\n\tRename %s to %s"%(strH5ad,h5adTmp))
      os.rename(strH5ad,h5adTmp)
    else:
      print("\tPrevious QC result is used: %s\n\tPlease remove/rename the above file if a new QC step is desired!"%strH5ad)
      return ""
  cluster_reso = 0.8 if config.get('clusterResolution') is None else config.get('clusterResolution')
  cluster_method = 'Louvain' if config.get('clusterMethod') is None else config.get('clusterMethod')
  return "python -u %s/processQC.py rawOBSM %s %s %.2f %s"%(strPipePath,postH5ad,strH5ad,cluster_reso,cluster_method)
  #cU.submit_cmd({"rawOBSM":"python -u %s/processQC.py rawOBSM %s %s %.2f %s"%(strPipePath,postH5ad,strH5ad,cluster_reso,cluster_method)},config)
  #if not os.path.isfile(strH5ad):
  #  msgError("Failed in obtaining the raw obsm step!")

if __name__ == "__main__":
  main()
