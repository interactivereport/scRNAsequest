import scPipe, os, time, sys, re
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
sns.set_style('whitegrid')

#strPipePath=""
UMIcol="h5path"
ANNcol="metapath"
CB_expCellNcol="expected_cells"
CB_dropletNcol="droplets_included"
CB_count="low_count_threshold"
CB_learningR="learning_rate"
sampleNameCol="Sample_Name"

def EXIT(msg):
  print(msg)
  exit()

def cellbender(strMeta):
  meta = pd.read_csv(strMeta)
  if not UMIcol in meta.columns or not CB_expCellNcol in meta.columns or not CB_dropletNcol in meta.columns:
    EXIT("Please make sure the following columns exist in %s\n:%s, %s, %s"%(strMeta,UMIcol,CB_expCellNcol,CB_dropletNcol))
  if not sampleNameCol in meta.columns:
    meta[sampleNameCol]=[re.sub(".h5$","",re.sub(".raw_feature_bc_matrix.h5$","",os.path.basename(one))) for one in meta[UMIcol]]
  strOut = os.path.join(os.path.dirname(strMeta),"cellbender")
  strH5out = os.path.join(strOut,"h5")
  os.makedirs(strH5out,exist_ok=True)
  
  H5pair=[]
  cmds = {}
  print("CellBender process starts ...\n\tFor more information, please visit https://cellbender.readthedocs.io/en/latest/usage/index.html")
  for i in range(meta.shape[0]):
    oneName=meta[sampleNameCol][i]
    if meta[CB_dropletNcol][i] < meta[CB_expCellNcol][i]:
      EXIT("Please make sure the %s column is larger than %s column for %s, details: \
        https://cellbender.readthedocs.io/en/latest/usage/index.html"%
        (CB_dropletNcol,CB_expCellNcol,oneName))
    oneH5 = os.path.join(strH5out,oneName+".cellbender.h5")
    H5pair += [{'old_path':meta[UMIcol][i],
      'new_path':re.sub(".h5$","_filtered.h5",oneH5),
      sampleNameCol:oneName}]
    if os.path.exists(oneH5):
      print("\t%s SKIP! CellBender H5 exists: \n\t\t%s"%(oneName,oneH5))
    else:
      cmds["CB_"+oneName] = "cellbender remove-background --input %s --output %s --cuda --expected-cells %d --total-droplets-included %d --fpr 0.01 --epochs 150 --low-count-threshold %d --learning-rate %f"%(
        meta[UMIcol][i],oneH5,meta[CB_expCellNcol][i],meta[CB_dropletNcol][i],meta[CB_count][i],meta[CB_learningR][i])
  if len(cmds)>0:
    scPipe.submit_cmd(cmds,{'parallel':'slurm','output':strOut,'gpu':True},core=1,memG=300000)
  cellbenderQC(H5pair,strOut)
  cellbenderInit(meta,H5pair,strOut)
  print("Before running the scAnalyzer, please check the log and pdf files in \n\t%s"%strH5out)

def cellbenderQC(H5pair,strOut):
  print("Cellbender QC ...")
  rmR={}
  with PdfPages("%s/cellbender_QC.pdf"%strOut) as pdf:
    cellN=[]
    for one in H5pair:
      print("\t"+one[sampleNameCol])
      old_adata = sc.read_10x_h5(one['old_path'])
      new_adata = sc.read_10x_h5(one['new_path'])
      cellN+=[{sampleNameCol:one[sampleNameCol],'cell_number':new_adata.shape[0]}]
      sc.pp.filter_cells(old_adata,min_counts=1)
      sc.pp.filter_cells(new_adata,min_counts=1)
      cInfo=pd.DataFrame({'CBrm':old_adata.obs.loc[new_adata.obs.index,'n_counts']-new_adata.obs['n_counts'],
        'CBkeepR':new_adata.obs['n_counts']/old_adata.obs.loc[new_adata.obs.index,'n_counts']})
      cInfo['CBrmR']=1-cInfo['CBkeepR']
      rmR[one[sampleNameCol]] = cInfo['CBrmR'].describe()
      one[ANNcol]=re.sub('.h5$','_qc.csv',one['new_path'])
      cInfo.to_csv(one[ANNcol])
      plotDensity(cInfo['CBrm'],
        "%s: removed UMI"%one[sampleNameCol],
        pdf)
      plotDensity(cInfo['CBkeepR'],
        "%s: UMI kept ratio"%one[sampleNameCol],
        pdf,bw=0.05)
    cellN=pd.DataFrame(cellN)
    plotBar(cellN,sampleNameCol,'cell_number',pdf,"CellBender filter cell number")
  pd.DataFrame(rmR).rename(index={'count':"cell_number"}).T.astype({'"cell_number"':'int'}).to_csv("%s/cellbender_rmRate.csv"%strOut,float_format='%.4f')

def plotDensity(dat,sTitle,pdf,bw='scott'):
  f, (ax_T, ax_B) = plt.subplots(2, sharex=True,figsize=(6,4),gridspec_kw={"height_ratios": (.8, .2)})
  sns.kdeplot(dat,bw_method=bw,ax=ax_T)
  sns.boxplot(dat,ax=ax_B)
  ax_T.set_title(sTitle)
  ax_B.set(xlabel='')
  pdf.savefig(bbox_inches="tight")
  plt.close(f)
def plotBar(dat,x,y,pdf,sTitle=None):
  h = 4.8
  w = max(6.4,dat.shape[0]*2/10)
  ax=dat.plot.bar(x=x,y=y,rot=90,legend=False,figsize=(w,h),grid=True)
  if not sTitle is None:
    ax.set_title(sTitle)
  pdf.savefig(bbox_inches="tight")
  plt.close()

def cellbenderInit(meta,H5pair,strOut):
  H5pair=pd.DataFrame(H5pair)
  H5pair.index=H5pair[sampleNameCol]
  meta[UMIcol]=list(H5pair.loc[meta[sampleNameCol],'new_path'])
  meta[ANNcol]=list(H5pair.loc[meta[sampleNameCol],ANNcol])
  meta.drop([CB_expCellNcol,CB_dropletNcol,CB_count,CB_learningR],axis=1,inplace=True,errors='ignore')
  scPipe.initSave(meta,strOut,False)

def main():
  #global strPipePath
  scPipe.strPipePath=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
  scPipe.checkInstallSetting()
  scPipe.MsgInit()
  
  strMeta = sys.argv[1]
  if os.path.isfile(strMeta):
    cellbender(os.path.realpath(strMeta))
  else:
    EXIT("The sample meta file is required, and %s doesn't exist!"%strPath)

if __name__ == "__main__":
  start_time = time.time()
  main()
  #print("--- total time passed %s seconds ---" % (time.time() - start_time))
