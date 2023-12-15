import scPipe, os, time, sys, re, shutil
import cmdUtility as cU
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from functools import reduce
#from PyPDF2 import PdfWriter,PdfReader
#from PyPDF2.generic import AnnotationBuilder
from pdf2image import convert_from_path
from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont
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

def cellbender(strMeta,nCore=0,useParallel=True):
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
  cupCMD=""
  print("CellBender process starts ...\n\tFor more information, please visit https://cellbender.readthedocs.io/en/latest/usage/index.html")
  if nCore==0:
    useGPU=True
    useCuda="--cuda "
    nCore=1
    mem=300000
  else:
    useGPU=False
    useCuda=""
    mem=16*int(nCore)
    cpuCMD=" --cpu-threads %d --estimator-multiple-cpu"%nCore
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
      cmds["CB_"+oneName] = "cellbender remove-background --input %s --output %s %s--expected-cells %d --total-droplets-included %d --fpr 0.01 --epochs 150 --low-count-threshold %d --learning-rate %f%s"%(
        meta[UMIcol][i],oneH5,useCuda,meta[CB_expCellNcol][i],meta[CB_dropletNcol][i],meta[CB_count][i],meta[CB_learningR][i],cpuCMD)
  if len(cmds)>0:
    print('\nrunning cellbeder for the following samples:\n%s'%'\n'.join([re.sub("^CB_","",i) for i in cmds.keys()]))
    cU.submit_cmd(cmds,{'parallel':'slurm' if useParallel else False,'output':strOut,'gpu':useGPU},core=nCore,memG=mem)
  cellbenderMergeLog(H5pair,strOut)
  cellbenderMergePdf(H5pair,strOut)
  cellbenderQC(H5pair,strOut)
  cellbenderInit(meta,H5pair,strOut)
  print("Before running the scAnalyzer, please check the log and pdf files in \n\t%s"%strH5out)

def cellbenderMergeLog(H5pair,strOut):
  print("\tmergeing log ...")
  with open(os.path.join(strOut,"all.log"),'wb') as outf:
    for one in H5pair:
      strF = re.sub("_filtered.h5$",".log",one['new_path'])
      if os.path.isfile(strF):
        outf.write(("\n\n***** %s *****\n"%one[sampleNameCol]).encode('utf-8'))
        with open(strF,'rb') as f:
          shutil.copyfileobj(f,outf)
def cellbenderMergePdf(H5pair,strOut):
  print("\tmergeing pdf ...")
  myFont = ImageFont.truetype('FreeSerifBold.ttf', 36)
  images=[]
  n=0
  for one in H5pair:
    n+=1
    strF = re.sub("_filtered.h5$",".pdf",one['new_path'])
    print("\t\t%d/%d: %s"%(n,len(H5pair),os.path.basename(strF)))
    if os.path.isfile(strF):
      img=convert_from_path(strF)[0]
      oneI = ImageDraw.Draw(img)
      oneI.text((400, 800),one[sampleNameCol],fill=(255, 0, 0),font=myFont)
      images = images+[img]
    else:
      print("\t\t\tmissing%s")
  images[0].save(os.path.join(strOut,'all_cellbender.pdf'),resolution=100,save_all=True,append_images=images[1:])
def cellbenderQC(H5pair,strOut):
  print("Cellbender QC ...")
  rmR={}
  with PdfPages("%s/cellbender_QC.pdf"%strOut) as pdf:
    cellN=[]
    geneRM={}
    geneRMr={}
    for one in H5pair:
      print("\t"+one[sampleNameCol])#,": %s vs %s"%(os.path.basename(one['old_path']),os.path.basename(one['new_path']))
      old_adata = sc.read_10x_h5(one['old_path'])
      old_adata.var_names_make_unique()
      new_adata = sc.read_10x_h5(one['new_path'])
      new_adata.var_names_make_unique()
      cellN+=[{sampleNameCol:one[sampleNameCol],'cell_number':new_adata.shape[0]}]
      sc.pp.filter_cells(old_adata,min_counts=1)
      sc.pp.filter_cells(new_adata,min_counts=1)
      sc.pp.filter_genes(old_adata, min_counts=1)
      sc.pp.filter_genes(new_adata, min_counts=0)
      # get cell removal information
      cInfo=pd.DataFrame({'CBrm':old_adata.obs.loc[new_adata.obs.index,'n_counts']-new_adata.obs['n_counts'],
        'CBkeepR':new_adata.obs['n_counts']/old_adata.obs.loc[new_adata.obs.index,'n_counts']})
      cInfo['CBrmR']=1-cInfo['CBkeepR']
      rmR[one[sampleNameCol]] = cInfo['CBrmR'].describe()
      one[ANNcol]=re.sub('.h5$','_qc.csv',one['new_path'])
      cInfo.to_csv(one[ANNcol])
      # get gene removal information
      gInfo=pd.DataFrame({'oriUMI':old_adata.var['n_counts'],
                          'newUMI':new_adata.var.loc[old_adata.var.index,'n_counts'],
                          'CBkeepR':new_adata.var.loc[old_adata.var.index,'n_counts']/old_adata.var['n_counts'],
                          'CBrmUMI':old_adata.var['n_counts']-new_adata.var.loc[old_adata.var.index,'n_counts']})
      gInfo['CBrmR']=1-gInfo['CBkeepR']
      gInfo.to_csv(re.sub('.h5$','_geneQC.csv',one['new_path']))
      geneRM[one[sampleNameCol]]=gInfo['CBrmUMI']
      geneRMr[one[sampleNameCol]]=gInfo['CBrmR']
      # plot
      plotDensity(cInfo['CBrm'],
        "%s: removed UMI"%one[sampleNameCol],
        pdf)
      plotDensity(cInfo['CBkeepR'],
        "%s: UMI kept ratio"%one[sampleNameCol],
        pdf,bw=0.05)
    cellN=pd.DataFrame(cellN)
    plotBar(cellN,sampleNameCol,'cell_number',pdf,"Cell number after CellBender")
    gRM=plotGeneQC(geneRM,geneRMr,pdf)
  pd.DataFrame(rmR).rename(index={'count':"cell_number"}).T.astype({'cell_number':'int'}).to_csv("%s/cellbender_rmRate.csv"%strOut,float_format='%.4f')
  gRM.to_csv("%s/cellbender_geneRM.csv"%strOut,float_format='%.4f')
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
def plotGeneQC(geneRM,geneRMr,pdf,topN=5):
  geneRM=pd.DataFrame(geneRM).fillna(0)
  geneRMr=pd.DataFrame(geneRMr).fillna(0)
  rmUMI = geneRM.apply(lambda x: x.max(),axis=1)
  rmUMIr = geneRMr.apply(lambda x: x.max(),axis=1)
  
  ## select the top removed genes
  gList = {'top %d genes removed by UMI'%topN:list(set(reduce(lambda x, y: x + y, geneRM.apply(lambda x: x.nlargest(topN).index).values.tolist(), []))),
    'top %d genes removed by rate'%topN:list(set(reduce(lambda x, y: x + y, geneRMr.loc[rmUMI.sort_values(ascending=False).index,].apply(lambda x: x.nlargest(topN).index).values.tolist(), [])))}
  selG = list(set(reduce(lambda x,y: x+y,gList.values(),[])))
  ## plot removed UMI against removed rate
  plt.scatter(rmUMI[rmUMI>0],rmUMIr[rmUMI>0],s=5,linewidths=0.5,marker='.',rasterized=True)
  plt.plot(rmUMI[selG],rmUMIr[selG],c='#ff7f0e',marker='*',ls='none',markersize=5)
  plt.xscale('log')
  plt.xlabel('Removed UMI')
  plt.ylabel('Removed UMI rate')
  plt.title('Largest removal across samples')
  pdf.savefig(bbox_inches="tight")
  plt.close()
  #print('Gain UMI in %d genes'%rmUMI[rmUMI<0].shape[0])
  if rmUMI[rmUMI<0].shape[0]>0:
    negG = rmUMI[rmUMI<0].index
    plt.scatter(rmUMI[negG],rmUMIr[negG],s=5,linewidths=0.5,marker='.',rasterized=True)
    for one in set(rmUMI.nsmallest(topN).index.tolist() + rmUMIr.nsmallest(topN).index.tolist()):
      plt.annotate(one, (rmUMI.abs()[one], rmUMIr.abs()[one]),xytext=(0,10),textcoords='offset pixels')
    plt.xlabel('Gained UMI')
    plt.ylabel('Gained UMI rate')
    plt.title('Largest UMI gained %d genes across samples'%negG.shape[0])
    pdf.savefig(bbox_inches="tight")
    plt.close()
  
  # plot top genes across samples
  for key in gList.keys():
    topRMgene=gList[key]
    A = geneRM.loc[topRMgene,:].melt(value_name="CBrmUMI",ignore_index=False).reset_index().merge(geneRMr.loc[topRMgene,:].melt(value_name="CBrmR",ignore_index=False).reset_index(),on=['index','variable'])
    fig = plt.figure(figsize=(max(5,len(topRMgene)/2),9))
    ax = plt.subplot(2,1,1)
    A.boxplot(by='index',column=['CBrmUMI'],ax=ax,rot=90)
    ax.set_yscale('log')
    ax.set_xlabel('')
    ax.set_ylabel('UMI removed')
    ax = plt.subplot(2,1,2)
    A.boxplot(by='index',column=['CBrmR'],ax=ax,rot=90)
    ax.set_xlabel('')
    ax.set_ylabel('UMI removed rate')
    fig.tight_layout()
    res=plt.suptitle(key)
    pdf.savefig(bbox_inches="tight")
    plt.close()
  
  return geneRM.add_suffix('_CBrmUMI').merge(geneRMr.add_suffix('_CBrmR'),left_index=True,right_index=True)
  
  
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
  nCore = int(sys.argv[2])
  useParallel= sys.argv[3]=="T"
  if os.path.isfile(strMeta):
    cellbender(os.path.realpath(strMeta),nCore,useParallel)
  else:
    EXIT("The sample meta file is required, and %s doesn't exist!"%strPath)

if __name__ == "__main__":
  start_time = time.time()
  main()
  #print("--- total time passed %s seconds ---" % (time.time() - start_time))
