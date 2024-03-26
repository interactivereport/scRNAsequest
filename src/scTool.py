import sys, argparse, time, os, warnings, yaml, gc, h5py, re
import anndata as ad
import pandas as pd
import numpy as np
from scipy.sparse import issparse, csc_matrix

maxG=500 # maximun number of genes can be exported
def MsgPower():
  strConfig = "%s/src/sys.yml"%strPipePath
  with open(strConfig,"r") as f:
    sysConfig = yaml.safe_load(f)
  print("\nPowered by %s"%sysConfig['powerby'])
  print("------------")
  exit()
  
def main():
  args = parseInput()
  distributeTask(args.tool)(args)

def parseInput():
  parser = argparse.ArgumentParser(description='Additional sc tools. WARNING: The input h5ad files will be over-written!',formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('tool',type=str,choices=["rm","add","export","pseudo"],help='Modify the h5ad by either "rm", "add" or "export" cell level annotations. \nCreate pseudo bulk as RNAsequest analysis input by "pseudo"')
  parser.add_argument('h5ad',type=str,help='Path to a h5ad file to be modified or extract from (raw UMI for pseudo)')
  parser.add_argument('changes',nargs="?",default="",help='Options:\n \
  1. (rm) A list of annotaions to be removed (separated by ",") \n \
  2. (add) A path to a csv file contains cell level annotations (first column is the cell ID)\n \
  3. (export) A list of genes (separated by ",", empty or max %d) to be exported along with all annotations.\n \
  4. (pseudo) A list of obs to group cells into pseudo bulk, separated by ","'%maxG)
  args = parser.parse_args()
  print(args)
  if args.tool in ["rm","add",'pseudo'] and len(args.changes)<3:
    print('"changes" are required for "rm", "add" and "pseudo" tool.')
    MsgPower()
  if not os.path.isfile(args.h5ad):
    print("%s does NOT exist!"%args.h5ad)
    MsgPower()
  args.h5ad = os.path.realpath(args.h5ad)
  return(args)

def errorTask(data):
  raise ValueError('Error task!')
def distributeTask(aTask):
  return {
    'rm':rmAnnotation,
    'add':addAnnotation,
    'export':exportAnnotation,
    'pseudo':pseudoBulk
  }.get(aTask,errorTask)
  
def rmAnnotation(args):
  try:# there is NO close for anndata 
    D = ad.read_h5ad(args.h5ad,backed="r+")
    args.changes = args.changes.split(",")
    lackAnno = [one for one in args.changes if not one in D.obs.columns]
    if len(lackAnno)>0:
      print("The following annotations are NOT available:")
      print("%s"%", ".join(lackAnno))
    rmAnno = [one for one in args.changes if one in D.obs.columns]
    if len(rmAnno)>0:
      D.obs.drop(rmAnno,axis=1,inplace=True)
      with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        D.write(args.h5ad)
      print("The following annotations were removed from %s:"%args.h5ad)
      print("%s"%", ".join(rmAnno))
  finally:
    for obj in gc.get_objects():   # Browse through ALL objects
      if isinstance(obj, h5py.File) and obj.__bool__():   # Just opened HDF5 files
        obj.close()

def addAnnotation(args):
  if not os.path.isfile(args.changes):
    print("The cell level annotation file (%s) does NOT exist!"%args.changes)
    MsgPower()
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    X = pd.read_csv(args.changes,sep=None,index_col=0,header=0)
  X.index=list(X.index)
  try:# there is NO close for anndata 
    D = ad.read_h5ad(args.h5ad,backed="r+")
    selAnn = [one for one in X.columns if not one in D.obs.columns]
    if len(selAnn)>0:
      print("The following annocation will be added for %d of %d cells of:"%(len(list(set(X.index)&set(D.obs_names))),D.shape[0]))
      print(", ".join(selAnn))
      for one in selAnn:
        if X[one].dtype=="category":
          X[one]=X[one].astype("str")
      with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tmp = D.obs.merge(X[selAnn],"left",left_index=True,right_index=True)[selAnn]
        D.obs = D.obs.merge(tmp.fillna("missing"),"left",left_index=True,right_index=True)
        D.write(args.h5ad)
  finally:
    for obj in gc.get_objects():   # Browse through ALL objects
      if isinstance(obj, h5py.File) and obj.__bool__():   # Just HDF5 files
        obj.close()

def exportAnnotation(args):
  D = ad.read_h5ad(args.h5ad,backed="r")
  X = D.obs
  # first two embedding
  for one in D.obsm_keys():
    emBed = pd.DataFrame(D.obsm[one],index=D.obs_names).iloc[:,range(2)]
    emBed.columns = ["%s_%d"%(one,i+1) for i in range(2)]
    X= X.merge(emBed,"left",left_index=True,right_index=True)
  
  if len(args.changes)>0:
    genes = args.changes.split(",")
    lackG = [one for one in genes if not one in D.var_names]
    if len(lackG)>0:
      print("The following genes are NOT available:")
      print("%s"%", ".join(lackG))
    genes = [one for one in genes if one in D.var_names]
    if len(genes)>0:
      if len(genes)>maxG:
        print("Only the first %n genes will be exported!"%maxG)
        genes = genes[0:maxG]
      with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        X = X.merge(D[:,genes].to_df(),'left',left_index=True,right_index=True)
  strCSV="%s.csv"%args.h5ad
  X.to_csv(strCSV)
  print("Exporting to: %s"%strCSV)

def pseudoBulk_init(args):
  strOut = os.path.join(os.path.dirname(args.h5ad),'pseudoBulk')
  if os.path.exists(strOut):
    raise Exception("%s exists!\n\tPlease rename/remove the above path!"%strOut)
  os.makedirs(strOut,exist_ok=True)
  return strOut
def pseudoBulk_check(D,col):
  if D.X is None and (D.raw is None or D.raw.X is None):
    raise AttributeError("Provided Matrix is None")
  if not D.raw is None and not D.raw.X is None:
    if D.raw.var is None:# or D.raw.obs is None
      raise AttributeError(".raw found in h5ad, but missing .raw.var")
    else:
      print(".raw exists, extract from .raw.X, .obs will be used instead of .raw.obs")
  if D.obs is None:
    raise AttributeError("Provided Obs are None")
  if D.var is None:
    raise AttributeError("Provided Var are None")
  if D.obs.columns.isin(col).sum() != len(col):
    print(col)
    raise AttributeError("At least one obs is not available in h5ad")
pseudoBulk_col="pseudo_sel"
def pseudoBulk_one(D,one):
  print("\t",one)
  # extract pseudo bulk sum
  if not D.raw is None and not D.raw.X is None:
    selC = D.raw.obs_names.isin(D.obs[D.obs[pseudoBulk_col]==one].index)
    ix = [i for i in range(len(selC)) if selC[i]]
    X = D.raw.X[ix,]
    gID = D.raw.var.index
  else:
    selC = D.obs[pseudoBulk_col]==one
    ix = [i for i in range(len(selC)) if selC[i]]
    X = D.X[ix,]
    gID = D.var.index
  if not issparse(X):
    X = csc_matrix(X.copy())
  exp = pd.DataFrame(X.sum(axis=0),columns=gID).iloc[0,:]
  exp.name=one
  # extract unique obs
  obs = D.obs[D.obs[pseudoBulk_col]==one].copy()
  obs = obs.loc[:,obs.nunique()==1].iloc[0,:]
  #obs['pseudo_cellN'] = X.shape[0]
  obs.name=one
  qc = pd.Series({'pseudo_cellN':int(X.shape[0]),
                  'pseudo_UMI':int(exp.sum())},name=one)
  return exp,obs,qc
def pseudoBulk_save(X,meta,qc,strOut):
  print("Saving ...")
  X.to_csv(os.path.join(strOut,"genes.UMI.tsv"),sep="\t",index_label="Gene")
  meta.to_csv(os.path.join(strOut,"samplesheet.tsv"),sep="\t",index_label="Sample_Name")
  qc.to_csv(os.path.join(strOut,"sampleQC.tsv"),sep="\t",index_label="Sample_Name")
  X1 = X.apply(lambda x: x/x.sum()*1e6)
  X1.to_csv(os.path.join(strOut,"genes.CPM.tsv"),sep="\t",index_label="Gene")
  X1 = X/qc['pseudo_cellN']*1000
  X1.to_csv(os.path.join(strOut,"genes.CPKcell.tsv"),sep="\t",index_label="Gene")
def pseudoBulk(args):
  strOut = pseudoBulk_init(args)
  D = ad.read_h5ad(args.h5ad,backed="r")
  metaCol = args.changes.split(",")
  pseudoBulk_check(D,metaCol)
  D.obs[pseudoBulk_col] = D.obs[metaCol].apply(lambda x: "_".join(x),axis=1)
  X = []
  meta = []
  qc=[]
  for one in D.obs[pseudoBulk_col].unique():
    x,m,q = pseudoBulk_one(D,one)
    X.append(x)
    meta.append(m)
    qc.append(q)
  X = pd.DataFrame(X).transpose()
  X1=X.astype(np.int32)
  meta = pd.DataFrame(meta)
  meta = meta.dropna(axis=1)
  qc = pd.DataFrame(qc)
  pseudoBulk_save(X,meta,qc,strOut)
  
if __name__ == "__main__":
  global strPipePath
  strPipePath=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
  main()
  MsgPower()
  #print("--- total time passed %s seconds ---" % (time.time() - start_time))

for one in selAnn:
  if X[one].dtype=="category":
    X[one]=X[one].astype("str")
