import sys, argparse, time, os, warnings, yaml, gc, h5py
import anndata as ad
import pandas as pd

maxG=500 # maximun number of genes can be exported
def MsgPower():
  strConfig = "%s/src/sys.yml"%strPipePath
  with open(strConfig,"r") as f:
    sysConfig = yaml.safe_load(f)
  print("\nPowered by %s"%sysConfig['powerby'])
  print("------------")
  exit()
  
def main():
  with open(sys.argv[1],"r") as f:
    args=parseInput(f.read())
  distributeTask(args.tool)(args)

def parseInput(strarg):
  parser = argparse.ArgumentParser(description='Additional sc tools. WARNING: The input h5ad files will be over-written!',formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('tool',type=str,choices=["rm","add","export"],help='Modify the h5ad by either "rm", "add" or "export" cell level annotations.')
  parser.add_argument('h5ad',type=str,help='Path to a h5ad file to be modified')
  parser.add_argument('changes',nargs="?",default="",help='Options:\n1. A list of annotaions to be removed (separated by ",") \n2. A path to a csv file contains cell level annotations (first column is the cell ID)\n3. A list of genes (separated by ",", empty or max %d) to be exported along with all annotations.'%maxG)

  if len(strarg)==0 or strarg=="-h" or strarg=="--help":
    parser.print_help()
    MsgPower()
  args = parser.parse_args(strarg.split())
  if args.tool in ["rm","add"] and len(args.changes)<3:
    print('"changes" are required for "rm" and "add" tool.')
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
    'export':exportAnnotation
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

if __name__ == "__main__":
  global strPipePath
  strPipePath=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
  main()
  MsgPower()
  #print("--- total time passed %s seconds ---" % (time.time() - start_time))

for one in selAnn:
  if X[one].dtype=="category":
    X[one]=X[one].astype("str")
