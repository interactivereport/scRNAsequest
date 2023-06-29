import random, re, warnings, os, glob, yaml, logging, functools
import anndata as ad
import scanpy as sc

print=functools.partial(print, flush=True)
warnings.simplefilter(action='ignore', category=FutureWarning)
logging.disable(level=logging.INFO)

def inputCheck(args):
  strH5ad = args[1]
  if not os.path.isfile(strH5ad):
    msgError("ERROR: %s does not exist!"%strH5ad)
  strConfig = args[2]
  if not os.path.isfile(strConfig):
    msgError("ERROR: %s does not exist!"%strConfig)
  with open(strConfig,"r") as f:
    config = yaml.safe_load(f)
  return config

def splitBatch(strH5ad,strOut,batchCell=None,batchKey=None,hvgN=None):
  os.makedirs(strOut,exist_ok=True)
  if batchCell is None:
    print("Batch process (batchCell) is not set in config, large amount of memory might be required")
    tmpH5ad=os.path.join(strOut,"tmp.h5ad")
    os.remove(tmpH5ad) if os.path.isfile(tmpH5ad) else None
    os.symlink(strH5ad,tmpH5ad)
    h5adList=[tmpH5ad]
  else:
    print("Initialize batch process")
    h5adList= glob.glob(os.path.join(strOut,"tmp*.h5ad"))
    if len(h5adList)==0:
      print("Reading ...")
      D=ad.read_h5ad(strH5ad)
      hvgN = 3000 if hvgN is None else hvgN
      print("\thighly_variable_genes seurat_v3 ...")
      span=0.3
      print("\t\thttps://github.com/scverse/scanpy/issues/1504")
      while span<=1:
        print("\t\ttry with span: %.2f"%span)
        try:
          hvg = sc.pp.highly_variable_genes(D,flavor='seurat_v3',inplace=False,batch_key=batchKey,n_top_genes=hvgN,span=span)
          break
        except:
          span += 0.1
      if span>1:
        msgError("Cannot find highly_variable_genes with span<%.2f, possible: too many low expressed genes in some batches!"%span)
      hvg[hvg.highly_variable].to_csv(os.path.join(strOut,"hvg.csv"))
      sampleCellN = D.obs[batchKey].value_counts()
      print(sampleCellN)
      sID=[]
      cellN=0
      batchN=0
      sName=list(sampleCellN.index)
      random.shuffle(sName)
      for one in sName:
        sID.append(one)
        cellN+=sampleCellN[one]
        #print(cellN)
        if cellN>batchCell or one==sName[-1]:
          print("batch %d: %d samples"%(batchN,len(sID)))
          strH5ad=os.path.join(strOut,"tmp_%d.h5ad"%batchN)
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
    else:
      print("--> Using previous batches in %s <--\nPlease remove/rename the folder, if new batches is desired."%strOut)
  return(h5adList)
