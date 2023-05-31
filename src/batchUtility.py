import random, re, warnings, os, glob
import anndata as ad
import scanpy as sc


def splitBatch(strH5ad,strOut,batchCell=None,batchKey=None,hvgN=None):
  if batchCell is None:
    print("Batch process (batchCell) is not set in config, large amount of memory might be required")
    h5adList=[strH5ad]
  else:
    print("Initialize batch process")
    os.makedirs(strOut,exist_ok=True)
    h5adList= glob.glob(strOut+"*.h5ad")
    if len(h5adList)==0:
      print("Reading ...")
      D=ad.read_h5ad(strH5ad)
      hvgN = 3000 if hvgN is None else hvgN
      hvg = sc.pp.highly_variable_genes(D,flavor='seurat_v3',inplace=False,batch_key=batchKey,n_top_genes=hvgN)
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