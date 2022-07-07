#!/usr/bin/env python

import sys, os, glob
strPath = sys.argv[1]
strH5ad = sys.argv[2]
print("working directory: %s"%strPath)

import pandas as pd
import anndata as annD
from scipy.sparse import csc_matrix
from scipy import sparse
import h5py

strH5 = '%s/X.h5'%strPath
strMTX = '%s/X.mtx'%strPath
strHdf5 ='%s/X.hdf5'%strPath
if os.path.isfile(strH5):
    D = annD.read_hdf(strH5,"Exp")
    D.obs_names = pd.Index([s.decode('utf-8') for s in D.obs_names])
    D.var_names = pd.Index([s.decode('utf-8') for s in D.var_names])
elif os.path.isfile(strMTX):
    D = annD.read_mtx(strMTX)
    D.var=pd.read_csv('%s/annData.var.csv'%strPath,index_col=0)
elif os.path.isfile(strHdf5):
    f = h5py.File(strHdf5,'r')
    X = sparse.csc_matrix((f['data'],f['indices'],f['indptr']),f['shape'])
    cID = [one.decode('utf-8') for one in f['row_names']]
    gID = [one.decode('utf-8') for one in f['col_names']]
    f.close()
    D = annD.AnnData(X)
    D.obs_names = pd.Index(cID)
    D.var_names = pd.Index(gID)
else:
    raise Exception("h5 or mtx file is missing in %s"%strPath)

X = pd.read_csv('%s/meta.csv'%strPath,index_col=0,low_memory=False)
for i in X.columns:
    if X[i].nunique()<100:
        X[i] = X[i].astype('category')
D.obs=X

for one in glob.glob("%s/layout*.csv"%strPath):
  X = pd.read_csv(one,index_col=0)
  strLayout = os.path.splitext(os.path.basename(one))[0].replace("layout.","")
  D.obsm["X_%s"%strLayout] = X.values

D.X = csc_matrix(D.X)

D.write(strH5ad)
print(strH5ad+" was created successfully!")


