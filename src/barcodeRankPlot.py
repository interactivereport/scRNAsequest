import os, math, warnings, logging
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import scanpy as sc
import pandas as pd
warnings.simplefilter(action='ignore', category=FutureWarning)
logging.disable()

UMIcol="h5path"
sampleNameCol="Sample_Name"
def EXIT(msg):
  print(msg)
  exit()

def plot(strRawMeta):
  meta = pd.read_csv(strRawMeta)#.iloc[0:2,:]
  if not sampleNameCol in meta.columns:
    meta[sampleNameCol]=[re.sub(".h5$","",re.sub(".raw_feature_bc_matrix.h5$","",os.path.basename(one))) for one in meta[UMIcol]]

  print("Barcode Rank Plot")
  cN = 3
  rN = math.ceil(meta.shape[0]/cN)
  fig = make_subplots(rows=rN,cols=cN,
    subplot_titles=list(meta[sampleNameCol]))
  for i in range(meta.shape[0]):
    fig.add_trace(plotOne(meta[UMIcol][i],meta[sampleNameCol][i]),
      row=int(i/cN)+1,col=i%cN+1)
  fig.update_xaxes(type='log',title_text="Cell Rank")
  fig.update_yaxes(type='log',title_text="UMI Count")
  fig.update_layout(height=400*rN, width=1500, title_text="Barcode Rank Plot",template='plotly_white')
  fig.write_html("%s/BarcodeRankPlot.html"%os.path.dirname(strRawMeta),
                full_html=False,
                include_plotlyjs='cdn')

def plotOne(oneH5,sName):
  print("\t"+sName)
  D = sc.read_10x_h5(oneH5)
  D.var_names_make_unique()
  sc.pp.filter_cells(D,min_counts=1)
  sc.pp.calculate_qc_metrics(D,inplace=True)
  x = list(range(1,D.shape[0]+1,1))
  y = sorted(D.obs['n_genes_by_counts'],reverse=True)
  if D.shape[0]>1000:
    x=list(range(100,len(y)+1,100))
    y=y[slice(99,len(y),100)]
  return go.Scatter(x=x, y=y,name="",showlegend=False)
