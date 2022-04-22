import yaml, io, os, sys, subprocess, errno, json, re, logging, warnings, shutil, time, random, pwd, math
import pandas as pd
from datetime import datetime
import scanpy as sc
import anndata as ad
import numpy as np
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

strPipePath=""
UMIcol="h5path"
ANNcol="metapath"
batchKey="library_id"
# maximum number of a parallel job can be re-submited is 5
maxJobSubmitRepN=5 
logging.disable()
def run_cmd(cmd):
  #print(cmd)
  try:
    cmdR = subprocess.run(cmd,shell=True,check=True,stdout=subprocess.PIPE)#capture_output=True
  except subprocess.CalledProcessError as e:
    if not e.returncode==1:
      print("%s process return above error"%cmd)#,e.stderr.decode("utf-8")
    cmdR = e
  return cmdR
def Exit(msg=""):
  if len(msg)>3:
    print(msg)
  MsgPower()
  exit()
def getConfig(strConfig="sys.yml",blines=False,bSys=True):
  if bSys:
    strConfig = "%s/src/%s"%(strPipePath,strConfig)
  with open(strConfig,"r") as f:
    if blines:
      config = f.readlines()
    else:
      config = yaml.safe_load(f)
  return(config)
def checkInstallSetting():
  ## check sys config file
  strConfig = "%s/src/sys.yml"%(strPipePath)
  if not os.path.isfile(strConfig):
    print("=====\nPlease contact admin to set the sys.yml in %s.\nAn Example is 'sys_example.yml'.\n====="%strPipePath)
    exit()

## parallel job management 
def submit_cmd(cmds,config,core=None):
  #cmds = {k:v for i, (k,v) in enumerate(cmds.items()) if not v is None}
  if len(cmds)==0:
    return
  if core is None:
    core=config['core']
  parallel = config["parallel"]
  if not parallel:
    for one in cmds.keys():
      if cmds[one] is None:
        continue
      print("submitting %s"%one)
      try:
        cmdR = subprocess.run(cmds[one],shell=True,check=True)#,capture_output=True
      except subprocess.CalledProcessError as e:
        print("%s process error return:\n %s"%(one,e.stderr.decode("utf-8")))
  elif parallel=="sge":
    jID = qsub(cmds,config['output'],core)
    print("----- Monitoring all submitted SGE jobs ...")
    ## in case of long waiting time to avoid Recursion (too deep)
    cmdN = {one:1 for one in cmds.keys()}
    while True:
      qstat(jID,config['output'],cmds,cmdN,core)
      if len(cmds)==0:
        break
      ## wait for 1 min if any running jobs
      time.sleep(60)
  elif parallel=="slurm":
    pass
  else:
    print("ERROR: unknown parallel setting: %s"%parallel)
    exit()

def qsub(cmds,strPath,core,jID=None):
  if jID is None:
    jID = "j%d"%random.randint(10,99)
  strWD = os.path.join(strPath,jID)
  try:
    os.makedirs(strWD)
  except FileExistsError:
    pass
  noRun = []
  for one in cmds.keys():
    oneScript = (qsubScript.replace('qsubCore',str(core))
                            .replace('jName',one)
                            .replace('wkPath',strWD)
                            .replace("jID",jID)
                            .replace("sysPath",strPipePath)
                            .replace("strCMD",cmds[one]))
    strF=os.path.join(strWD,"%s.sh"%one)
    with open(strF,"w") as f:
      f.write(oneScript)
    run_cmd("qsub %s"%strF)
  for one in noRun:
    del cmds[one]
  return jID
def qstat(jID,strPath,cmds,cmdN,core):
  print(".",end="")
  strWD = os.path.join(strPath,jID)
  qstateCol=4 # make sure the state of qstat is the 4th column (0-index)
  qJobCol=0 # make sure the job-id of qstat is the 0th column (0-index)
  # remove the job in wrong stats: Eqw and obtain the running and waiting job names
  qs = run_cmd("qstat | grep '%s_'"%jID).stdout.decode("utf-8").split("\n")
  jNames = []
  for one in qs:
    if len(one)<5:
      continue
    oneJob = [a.strip() for a in one.split(" ") if len(a.strip())>0 ]
    if oneJob[qstateCol]=='Eqw':
      run_cmd("qdel %s"%oneJob[qJobCol])
      continue
    jobDetail = [' '.join(one.split()) for one in run_cmd("qstat -j %s"%oneJob[qJobCol])
                                                    .stdout
                                                    .decode("utf-8")
                                                    .split("\n")]
    aProperty="job_name"
    jName=[a.replace("%s: %s_"%(aProperty,jID),"").strip() for a in jobDetail if aProperty in a]
    if len(jName)==0:
      continue
    jNames.append(jName[0])
  # check the finished job, resubmit if not sucessfully finished
  resub = {}
  finishedJob = []
  for one in cmds.keys():
    if not one in jNames:
      strLog = os.path.join(strWD,"%s.log"%one)
      if not "DONE" in run_cmd("tail -n 1 %s"%strLog).stdout.decode("utf-8"):
        cmdN[one] += 1
        if cmdN[one]>maxJobSubmitRepN:
          print("\n--->ERROR failed with %d times qsub: %s"%(maxJobSubmitRepN,one))
          finishedJob.append(one)
          continue
        resub[one] = cmds[one]
      else:
        print("\n\t===== %s ====="%one)
        with open(strLog,"r") as f:
          print(f.read())
        print("\tFinished: %s"%one)
        finishedJob.append(one)
  if len(resub)>0:
    re1=qsub(resub,strPath,core,jID)
    time.sleep(5) #might not needed for the qsub to get in
  for one in finishedJob:
    del cmds[one]
  ### in case of long waiting time to avoid Recursion (too deep)
  #if len(cmds)>0:
  #  ## wait for 1 min if any running jobs
  #  time.sleep(60)
  #  qstat(jID,strPath,cmds,core,iN)
  
## msg
def MsgPower():
  sysConfig = getConfig()
  print("\nPowered by %s"%sysConfig['powerby'])
  print("------------")
def MsgHelp():
  MsgInit()
  print("\nscAnalyzer /path/to/a/DNAnexus/download/folder === or === scAnalyzer /path/to/a/config/file\n")
  print("The config file will be generated automatically when a DNAnexus download folder is provided")
  print("Available reference data:")
  sysConfig = getConfig()
  for i in sysConfig['ref']:
    print("\t%s: more information @ %s"%(i,sysConfig[i]["ref_link"]))
  print("If one of the above can be used as a reference for your datasets, please update the config file with the name in 'ref_name'.\n")
  Exit()
def MsgInit():
  cmdURL="cd %s;git config --get remote.origin.url"%strPipePath
  cmdDate="cd %s;git show -s --format=%%ci"%strPipePath
  cmdHEAD="cd %s;git rev-parse HEAD"%strPipePath
  print("###########\n## scAnalyzer: %s"%run_cmd(cmdURL).stdout.decode("utf-8").replace("\n",""))
  print("## Pipeline Path: %s"%strPipePath)
  print("## Pipeline Date: %s"%run_cmd(cmdDate).stdout.decode("utf-8").replace("\n",""))
  print("## git HEAD: %s###########\n"%run_cmd(cmdHEAD).stdout.decode("utf-8"))
  #print("\nLoading resources")

## init projects
def initProject(strDNAnexus):
  if os.path.isdir(strDNAnexus):
    strDNAnexus = os.path.realpath(strDNAnexus)
    meta = initMeta(strDNAnexus)
    if meta is None:
      strConfig = initExternal(strDNAnexus)
    else:
      strConfig = initSave(meta,strDNAnexus)
    initMsg(strConfig)
    exit()
def initMeta(strInput):
  print("Procing sample meta information ...")
  strMeta = "%s/samplesheet.tsv"%strInput
  if not os.path.isfile(strMeta):
    return None
  meta = pd.read_csv(strMeta,sep="\t")
  sysConfig = getConfig()
  config = getConfig("template.yml")
  
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    if not 'notMeta' in sysConfig.keys():
      sysConfig['notMeta'] = []
    meta = meta[[one for one in meta.columns if not one in sysConfig['notMeta']]].dropna(1,'all')
  meta.insert(0,UMIcol,["%s/%s.filtered_feature_bc_matrix.h5"%(strInput,one) for one in meta[config["sample_name"]]])
  for oneH5 in meta[UMIcol]:
    if not os.path.isfile(oneH5):
      Exit("The UMI h5 file %s does not exist, please correct the sample sheet"%oneH5)
  return meta
def initSave(meta,strInput):
  print("Saving the initialized project ...")
  ix = 0
  strOut="%s/sc%s_%d/"%(strInput,datetime.today().strftime('%Y%m%d'),ix)
  while os.path.isdir(strOut):
    ix +=1
    strOut="%s/sc%s_%d"%(strInput,datetime.today().strftime('%Y%m%d'),ix)
  os.makedirs(strOut)
  
  # save meta file
  strMeta = os.path.join(strOut,"sampleMeta.csv")
  meta.to_csv(strMeta,index=False)
  
  # save empty DE csv table
  strDEG = os.path.join(strOut,"DEGinfo.csv")
  with open(strDEG,"w") as f:
    f.write("sample,cluster,group,alt,ref,covars[+ separated],method[default NEBULA],model[default HL]\n")

  # save config file
  config = getConfig("template.yml",True)
  strJson = "%s/samplesheet.json"%strInput
  if os.path.isfile(strJson):
    with open(strJson,"r") as f:
      prjInfo = json.load(f)
    config = [one.replace("initPrjName",prjInfo['Project']["TSTID"])
                  .replace("initPrjTitle",prjInfo['Project']["Study_Title"]) for one in config]
  else:
    print("Project information (%s) is not avaiable."%strJson)
    print("===> Please update the prj_name and prj_title manually in config file <===")
  config = [one.replace("initOutput",strOut)
                .replace("initPrjMeta",strMeta)
                .replace("initDEG",strDEG)for one in config]
  strConfig = os.path.join(strOut,"config.yml")
  with open(strConfig,"w") as f:
    f.writelines(config)
  return strConfig
def initExternal(strInput):
  config = getConfig("template.yml")
  print("*****\nExternal data, please fill the config file and provide sample sheet with two mandatory columns:")
  print("Two columns with names '%s' and '%s' indicating sample name and UMI path"%(config['sample_name'],UMIcol))
  print("Option columns for sample annotation")
  print("Option column '%s' can be provided for cell level annotation first column cell bar code\n*****"%ANNcol)
  # save empty sample meta file with columns
  strMeta = os.path.join(strInput,"sampleMeta.csv")
  with open(strMeta,"w") as f:
    f.writelines(["%s,%s,%s\n"%(config['sample_name'],UMIcol,ANNcol)])
  # save empty DE csv table
  strDEG = os.path.join(strInput,"DEGinfo.csv")
  with open(strDEG,"w") as f:
    f.write("sample,cluster,group,alt,ref,covars[+ separated],method[default NEBULA],model[default HL]\n")

  configL = getConfig("template.yml",blines=True)
  configL = [one.replace("initOutput",strInput) for one in configL]
  configL = [one.replace("initPrjMeta",strMeta) for one in configL]
  configL = [one.replace("initDEG",strDEG) for one in configL]
  configL = [re.sub(r'\b(init\w*)',' #required',one) for one in configL]
  strConfig = os.path.join(strInput,"config.yml")
  with open(strConfig,"w") as f:
    f.writelines(configL)
  return strConfig
def initMsg(strConfig):
    print("\nSingle cell/nuclei RNAseq project was created @%s"%os.path.dirname(strConfig))
    print("Please update the config file, if any reference can be used.")
    print("\tReference information can be access by 'scAnalyzer' without any input argument.")
    print("Additional columns can be provided into the sample meta table (sampleMeta.csv).")
    print("DEG table (DEGinfo.csv) can be filled later and rerun the following command.")
    print("Please run the following command to use the pipeline for the input dataset.")
    print("\n===> scAnalyzer %s"%strConfig)
    MsgPower()

## pipeline run 
def runPipe(strConfig):
  MsgInit()
  sc.settings.n_jobs=1
  config = getConfig(strConfig,bSys=False)
  checkConfig(config)
  meta = getSampleMeta(config["sample_meta"])
  prefix = runQC(config,meta)
  
  methods = runMethods(prefix,strConfig)
  combine(methods,prefix,config)
  moveCellDepot(prefix)
  runDEG(strConfig,prefix,config)
  MsgPower()
def checkConfig(config):
  if config["prj_name"] is None:
    Exit("'prj_name' is required in config file!")
  if config['sample_meta'] is None:
    Exit("'sample_meta' is required in config file!")
def runQC(config,meta):
  plotSeqQC(meta,config["sample_name"],config["output"],config["group"])
  prefix = os.path.join(config["output"],config["prj_name"])
  if not config["runAnalysis"] or not os.path.isfile("%s_raw.h5ad"%prefix) or config["newProcess"]:
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      if not os.path.isfile("%s_raw_prefilter.h5ad"%prefix) or config["newProcess"]:
        adata = getData(meta,config["sample_name"])
        adata.write("%s_raw_prefilter.h5ad"%prefix)
      else:
        adata = sc.read_h5ad("%s_raw_prefilter.h5ad"%prefix)
      adata = preprocess(adata,config)
      adata.obs,adata.obsm=obtainRAWobsm(adata.copy())
      plotQC(adata,config["output"],config["group"])
      adata.write("%s_raw.h5ad"%prefix)
    if not config["runAnalysis"]:
      runQCmsg(config)
      exit()
  return prefix
def runQCmsg(config):
  print("Please check the following QC files @%s:\n\tsequencingQC.csv\n\tsequencingQC.pdf\n\tumiQC.pdf"%config['output'])
  print("And then:")
  print("\t1. Remove any outlier sample from the sample meta table %s")
  print("\t2. Update config file on cutoff values for cell filtering")
  print("\t3. After making sure the cell filtering setting (might several iteration), set 'runAnalysis: True' in the config file.")
  print("\t4. (Optional) consider to enable parallel by setting: 'parallel: sge' for CAMHPC or 'parallel: slurm' for EdgeHPC.")
  MsgPower()
def plotSeqQC(meta,sID,strOut,grp=None):
  print("plotting sequence QC ...")
  seqQC = []
  for i in range(meta.shape[0]):
    strF = os.path.join(os.path.dirname(meta[UMIcol][i]),"%s.metrics_summary.csv"%meta[sID][i])
    if os.path.isfile(strF):
      one = pd.read_csv(strF,thousands=",")
      one.index=[meta[sID][i]]
      seqQC.append(one)
    else:
      return
  QC = pd.concat(seqQC)
  k=list(QC.columns)
  for i,one in enumerate(k):
    if pd.api.types.is_string_dtype(QC[one]) and QC[one][0].endswith("%"):
      QC[one] = QC[one].str.rstrip('%').astype('float')
      k[i]=one+"%"
  QC.columns=k
  QC.to_csv("%s/sequencingQC.csv"%strOut)
  QC['sample'] = QC.index
  with PdfPages("%s/sequencingQC.pdf"%strOut) as pdf:
    for one in k:
      ax = QC.plot.bar(x='sample',y=one,rot=90,legend=False)
      ax.set_title(one)
      plt.grid()
      pdf.savefig(bbox_inches="tight")
      plt.close()
    if not grp==None:
      for oneG in grp:
        if oneG in meta.columns:
          QC[oneG] = list(meta[oneG])
          for one in k:
            ax = QC[[one,oneG]].boxplot(by=oneG,rot=90)
            ax.set_title(one)
            plt.grid()
            pdf.savefig(bbox_inches="tight")
            plt.close()
def getSampleMeta(strMeta):
  print("processing sample meta information ...")
  if not os.path.isfile(strMeta):
    Exit("Sample meta information (%s) specified in config file does not exist!"%strMeta)
  meta = pd.read_csv(strMeta)
  if not UMIcol in meta.columns:
    Exit("'%s' columns is required in meta file (%s)"%(UMIcol,strMeta))
  for oneH5 in meta[UMIcol]:
    if not os.path.isdir(oneH5) and not (os.path.isfile(oneH5) and oneH5.endswith(".h5")):
      Exit("The UMI file %s does not exist, please correct the sample sheet"%oneH5)
  return(meta)
def getData(meta,sID):
  print("processing sample UMI ...")
  adatals = []
  for i in range(meta.shape[0]):
    print("\t%s"%meta[sID][i])
    if os.path.isdir(meta[UMIcol][i]):
      adata = sc.read_10x_mtx(meta[UMIcol][i])
    elif meta[UMIcol][i].endswith('.h5'):
      adata = sc.read_10x_h5(meta[UMIcol][i])
    else:
      Exit("Unsupported UMI format: %s"%meta[UMIcol][i])
    adata.var_names_make_unique()

    if ANNcol in meta.columns:
      annMeta = pd.read_csv(meta[ANNcol][i],index_col=0)
      adata = adata[adata.obs.index.isin(list(annMeta.index))]
      adata.obs = pd.merge(adata.obs,annMeta,left_index=True,right_index=True)
      print("\t\tCell level meta available, cell number: %d"%adata.shape[0])
      
    for one in meta.columns:
      if not 'path' in one and not one==sID:
        adata.obs[one]=meta[one][i]
    adatals.append(adata)
  
  adata = sc.AnnData.concatenate(*adatals,
    batch_categories=meta[sID],
    batch_key=batchKey)
  ## remove duplicated columns in var
  varCol = [one.split("-")[0] for one in adata.var.columns]
  varInx = [i for i,v in enumerate(varCol) if not v in varCol[:i]]
  adata.var = adata.var.iloc[:,varInx]
  adata.var.columns=[varCol[i] for i in varInx]
  return(adata)
def preprocess(adata,config):
  MTstring=config["MTstring"]
  min_cells=config["min.cells"]
  min_features=config["min.features"]
  mt_cutoff=config["mt.cutoff"]
  highCount_cutoff=config["highCount.cutoff"]
  highGene_cutoff=config["highGene.cutoff"]

  print("preprocessing ...")
  sc.pp.calculate_qc_metrics(adata, inplace=True)
  #get MT genes
  if MTstring is None or len(MTstring)<2:
    for one in ["MT-","Mt-","mt-"]:
      mito_genes = adata.var_names.str.startswith(one)
      if (mito_genes).sum()>0:
        break
  else:
    mito_genes = adata.var_names.str.startswith(MTstring)
  adata.obs['pct_counts_mt'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1 * 100
  print("\tall mitochondrial genes removed")
  adata = adata[:, np.invert(mito_genes)]
  ## filtering low content cells and low genes
  print("\tfiltering cells and genes")
  sc.pp.filter_genes(adata,min_cells=min_cells)
  print("\t\tfiltered genes with min.cells %d left %d genes"%(min_cells,adata.shape[1]))
  sc.pp.filter_cells(adata,min_genes=min_features)
  print("\t\tfiltered cells with min.features %d left %d cells"%(min_features,adata.shape[0]))
  
  adata = adata[adata.obs.pct_counts_mt<mt_cutoff,:]
  print("\t\tfiltered cells with mt.cutoff %d left %d cells"%(mt_cutoff,adata.shape[0]))
  adata = adata[adata.obs.n_genes_by_counts<highGene_cutoff,:]
  print("\t\tfiltered cells with highGene.cutoff %d left %d cells"%(highGene_cutoff,adata.shape[0]))
  adata = adata[adata.obs.total_counts<highCount_cutoff,:]
  print("\t\tfiltered cells with highCount.cutoff %d left %d cells"%(highCount_cutoff,adata.shape[0]))
  if adata.shape[0]<10:
    Exit("Few cells (%d<10) left after filtering, please check the filtering setting in config to contitue!"%adata.shape[0])
  return adata
def obtainRAWobsm(D):
  # 95 percentile to normalize
  print("\tinitializing layout")
  sc.pp.normalize_total(D,target_sum=math.ceil(np.percentile(D.X.sum(axis=1).transpose().tolist()[0],95)))
  sc.pp.log1p(D)
  sc.pp.highly_variable_genes(D, min_mean=0.01, max_mean=3, min_disp=0.5)
  D = D[:, D.var.highly_variable]
  sc.pp.regress_out(D, ['total_counts', 'pct_counts_mt'])
  sc.pp.scale(D, max_value=10)
  sc.tl.pca(D, svd_solver='arpack',n_comps = 100)
  npcs = 50
  sc.pp.neighbors(D, n_neighbors=10, n_pcs=npcs)
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    sc.tl.louvain(D,key_added="raw.louvain")
    ## umap embedding 
    print("\tembedding ...")
    sc.tl.umap(D, init_pos='spectral')
    sc.tl.rank_genes_groups(D, 'raw.louvain')
    sc.tl.tsne(D, n_pcs=npcs)
  return D.obs,D.obsm #, D.var.highly_variable
def plotQC(adata,strOut,grp=None):
  print("plotting UMI QC ...")
  with PdfPages("%s/umiQC.pdf"%strOut) as pdf:
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',color=batchKey,alpha=0.5)
    pdf.savefig(bbox_inches="tight")
    sc.pl.highest_expr_genes(adata, n_top=20)
    pdf.savefig(bbox_inches="tight")
    sc.pl.violin(adata, keys = 'n_genes_by_counts',groupby=batchKey,rotation=90)
    pdf.savefig(bbox_inches="tight")
    sc.pl.violin(adata, keys ='pct_counts_mt',groupby=batchKey,rotation=90)
    pdf.savefig(bbox_inches="tight")
    if not grp==None:
      for oneG in grp:
        if oneG in adata.obs.columns:
          sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',color=oneG,alpha=0.5)
          pdf.savefig(bbox_inches="tight")
          sc.pl.violin(adata, keys = 'n_genes_by_counts',groupby=oneG,rotation=90)
          pdf.savefig(bbox_inches="tight")
          sc.pl.violin(adata, keys ='pct_counts_mt',groupby=oneG,rotation=90)
          pdf.savefig(bbox_inches="tight")

def runMethods(prefix,strConfig):
  print("starting the process by each method ...")
  sysConfig = getConfig()
  config = getConfig(strConfig,bSys=False)
  checkLock(config,sysConfig)
  
  cmds = {}
  allM = sysConfig['methods']
  for m in allM.keys():
    if not config["newProcess"] and os.path.isfile("%s_%s.h5ad"%(prefix,m)):
      continue
    # /home/zouyang/projects/scRNAsequest/src/SCT.py 
    # /home/zouyang/projects/scRNAsequest/src/harmony.py
    # /home/zouyang/projects/scRNAsequest/src/seurat.py
    # /home/zouyang/projects/scRNAsequest/src/seuratRef.py
    
    # /camhpc/ngs/projects/TST11837/dnanexus/20220311155416_maria.zavodszky/sc20220403_0/TST11837_raw.h5ad config.yml
    cmd="%s/src/%s %s_%s.h5ad %s"%(strPipePath,
                                    allM[m][0],
                                    prefix,
                                    allM[m][1],
                                    strConfig)
    cmds[m]=cmd
  submit_cmd(cmds,config)
  return allM.keys()
def checkLock(config,sysConfig):
  strH5ad = os.path.join(sysConfig['celldepotDir'],"%s.h5ad"%config['prj_name'])
  strLock = "%s.lock"%strH5ad
  if not config['overwrite']:
    if os.path.isfile(strH5ad):
      uName = pwd.getpwuid(os.stat(strH5ad).st_uid).pw_name
      Exit("ERROR: The project %s (owned by %s) does already exist! If you are sure to overwrite, please update config with 'overwrite: True'!"%(config['prj_name'],uName))
    if os.path.isfile(strLock):
      uName = pwd.getpwuid(os.stat(strLock).st_uid).pw_name
      Exit("ERROR: User %s is in the process of project %s. If you are sure to overwrite, please update config with 'overwrite: True'!"%(uName,config['prj_name']))
  if os.path.isfile(strLock):
    os.remove(strLock)
  run_cmd("touch %s"%strLock)
def combine(mIDs,prefix,config):
  print("Evaluating all methods ...")
  CKmethods = [one for one in mIDs if os.path.isfile("%s_%s.h5ad"%(prefix,one))]+['raw']
  cmd="Rscript %s/src/kBET.R %s %s"%(strPipePath,prefix,",".join(CKmethods))
  submit_cmd({"kBET":cmd,
              "silhouette":"%s/src/silhouette.py %s %s"%(strPipePath,prefix,",".join(CKmethods))},
            config)
  #run_cmd("%s/src/silhouette.py %s %s"%(strPipePath,prefix,",".join(CKmethods)))
  if not "SCT" in CKmethods:
    Exit("SCT is missing! and it is required expression for visualization!")
  
  print("combining all methods results ...")
  D = integrateH5ad("%s_SCT.h5ad"%prefix,CKmethods,prefix)
  Draw = integrateH5ad("%s_raw.h5ad"%prefix,CKmethods,prefix)
  print("saving combined results ...")
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    D.write("%s.h5ad"%prefix)
    Draw.write("%s_raw_added.h5ad"%prefix)
    if os.path.isfile("%s.h5ad.lock"%prefix):
      os.remove("%s.h5ad.lock"%prefix)
def integrateH5ad(strH5ad,methods,prefix):
  D = ad.read_h5ad(strH5ad)
  obsm = pd.DataFrame(index=D.obs.index)
  for one in methods:
    if one in ["SCT","raw"]:
      continue
    if not os.path.isfile("%s_%s.h5ad"%(prefix,one)):
      print("Warning: ignore missing h5ad for method %s"%one)
      continue
    D1 = ad.read_h5ad("%s_%s.h5ad"%(prefix,one),backed=True)
    obs = D1.obs.copy()
    addObs = [one for one in obs.columns if not one in D.obs.columns]
    if len(addObs)>0:
      D.obs=D.obs.merge(obs[addObs],how="left",left_index=True,right_index=True)
      
    # check the order of cells
    for k in D1.obsm.keys():
      kname = k.replace("X_","X_%s_"%one)
      obsm1 = obsm.merge(pd.DataFrame(D1.obsm[k],index=D1.obs.index),how="left",left_index=True,right_index=True)
      D.obsm[kname] = obsm1.fillna(0).to_numpy()
  ## remove NA/nan
  for one in D.obs.columns:
    if D.obs[one].isna().sum()>0:
      print("Fix NA: %s"%one)
      if D.obs[one].dtype == 'category':
        D.obs[one] = list(D.obs[one])
        D.obs[one].fillna("NA")
  return D
def moveCellDepot(prefix):
  sysConfig = getConfig()
  shutil.copy("%s.h5ad"%prefix, sysConfig['celldepotDir'])
  print("\nTo check the result, please visit: %s%s.h5ad/"%(sysConfig['celldepotHttp'],os.path.basename(prefix)))
  print("\nAfter confirm the results, please update the publish section of config file before running scAnalyzer with 'publish: True' ")
  print("=== scAnalyzer is completed ===")

def runDEG(strConfig,prefix,config):
  #Rscript /home/zouyang/projects/scRNAsequest/src/scRNAseq_DE.R
  if config['DEG_desp'] is None or not os.path.isfile(config['DEG_desp']):
    return
  D = pd.read_csv(config['DEG_desp'],header=0)
  if D.shape[0]==0:
    return
  cmd = "Rscript %s/src/scRNAseq_DE.R %s"%(strPipePath,strConfig)
  msg = run_cmd(cmd).stdout.decode("utf-8")
  if "scDEG task creation completed" in msg:
    with open("%_scDEG.cmd.json"%prefix,"r") as f:
      scDEGtask = json.load(f)
    submit_cmd(scDEGtask,config,1)

def main():
  global strPipePath
  strPipePath=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
  checkInstallSetting()
  
  if len(sys.argv)<2:
    MsgHelp()
    
  strPath = sys.argv[1]
  if os.path.isdir(strPath):
    initProject(strPath)
  elif os.path.isfile(strPath):
    runPipe(os.path.realpath(strPath))
  else:
    print("The config file is required, and %s doesn't exist!"%strPath)
  
qsubScript='''#!/bin/bash
#$ -N jID_jName
#$ -wd wkPath
#$ -pe node qsubCore
#$ -o jName.log
#$ -e jName.log
#- End UGE embedded arguments
: > $SGE_STDOUT_PATH
cat $PE_HOSTFILE
echo 'end of HOST'

# exit 
set -e

env -i bash -c 'set -o allexport;source sysPath/src/.env;set +o allexport;eval $condaEnv;strCMD'
echo 'DONE'
'''
sbatch='''#!/bin/bash
'''
if __name__ == "__main__":
  start_time = time.time()
  main()
  print("--- total time passed %s seconds ---" % (time.time() - start_time))
