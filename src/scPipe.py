import yaml, io, os, sys, subprocess, errno, json, re, logging, warnings, shutil, time, random, pwd, math, configparser, sqlite3, glob, gzip, tracemalloc
import pandas as pd
from datetime import datetime
import scanpy as sc
import anndata as ad
import numpy as np
from scipy.sparse import csc_matrix
import barcodeRankPlot as BRP
import processQC as pQ
import cmdUtility as cU
import mergeH5ad as mH

logging.disable(level=logging.INFO)
warnings.simplefilter(action='ignore', category=FutureWarning)
sc.set_figure_params(vector_friendly=True, dpi_save=300)
strPipePath=""
UMIcol="h5path"
ANNcol="metapath"
IntronExon="intron_exon_count_path"
batchKey="library_id"
CB_expCellNcol="expected_cells"
CB_dropletNcol="droplets_included"
CB_count="low_count_threshold"
CB_learningR="learning_rate"
beRaster=True
tempDir="raw"
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
    print("=====\nPlease set the sys.yml in %s.\n"%strPipePath)
    print("An example is '%s/src/sys_example.yml.\n====="%strPipePath)
    exit()

## msg
def MsgPower():
  sysConfig = getConfig()
  print("\nPowered by %s"%sysConfig['powerby'])
  print("------------")
def MsgHelp():
  sysConfig = getConfig()
  if sysConfig.get("prehelp") is not None:
    print(sysConfig.get("prehelp"))
  #print("The config file will be generated automatically when a DNAnexus download folder is provided")
  strSeuratRef = os.path.join(strPipePath,"seuratdata.csv")
  subprocess.run("R -q -e 'data.table::fwrite(SeuratData::AvailableData(),\"%s\")'"%strSeuratRef,
    shell=True,check=True,stdout=subprocess.PIPE)
  refInfo = pd.read_csv(strSeuratRef).iloc[:,0:7]
  strSysRef = os.path.join(sysConfig['refDir'],"scRNAsequest_ref.csv")
  if os.path.isfile(strSysRef):
    refInfo = pd.concat([refInfo,pd.read_csv(strSysRef)]).reset_index(drop=True)
  print("Available reference data (use 'Dataset' column in config):")
  print(refInfo.iloc[:,0:3].to_string(index=False))
  
  print("\n\nusage: scAnalyzer </path/to/the/config/file> or scAnalyzer </path/to/a/folder>")
  print("\tAn template config file will be created if </path/to/a/folder> is provided")
  Exit()
def MsgInit():
  print("\n\n*****",datetime.now().strftime("%Y-%m-%d %H:%M:%S"),"*****")
  if os.path.isdir(os.path.join(strPipePath,".git")):
    gitConfig = configparser.ConfigParser()
    tmp = gitConfig.read(os.path.join(strPipePath,".git","config"))
    url = gitConfig['remote "origin"']['url']

    gitLog = pd.read_csv(os.path.join(strPipePath,".git","logs","HEAD"),sep="\t",header=None)
    gitLog = gitLog.iloc[-1,0].split(" ")
    print("###########\n## scRNAsequest: %s"%url)
    print("## Pipeline Path: %s"%strPipePath)
    print("## Pipeline Date: %s %s"%(datetime.fromtimestamp(int(gitLog[-2])).strftime('%Y-%m-%d %H:%M:%S'),gitLog[-1]))
    print("## git HEAD: %s\n###########\n"%gitLog[1])

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
  print("Processing sample meta information ...")
  strMeta = "%s/samplesheet.tsv"%strInput
  config = getConfig("template.yml")
  if not os.path.isfile(strMeta):
    print("\tNo samplesheet.tsv! Scanning *filtered_feature_bc_matrix.h5")
    sName = [re.sub(".filtered_feature_bc_matrix.h5","",os.path.basename(one)) for one in glob.glob(os.path.join(strInput,"*filtered_feature_bc_matrix.h5"))]
    if len(sName)==0:
      print("\tSkip meta information: No h5 files found!")
      return None
    meta = pd.DataFrame({config["sample_name"]:sName})
  else:
    meta = pd.read_csv(strMeta,sep="\t")
  sysConfig = getConfig()
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    if not 'notMeta' in sysConfig.keys():
      sysConfig['notMeta'] = []
    meta = meta[[one for one in meta.columns if not one in sysConfig['notMeta']]].dropna(1,'all')
  # if intron/exon count files are available
  strInEx = findIntronExon(list(meta[config["sample_name"]]),strInput)
  if strInEx is not None:
    meta.insert(0,IntronExon,strInEx)
  meta.insert(0,UMIcol,["%s/%s.filtered_feature_bc_matrix.h5"%(strInput,one) for one in meta[config["sample_name"]]])
  for oneH5 in meta[UMIcol]:
    if not os.path.isfile(oneH5):
      Exit("The UMI h5 file %s does not exist, please correct the sample sheet"%oneH5)
  return meta
def initRawMeta(meta):
  print("Create raw h5 sample meta ...")
  metaRaw = meta.copy()
  h5raw = []
  expCellN=[]
  for oneH5 in metaRaw[UMIcol]:
    oneH5raw = re.sub("filtered_feature_bc_matrix","raw_feature_bc_matrix",oneH5)
    h5raw += [oneH5raw if os.path.isfile(oneH5raw) else ""]
    oneMetrics = re.sub("filtered_feature_bc_matrix.h5","metrics_summary.csv",oneH5)
    oneExpN = '0'
    if os.path.isfile(oneMetrics):
      oneMe = pd.read_csv(oneMetrics,header=0)
      oneExpN = oneMe['Estimated Number of Cells'] if 'Estimated Number of Cells' in oneMe.columns else '0'
    expCellN += [str(list(oneExpN)[0])]
  metaRaw[UMIcol] = h5raw
  metaRaw[CB_expCellNcol] = [int(i.replace(',','')) for i in expCellN]
  metaRaw[CB_dropletNcol] = [0]*meta.shape[0]
  metaRaw[CB_count] = [15]*meta.shape[0]
  metaRaw[CB_learningR] = [0.0001]*meta.shape[0]
  return metaRaw[metaRaw[UMIcol].str.len()>0]
def findIntronExon(sNames,strInput):
  sFile = []
  n = 0
  for i in sNames:
    oneF = glob.glob("%s/%s.intron_exon_UMI*"%(strInput,i))
    if len(oneF)>0:
      n +=1
      sFile.append(oneF[0])
    else:
      sFile.append("")
  if n>0:
    return(sFile)
  return(None)
def initSave(meta,strInput,saveRaw=True):
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
  if saveRaw:
    rawMeta = initRawMeta(meta)
    print("*** Total %d samples with %d raw h5"%(meta.shape[0],rawMeta.shape[0]))
    if rawMeta.shape[0]>0:
      rawMeta.to_csv(re.sub(".csv$","_raw.csv",strMeta),index=False)
      BRP.plot(re.sub(".csv$","_raw.csv",strMeta))
    else:
      print("\tSkip raw h5: No raw h5 for any sample!")

  # save empty DE csv table
  strDEG = os.path.join(strOut,"DEGinfo.csv")
  with open(strDEG,"w") as f:
    f.write("comparisonName,sample,cluster,group,alt,ref,covars[+ separated],method[default NEBULA],model[default HL]\n")

  # save config file
  sysConfig = getConfig()
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
                .replace("initDEG",strDEG)
                .replace("initMethods",'[%s]'%','.join(sysConfig['methods'].keys()))
                .replace("initJob","j%d"%random.randint(10,99)) for one in config]
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
    f.write("comparisonName,sample,cluster,group,alt,ref,covars[+ separated],method[default NEBULA],model[default HL]\n")

  sysConfig = getConfig()
  configL = getConfig("template.yml",blines=True)
  configL = [one.replace("initOutput",strInput)
                .replace("initPrjMeta",strMeta)
                .replace("initDEG",strDEG)
                .replace("initPrjName",' #required')
                .replace("initMethods",'[%s]'%','.join(sysConfig['methods'].keys())) for one in configL]

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
    print("\n===> scRMambient %s"%os.path.join(os.path.dirname(strConfig),"sampleMeta_raw.csv"))
    print("Please make sure to update "+CB_dropletNcol+" column in above'sampleMeta_raw.csv'")
    print("\n===> scAnalyzer %s"%strConfig)
    MsgPower()

## pipeline run
def runPipe(strConfig):
  sc.settings.n_jobs=1
  config = getConfig(strConfig,bSys=False)
  config = checkConfig(config)
  meta = getSampleMeta(config["sample_meta"])
  rawCMD = pQ.runQC(config,meta,strConfig)
  if config["runAnalysis"]:
    prefix = os.path.join(config["output"],config["prj_name"])
    methods = runMethods(strConfig,rawCMD)
    scaleF = mH.combine(methods,prefix,config)
    runDEG(strConfig,prefix,config)
    moveCellDepot(prefix,config,scaleF)
    print("=== scAnalyzer is completed ===")
  
  MsgPower()
def setupDir(strOut):
  try:
    os.makedirs(strOut)
  except FileExistsError:
    pass
def checkConfig(config):
  global beRaster
  if config["prj_name"] is None:
    Exit("'prj_name' is required in config file!")
  if config['sample_meta'] is None:
    Exit("'sample_meta' is required in config file!")
  if 'initPrjName' in config["prj_name"]:
    Exit("Please update the 'prj_name' in the config file")
  if 'initPrjTitle' in config["prj_title"]:
    Exit("Please update the 'prj_title' in the config file")
  if not "rasterizeFig" in config.keys():
    config["rasterizeFig"] = True
  beRaster = config["rasterizeFig"]
  return config
def getSampleMeta(strMeta):
  print("processing sample meta information ...")
  if not os.path.isfile(strMeta):
    Exit("Sample meta information (%s) specified in config file does not exist!"%strMeta)
  meta = pd.read_csv(strMeta)
  if meta.shape[0]==0:
    Exit("Sample definition file (%s) is empty"%strMeta)
  if not UMIcol in meta.columns:
    Exit("'%s' columns is required in meta file (%s)"%(UMIcol,strMeta))
  for oneH5 in meta[UMIcol]:
    if os.path.isdir(oneH5):
      continue
    elif os.path.isfile(oneH5) and (oneH5.endswith(".h5") or oneH5.endswith(".h5ad") or oneH5.endswith(".csv") or oneH5.endswith(".tsv")):
      continue
    Exit("The UMI file %s does not exist/not supported (only supports .h5/h5ad/csv/tsv matrix and mtx folder)"%oneH5)
  return(meta)

def runMethods(strConfig,rawCMD):
  print("starting the process by each method ...")
  sysConfig = getConfig()
  config = getConfig(strConfig,bSys=False)
  checkLock(config,sysConfig)

  allM = sysConfig['methods']
  if 'methods' in config.keys():
    selM = list(set(config['methods']) & set(allM.keys()))
    if not 'SCT' in selM:
      Exit("The method 'SCT' including log1p is required but missing in 'methods' setting in config file:\n\t %s"%strConfig)
    allM = {a:allM[a] for a in selM}

  prefix = os.path.join(config["output"],tempDir,config["prj_name"])
  cmds = {}
  if len(rawCMD)>5:
    cmds = {'raw':rawCMD}
  for m in allM.keys():
    strH5ad = "%s.h5ad"%os.path.join(config["output"],m,config["prj_name"])
    if not config["newProcess"] and os.path.isfile(strH5ad):
      print("\tUsing previous %s results: %s\n\t\tPlease rename/remove the above file to rerun!"%(m,strH5ad))
      continue
    setupDir(os.path.dirname(strH5ad))
    cmds[m]="%s/src/%s %s_raw_postfilter.h5ad %s"%(strPipePath,#_%s
                                    allM[m][0],
                                    prefix,
                                    #allM[m][1],
                                    strConfig)
    
  if len(cmds)>0:
    cU.submit_cmd(cmds,config)
  return allM.keys()
def checkLock(config,sysConfig):
  if sysConfig['celldepotDir'] is not None:
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
    cU.run_cmd("touch %s"%strLock)
def rmLock(config):
  sysConfig = getConfig()
  strH5ad = os.path.join(sysConfig['celldepotDir'],"%s.h5ad"%config['prj_name'])
  strLock = "%s.lock"%strH5ad
  if os.path.isfile(strLock):
      os.remove(strLock)
  
def description(strF,strDesc):
  #print(strF)
  try:
    with open(strF,"w") as f:
      f.write(strDesc+"\n")
  except:
    print("Error: fail to write/create %s"%strF)
    return
  
  fstat = os.stat(strF)
  if fstat.st_uid==os.getuid():
    os.chmod(strF,0o664)
  else:
    if bool(re.search('6.$',oct(os.stat(strF).st_mode))):
      print("Warning: no group written permission for ",strF)

def moveCellDepot(prefix,config,scaleF=None):
  sysConfig = getConfig()
  if sysConfig['celldepotDir'] is None:
    print("*** CellDeport is NOT setup ***")
    return()
  strCDfile = os.path.join(sysConfig['celldepotDir'],os.path.basename("%s.h5ad"%prefix))
  strRawfile = os.path.join(sysConfig['celldepotDir'],os.path.basename("%s_raw_obsAdd.h5ad"%prefix))
  if os.path.isfile(strCDfile) and os.access(strCDfile,os.W_OK):
    os.remove(strCDfile)
  if os.path.isfile(strRawfile) and os.access(strRawfile,os.W_OK):
    os.remove(strRawfile)
  if not os.path.isfile(strCDfile):
    shutil.copy("%s.h5ad"%prefix, sysConfig['celldepotDir'])
  else:
    print("Error: cannot copy the file to celldepot folder, permission issue!")
    return()
  if not os.path.isfile(strRawfile):
    shutil.copy("%s_raw_obsAdd.h5ad"%prefix, sysConfig['celldepotDir'])
  rmLock(config)
  
  # create description file
  expScaler = config.get('expScaler')
  if expScaler is None or expScaler==0:
    descTxt = "Description: %s\nData: SCT normalized in 'ln' scale"%config['prj_title']
  else:
    descTxt = "Description: %s\nData: LogNormalize in 'ln' scale"%config['prj_title']
  if scaleF is not None:
    descTxt += "\nScale Factor:%d"%round(scaleF)
  description("%s/%s.txt"%(sysConfig['celldepotDir'],os.path.basename(prefix)),descTxt)

  description("%s/%s_raw_added.txt"%(sysConfig['celldepotDir'],os.path.basename(prefix)),
    "Description: %s\nData: UMI"%config['prj_title'])

  if os.path.isfile("%s.db"%prefix):
    shutil.copy("%s.db"%prefix, sysConfig['celldepotDir'])
    print("scDEG is available in VIP")
  #print("\nTo check the result, please visit: %s%s.h5ad/"%(sysConfig['celldepotHttp'],os.path.basename(prefix)))
  print("\nAfter confirm the results, please use 'Create Project' on CellDepot to add this project with 'File Name' of %s.h5ad."%os.path.basename(prefix))
  

def runDEG(strConfig,prefix,config):
  if os.path.isfile("%s.db"%prefix) and not config["newProcess"]:
    print("Skip scDEG, db file exists: %s.db"%prefix)
    return
  if config['DEG_desp'] is None or not os.path.isfile(config['DEG_desp']):
    print("Skip scDEG: Missing DEG description file!")
    return
  D = pd.read_csv(config['DEG_desp'],header=0)
  if D.shape[0]==0:
    return
  
  cmd = "Rscript %s/src/scRNAseq_DE.R %s"%(strPipePath,strConfig)
  msg = cU.run_cmd(cmd).stdout.decode("utf-8")
  #msg="scDEG task creation completed"
  if "scDEG task creation completed" in msg:
    with open("%s_scDEG.cmd.json"%prefix,"r") as f:
      scDEGtask = json.load(f)
    umiF = config.get('UMI')
    if not config.get('memory') is None:
      memG=int(re.sub("G$","",config.get('memory')))
    elif umiF is None:
      memG = math.ceil(os.path.getsize("%s_raw_added.h5ad"%prefix)*50/1e9)
    else:
      memG = math.ceil(os.path.getsize(umiF)*50/1e9)
    cU.submit_cmd(scDEGtask,config,math.ceil(memG/16),memG)
    formatDEG(prefix)
def formatDEG(prefix):
  print("=== Formating scDEG results to create the db file ===")
  with open("%s_scDEG.cmd.json"%prefix,"r") as f:
      DEGcmds = json.load(f)
  DEGpaths = list(set([one.split(";")[0].replace("cd ","") for k,one in DEGcmds.items()]))
  csv = []
  for onePath in DEGpaths:
    for f in os.listdir(onePath):
      strCSV = os.path.join(onePath,f)
      if not os.path.isfile(strCSV) or not f.endswith("csv"):
        continue
      print("\tprocessing: ",f)
      tab = pd.read_csv(strCSV).iloc[:,0:4]
      tab.columns = ["gene","log2fc","pval","qval"]
      tags = f[:-4].split("__")
      tab["contrast"] = tags[0]
      tab["tags"] = ";".join(tags[1:]+[os.path.basename(onePath)])
      csv += [tab]
  if len(csv)==0:
    return
  data = pd.concat(csv,ignore_index=True)
  data = data.dropna()
  D = data[["log2fc","pval","qval"]]
  D.index = pd.MultiIndex.from_frame(data[["gene","contrast","tags"]])
  conn = sqlite3.connect('%s.db'%prefix)
  D.to_sql("DEG",conn,if_exists="replace")
  conn.close()

def main():
  global strPipePath
  strPipePath=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
  checkInstallSetting()
  MsgInit()

  if len(sys.argv)<2:
    MsgHelp()

  strPath = sys.argv[1]
  if os.path.isdir(strPath):
    initProject(strPath)
  elif os.path.isfile(strPath):
    runPipe(os.path.realpath(strPath))
  else:
    print("The config file is required, and %s doesn't exist!"%strPath)

if __name__ == "__main__":
  start_time = time.time()
  main()
  print("--- total time passed %s seconds ---" % (time.time() - start_time))
