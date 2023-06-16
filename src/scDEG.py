import scPipe, os, time, sys

def MsgHelp():
  print("\nscDEG /path/to/a/folder === or === scDEG /path/to/a/config/file\n")
  print("An empty config file will be generated automatically when a folder is provided")
  scPipe.Exit()
  
def initProject(strInput):
  strInput=os.path.realpath(strInput)
  print("*****\nCreating an empty config file and a empty comparison definition file in ",strInput)
  print("\t Please Fill in the required information in the config and comparison files")

  # save empty DE csv table
  strDEG = os.path.join(strInput,"DEGinfo.csv")
  with open(strDEG,"w") as f:
    f.write("comparisonName,sample,cluster,group,alt,ref,covars[+ separated],method[default NEBULA],model[default HL]\n")
  
  configL = scPipe.getConfig("template.yml",blines=True)
  pos = [i for i in range(len(configL)) if configL[i].startswith("## DEG analysis") or configL[i].startswith("## publish")]
  DEGconfig = configL[pos[0]:pos[1]]
  DEGconfig = [one.replace("initDEG",strDEG) for one in DEGconfig]
  DEGconfig = ["UMI: #required, can be a matrix rds or a h5ad file\n",
    "meta: #required, can be a cell annotation data.frame rds or a h5ad file\n",
      "output: "+strInput+"\n",
      "DBname: cellxgeneVIP\n",
      'parallel: slurm # False or "sge" or "slurm"\n',
      'memory: 16G\n',
      'newProcess: False #False use existing DEG results']+DEGconfig

  strConfig = os.path.join(strInput,"config.yml")
  with open(strConfig,"w") as f:
    f.writelines(DEGconfig)
  
  scPipe.Exit()

def runDEG(strConfig):
  config = scPipe.getConfig(strConfig,bSys=False)
  prefix = config['output']+"/"+config['DBname']
  scPipe.runDEG(strConfig,prefix,config)

def main():
  scPipe.strPipePath=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
  if len(sys.argv)<2:
    MsgHelp()
  
  strPath = sys.argv[1]
  if os.path.isdir(strPath):
    initProject(strPath)
  elif os.path.isfile(strPath):
    runDEG(os.path.realpath(strPath))
  else:
    print("The config file is required, and %s doesn't exist!"%strPath)


if __name__ == "__main__":
  start_time = time.time()
  main()
  print("--- total time passed %s seconds ---" % (time.time() - start_time))
