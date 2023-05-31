import os, subprocess, time, re, random, glob
from natsort import natsorted
import pandas as pd

strPipePath=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
maxJobSubmitRepN = 2

def run_cmd(cmd):
  #print(cmd)
  try:
    cmdR = subprocess.run(cmd,shell=True,check=True,stdout=subprocess.PIPE)#capture_output=True,
  except subprocess.CalledProcessError as e:
    if not e.returncode==1:
      print("%s process return above error"%cmd)#,e.stderr.decode("utf-8")
    cmdR = e
  return cmdR
## parallel job management
def submit_cmd(cmds,config,core=None,memG=0):
  #cmds = {k:v for i, (k,v) in enumerate(cmds.items()) if not v is None}
  if len(cmds)==0:
    return
  if core is None:
    core=config['core']
  parallel = config["parallel"]
  if not parallel:
    os.makedirs(os.path.join(config["output"],"log"),exist_ok=True)
    for one in cmds.keys():
      if cmds[one] is None:
        continue
      print("\n\n\nsubmitting %s"%one)
      oneCMD=cmds[one] #+" 2>&1 | tee "+ strLog
      try:
        subprocess.run(oneCMD,shell=True,check=True)
      except:
        print("%s process return error!"%one)
  elif parallel=="sge":
    jID = qsub(cmds,config['output'],core,memG=memG)
    print("----- Monitoring all submitted SGE jobs: %s ..."%jID)
    ## in case of long waiting time to avoid Recursion (too deep)
    cmdN = {one:1 for one in cmds.keys()}
    failedJobs = {}
    while True:
      qstat(jID,config['output'],cmds,cmdN,core,memG,failedJobs)
      if len(cmds)==0:
        break
      ## wait for 1 min if any running jobs
      time.sleep(60)
    if len(failedJobs)>0:
      print("\n\n*** The following jobs failed:")
      for k,one in failedJobs.items():
        print("\t--->ERROR failed with %d times qsub: %s"%(maxJobSubmitRepN,one))
      print("\n*** Please check log files in %s folder and consider rerun the analysis!"%jID)
  elif parallel=="slurm":
    if memG==0 and config.get('memory') is not None and config.get('memory').endswith("G"):
      print("\tAssign memory according to config: ",config.get('memory'))
      memG=int(re.sub("G$","",config.get('memory')))
    jID = sbatch(cmds,config['output'],core,memG=memG,jID=config.get("jobID"),gpu=config.get('gpu'))
    print("\t----- Monitoring all submitted SLURM jobs: %s ..."%jID)
    ## in case of long waiting time to avoid Recursion (too deep)
    cmdN = {one:1 for one in cmds.keys()}
    failedJobs = {}
    while True:
      squeue(jID,config['output'],cmds,cmdN,core,memG,failedJobs,gpu=config.get('gpu'))
      if len(cmds)==0:
        break
      ## wait for 1 min if any running jobs
      time.sleep(60)
    if len(failedJobs)>0:
      print("\n\n*** The following jobs failed:")
      for k,one in failedJobs.items():
        print("\t--->ERROR failed with %d times qsub: %s"%(maxJobSubmitRepN,one))
      print("\n*** Please check log files in %s folder and consider rerun the analysis!"%jID)
  else:
    print("ERROR: unknown parallel setting: %s"%parallel)
    exit()

def qsub(cmds,strPath,core,memG=0,jID=None):
  if jID is None:
    jID = "j%d"%random.randint(10,99)
  strWD = os.path.join(strPath,jID)
  try:
    os.makedirs(strWD)
  except FileExistsError:
    pass
  noRun = []
  if memG==0:
    qsubTmp = "\n".join([i for i in qsubScript.split("\n") if not "MEMFREE" in i])
  else:
    qsubTmp = qsubScript.replace('MEMFREE',str(memG))
  for one in cmds.keys():
    oneScript = (qsubTmp.replace('qsubCore',str(core))
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
def qstat(jID,strPath,cmds,cmdN,core,memG,failedJobs):
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
          #print("\n--->ERROR failed with %d times qsub: %s"%(maxJobSubmitRepN,one))
          failedJobs[len(failedJobs)+1] = one
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
    re1=qsub(resub,strPath,core,memG,jID)
    time.sleep(5) #might not needed for the qsub to get in
  for one in finishedJob:
    del cmds[one]
  ### in case of long waiting time to avoid Recursion (too deep)
  #if len(cmds)>0:
  #  ## wait for 1 min if any running jobs
  #  time.sleep(60)
  #  qstat(jID,strPath,cmds,core,iN)

def sbatch(cmds,strPath,core,memG=0,jID=None,gpu=False):
  gpu=False if gpu is None else gpu
  sbatchScript=sbatchScriptCPU
  if gpu:
    sbatchScript=sbatchScriptGPU
  if jID is None:
    jID = "j%d"%random.randint(10,99)
  strWD = os.path.join(strPath,jID)
  os.makedirs(strWD,exist_ok=True)
  if memG==0:
    pScript = "\n".join([i for i in sbatchScript.split("\n") if not "MEMFREE" in i])
  else:
    pScript = sbatchScript.replace('MEMFREE',str(memG))
  for one in cmds.keys():
    oneScript = (pScript.replace('CoreNum',str(core))
                            .replace('jName',one)
                            .replace('wkPath',strWD)
                            .replace("jID",jID)
                            .replace("sysPath",strPipePath)
                            .replace("strCMD",cmds[one]))
    strF=os.path.join(strWD,"%s.sh"%one)
    with open(strF,"w") as f:
      f.write(oneScript)
    run_cmd("sbatch %s"%strF)
  return jID
def squeue(jID,strPath,cmds,cmdN,core,memG,failedJobs,gpu=False):
  gpu=False if gpu is None else gpu
  print(".",end="")
  strWD = os.path.join(strPath,jID)
  qstateCol=4 # make sure the state of qstat is the 4th column (0-index)
  qJobCol=0 # make sure the job-id of qstat is the 0th column (0-index)
  # remove the job in wrong stats: Eqw and obtain the running and waiting job names
  qs = run_cmd("squeue")
  if qs.returncode!=0:
    return
  qs = pd.DataFrame([[a.strip() for a in one.split(" ") if len(a.strip())>0 ] for one in qs.stdout.decode("utf-8").split("\n") if len(one)>6 and '%s_'%jID in one])
  errID=[]
  runID=[]
  if qs.shape[0]>0:
    qs[0] = [a.split("_")[0] for a in qs[0]]
    errID = set(qs[qs[qstateCol].isin(['S','ST'])][qJobCol])
    runID = set(qs[~qs[qstateCol].isin(['S','ST'])][qJobCol])

  for one in errID:
    run_cmd("scancel %s"%one)
    runID.discard(one)
  jNames = []
  for oneJob in runID:
    if len(oneJob)<5:
      continue
    jobDetail = parseSlurmJob(oneJob)
    jName=jobDetail.get('JobName')
    if jName is None:
      continue
    jNames.append(jName.replace(jID+"_",""))
  # check the finished job, resubmit if not sucessfully finished
  jNames=list(set(jNames))
  resub = {}
  finishedJob = []
  for one in cmds.keys():
    if not one in jNames:
      strLog = glob.glob(os.path.join(strWD,"%s.log"%one))
      if len(strLog)==0 or not "DONE" in run_cmd("tail -n 1 %s"%natsorted(strLog,reverse=True)[0]).stdout.decode("utf-8"):
        cmdN[one] += 1
        if cmdN[one]>maxJobSubmitRepN:
          #print("\n--->ERROR failed with %d times qsub: %s"%(maxJobSubmitRepN,one))
          failedJobs[len(failedJobs)+1] = one
          finishedJob.append(one)
          continue
        resub[one] = cmds[one]
      else:
        print("\n\t===== %s ====="%one)
        with open(strLog[0],"r") as f:
          print(f.read())
        print("\tFinished: %s"%one)
        finishedJob.append(one)
  if len(resub)>0:
    re1=sbatch(resub,strPath,core,memG,jID,gpu=gpu)
    time.sleep(5) #might not needed for the qsub to get in
  for one in finishedJob:
    del cmds[one]
def parseSlurmJob(jobID):
  jInfo = {}
  for one in (run_cmd("scontrol show job %s"%jobID)
              .stdout
              .decode("utf-8")
              .split("\n")):
    tmp = one.split("=")
    if len(tmp)<2:
      continue
    if len(tmp)==2:
      jInfo[tmp[0].strip()]=tmp[1].strip()
      continue
    k=tmp[0].strip()
    if k in jInfo.keys():
      break
    v=""
    for i in range(1,len(tmp)-1):
      vk=tmp[i].rsplit(" ",1)
      if len(vk)<2:
        v += "="+vk[0]
        continue
      if len(v)>0:
        vk[0]= v.strip("=")+"="+vk[0]
      jInfo[k] = vk[0]
      k=vk[1].strip()
    if len(v)>0:
      tmp[-1] = v.strip("=")+"="+tmp[-1]
    jInfo[k]=tmp[-1].strip()
  return(jInfo)

qsubScript='''#!/bin/bash
#$ -N jID_jName
#$ -wd wkPath
#$ -pe node qsubCore
#$ -l m_mem_free=MEMFREEG
#$ -o jName.log
#$ -e jName.log
#- End UGE embedded arguments
: > $SGE_STDOUT_PATH
cat $PE_HOSTFILE
echo 'end of HOST'

# exit
set -e

env -i bash -c 'source sysPath/src/.env;eval $condaEnv;strCMD'
echo 'DONE'
'''
sbatchScriptCPU='''#!/bin/bash
#SBATCH -J jID_jName
#SBATCH -D wkPath
#SBATCH -n CoreNum
#SBATCH -t 72:0:0
#SBATCH --mem=MEMFREEG
#SBATCH -o jName.log
#SBATCH -e jName.log
#- End embedded arguments
echo $SLURM_JOB_NODELIST
echo 'end of HOST'
# exit
set -e
env -i bash -c 'source sysPath/src/.env;eval $condaEnv;strCMD'
echo 'DONE'
'''
sbatchScriptGPU='''#!/bin/bash
#SBATCH --job-name=jID_jName
#SBATCH -D wkPath
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --requeue
#SBATCH -p gpu
#SBATCH --mem=MEMFREE
#SBATCH -o jName.log
#SBATCH -e jName.log
#- End embedded arguments
echo $SLURM_JOB_NODELIST
echo 'end of HOST'
# exit
set -e
env -i bash -c 'source sysPath/src/.env;eval $condaEnv;strCMD'
echo 'DONE'
'''
