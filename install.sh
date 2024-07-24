#!/usr/bin/env bash
appEnvPath="${1:-~/.conda/envs/scRNAsequest}"

appEnvPath=$(realpath ${appEnvPath/\~/$HOME})

#echo -e "\n*** If GPU is avaiable, please run this install on a GPU available node"
echo -e "\n***CONDA_SSL If SSL certificate (../...crt) needs to be added into this conda env, please export environment variabble 'CONDA_SSL' with the path to the certificate file"
echo "***Warning if the following conda env exists, it will be REMOVED!"
echo -e "\nInstallation location (first position parameter): $appEnvPath"
echo "Wait for 10 seconds, use Ctrl+C to terminate"

for i in {10..1}; do echo -ne "$i\033[0K\r"; sleep 1; done;

appEnvPathCB=$appEnvPath"_CB"

set -e
# get conda path
condaPath=$(which conda)
if [[ ${#condaPath} -lt 3 ]]; then
    echo "Missing conda"
    echo "Please install conda and add it into PATH"
    exit
else
    echo "conda in $condaPath"
fi
src="$(dirname $0)/src"
condaPath=$(dirname $(dirname $condaPath))

# remove old conda env
conda env remove -p $appEnvPath
conda env remove -p $appEnvPathCB

# Create conda env, mamba is not in the base conda
conda create -y -p $appEnvPath python=3.10 mamba -c conda-forge
conda create -y -p $appEnvPathCB python=3.7 -c conda-forge

# set SSL if available
if [[ -n "$CONDA_SSL" ]] &&  [[ -f "$CONDA_SSL" ]]; then
    cat $CONDA_SSL >> $appEnvPath/ssl/cacert.pem
    cat $CONDA_SSL >> $appEnvPathCB/ssl/cacert.pem
fi

# setup needed codna env variables
cleanPath () {
  IFS="$2" read -ra path <<< "$1"
  selPath=""
  for i in "${path[@]}"; do
    [[ $i == "/home"* ]] && continue
    selPath="$selPath:$i"
  done
  echo $(echo ${selPath:1} | sed 's|:*$||')
}
getEnvPrefix () {
  a=($(env))
  b=""
  for i in "${a[@]}"; do
    if [[ $i == "$1"* ]]; then
      b+="  ${i/=/: }\n"
    fi
  done
  echo "$b"
}
getSSL () {
	b=""
	if [[ -f "$1" ]]; then
		b+="  GIT_SSL_CAPATH: $1\n"
		b+="  REQUESTS_CA_BUNDLE: $1\n"
		b+="  SSL_CERT_FILE: $1\n"
		b+="  CURL_CA_BUNDLE: $1\n"
	fi
  echo "$b"
}
cp install/variables.yml install/install_variables.yml
cp install/variables.yml install/CB_variables.yml
echo -e "  PATH: $appEnvPath/bin:$(cleanPath $PATH ':')" >> install/install_variables.yml
echo -e "  PATH: $appEnvPathCB/bin:$(cleanPath $PATH ':')" >> install/CB_variables.yml
echo -en "$(getEnvPrefix 'SGE')" >> install/install_variables.yml
echo -en "$(getEnvPrefix 'SGE')" >> install/CB_variables.yml
echo -en "$(getEnvPrefix 'SLURM')" >> install/install_variables.yml
echo -en "$(getEnvPrefix 'SLURM')" >> install/CB_variables.yml
echo -en "$(getSSL $CONDA_SSL)">> install/install_variables.yml
echo -en "$(getSSL $CONDA_SSL)" >> install/CB_variables.yml
echo -e "  RETICULATE_PYTHON: $appEnvPath/bin/python" >> install/install_variables.yml
echo -e "  RETICULATE_PYTHON: $appEnvPathCB/bin/python" >> install/CB_variables.yml
#echo -e "  LD_LIBRARY_PATH: $appEnvPath/lib:$(cleanPath $LD_LIBRARY_PATH ':')" >> install/install_variables.yml
#echo -e "  LD_LIBRARY_PATH: $appEnvPathCB/lib:$(cleanPath $LD_LIBRARY_PATH ':')" >> install/CB_variables.yml
echo -e "  PKG_CONFIG_PATH: $appEnvPath/lib/pkgconfig" >> install/install_variables.yml
echo -e "  PKG_CONFIG_PATH: $appEnvPathCB/lib/pkgconfig" >> install/CB_variables.yml

# setup .env
echo "export condaEnv='source $condaPath/etc/profile.d/conda.sh;conda activate $appEnvPath'" > $src/.env
echo "export condaEnvCB='source $condaPath/etc/profile.d/conda.sh;conda activate $appEnvPathCB'" >> $src/.env

# install dependencies for conda env
env -i src="$src" bash -c 'source $src/.env;$src/../install/install.sh'
env -i src="$src" bash -c 'source $src/.env;$src/../install/CB.sh'

echo ""
echo "Please add $(dirname $src) into your 'PATH'!"
echo "If no errors above, scRNAsequest installation is successful!"
