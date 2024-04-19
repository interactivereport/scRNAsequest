#!/usr/bin/env bash
appEnvPath="${1:-~/.conda/envs/scRNAsequest}" #"/home/zouyang/.conda/envs/scRNAsequest_test" #"/edgehpc/dept/compbio/edge_condaEnv/scRNAsequest" #

appEnvPath=$(realpath ${appEnvPath/\~/$HOME})
echo -e "\n*** If GPU is avaiable, please run this install on a GPU available node"
echo "*** If SSL certificate (../...crt) needs to be added into this conda env, please export environment variabble 'CONDA_SSL' with the path to the certificate file"
echo -e "\nInstallation location (first position parameter): $appEnvPath"
echo "Wait for 10 seconds, use Ctrl+C to terminate"

for i in {10..1}; do echo -ne "$i\033[0K\r"; sleep 1; done;

set -e
condaPath=$(which conda)
if [[ ${#condaPath} -lt 3 ]]; then
    echo "Missing conda"
    echo "Please install conda and add it into PATH"
    exit
else
    echo "conda in $condaPath"
fi

src="$(dirname $0)/src"
conda env remove -p $appEnvPath
# mamba is not in the base conda
conda create -y -p $appEnvPath python=3.10 mamba -c conda-forge
condaPath=$(dirname $(dirname $condaPath))
# setup needed env variables
if [[ -n "$CONDA_SSL" ]] &&  [[ -f "$CONDA_SSL" ]]; then
    cat $CONDA_SSL >> $appEnvPath/ssl/cacert.pem
fi
source $condaPath/etc/profile.d/conda.sh
conda activate $appEnvPath

mamba env update -f install/install.yml

echo "export condaEnv='source $condaPath/etc/profile.d/conda.sh;conda activate $appEnvPath'" > $src/.env
echo "export PATH=$PATH" >> $src/.env
echo "export OPENBLAS_NUM_THREADS=1" >> $src/.env
echo "export SGE_EXECD_PORT=$SGE_EXECD_PORT" >> $src/.env
echo "export SGE_QMASTER_PORT=$SGE_QMASTER_PORT" >> $src/.env
echo "export SGE_ROOT=$SGE_ROOT" >> $src/.env
echo "export SLURM_CONF=$SLURM_CONF" >> $src/.env
echo "export RETICULATE_PYTHON=$(which python)"
conda deactivate

## additional packages which are not available on anaconda
env -i src="$src" bash -c 'source $src/.env;eval $condaEnv;$src/../install/install.extra'

echo ""
echo "If no errors above, scRNAsequest installation is successful!"
