#!/usr/bin/env bash
appEnvName="scRNAsequest"

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
if { conda env list | grep "^$appEnvName"; } >/dev/null 2>/dev/null; then conda env remove -n $appEnvName; fi
# mamba is not in the base conda=h582c2e5_0_cpython
conda create -y -n $appEnvName "python=3.8.13" "mamba=1.1.0" -c conda-forge
#conda env create -f install.yml
condaPath=$(dirname $(dirname $condaPath))
# setup needed env variables
source $condaPath/etc/profile.d/conda.sh
conda activate $appEnvName
mamba env update -f install/install.yml

echo "export condaEnv='source $condaPath/etc/profile.d/conda.sh;conda activate $appEnvName'" > $src/.env
echo "export PATH=$PATH" >> $src/.env
echo "export OPENBLAS_NUM_THREADS=1" >> $src/.env
echo "export MKL_NUM_THREADS=1" >> $src/.env
echo "export SGE_EXECD_PORT=$SGE_EXECD_PORT" >> $src/.env
echo "export SGE_QMASTER_PORT=$SGE_QMASTER_PORT" >> $src/.env
echo "export SGE_ROOT=$SGE_ROOT" >> $src/.env
echo "export SLURM_CONF=$SLURM_CONF" >> $src/.env
echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH" >> $src/.env
conda deactivate

## additional packages which are not available on anaconda
env -i src="$src" bash -c 'source $src/.env;eval $condaEnv;$src/../install/install.extra'

echo "If no errors above, scRNAsequest installation is successful!"
