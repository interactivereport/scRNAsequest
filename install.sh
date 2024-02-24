#!/usr/bin/env bash
# Please update the "appEnvPath" below for the location of the conda env
# if SSL certificate (../...crt) needs to be added into this conda env,
# please specify environment variabble "CONDA_SSL" with the path to the certificate file
appEnvPath="~/.conda/envs/scRNAsequest"

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
# mamba is not in the base conda=h582c2e5_0_cpython
conda create -y -p $appEnvPath "python=3.8.13" "mamba=1.1.0" -c conda-forge
if [[ -n "$CONDA_SSL" ]] &&  [[ -f "$CONDA_SSL" ]]; then
    cat $CONDA_SSL >> $appPATH/ssl/cacert.pem
fi

#conda env create -f install.yml
condaPath=$(dirname $(dirname $condaPath))
# setup needed env variables
source $condaPath/etc/profile.d/conda.sh
conda activate $appEnvPath
mamba env update -f install/install.yml

#avoid user local python env for reticulate
echo "RETICULATE_PYTHON=$appEnvPath/bin/python" >> $appEnvPath/lib/R/etc/Renviron

echo "export condaEnv='source $condaPath/etc/profile.d/conda.sh;conda activate $appEnvPath'" > $src/.env
echo "export PATH=$PATH" >> $src/.env
echo "export PYTHONNOUSERSITE=1" >> $src/.env
echo "export OPENBLAS_NUM_THREADS=1" >> $src/.env
echo "export MKL_NUM_THREADS=1" >> $src/.env
echo "export SGE_EXECD_PORT=$SGE_EXECD_PORT" >> $src/.env
echo "export SGE_QMASTER_PORT=$SGE_QMASTER_PORT" >> $src/.env
echo "export SGE_ROOT=$SGE_ROOT" >> $src/.env
echo "export SLURM_CONF=$SLURM_CONF" >> $src/.env
#echo "export LD_LIBRARY_PATH=$appEnvPath/lib:$LD_LIBRARY_PATH" >> $src/.env
conda deactivate

## additional packages which are not available on anaconda
env -i src="$src" bash -c 'source $src/.env;eval $condaEnv;$src/../install/install.extra'

echo "If no errors above, scRNAsequest installation is successful!"
