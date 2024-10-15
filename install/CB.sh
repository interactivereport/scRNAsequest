#!/usr/bin/env bash

set -e
#Setting environment variables: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#setting-environment-variables
eval $condaEnvCB
conda env update -f $(dirname $0)/CB_variables.yml
conda deactivate
eval $condaEnvCB
pip install tables==3.7.0
#conda install -y -c conda-forge python=3.7 pytables=3.6.1
pip3 install lxml==5.1.0 torch==1.13.1+cu117 torchvision==0.14.1+cu117 torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cu117
#cellbender v0.3.0
# issue with v0.3.1: https://github.com/broadinstitute/CellBender/issues/306
pip3 install -U git+https://github.com/broadinstitute/CellBender@4990df713f296256577c92cab3314daeeca0f3d7



