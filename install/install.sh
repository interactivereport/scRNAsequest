#!/usr/bin/env bash

set -e
#Setting environment variables: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#setting-environment-variables
eval $condaEnv
mamba env update -f $(dirname $0)/install_variables.yml
conda deactivate
eval $condaEnv
echo -e "\t*** install packages ***"
mamba env update -f $(dirname $0)/install.yml

#export CBLAS_H_DIR=/edgehpc/apps/gb/anaconda3/4.9.2/include
#export PKG_CONFIG_PATH=/mnt/depts/dept04/compbio/edge_condaEnv/scRNAsequest_seurat5/lib/pkgconfig
R -q -e 'if(!require(peakRAM)) install.packages("peakRAM",repos="https://cran.rstudio.com/")'
R -q -e 'if(!require(revealjs)) install.packages("revealjs",repos="https://cran.rstudio.com/")'
R -q -e 'if(!require(kBET)) devtools::install_github("theislab/kBET",upgrade="never",dependencies=F)'
R -q -e 'if(!require(nebula)) devtools::install_github("lhe17/nebula",ref="v.1.5.3",upgrade="never",dependencies=F)'
R -q -e 'if(!require(SeuratData)) devtools::install_github("satijalab/seurat-data",upgrade="never",dependencies=F)'
R -q -e 'if(!require(Azimuth)) devtools::install_github("satijalab/azimuth",ref="v0.5.0",upgrade="never",dependencies=F)'
R -q -e 'if(!require(RcppPlanc)) devtools::install_github("welch-lab/RcppPlanc",upgrade="never",dependencies=F)' #required by Liger
# for all UMAP: https://github.com/bwlewis/irlba/issues/70
R -q -e 'install.packages("irlba", type = "source",repos="https://cran.rstudio.com/")'
