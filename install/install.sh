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
# A failed R install only prints a warning, so each line ends with library(): if the package can't be loaded,
# R exits with an error and `set -e` stops the installation instead of finishing with the package missing.
R -q -e 'if(!require(peakRAM)) install.packages("peakRAM",repos="https://cran.rstudio.com/"); library(peakRAM)'
R -q -e 'if(!require(revealjs)) install.packages("revealjs",repos="https://cran.rstudio.com/"); library(revealjs)'
R -q -e 'if(!require(kBET)) devtools::install_github("theislab/kBET",upgrade="never",dependencies=F); library(kBET)'
R -q -e 'if(!require(nebula)) devtools::install_github("lhe17/nebula",ref="v.1.5.3",upgrade="never",dependencies=F); library(nebula)'
R -q -e 'if(!require(SeuratData)) devtools::install_github("satijalab/seurat-data",upgrade="never",dependencies=F); library(SeuratData)'
R -q -e 'if(!require(Azimuth)) devtools::install_github("satijalab/azimuth",ref="v0.5.0",upgrade="never",dependencies=F); library(Azimuth)'
# RcppPlanc builds a bundled hwloc with its own ./configure, which looks for gcc/cc on PATH; the conda env only has the
# prefixed compiler, so pass R's compiler explicitly (otherwise this relies on a system gcc, which Docker images lack)
CC="$(R CMD config CC)" R -q -e 'if(!require(RcppPlanc)) devtools::install_github("welch-lab/RcppPlanc",ref="e1989e675711c70147c10f1f7f10f21bdee3b026",upgrade="never",dependencies=F); library(RcppPlanc)' #required by Liger
# for all UMAP: https://github.com/bwlewis/irlba/issues/70
# irlba is already installed (conda binary) and is reinstalled from source, so library() can't detect a failure; treat the failed-install warning as an error instead.
R -q -e 'withCallingHandlers(install.packages("irlba", type = "source",repos="https://cran.rstudio.com/"), warning=function(w) if(grepl("non-zero exit status",conditionMessage(w))) stop(w)); library(irlba)'
