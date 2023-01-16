# scRNASequest: an ecosystem of scRNA-seq analysis, visualization and publishing

Tutorial: https://interactivereport.github.io/scRNAsequest/tutorial/docs/index.html

![scRNASequest](https://interactivereport.github.io/scRNAsequest/images/Cover.png?raw=true "scRNASequest")

**Overview of the RNASequest workflow.** User provides gene expression matrix (h5 or MTX) from Cell Ranger and sample metadata to the semi-automated workflow, scRNASequest. It generates basic quality control (QC) reports and allows users to choose from popular data harmonization tools such as [Harmony](https://www.nature.com/articles/s41592-019-0619-0), [LIGER](https://www.cell.com/cell/fulltext/S0092-8674(19)30504-5), and [Seurat](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8) to remove batch effects. [Azimuth](https://www.cell.com/cell/fulltext/S0092-8674(21)00583-3) reference-based cell label transfer is enabled as optional to per-form cell type identification. Cluster- or cell-type-specific multi-sample multi-condition single cell level DE analysis is by default performed through [NEBULA](https://www.nature.com/articles/s42003-021-02146-6). Finally, an h5ad file will be generated to be loaded into the [cellxgene VIP](https://www.biorxiv.org/content/10.1101/2020.08.28.270652v2.full) framework or [CellDepot](https://www.sciencedirect.com/science/article/pii/S0022283621006665) single-cell data management system for interactive visualization and analysis.

## 1. Installation

First we install scRNASequest by downloading the scripts from GitHub:

```
git clone https://github.com/interactivereport/scRNASequest.git
cd scRNASequest

# Install scRNASequest conda environment
# Please make sure you have conda installed before, and this step may take a while
bash install

# The .env will be created under the src directory
ls ~/scRNASequest/src/.env

# Check the path of current directory and add it to $PATH:
CurrentDir=`pwd`
export PATH="$CurrentDir:$PATH"

# However, the above command only adds the RNASequest directory to $PATH temporarily
# To add it to the environment permanently, edit ~/.bash_profile or ~/.bashrc:
vim ~/.bash_profile
# Add the full path of the RNASequest directory to $PATH, for example, $HOME/scRNASequest
PATH=$PATH:$HOME/scRNASequest
# Source the file
source ~/.bash_profile

#To verify the installation, typing the main program name, and the manual page will show up:
scAnalyzer
```

## 2. Quick start

This is a quick start guide of the pipeline. Please refer to the [**full tutorial**](https://interactivereport.github.io/scRNAsequest/tutorial/docs/index.html) for more details.

### 2.1 File preparation

This pipeline accepts **.h5** or **MTX** (an mtx file and associated barcodes file as well as a features file) containing cell count information after running Cell Ranger. You can also include the Cell Ranger QC results: metrics_summary.csv, but this is optional. When downloading data from NCBI/GEO, it may be necessary to rename the files.

An example of .h5 input file hierarchy, using data from [E-MTAB-11115](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11115):

```
5705STDY8058280.raw_feature_bc_matrix.h5
5705STDY8058280.metrics_summary.csv
5705STDY8058281.raw_feature_bc_matrix.h5
5705STDY8058281.metrics_summary.csv
5705STDY8058282.raw_feature_bc_matrix.h5
5705STDY8058282.metrics_summary.csv
5705STDY8058283.raw_feature_bc_matrix.h5
5705STDY8058283.metrics_summary.csv
5705STDY8058284.raw_feature_bc_matrix.h5
5705STDY8058284.metrics_summary.csv
5705STDY8058285.raw_feature_bc_matrix.h5
5705STDY8058285.metrics_summary.csv
```

An example of MTX input file hierarchy, using data from [GSE185538](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185538). Files under each separate directory must follow the naming criteria: barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz:

```
GSE185538/
    ├── GSM5617891_snRNA_FCtr
        ├── barcodes.tsv.gz
        ├── features.tsv.gz
        └── matrix.mtx.gz
    ├── GSM5617892_snRNA_FEcig
        ├── barcodes.tsv.gz
        ├── features.tsv.gz
        └── matrix.mtx.gz
    ├── GSM5617893_snRNA_MCtr
        ├── barcodes.tsv.gz
        ├── features.tsv.gz
        └── matrix.mtx.gz
    └── GSM5617894_snRNA_MEcig
        ├── barcodes.tsv.gz
        ├── features.tsv.gz
        └── matrix.mtx.gz
```

### 2.2 Config file



## Output

## Cellxgene VIP visualization
