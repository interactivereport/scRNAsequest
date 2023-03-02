# scRNASequest: an ecosystem of scRNA-seq analysis, visualization and publishing

Tutorial: https://interactivereport.github.io/scRNAsequest/tutorial/docs/index.html

![scRNASequest](https://interactivereport.github.io/scRNAsequest/images/Cover.png?raw=true "scRNASequest")

**Fig. 1. Overview of the scRNASequest workflow.** User provides single-cell or single-nucleus RNA-seq gene expression matrix (h5 or MTX) from Cell Ranger and sample metadata to the semi-automated workflow, scRNASequest. It generates basic quality control (QC) [Bookdown]() reports and allows users to choose from popular data harmonization tools such as [Harmony](https://www.nature.com/articles/s41592-019-0619-0), [LIGER](https://www.cell.com/cell/fulltext/S0092-8674(19)30504-5), and [Seurat](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8) to remove batch effects. [Azimuth](https://www.cell.com/cell/fulltext/S0092-8674(21)00583-3) reference-based cell label transfer is enabled as optional to perform cell type identification. Cluster- or cell-type-specific multi-sample multi-condition single cell level DE analysis is by default performed through [NEBULA](https://www.nature.com/articles/s42003-021-02146-6). Finally, an h5ad file will be generated to be loaded into the [cellxgene VIP](https://www.biorxiv.org/content/10.1101/2020.08.28.270652v2.full) framework or [CellDepot](https://www.sciencedirect.com/science/article/pii/S0022283621006665) single-cell data management system for interactive visualization and analysis.

## 1. Installation

First, please make sure you have [Conda](https://docs.conda.io/en/latest/) installed, or, [Anaconda](https://docs.anaconda.com/anaconda/install/index.html)/[Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed.

```
which conda
# Your conda path will be returned
```

Then, we install scRNASequest by downloading the source code from GitHub:

```
git clone https://github.com/interactivereport/scRNASequest.git
cd scRNASequest

# Install scRNASequest conda environment
# Please make sure you have conda installed before, and this step may take a while
bash install

# The .env will be created under the src directory
ls ~/scRNASequest/src/.env

# Check the path of the current directory and add it to $PATH:
CurrentDir=`pwd`
export PATH="$CurrentDir:$PATH"

# However, the above command only adds the scRNASequest directory to $PATH temporarily
# To add it to the environment permanently, edit ~/.bash_profile or ~/.bashrc:
vim ~/.bash_profile
# Add the full path of the scRNASequest directory to $PATH, for example, $HOME/scRNASequest
PATH=$PATH:$HOME/scRNASequest
# Source the file
source ~/.bash_profile

#To verify the installation, type the main program name, and the manual page will show up:
scAnalyzer
```

## 2. Quick start

This is a quick start guide to the pipeline. Please refer to the [**full tutorial**](https://interactivereport.github.io/scRNAsequest/tutorial/docs/index.html) for more details.

### 2.1 Data preparation

This pipeline accepts [**h5**](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices) or [**MTX**](https://kb.10xgenomics.com/hc/en-us/articles/115000794686-How-is-the-MEX-format-used-for-the-gene-barcode-matrices) (an mtx file and associated barcodes file as well as a features file) containing cell count information [after running Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/overview). You can also include the Cell Ranger QC results: DataPrefix.metrics_summary.csv, but this is optional. When you use downloaded data from NCBI/GEO, it may be necessary to rename the files.

An example of **h5** input file hierarchy, using data from [E-MTAB-11115](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11115):

```
E-MTAB-11115/
    ├── 5705STDY8058280.raw_feature_bc_matrix.h5
    ├── 5705STDY8058280.metrics_summary.csv
    ├── 5705STDY8058281.raw_feature_bc_matrix.h5
    ├── 5705STDY8058281.metrics_summary.csv
    ├── 5705STDY8058282.raw_feature_bc_matrix.h5
    ├── 5705STDY8058282.metrics_summary.csv
    ├── 5705STDY8058283.raw_feature_bc_matrix.h5
    ├── 5705STDY8058283.metrics_summary.csv
    ├── 5705STDY8058284.raw_feature_bc_matrix.h5
    ├── 5705STDY8058284.metrics_summary.csv
    ├── 5705STDY8058285.raw_feature_bc_matrix.h5
    └── 5705STDY8058285.metrics_summary.csv
```

An example of **MTX** input file hierarchy, using data from [GSE185538](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185538). Files under each separate directory must follow the naming criteria: **barcodes.tsv.gz**, **features.tsv.gz**, and **matrix.mtx.gz**:

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

### 2.2 Generate templates of configuration files

Users can initiate the pipeline by running the `scAnalyzer` script with a working directory, where future outputs will be generated.

```
#Example:
scAnalyzer ~/Working_dir
```

The output files include: **config.yml** (a template config file), **DEGinfo.csv** (a template for differential expression analysis), and an empty **sampleMeta.csv** file. The config.yml and sampleMeta.csv are required to run the pipeline, and the DEGinfo.csv is only required for the DE analysis. You can leave DEGinfo.csv empty (by default, just a header line there) currently.

### 2.3 Prepare the sampleMeta.csv file

The empty **sampleMeta.csv** file contains three preset column headers that can be recognized by the pipeline: 

- **Sample_Name** is used to provide a simplified sample name.
- **h5path** is for the full data path. If the inputs are h5 files, this would be the full path of each file. However, for MTX inputs, this should be the path to the data directory, where three data files are stored.
- **metapath** is for cell annotation information, which is optional [(See more details about the annotation file)](https://interactivereport.github.io/scRNAsequest/tutorial/docs/data-preparation.html#public-data-in-h5-format). If you don't have this information, you can delete it from the sampleMeta.csv

Users can add more meta information to this sampleMeta.csv. Here are two examples:

```
#Example 1, for h5 data
Sample_Name,h5path,Sex,Age
5705STDY8058280,~/E-MTAB-11115/5705STDY8058280_filtered_feature_bc_matrix.h5,Female,56d
5705STDY8058281,~/E-MTAB-11115/5705STDY8058281_filtered_feature_bc_matrix.h5,Female,56d
5705STDY8058282,~/E-MTAB-11115/5705STDY8058282_filtered_feature_bc_matrix.h5,Female,56d
5705STDY8058283,~/E-MTAB-11115/5705STDY8058283_filtered_feature_bc_matrix.h5,Male,56d
5705STDY8058284,~/E-MTAB-11115/5705STDY8058284_filtered_feature_bc_matrix.h5,Male,56d
5705STDY8058285,~/E-MTAB-11115/5705STDY8058285_filtered_feature_bc_matrix.h5,Male,56d

#Example 2, for MTX data
Sample_Name,h5path,Treatment,Sex
FCtr,~/GSE185538/GSM5617891_snRNA_FCtr,Control,Female
FEcig,~/GSE185538/GSM5617892_snRNA_FEcig,EcTreated,Female
MCtr,~/GSE185538/GSM5617893_snRNA_MCtr,Control,Male
MEcig,~/GSE185538/GSM5617894_snRNA_MEcig,EcTreated,Male
```

### 2.4 Prepare the config.yml file

This **config.yml** file contains critical configuration parameters to run the pipeline. Please use the following template as an example to prepare this file: [**config.yml**](https://github.com/interactivereport/scRNAsequest/blob/main/src/template.yml).

Some tips:

- **ref_name**: If left blank, the pipeline won't run reference-based cell type annotation. Users can provide an RDS or H5ad file following [Azimuth reference file format](https://github.com/satijalab/azimuth/wiki/Azimuth-Reference-Format) or download from [Azimuth references website](https://azimuth.hubmapconsortium.org/references/). To build a reference, use the `scRef` script from our pipeline [(See more details)](https://interactivereport.github.io/scRNAsequest/tutorial/docs/reference-building.html).
- **gene_group**: You can define your own gene groups to run QC. If the "**rm: False**" is set to False, gene counts will be checked, and cells will be filtered out based on the "cutoff" percentage. To completely eliminate contaminating genes, such as mitochondria genes, set "**rm: True**".
- **runAnalysis: False**: Run analysis means performing data integration and DE analysis. If set to False, the pipeline will only run the initial QC step, which allows users to examine whether the default cell filtering cutoffs are adequate. If they look good, set it to True to run the whole pipeline.
- **overwrite: False**: Set to True if you rerun the pipeline and want to overwrite the previous results.
- **DEG_desp**: Path to the DE configuration file. If the file is empty, it won't perform any analysis. See section 3 about how to fill in this file.

### 2.5 Start the pipeline

Now we have prepared the minimal files (Data files, sampleMeta.csv, and config.yml) to start the pipeline.

- a. Initiate the analysis with `runAnalysis` in config.yml set to False

Here, we pass the path to the config.yml file to run the pipeline:

```
scAnalyzer ~/Working_dir/config.yml
```

This will only run the QC step and generate a [Bookdown report here](https://interactivereport.github.io/scRNAsequest/examples/E-MTAB-11115/bookdown/index.html).

- b. Examine the QC parameters

Check the Bookdown report and adjust the filtering parameters if needed. Repeat this step until the filtering criteria are satisfied.

- c. Run the full analysis

```
#Set config.yml:
#runAnalysis: True
#overwrite: True

#Run the full analysis
scAnalyzer ~/Working_dir/config.yml
```

## 3. Output

After running the above steps, you will see a series of files generated. The main results include:

```
outputdir
    ├── QC/                                  # QC plots
    ├── raw/                                 # Raw h5ad files before and after cell filtering
    ├── Liger/                               # Liger results
    ├── sctHarmony/                          # Harmony results
    ├── SeuratRef/                           # Seurat results
    ├── SeuratRPCA/                          # SeuratRPCA results
    ├── evaluation/                          # kBET and Silhouette plots
    ├── project_name_BookdownReport.tar.gz   # Bookdown report
    ├── project_name.h5ad                    # Final h5ad file including all integration results
    ├── project_name_raw_added.h5ad          # This h5ad also contains raw UMI counts
    ...
```

The `project_name.h5ad` file can be updated to [**Cellxgene VIP**](https://github.com/interactivereport/cellxgene_VIP) platform for visualization (See section 5).

A full description of output files can be seen [here](https://github.com/interactivereport/scRNAsequest/blob/main/src/file.description.yml).

## 4. Differential expression (DE) analysis

The DE analysis in this pipeline is designed to compare **“alt”** and **“ref”** cells from **“group”** within each entry of **“cluster”** considering **“sample”** variations. The **“group”** variable should contain conditions to compare, such as Mutant v.s. Control. Thus, this pipeline is designed to loop through each cluster, and perform DE analysis between **“alt”** v.s. **“ref”**.

We strongly suggest that the reference-based label transfer has been run before DE analysis, so we can run DE based on meaningful cell types, rather than the cluster numbers automatically assigned by the pipeline.

Here, we demo the DE analysis using one real dataset. The screenshot displays the Cellxgene VIP platform:

![DEGinfo](https://interactivereport.github.io/scRNAsequest/tutorial/images/UMAP.png?raw=true "DEGinfo")

In the first example, we would like to run DE analysis between **‘Female’** and **‘Male’** for each cluster annotated by **predicted.celltype1**(Created by reference-based label transfer, setting ref_name in the config.yml file). The first column assigns a name to the current comparison. In the second column, we input the header name **library_id**, which annotates the data sources. Then we add **predicted.celltype1** in the *cluster* column, which allows the pipeline to loop through each cluster in **predicted.celltype1**. The *group* column contains the header name storing the comparison groups, and here we use the **Sex** annotation. Each time, the pipeline can only compare two conditions, such as **‘Female’** and **‘Male’**. If the group column contains more groups, please list them as multiple lines in the DEGinfo.csv file.

We can also add covariates if needed, but this is optional.

Here is the DEGinfo.csv we described above:

```
comparisonName,sample,cluster,group,alt,ref,covars[+ separated],method[default NEBULA],model[default HL]
Compare_Female_vs_Male,library_id,predicted.celltype1,Sex,Female,Male,,NEBULA,HL
```

After preparing this DEGinfo.csv file, we can simply rerun the pipeline (previous steps before DE will be skipped):

```
scAnalyzer ~/Working_dir/config.yml
```

Please refer to the [full tutorial](https://interactivereport.github.io/scRNAsequest/tutorial/docs/differential-expression-de-analysis.html) for more details related to DE analysis.

## 5. Cellxgene VIP visualization

scRNASequest pipeline generates an h5ad file that is fully compatible with Cellxgene VIP for data analysis and visualization. 

We processed an example dataset here:

https://apps.bxgenomics.com/scrnaview/app/core/app_project_launcher.php?ID=422

You can 'Sign in as Guest' and navigate this demo data. First, make sure you select the optimal embedding using the button at the botton of the plot. We recommend **'Liger_umap'** as the embedding for visualization. We can also see that categorical variables are shown on the left, while numeric variables are listed on the right. You can use the drop-shape button to highlight each variable on the plot.

The Cellxgene VIP window is close to the 'Cellxgene' logo on the top left corner. You can maximize the Cellxgene VIP window to use its functions. The *annotation_1* and *annotation_1_print* are the labels provided by the author, see setting [here](https://interactivereport.github.io/scRNAsequest/tutorial/docs/data-preparation.html#public-data-in-h5-format). The *libary_id* indicates data names. The following *liger_cluster* and *liger_cluster_predicted.celltype1* show default cluster assignment and cell type label transfer results. The last item *sex* shows the sex information acquired from the sample meta file.

Please refer to the [GitHub website](https://github.com/interactivereport/cellxgene_VIP) and [Online tutorial](https://interactivereport.github.io/cellxgene_VIP/tutorial/docs/) for more details related to the Cellxgene VIP platform.

## 6. CellDepot data management and publishing

Please refer to this tutorial to manage and publish the project using CellDepot: https://interactivereport.github.io/CellDepot/bookdown/docs/.

