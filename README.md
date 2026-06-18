# T1D_Publication_Steenbuck-Damond
This repository contains all code scripts for the type 1 diabetes study from Steenbuck, Damond et al., 2025.  Read the preprint here: https://www.biorxiv.org/content/10.1101/2025.03.05.641526v1.

# Abstract
The pathogenesis of type 1 diabetes, particularly at autoantibody-positive preclinical stages, remains poorly understood, in part due to limited sample availability. Here, we show imaging mass cytometry data of pancreas samples from 88 organ donors, including 28 single and 10 multiple autoantibody-positive donors. We imaged 16 million single-cells using 79 antibodies to characterize β-cell states and the islet-immune interface, correcting for relevant covariates. We identified interactions between pro-inflammatory macrophages and exhausted-like (PD1+TIM3+) T-cells. These interactions were characteristic of early disease and of insulitic islets, indicating a key role of macrophages in disease development. β-cells showed loss of IAPP in pre-clinical disease, insulitic Interferon signatures, and no increase in three measured ER stress markers in disease samples relative to control. Multiple immune cell subtypes were associated with young age and insulitis, potentially contributing to greater disease severity in younger patients. Our data present a resource describing early type 1 diabetes progression and reveal potentially clinically actionable features before extensive β-cells loss. 

# Repository structure
This repository contains data analysis scripts for re-creating all results and manuscript figures based on provided data objects from zenodo (https://zenodo.org/records/14968076). This repository contains the `/T1D_processing` and `/T1D_analysis` directories. For reproducing the analysis: `/T1D_processing/processing` and the `/T1D_analysis` are the most important directories. All directories contain extensive .README files giving further information. See below for details.

# Data access

Raw and processed data is provided on zenodo (IMC: https://zenodo.org/records/14968076). 
Please follow the respective linked directories containing the raw image data.
An informative subset of the data is also available on Pancreatlas (https://www.pancreatlas.org/) as a browsable set of annotated images. 
The processed dataset is currently being uploaded to the *imcdatasets* R/Bioconductor package (https://github.com/BodenmillerGroup/imcdatasets). 
We are currently adding a vignette to *imcdatasets* showcasing how 

# How can you use this dataset in your work?

## T1D scientist:
The dataset contains imaging data of 79 unique antibodies across 88 donors (incl. 28 sAAb+ and 10 mAAb+) with rich metadata and image annotations. 
It therefore presents a rich resource for biological discovery and for validating your targets of interest.
For example, you aim to knock-out a given target or found a potentially interesting transcript of a protein included in this study?
Use this dataset, to check its protein levels across biological contexts (disease stages, cell types, insulitis, ICIs vs IDIs).

## Bioinformatician:
The dataset is highly applicable to benchmarking studies or building methods. Tasks could include:
- Cell type annotation & subclustering.
- Cell type domains.
- Methods to extract biological insights from ordinal cross-sectional spatial atlases
- Object detection and extraction of biologically informative geometric features (based on islets)
- Joint Multi-modal analysis of spatial data acquired on two sections, for example in matching similar cells.
- Scalable spatial methods as the dataset is large (16M cells; > 7,000 paired images (=14,050 IMC images).

-> Check the vignette to be added to imcdatasets; it will cover details of key columns (e.g. cell-type, cell-id, spatial domain etc.).
-> For any ideas or queries, feel reach out to me (nathan.steenbuck@uzh.ch) or open an Issue.

# Repository details:
Pre-processing
The `/T1D_processing/` directory contains all **pre-processing** scripts. This includes subdirectories containing scripts for:  
1. `/titration`: Antibody titration with IMC.
2. `/region_selection`: ImageJ and Python scripts for ROI selection and registration across the Immune and Islet panel.
3. `/islet_segmentation`: Code and files for islet segmentation of the Pilot Study. Islet masks are used downstream as training set. 
4. `/samples`: Information about donors imaged in study. Used for variance minimization across batches.
5. `/processing`: Main folder for pre-processing. Includes scripts for **cell segmentation**, **islet segmentation** and **feature extraction** from the **final study.** 
6. `/ext`: contains misc. required files for pre-processing.

## Analysis
The `/T1D_analysis/` directory contains all **analysis** scripts. This includes all cell type annotation and downstream analysis scripts (e.g. spatial + differential abundance analysis).

The `/Snake_pipe/` directory contains Snakefiles and configuration files to reproducibly run the pre-processing (`/T1D_processing/processing/`) including cell type annotation (01-08 scripts in `/T1D_analysis/`). 
Please note: we provide the full SCE_object generated after this step on Zenodo.

# Software note
Most scripts are written in the statistical programming language R and Python. 
Dependencies for pre-processing and running 01-08 script of T1D_analysis are at dockerhub (**https://hub.docker.com/r/nathanste/t1d_analysis**).

For the further downstream analysis, all **R** code was run with R version 4.5.1 on *Ubuntu 24.04*. Each .Rmd file contains a Software requirements section that lists all relevant packages for data analysis and processing. File paths have to be adjusted by the user.
