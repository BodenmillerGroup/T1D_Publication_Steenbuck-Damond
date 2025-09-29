# T1D_Publication_Steenbuck-Damond
This repository contains all code scripts for the type 1 diabetes study from Steenbuck, Damond et al., 2025.  Read the preprint here: https://www.biorxiv.org/content/10.1101/2025.03.05.641526v1

# Abstract
The natural history and pathogenesis of type 1 diabetes, particularly during the autoantibody-positive stages preceding clinical onset, are not well understood, in part, due to limited availability of human pancreatic samples. Here, we studied 88 pancreas samples of organ donors, including 28 single autoantibody-positive and 10 multiple autoantibody-positive donors, by imaging mass cytometry. We imaged over 10,000 islets and 16 million single-cells using 79 antibodies revealing both β-cell states and the islet-immune interface. Analysis of the adaptive and innate immune cell states identified a novel pro-inflammatory macrophage – exhausted-like T-cell axis as characteristic for early disease stages, and as a central component of (peri-)insulitic lesions, indicating a key role of macrophages in type 1 diabetes development. In beta cells, alterations in Interferon signatures and downregulation across lineage and functional markers, including markers of endoplasmic reticulum stress, were characteristic of recent-onset disease. We also identified IAPP loss, prior to Insulin loss, from β-cells as a major indicator of pre-clinical disease. Multiple immune cell subtypes were associated with young age and insulitis, thus potentially modulating the higher disease severity observed in younger type 1 diabetics. Together, our results provide understanding of earliest stages of type 1 diabetes progression and reveals multiple novel and potentially clinically actionable disease features when most beta cells are still alive.  

# Repository structure
This repository contains data analysis scripts for re-creating all results and manuscript figures based on provided data objects from zenodo (https://zenodo.org/records/14968076). This repository contains the `/T1D_processing` and `/T1D_analysis` directories. For reproducing the analysis: `/T1D_processing/processing` and the `/T1D_analysis` are the most important directories. All directories contain extensive .README files giving further information.

The `/T1D_processing/` directory contains all **pre-processing** scripts. This includes subdirectories containing scripts for:  
1. `/titration`: Antibody titration with IMC.
2. `/region_selection`: ImageJ and Python scripts for ROI selection and registration across the Immune and Islet panel.
3. `/islet_segmentation`: Code and files for islet segmentation of the Pilot Study. Islet masks are used downstream as training set. 
4. `/samples`: Information about donors imaged in study. Used for variance minimization across batches.
5. `/processing`: Main folder for pre-processing. Includes scripts for **cell segmentation**, **islet segmentation** and **feature extraction** from the **final study.** 
6. `/ext`: contains misc. required files for pre-processing.

The `/T1D_analysis/` directory contains all **analysis** scripts. This includes all cell type annotation and downstream analysis scripts (e.g. spatial + differential abundance analysis).

The `/Snake_pipe/` directory contains Snakefiles and configuration files to reproducibly run the pre-processing (`/T1D_processing/processing/`) including cell type annotation (01-08 scripts in `/T1D_analysis/`). 
Please note: we provide the full SCE_object generated after this step on Zenodo.

# Software note
Most scripts are written in the statistical programming language R and Python. 
Dependencies for pre-processing and running 01-08 script of T1D_analysis are at dockerhub (**https://hub.docker.com/r/nathanste/t1d_analysis**).

For the further downstream analysis, all **R** code was run with R version 4.5.1 on *Ubuntu 24.04*. Each .Rmd file contains a Software requirements section that lists all relevant packages for data analysis and processing. File paths have to be adjusted by the user.

# Data access

Raw and processed data is provided on zenodo (IMC: https://zenodo.org/records/14968076). Please follow the respective linked directories containing the raw image data.

NOTE: We plan to add the processed dataset to the imcdatasets R/Bioconductor package (https://github.com/BodenmillerGroup/imcdatasets), which will facilitate easy community access and direct integration into computational workflows.