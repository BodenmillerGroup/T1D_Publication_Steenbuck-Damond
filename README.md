# T1D_Publication_Steenbuck-Damond
This repository contains all code scripts for the type 1 diabetes study from Steenbuck, Damond et al., 2026 published in Nature Metabolism.  
Read the study here: https://www.nature.com/articles/s42255-026-01559-z
and an associated commentary here: https://www.nature.com/articles/s42255-026-01547-3.

# Abstract
The pathogenesis of type 1 diabetes, particularly at autoantibody-positive preclinical stages, remains poorly understood, in part due to limited sample availability. Here, we show imaging mass cytometry data of pancreas samples from 88 organ donors, including 28 single and 10 multiple autoantibody-positive donors. We imaged 16 million single-cells using 79 antibodies to characterize β-cell states and the islet-immune interface, correcting for relevant covariates. We identified interactions between pro-inflammatory macrophages and exhausted-like (PD1+TIM3+) T-cells. These interactions were characteristic of early disease and of insulitic islets, indicating a key role of macrophages in disease development. β-cells showed loss of IAPP in pre-clinical disease, insulitic Interferon signatures, and no increase in three measured ER stress markers in disease samples relative to control. Multiple immune cell subtypes were associated with young age and insulitis, potentially contributing to greater disease severity in younger patients. Our data present a resource describing early type 1 diabetes progression and reveal potentially clinically actionable features before extensive β-cells loss. 

# Data access

Raw and processed data is provided on zenodo (IMC: https://zenodo.org/records/14968076). 
Please follow the respective linked directories containing the raw image data.
An informative subset of the data is also available on Pancreatlas (https://www.pancreatlas.org/) as a browsable set of annotated images. 
The processed dataset is currently being uploaded to the *imcdatasets* R/Bioconductor package (https://github.com/BodenmillerGroup/imcdatasets). 
Until then, find the data to be uploaded here: Zenodo; ().

# 🚨 How can you use this dataset in your work?

## T1D scientist:
The dataset contains imaging data of 79 unique antibodies across 88 donors (incl. 28 sAAb+ and 10 mAAb+) with rich metadata and image annotations. 
It therefore presents a rich resource for biological discovery and for validating your targets of interest.
For example, you aim to knock-out a given target, found a potentially interesting transcript of a protein included in this study, or have identified a new discrete islet type?
Use this dataset, to check its protein levels across biological contexts (disease stages, cell types, insulitis, ICIs vs IDIs etc.).

Furthermore, it will be interesting to re-use the developed *Infiltration Score*. 
This scores takes islet area, immune islet cell infiltration and distance to the islet into account. It helps to distinguish islets with strong islet infiltration,
peri-islet infiltration and islets with more distal immune cell accumulation (see Manuscript).

We also strongly encourage to browse the annotated images via Pancreatlas (https://www.pancreatlas.org/). Both Panels were acquired on directly consecutive sections
and all 79 unique markers can be jointly visualized and analyzed via Pancreatlas.

## Bioinformatician:
The dataset is highly applicable to benchmarking studies or method development.
Please note, we host all data via Zenodo, but suggest to access SCE, image (subset only), and mask data via `imcdatasets`. 
Scripts typically analyze data by panel (-> Immune panel, Islet panel).

Tasks could include:
- Cell type annotation & subclustering.
- Detection of Tissue domains / niches. 
- Methods to extract biological insights from *ordinal* cross-sectional spatial atlases.
- Extraction of biologically informative geometric features from tissue domains (-> based on islets).
- Image Registration as data is acquired on two directly consecutive 4 µm sections
- Scalable spatial methods as the dataset is large (16M cells; > 7,000 paired images (=14,050 IMC images)).
- Training of foundation models.
- Human T1D or human pancreas atlas building.
- Infiltration Scores.

🚨 Please also see our brief tutorial (`T1D_analysis/data_tutorial.Rmd`) in which we cover details of the data and explain key metadata annotations (e.g. cell-type, cell-id, ICI etc.), 
which are necessary for the tasks outlined above.
-> For any ideas or queries, feel free to reach out to me anytime (nathan.steenbuck@uzh.ch) or open an Issue in this GitHub repository.

# Repository structure
This repository contains data analysis scripts for re-creating all results and manuscript figures based on provided data objects from zenodo (https://zenodo.org/records/14968076). This repository contains the `/T1D_processing` and `/T1D_analysis` directories. For reproducing the analysis: `/T1D_processing/processing` and the `/T1D_analysis` are the most important directories. All directories contain extensive .README files giving further information. See below for details.

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

The `/Snake_pipe/` directory contains Snakefiles and configuration files to reproducibly run the pre-processing (`/T1D_processing/processing/`) including cell type annotation (01-05 scripts in `/T1D_analysis/`). 
Please note: we provide the full SCE_object generated after this steps on Zenodo.

# Software note
Most scripts are written in the statistical programming language R and Python. 
Dependencies for pre-processing and running 01-08 script of T1D_analysis are at dockerhub (**https://hub.docker.com/r/nathanste/t1d_analysis**).

For the further downstream analysis, all **R** code was run with R version 4.5.1 on *Ubuntu 24.04*. Each .Rmd file contains a Software requirements section that lists all relevant packages for data analysis and processing. File paths have to be adjusted by the user.
