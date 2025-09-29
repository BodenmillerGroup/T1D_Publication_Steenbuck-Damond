# Samples

This folder contains various information about the pancreas donors (cases) that
were used in this study.
All cases and clinical data were collected by nPOD (https://www.jdrfnpod.org).

## Case List
The `CaseList.csv` file is a list of all cases used in this study, with
available clinical data. Full information about the cases can be found at
https://www.jdrfnpod.org/for-investigators/online-pathology-information/.
For simplicity, cases HLD011 and CV02 were renamed as 8011 and 8002,
respectively, in the data processing and analysis sections.

## Batch allocation
Cases were split into four batches that were sectioned, stained, and imaged
separately.  
Variance minimization was used to allocate cases from the single
autoantibody-positive (sAAb+), multiple autoantibody-positive (mAAb+), and T1D
recent-onset groups to the four batches. Cases 6553, 6558, 6562, and 6563 were
not available at the beginning of the study and were added later on to batch 4.
The variance minimization procedure is described in `Sella, F., Raz, G., &
Cohen Kadosh, R. (2021). When randomisation is not good enough: Matching groups
in intervention studies. Psychonomic bulletin & review, 28(6), 2085–2093`.  

Batch allocation was performed using the R function `VarMin.R` provided with
that article and available here:
https://osf.io/6jfvk/?view_only=1998adab67f64f6895168be4684388ac.
The `VarianceMinimization_Template.csv` file, which containing the
characteristics of pancreas donors used in this study, was provided as input to
the VarMin function.

## Cohort statistics
The `CohortAnalysis.Rmd` file is an R script for plotting the major
characteristics of the cohort, such as gender repartition and age distribution.

## Installation

To run the `CohortAnalysis.Rmd` script, a anaconda environment containing all required packages
can be installed using the provided `environment.yml` file.

The environment can be created using the following commands:
- `conda env create -f environment.yml` # Adapt the path to the environment.yml file if needed.  
- `conda activate samples_env`

## Optimal ROI number
The optimal number of region of interest (ROIs) to measure per image for this
study was based on a previous pilot study, performed on a smaller number of
cases (16) but with a large number of regions measured per case (>100).  

We calculated that subsetting the number of ROIs to 75 per case still allowed
to recover at least 95% of the immune clusters identified in the pilot study.
The results are available in the `optimal_roi_number_*.csv` files, with
`Immune`, `Lympho`, and `Myelo` corresponding to the three immune panels of the
Pilot study.
