# T1D_preprocessing

The code in this repository is linked to the project ___Mapping type 1 diabetes progression by imaging mass cytometry___.

It contains different folders with scripts for acquisition and processing of imaging mass cytometry (IMC) data:
* [islet_segmentation](./islet_segmentation): jupyter notebook to train a convolutional neural network for islet segmentation.
* [processing](./processing): jupyter notebooks to extract images from IMC acquisitions, segment cell and islets, and export data for downstream analysis in R.
* [region_selection](./region_selection): code to align immunofluorescence (IF) images from consecutive sections and to select regions of interest to measure by IMC.
* [samples](./samples): list of cases (patients), code for batch allocation, and code to plot cohort statistics.
* [titration](./titration): notebook to calculate optimal antibody concentrations.
* [ext](./ext): additional files for documentation.

Detailed documentation is available in each folder.


# zenodo folder

Raw and processed data should be uploaded to zenodo and made publicly available upon publication.

## Things to update before publication of this repository

Add link to the publication to this `README` file.

Add link to the zenodo folder to this `README` file.

*islet_segmentation*:  
- Add zenodo link to the processed images and masks, and to the model in the `README` file.  
- Add zenodo link to the processed images and masks that are used for training in `islet_segmentation/islet_segmentation.ipynb`.
- Add zenodo link to the trained model in `islet_segmentation/islet_segmentation.ipynb`.  

*processing*:
- Add zenodo link to the raw data files and to panel csv files to `processing/01_Preprocessing.ipynb`.
- Add zenodo link to the islet segmentation model in `processing/02_IsletSegmentation.ipynb` (`Load the model` paragraph).  

*region_selection*:
- Add zenodo link the the islet classifier and to the czi bf and if files to the `README`.

*samples*: OK

*titration*:
- Add link to the zenodo titration repository in the `README.md` file.
- Add link to the zenodo titration repository in the `titration/titration.ipynb` notebook.
