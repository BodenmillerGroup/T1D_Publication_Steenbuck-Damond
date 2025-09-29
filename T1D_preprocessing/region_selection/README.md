# **REGION SELECTION**

## **Introduction**

The goal of this workflow is to define the regions (on tissue sections) that will be measured by Imaging Mass Cytometry (IMC). This selection is based on immunofluorescence (IF) images of the same section (combined IF-IMC protocol, see the "Methods" section of the associated publication).  

This workflow works with two or more antibody panels acquired on consecutive sections. The panel used to define the regions of interest (ROIs) to be measured by IMC is thereafter referred to as the "base" panel. In the context of this workflow, the `Islet` panel is generally used to define ROIs.  

The main steps of this workflow are:
- Acquire immunofluorescence (IF) and brightfield (BF) images on a slide scanner.
- Convert the slide scanner images to `.tiff` files and process them to be compatible with this workflow.
- Define ROIs on the "base" panel IF images.
- Automatically register (align) the IF images from the "base" panel and from the other panel(s) that were acquired on consecutive section(s).
- Apply the transformation learned during registration to adapt the ROIs defined on the "base" panel to the consecutive sections.
- Record landmarks (identical points) on the IF images and on Fluidigm's CyTOF software to perform IF to IMC registration.
- Based on the landmarks, learn the IF to IMC transformation
- Convert the coordinates of the ROIs selected on IF images to IMC coordinates.

A constant naming scheme must be used throughout this workflow. The following elements are relevant:
- `CaseID`: 4-digit identification number of the pancreas donor, corresponding to nPOD numbering system (https://www.jdrfnpod.org). Below, we use `1234` as CaseID example.
- `Panel`: name of the antibody panel, e.g. `Islet` or `Immune`.  

All files generated during execution of this workflow are stored in a folder of your choice on your storage disk. This directory is herein referred to as the `BASE` folder. This workflow creates a separate subfolder for each case (patient) in the dataset to organize the generated files.

## **Requirements**
This workflow has been tested on Windows 10. It requires the following tools:
- **Fiji**: https://imagej.net/software/fiji. All `.ijm` and `.py` scripts in this repository are run with FiJi. This workflow has been tested with FiJi based on ImageJ v1.53p.
  - Image registration requires version 3.0.7 of the "*register virtual stack slices plugin*" not the version 3.0.5 that is installed by default in FiJi => Remove the version 3.0.5 in the FiJi/plugins folder and replace it by the "*register_virtual_stack_slices-3.0.7.jar*" file which is present in the `ext` folder of this repository, then restart FiJi.
- **Jupyter notebook**: ideally with the anaconda environment defined in the `../environment.yml` file.
- **CZItoTIFF batch converter**: a copy of this software (Windows only) and of the associated documentation can be found in the `ext` folder. The `.tiff` files can be generated using other methods but the image importation scrips (`01_ImportSlideScanner_IF.ijm` and  `02_ImportSlideScanner_BF.ijm`) should be modified accordingly.
- **CyTOF software** Is installed on all Imaging Mass Cytometry machines. It can also be downloaded for local use from https://www.fluidigm.com/products-services/software (Hyperion imaging system - Non-acquisition workstations only). Tested with version 7.0.8493
- **CellProfiler** *(required for automated ROI selection only)*: https://cellprofiler.org. Tested with CellProfiler 4.2.1.
- **ilastik** *(required for automated ROI selection only)*: https://www.ilastik.org. Tested with ilastik 1.3.3post3.

## **Tutorial**

### A. Slide scanner image acquisition
The immunostaining protocol for combined IF and IMC imaging can be found in the "Methods" section of the publication.   
After immunostaining, brightfield and fluorescence imaging were acquired with a ZEISS AxioScan.Z1 slide scanner.
Raw brightfield and immunofluorescence `*.czi` files can be downloaded from the associated zenodo repository ***ADD LINK TO CZI FILES***

#### A.1. Brightfield imaging
- Immediately after the end of the immunostaining protocol, slides were inserted in the slide scanner for brightfield imaging.
- The imaged region corresponds to the entire glass slide area (excluding the label).
- All slides are imaged using the `T1D_BF_2.5X.czspf` scan profile (a copy of this scan profile can be found in the `ext` folder).
- The output images are named as `CaseID_BF_Panel.czi` (e.g. "1234_BF_Islet.czi").

#### A.2. Fluorescence imaging
- Fluorescence images are acquired immediately after brightfield images (within 24 hours after the end of the immunostaining protocol).
- The imaged region corresponds to the entire tissue area (selected manually).
- Image all tissue sections using the `T1D_IF_10X.czspf` scan profile (a copy of this scan profile can be found in the `ext` folder).
- The output images are named as `CaseID_Panel.czi` (e.g. "1234_Islet.czi").
- After the acquisition, the quality of the images is checked. If some images or image regions are out of focus, the corresponding slides should be imaged again.

### B. Image file conversion
The images generated by the slide scanner (`.czi` format) are converted to `.tiff` files using the CZItoTIFF batch converter (Windows only). A copy of this software and the associated documentation can be found in the `ext` folder.  

#### B.1. Fluorescence images
- In the `BASE` folder (a folder of your choice on your storage disk where all the files generated by this workflow will be stored), create a subfolder named `czi_if` to store the fluorescence images.
- Place all the fluorescence CZI files generated by the slide scanner in the `BASE/czi_if` folder.
- Convert the CZI files using the CZI-TIFF converter using the following settings:
![Conversion settings for brightfield images](../ext/img/img_czi-to-tiff_IF.png)

#### B.2. Brightfield images
- In the `BASE` folder create a subfolder named `czi_bf` to store the brightfield images.
- Place all the brightfield CZI files generated by the slide scanner in the `BASE/czi_bf` folder.
- Convert the CZI files using the CZI-TIFF converter using the following settings:
![Conversion settings for brightfield images](../ext/img/img_czi-to-tiff_BF.png)

### C. Slide scanner images import
In this section the `.tiff` images generated at the previous step are processed (rotated, auto-contrasted) so that they can be used in this workflow.  
A folder tree is created in the `BASE` folder that will be used to store newly generated files: one subfolder is created per Case ID. In each "case" folder, the following subfolders are created and will be filled up by running the different scripts in this workflow:
- *mask*: masks representing the selected regions of interest (ROIs).
- *mcd*: images that can be loaded in `.mcd` files (proprietary file format used on Imaging Mass Cytometry machines).
- *merged*: merged IF images that combine the (three) IF channels acquired on the slide scanner.
- *original*: DAPI images that are used for registration of consecutive sections and the related `.xml` files containing the registration information.
- *registered*: registered IF images (obtained after applying the transformation to the merged IF images).

#### C.1. Fluorescence images
In FiJi, open the script `01_ImportSlideScanner_IF.ijm`. This script processes DAPI images so that they can be used for registration of consecutive tissue sections and combines the individual channel fluorescent images to generate merged IF images that will be used for region selection.
- Adapt the settings if needed (see the instructions inside the script).
- Run the script. Enter the `BASE` folder when prompted to select the input folder.

#### C.2. Brightfield images
In FiJi, open the script `02_ImportSlideScanner_BF.ijm`. This script processes brightfield images so that they can be loaded in `.mcd` files.
- Adapt the settings if needed (see the instructions inside the script).
- Run the script. Enter the `BASE` folder when prompted to select the input folder.

### D. Image registration
Here, DAPI images from consecutive tissue sections are registered. The transformation is learned and stored in `.xml` files. This will allow to apply the exact same transformation to merged IF images and to ROI masks.   
Requires version 3.0.7 of the virtual stack slices registration plugin (see requirements above) !

#### D.1. Automated image registration
In FiJi, open the `03_RegisterImages.py` script.
- Enter the input directory (your `BASE` folder).
- Enter the list of case IDs and the panels as indicated.
- Run the script.
- After registration, an image stack is automatically opened. Check that the two images in the stack are well aligned and close the stack.

#### D.2. Manual image registration
Usually, the automated registration works well. In case it does not the images can be registered manually (one case at a time):
- In FiJi, open Plugins > Register > Register Virtual Stack Slices.
- Select the `original/Panel` (e.g. `original/Immune`) folder as the source directory.
- Select the `registered/Panel` (e.g. `registered/Immune`) folder as the output directory.
- Select "*Affine*" for feature extraction model and registration model.
- Tick the "*Save transforms*" box and OK.
- Select the `original/Panel` folder as the folder where to save transformations.
- Select the merged IF image corresponding to the current panel (e.g. `1234_Merge_Immune.tif`) in the `original` as the reference image. Do *not* select the base panel image as the reference image.
- Check that the registration worked on the virtual stack, then close it.

If this does not work, the images in the `original/Panel` should be made more similar. This can for instance be done by adjusting the contrast/brightness to similar values. Another option is to apply Process > CLAHE in FiJi. If some structure is present on one image but not the other, that structure can be cropped out. The merged IF images in the `merged` folder can also be used for as source images if the DAPI images in the `original/Panel` folder do not work (in that case, copy the output `.xml` files to the corresponding `original/Panel` folder).

### E. Creation of MCD files
At this stage, `.mcd` files for CyTOF acquisition can be created. Brightfield images are imported and panorama regions are selected. The panoramas are then acquired and will be used for registration of the IF and the IMC images.

- Open the CyTOF software and create a new file for each slide to measure.
- *OPTIONAL* Import the low resolution slide image from the `mcd` folder (named as `CaseID_BF_Panel_2_pt2.tif` e.g. 1234_BF_Islet_2_pt2.tif). In "Image properties" (right of the Import Slide Image window), enter the image size in micrometers. The image size can be viewed by opening the image in FiJi.
- Import the high resolution tissue image from the `mcd` folder (named as `CaseID_BF_Panel.tif` e.g. "1234_BF_Islet.tif"). Enter the image size as above.
- Place the imported image(s) at their (approximately) correct position on the canvas.
- Draw panoramas. Select 4 to 6 regions at the periphery of the tissue section.
- Save the file as `CaseID_Panel.mcd` (e.g. "1234_Islet.mcd").

Once all the `.mcd` files have been created, go to the CyTOF/Hyperion, open the `.mcd` file and load the corresponding slide in the machine. Create the selected panoramas and save the file.

### F. Automated region selection
In this section, regions of interest (ROIs) are automatically selected based on islets. To this aim, an ilastik classifier is trained to recognize islets and CellProfiler is used to segment the islets, which are then imported in FiJi's ROI manager for downstream processing.  
*Note*: The ROI selection is performed on the "base" panel only. For other panels (acquired on consecutive sections) the ROIs defined here will be transformed.

#### F.1. Train an ilastik classifier for islet segmentation
First, random crops of the IF image corresponding to the islet channel are generated. These crops are then loaded in ilastik and a classifier is trained. The probabilities for islet and non-islet regions are then exported.
A pre-trained classifier can be downloaded from ***ADD LINK TO CLASSIFIER***.
- In FiJi, run the `04a_Autoselect_BatchCrop.ijm` script to generate random crops of the slide scanner IF images corresponding to the islet marker (usually CD99, SYP or CHGA).
- Open ilastik and create a new classifier or open the pre-trained `islet_selection.ilp` classifier (in this case, go directly to *Prediction Import* below).
  - In *Input Data* load the random crops from the `BASE/islet_selection/crops` folder.
  - In *Feature Selection* select all the features.
  - In *Training* create two labels ("Islet" and "Exo_bg"). Train the classifier to recognize islets (label "Islet") and regions corresponding to non-islet tissue and to tissue-less regions (label "Exo_bg").
  - In *Prediction Export* choose *Probabilities* as the source and modify the export settings as follows:
    - [X] Convert to Data Type unsigned 16-bit.
    - [X] Renormalize.
    - Output format = `.tiff`.
    - File = `BASE/islet_selection/probabilities/{nickname}_{result_type}.tiff`.
  - In *Batch processing*, select raw data files, import all the `CaseID_Panel-Channel_BIN.tif` from the `BASE/CaseID/mcd` folders (e.g., "1234_Islet-Cy5_BIN.tif" if *Cy5* is the islet marker channel and *Islet* is the base panel used for ROI selection).
  - Click on *Process all files*.

#### F.2. Segment islets with CellProfiler
Probabilities generated from ilastik are loaded into CellProfiler, islets are segmented and islet masks are exported.
- In CellProfiler, load the `04b_Autoselect_SegmentIslets.cppipe` pipeline (File > Import > Pipeline from file).
  - In *Images* load all files in the `BASE/islet_selection/probabilities` as input.
  - In *Metadata* make sure that the file names are correctly parsed (the CaseID, Channel and Panel columns should be correctly populated).
  - In *Output Settings* enter the `BASE` folder as the default output folder.
  - Click on *Analyze Images*

#### F.3. Import segmented islets as FiJi ROIs
In this section islets are identified from islet masks exported from CellProfiler. The islets masks are thresholded and loaded as a set of ROIs in FiJi. The bounding boxes of these islets are calculated to obtain rectangular ROIs. Finally, ROIs are randomly subsetted.

- Open the `04c_Autoselect_ImportROIs.ijm` script in FiJi.
- Adapt the parameters as indicated in the script. In particular, the *nb_rois* variable indicates the maximal number of ROIs to return.

### G. ROI selection
The ROIs automatically selected at the previous step are interactively adapted to define a final selection of regions of interest to measure by IMC. This initial selection is performed on the "base" panel. A color mask indicating the ROI number, position and size is created. This color mask is then transformed to adapt the ROI selection to other panels acquired on consecutive sections. The coordinates and size of all ROIs are then recorded, saved as FiJi ROI sets and exported as `.csv` files.

#### G.1. Interactively modify ROIs
Automatically generated ROIs can here be manually adapted using FiJi's ROI manager.
- Open the `05_DefineROIs.ijm` script in FiJi and adapt the parameters if needed.
- The script automatically opens merged IF images and displays the automatically selected ROIs.
- Adapt the ROIs if needed. ROIs can be deleted, adjusted or repositioned using FiJi rectangle tool (click on "update" or type "u" after each ROI modification). New ROIs can also be created (type "t" to add the selected rectangle to the ROI manager).
  - If needed, use *Image > Adjust > Color Balance* in FiJi to increase the intensity of individual IF channels.
  - Delete ROIs that are too close to each other (to avoid overlapping during IMC acquisition).
  - Delete ROIs that are outside the IMC acquisition area (shown as a green rectangle in `.mcd` files).
- Once the ROI selection is satisfactory, press "OK" to go to the next slide.
- All ROI sets are saved in case-specific folders as `CaseID_iBoxes_Original.zip`.

#### G.2. Create ROI color masks
Color masks containing the ROI number, size and position are created.
- Open the `06_CreateMaskFromROI.ijm` script in FiJi and adapt the parameters if needed.
- The ROIs defined at the previous step are loaded and enlarged by a fixed value defined by the "enlargement" parameter (pixel units).
- Each ROI is measured (color intensity, x and y coordinates, w and h size) and these measurements are exported to `.csv` files named `TEMPBoxes_Panel.csv`.
- A color mask representing all ROIs is created for each slide. Each ROI has a different color that will later be used to identify and match ROIs from different panels.
- A dummy (blank) image of the same size as the merged IF images is created for each non-base panel (panels that were *not* used to select ROIs). This dummy image will be used for mask transformation at the next step.

#### G.3. Transform ROI masks
ROIs defined on the "base" panel are transformed to fit images from other panels acquired on consecutive sections. For this, the transforms learned in section D are applied to the ROI color masks.
- In FiJi, open *Analyze > Set Measurements* and make sure that *Min & Max gray value*, *Bounding rectangle*, and *Add to overlay* are selected.
- Open the `07_TransformROIMask.py` script in FiJi and adapt the input parameters if needed (the same parameters as in `03_RegisterImages.py` can be used).
- The script automatically transforms the color masks (output stored in the `mask` subfolder).
- The merged IF images created in section C are transformed in the same way (output stored in the `registered` subfolder).
- In addition, the DAPI images used for registration of consecutive sections in section D are deleted to save space.

#### G.4. Interactively modify ROIs
Here, transformed ROIs can be interactively modified, exactly like in paragraph G.1.
- In FiJi, open the `08_AdaptTransformedROIs.ijm` script and adapt the parameters if needed.
- This script loads ROIs from the transformed color masks and displays the ROIs on the corresponding merged IF image.
- The ROIs can be repositioned or resized as in paragraph G.1.
- ***Do not add or delete ROIs!*** Adding or deleting ROIs will break the matching of ROIs across different panels.
- In the end, each ROI is measured (color intensity, x and y coordinates, w and h size) and these measurements are exported to `.csv` files named `TEMPBoxes_Panel.csv`.

#### G.5. Collect all ROI information
All ROIs are measured and the measurements are exported to `.csv` files.
- Open Juypter notebook and load the `09_SummarizeROIs.ipynb` pipeline.
- Modify the input directory, and lists of cases and panels if needed. Run all notebook cells.
- All ROI measurements are exported as `CaseID_Boxes.csv` and `CaseID_Boxes_Panel.csv`.
- A per slide summary of ROI numbers and total area is exported in the `BASE` folder as `acquisition_summary.csv`.

### H. Landmark selection for IF to IMC registration
In this section, landmarks found on both the merged IF image and on the brightfield image in the CyTOF software are identified. These landmarks represent identical points on the two images of the same section and will be used at the next step to convert the ROI coordinates defined on IF images to IMC coordinates. This is first done on the "base panel" images, then the coordinates are transformed to facilitate the selection of similar landmarks on the "non-base panel" images.

#### H.1. Select the landmarks
Repeat these operations for all "base panel" images in the `BASE` directory (one image per case):
- Open the "base" panel merged IF image (folder `merged`) in FiJi.
- If needed, use *Image > Adjust > Color Balance* in FiJi to increase the intensity of individual IF channels.
- Open the corresponding `.mcd` file in the CyTOF software.
- Open the `CaseID_Coordinates_Panel.csv` file for the "base" panel.
- Identify 4 to 6 landmarks that are recognizable on the IF image in FiJi and on the brightfield image in the CyTOF software.
- Add each IF landmark to FiJi's ROI manager by clicking on the landmark with the "Point" tool and typing "t".
- Record the coordinates of each landmark in the acquired panoramas of the `.mcd` file by moving the mouse cursor on the landmark and writing down its X and Y coordinates in the opened `.csv` file (columns *IMCX* and *IMCY*).
- Do not record landmarks in areas where no panorama has been acquired!
- When all landmarks have been defined, select all points in the ROI manager and click on *More > List* to display the landmark coordinates.
- Copy the content of the *X* and *Y* columns to the *SlideX* and *SlideY* columns of the `.csv` file.
- Number the landmarks from 1 to N in the *CoordNb* column and save the `.csv` file.
- In the end the `.csv` file should look as follows:  
![Landmarks](../ext/img/landmarks.png)
- In the ROI manager, click on *More > Save* and save the ROI set as `CaseID_Coordinates_Panel.zip` (e.g. "1234_Coordinates_Islet.zip").

#### H.2. Transform the landmarks
After landmark coordinates have been recorded for the "base panel" images of all cases in the `BASE` folder, the landmarks are transformed, similar to what was done for ROIs in paragraphs G.2, G.3 and G.4.
- Run the `10_CreateMaskFromCoords.ijm` script in FiJi.
  - This script creates a color mask from the coordinates that were recorded in H.1.
  - It also creates a dummy (blank) image that will be used for transformation of the color mask.
- Run the `11_TransformCoordMask.py` python script. The same input parameters can be used as in scripts 03 and 07.
	- Similar to the script "07_TransformROIMask.py", this script applies the transformation learned in section D to the coordinates mask to match the images obtained from "non-base panels".

#### H.3. Adapt the landmarks
Here, landmarks for the non-base panels are defined. The transformed landmarks obtained at the previous step are used as a starting point but they must be adapted and the matching IMC coordinates have to be recorded.
- Open the `12_AdaptTransformedCoords.ijm` script in FiJi. Modify the input parameters if needed and run the script.
- The script will automatically open the registered merged IF images (from the `registered` folder) and display the transformed landmarks.
- The next steps are similar to paragraph H.1:
  - Open the corresponding `.mcd` file in the CyTOF software.
  - Open the `CaseID_Coordinates_Panel.csv` file of the corresponding panel.
  - Landmarks can be modified using the "Point" tool in FiJi (press "update" in the ROI manager), added (type "t") or deleted.
  - Record the matching coordinates in the acquired panoramas of the `.mcd` file.
  - Do *not* record landmarks in areas where no panorama has been acquired!
  - Write down the CyTOF coordinates in the `.csv` file as above (columns *IMCX* and *IMCY*).
  - *Warning:* make sure that the landmark number in the *CoordNb* column matches the landmark number in the ROI manager!
  - When all landmarks have been defined, select all points in the ROI manager and click on *More > List* to display the landmark coordinates. Copy the content of the *X* and *Y* columns to the *SlideX* and *SlideY* columns of the `.csv` file.
  - The `.csv` file content should look like the one shown in paragraph H.1.
  - Press *OK* to continue.
  - The landmark ROI set is automatically saved as `CaseID_Coordinates_Panel.zip` (e.g. "1234_Coordinates_Immune.zip").
- Repeat this operation for all "non-base" panel images in the dataset.

### I. IF to IMC registration
The matched landmarks identified on the IF images and in the CyTOF software are used to transform the ROIs defined on the IF images. These transformed ROIs can then be directly imported in the CyTOF software.

#### I.1. Register landmarks and transform ROIs
- Open Jupter notebook and load the `13_Register_IF_IMC.ipynb` pipeline.
- Adapt the parameters. One slide (case and panel combination) is performed at a time.
- Run all cells and make sure that the transformed landmarks do overlap in the quality control step.
- For every slide, a file named `CaseID_IMCBoxes_Panel.csv` (e.g. "1234_IMCBoxes_Islet.csv") is generated.


#### I.2. Finalize MCD files for IMC acquisition
Finally, the transformed ROIs are imported in the CyTOF software for Imaging Mass Cytometry acquisition.
- In the CyTOF software, open an `.mcd` file.
- Click on *Import* and load the corresponding `CaseID_IMCBoxes_Panel.csv` file that was generated at the previous step.
- Check that the imported ROIs are located in the right place with regard to the acquired panoramas.
- Click on *Template* to load the corresponding antibody panel.
- The two panel files for the CyTOF software can be found in the `ext` subfolder of this repository (`Panel_Immune.conf` and `Panel_Islet.conf`).
- Adapt the ablation energy and frequency if needed.
- Save the `.mcd` file, which is now ready for IMC acquisition.
- The two IMC panels (to use with the CyTOF)
