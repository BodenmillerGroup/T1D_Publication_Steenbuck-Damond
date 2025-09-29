
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/mnt/snakemake', '/scratch/nsteen/T1D_preprocessing/processing']); import pickle; snakemake = pickle.loads(b'\x80\x04\x95\x86\x05\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94(\x8cD/scratch/nsteen/T1D_preprocessing/processing/05_ImageRegistration.py\x94\x8c9/scratch/nsteen/processing/txt_output/04_Measurements.out\x94\x8c@/scratch/nsteen/T1D_preprocessing/processing/helper_functions.py\x94e}\x94(\x8c\x06_names\x94}\x94\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x12\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x18)}\x94\x8c\x05_name\x94h\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94\x8c>/scratch/nsteen/processing/txt_output/05_ImageRegistration.out\x94a}\x94(h\x0e}\x94h\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94}\x94(h\x0e}\x94h\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94}\x94(h\x0e}\x94h\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bub\x8c\x07threads\x94K\x08\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x08K\x01M\x00\xfaM\xe8\x03\x8c\r/sctmp/nsteen\x94\x8c\x073:00:00\x94e}\x94(h\x0e}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06mem_mb\x94K\x02N\x86\x94\x8c\x07disk_mb\x94K\x03N\x86\x94\x8c\x06tmpdir\x94K\x04N\x86\x94\x8c\x04time\x94K\x05N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bhYK\x08h[K\x01h]M\x00\xfah_M\xe8\x03hahU\x8c\x04time\x94hVub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x0e}\x94h\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bub\x8c\x06config\x94}\x94(\x8c\x08data_dir\x94\x8c\x10/scratch/nsteen/\x94\x8c\x08home_dir\x94\x8c"/home/nsteen/scripts/T1D_analysis/\x94\x8c\rcontainer_dir\x94\x8c&/home/nsteen/data/singularity/sandbox/\x94u\x8c\x04rule\x94\x8c\x12Registration_05_py\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8c,/scratch/nsteen/T1D_preprocessing/processing\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = True; __real_file__ = __file__; __file__ = '/scratch/nsteen/T1D_preprocessing/processing/05_ImageRegistration.py';
######## snakemake preamble end #########
#!/usr/bin/env python
# coding: utf-8

# # Image Registration
# 
# This is the fifth and final script in the processing pipeline for IMC data.
# 
# The goal is to register the islet panel and immune panel IMC images. 
# This is performed using functions from scikit-image and opencv.
# 
# By registration, we will be able to estimate distances between panels in downstream analyses.
# For example, the distance of a beta-cell subpopulation on the Islet panel and an immune-cell subpopulation on the Immune panel.
# Note: The IMC images are taken from 4 µm consecutive sections.
# Hence, cells between panels are expected to be the same cell in around 50 % of the time.
# 
# The following steps are performed:
# 
# **Global registration**
# - Global registration using Phase Cross-Correlation
# 
# **Local registraiton**
# - Local registration using ICP (Iterative Closest Point).
# 
# **QC**:
# - The registration is used to perform some QC to delete/edit images.
# 

# ## **Configuration**
# 
# ### **Import packages**

import numpy as np
from pathlib import Path
import pickle
import sys
import pandas as pd

import matplotlib.pyplot as plt
from steinbock import io
import seaborn as sns
import cv2 as cv

from skimage.registration import phase_cross_correlation

sys.path.append("/home/T1D_preprocessing/processing/")
from helper_functions import pad_images
from helper_functions import icp
from helper_functions import plot_centroids
from helper_functions import write_out_dfs


# ### **Load directories and panels**
# 
# Paths to input and output folders as well as antibody panels were exported by the first script (`01_Preprocessing.ipynb`). Here they are imported again.
print(sys.path)
print(sys.executable)
base_dir = Path("/home/processing/")

with open(base_dir / "variables/folders.txt", "rb") as handle:
    folders = pickle.loads(handle.read())
print(folders)

# FIXME: uncomment this for development; delete for actual run.
#folders["img"] = Path("/mnt/T1D_Vol/processing/img")
#folders["data_cells"] = Path("/mnt/T1D_Vol/processing/data_cells")
#folders["seg_cells"] = Path("/mnt/T1D_Vol/processing/seg_cells")
#folders["seg_islets"] = Path("/mnt/T1D_Vol/processing/seg_islets")
#folders["masks_islets"] = Path("/mnt/T1D_Vol/processing/masks_islets")


with open(base_dir / "variables/panels.txt", "rb") as handle:
    panels = pd.read_pickle(handle)

for panel_name, panel in panels.items():
    print(panel_name, "\n", panel.head())
panel_names = list(panels.keys())


# # Prepare Data
# 
# First, a dataframe (`image_dfs`) containing paths to all three types of input data is prepared.
# 
# Load in:
# 1) Merged islet segmentation channels (CHGA, SYP)
# 2) Binary segmentation masks
# 3) Regionprops-files: contain coordinates, metadata


image_dfs = dict.fromkeys(panel_names)

for pan_idx, panel_name in enumerate(panels):
    image_dfs[panel_name] = pd.DataFrame(columns=["image_name", "img_path", "seg_subdir"])

print(image_dfs)

measurement = "regionprops"
for pan_idx, panel_name in enumerate(panels):
    print(pan_idx, panel_name)

    #------------------------------------------------------------------------------
    ## Load combined Islet channels (CHGA, SYP) images:
    ## These were used to segment the islets.
    seg_subdir = folders["seg_islets"] / panel_name
    input_images = sorted([x.name for x in Path.iterdir(seg_subdir) if x.name.endswith(".tiff")])
    print("Number of images in", panel_name, "panel:", len(input_images))
    
    for img_path in input_images:
        # Add image metadata to the data frame
        image_name = img_path.replace(".tiff", "")
        image_dfs[panel_name] = pd.concat(
            [image_dfs[panel_name], pd.DataFrame({
                "seg_subdir": seg_subdir,
                "image_name": image_name,
                "img_path": img_path,
                "measurement": "image"},
                index = [len(image_dfs[panel_name].index)])])
    
    #------------------------------------------------------------------------------
    ## Load binary Islet masks:
    seg_subdir = folders["masks_islets"] / panel_name
    input_images = sorted([x.name for x in Path.iterdir(seg_subdir) if x.name.endswith(".tiff")])
    print("Number of images in", panel_name, "panel:", len(input_images))
    
    for img_path in input_images:
        # Add image metadata to the data frame
        image_name = img_path.replace(".tiff", "")
        image_dfs[panel_name] = pd.concat(
            [image_dfs[panel_name], pd.DataFrame({
                "seg_subdir": seg_subdir,
                "image_name": image_name,
                "img_path": img_path,
                "measurement": "bin_mask"},
                index = [len(image_dfs[panel_name].index)])])
    
    #------------------------------------------------------------------------------
    ## Load regionprops measurements:
    ## These contains the centroids of the islets.
    seg_subdir = folders["data_cells"] / panel_name / measurement
    input_csvs = sorted([x.name for x in Path.iterdir(seg_subdir) if x.name.endswith(".csv")])
    print("Number of images in", panel_name, "panel:", len(input_csvs))

    for img_path in input_csvs:
        # Add image metadata to the data frame
        image_name = img_path.replace(".csv", "")
        image_dfs[panel_name] = pd.concat(
            [image_dfs[panel_name], pd.DataFrame({
                "seg_subdir": seg_subdir,
                "image_name": image_name,
                "img_path": img_path,
                "measurement": measurement},
                index = [len(image_dfs[panel_name].index)])])
        
image_dfs["Immune"].head()


# ## **Process Image DF**
# 
# Further process the image DF. 
# 
# Remove 2 images that are missing in the data_cells folder. 
# 
# Finally, merge the 2 dataframes (Immune/Islet) using the image identifier column.

for pan_idx, panel_name in enumerate(panels):
    image_dfs[panel_name].pivot(index = "image_name", columns="measurement", values="img_path")
    image_dfs[panel_name]["case_id"] = image_dfs[panel_name]["image_name"].apply(lambda x: x.split("_")[0]).astype(str)
    image_dfs[panel_name]["ROI"] = image_dfs[panel_name]["image_name"].apply(lambda x: x.split("_")[3]).astype(str) 
    image_dfs[panel_name]["img_id"] = "roi_" + image_dfs[panel_name]["case_id"] + "_" + image_dfs[panel_name]["ROI"]
    image_dfs[panel_name] = image_dfs[panel_name].drop(columns=["ROI", "case_id"])
    # Remove rows with NA values -i.e. missing images.
    image_dfs[panel_name] = image_dfs[panel_name].dropna()

image_dfs["Immune"].head()


# Pivot the dataframe, which introduces hierarchical indexing.

image_dfs["Islet"] = image_dfs["Islet"].pivot(index = "img_id", columns="measurement", values=["img_path", "seg_subdir"])
image_dfs["Immune"] = image_dfs["Immune"].pivot(index = "img_id", columns="measurement", values=["img_path", "seg_subdir"])

df = image_dfs["Islet"].merge(image_dfs["Immune"], on = "img_id", how = "inner", suffixes = ("_Islet", "_Immune"))
df.head()


# Further add the output directories to the image_df.
# The following directories are created:
# - registered_regionprops -> the transformed coordinates
# - inter_panel_neighbors -> nearest neighbors across panels (can be same cell)

# Add output directories for registration and inter-panel neighbors.
for panel_name in panel_names:
    seg_subdir = "seg_subdir_" + panel_name
    img_path = "img_path_" + panel_name

    # Registration subdirectories and paths (output):
    df.loc[:, (img_path, 'registered_regionprops')] = df.loc[:, (img_path, 'regionprops')].values
    df.loc[:, (seg_subdir, 'registered_regionprops')] = folders["data_cells"] / panel_name / "registered_regionprops"

    # Inter-panel neighbor Subdirectories and paths (output):
    df.loc[:, (seg_subdir, 'inter_panel_neighbors')] = folders["data_cells"] / panel_name / "inter_panel_neighbors"
    df.loc[:,(img_path, 'inter_panel_neighbors')] = df.loc[:, (img_path, 'regionprops')].values
    

## Final processing of the dataframe.
df = df.dropna()
df["case_id"] = df["img_path_Islet"]["bin_mask"].str.split("_", expand = True)[0]
df["islet_roi_id"] = df.loc[:, ("img_path_Islet", "bin_mask")].apply(lambda x: x.replace(".tiff", "").replace("ROI_", ""))
df["immune_roi_id"] = df["img_path_Immune"]["bin_mask"].apply(lambda x: x.replace(".tiff", "").replace("ROI_", ""))
df.head()


# ## Quality Control:
# 
# Acquisitions from case_id `6526` are rotated. Here, helper functions are defined to rotate these images. `
# 
# Another QCed acquisition `6209_Islet_ROI_013` will be deleted in the first script of the analysis pipeline (`01_ImportData_cells_{panel}.Rmd`). Here, Islet and Immune panel images do not match.
# 
# `6362_Islet_ROI_052` + `6362_Islet_ROI_056` are switched. -> has been fixed upstream (`01_Processing.ipynb`)
# 
# ### Helper Functions 
# Define Helper Functions to rotate the affected images.

def flip_coords(coordinates, y_max):
    rotated_coordinates = [(-y+y_max, x) for x, y in coordinates]
    df_rot = pd.DataFrame(rotated_coordinates)
    return df_rot


# # Perform Registration
# 
# The Registration between Immune and Islet panel images is performed by first performing global registration and then local registration. 
# 
# ## Global Registration
# 
# Use **phase cross-correlation** to perform rigid transformation (translation only) for global registration.
# 
# Given that **a)** the panel-images have been pre-aligned and **b)** that IMC ROIs are small we make the assumption that only translation is sufficient for accurate registration (=`checked`).
# 
# For the global registration we use both the **binary islet segmentation images (0/1)** and the **merged islet channel images (CHGA+SYP)** used for the islet segmentation (see `02_IsletSegmentation.ipynb`). 
# 
# ## Local Registration:
# 
# Use **ICP (Iterative Closest Point)** to perform local registration. ICP aligns two point clouds, which are the globally registered cell-centroid coordinates.
# 
# We compute the **residual error** (`error_max`) between source and destination point cloud (Islet and Immune coordinates) for both inputs (merged image and binary segmentation mask) and choose the registration with lower residual error.
# 
# Note: 
# - ICP is used for cloudpoints with minor transformations. Thus, the overall registration depents greatly on the initial pose estimation, i.e. global registration.
# 

# ### Registration Function:
# 
# `df` = Input Image_df
# 
# `idx` = Index in the Image_df. One Acquisition is treated at a time.
# 
# `sus` = Array to store indices of potentially dangerous acquisition. Used for QC.
# 
# `upsample_factor` = Value above 1 to use upsampled matrix-multiple DFT for phase cross-correlation
# see (https://scikit-image.org/docs/stable/auto_examples/registration/plot_register_translation.html).
# 
# `binary` = Use either the binary islet segmentation mask (=True) or the merged islet channels. 
# 


def register_images(df, idx, sus = [], upsample_factor=1, binary = False, verbose = False):
    #------------------------------------------------------------------------------------------------------------------
    ## Read in the Regionprop files (coordinates)
    case_id = df["case_id"].iloc[idx]
    islet_roi_id = df["islet_roi_id"].iloc[idx]
    immune_roi_id = df["immune_roi_id"].iloc[idx]

    isl_regionprops_file = df["seg_subdir_Islet"]["regionprops"].iloc[idx] / df["img_path_Islet"]["regionprops"].iloc[idx]
    isl_regionprops = pd.read_csv(str(isl_regionprops_file))
    
    imm_regionprops_file =  df["seg_subdir_Immune"]["regionprops"].iloc[idx] / df["img_path_Immune"]["regionprops"].iloc[idx]
    imm_regionprops = pd.read_csv(str(imm_regionprops_file))

    ## Read in the binary masks or the merged Islet channel images (CHGA+SYP).
    if binary:
        isl_img_file = df["seg_subdir_Islet"]["bin_mask"].iloc[idx] / df["img_path_Islet"]["bin_mask"].iloc[idx]
        imm_img_file = df["seg_subdir_Immune"]["bin_mask"].iloc[idx] / df["img_path_Immune"]["bin_mask"].iloc[idx]
    else:
        isl_img_file = df["seg_subdir_Islet"]["image"].iloc[idx] / df["img_path_Islet"]["image"].iloc[idx]
        imm_img_file = df["seg_subdir_Immune"]["image"].iloc[idx] / df["img_path_Immune"]["image"].iloc[idx]
    isl_img = io.read_image(str(isl_img_file))
    imm_img = io.read_image(str(imm_img_file))

    ## Get the original centroids for later plotting.
    isl_regionprops["Panel"] = "Islet"
    imm_regionprops["Panel"] = "Immune"

    isl_regionprops["centroid_1_dx"] = isl_regionprops["centroid-1"]
    isl_regionprops["centroid_0_dy"] = isl_regionprops["centroid-0"]

    ## Flip the coordinates for case_id 6526
    if case_id == "6526":
        df_temp = isl_regionprops[["centroid-1", "centroid-0"]].to_numpy()
        y_max = max(df_temp[:,1])
        df_out = flip_coords(df_temp, y_max)
        isl_regionprops["centroid_1_dx"] = df_out[0]
        isl_regionprops["centroid_0_dy"] = df_out[1]

    imm_regionprops["centroid_1_dx"] = imm_regionprops["centroid-1"]
    imm_regionprops["centroid_0_dy"] = imm_regionprops["centroid-0"]

    df_regionprops = pd.concat([isl_regionprops, imm_regionprops], axis=0)
    df_original = df_regionprops.assign(panel_islet_parent=df_regionprops["Panel"] + "_" + df_regionprops["islet_parent"].astype(str))
    #------------------------------------------------------------------------------------------------------------------
    ## Pad the images.
    isl_img = isl_img[0,:,:]
    imm_img = imm_img[0,:,:]

    ## Flip the images by 270° for case_id 6526
    if case_id == "6526":
        isl_img = np.rot90(isl_img, k = 3)
        if verbose:
            print("Flipping images")

    ## uncomment for QC:
    #if (isl_img.shape[0] > imm_img.shape[0]) or (isl_img.shape[1] > imm_img.shape[1]):
    #    print("Islet image is larger than the Immune image in at least 1 dimension")
    #    print(f"Isl image: {isl_img.shape}")
    #    print(f"Imm image: {imm_img.shape}")
    #    sus.append(idx)

    isl_img = (isl_img * 255).astype('uint8')
    imm_img = (imm_img * 255).astype('uint8')

    # Pad the smaller image and mask, to make sure both images/masks have the same size
    isl_padded, imm_padded, isl_off, imm_off = pad_images(isl_img, imm_img, (0))
    
    # isl_off = (offset_y1, offset_y1+y1, offset_x1, offset_x1+x1)
    # Translate the x,y-coordinates by the image offsets introduced by the padding
    isl_regionprops["centroid_1_dx"] = isl_regionprops["centroid_1_dx"] + isl_off[2]
    isl_regionprops["centroid_0_dy"] = isl_regionprops["centroid_0_dy"] + isl_off[0]

    imm_regionprops["centroid_1_dx"] = imm_regionprops["centroid_1_dx"] + imm_off[2]
    imm_regionprops["centroid_0_dy"] = imm_regionprops["centroid_0_dy"] + imm_off[0]

    df_padded = pd.concat([isl_regionprops, imm_regionprops], axis=0)
    df_padded = df_padded.assign(panel_islet_parent=df_padded["Panel"] + "_" + df_padded["islet_parent"].astype(str))
    #------------------------------------------------------------------------------------------------------------------
    ## Perform Global Registration using phase cross-correlation with subpixel accuracy (upsample-factor)
    shift, glob_error, diffphase = phase_cross_correlation(isl_padded, imm_padded, upsample_factor=upsample_factor)  

    ## Warp
    dx = -shift[1]
    dy = -shift[0]
    affine_matrix = np.array([[1, 0, dx], [0, 1, dy]])
    if verbose:
        print("Global Registration: ")
        print(f"Translate the padded Islet image by dx = {dx}, dy = {dy} to align it with the padded Immune image")

    # Apply the translation to the Islet image -> This is the globally registered image!
    isl_img_reg_global = cv.warpAffine(isl_padded, affine_matrix, (isl_padded.shape[1], isl_padded.shape[0]))

    # Apply the registration to the padded coordinates. 
    isl_regionprops["centroid_1_dx"] = isl_regionprops["centroid_1_dx"] + dx 
    isl_regionprops["centroid_0_dy"] = isl_regionprops["centroid_0_dy"] + dy 

    df_pre_icp = pd.concat([isl_regionprops, imm_regionprops], axis=0)
    df_pre_icp = df_pre_icp.assign(panel_islet_parent=df_pre_icp["Panel"] + "_" + df_pre_icp["islet_parent"].astype(str))
    #------------------------------------------------------------------------------------------------------------------
    ## Perform ICP
    a = isl_regionprops[['centroid_1_dx', 'centroid_0_dy']].T.to_numpy()
    b = imm_regionprops[['centroid_1_dx', 'centroid_0_dy']].T.to_numpy()

    # ICP converges very rapidly, usually in first 5 iterations.
    T_opt, error_max, indices_opt = icp(a, b, max_iter = 20)

    ## Get the indices of the Immune cells that are closest to the Islet cells.
    opt_ind = indices_opt.squeeze().astype(int)
    imm_keys = imm_regionprops.iloc[opt_ind, :]["Object"].to_list()
    d = {'Object_Islet': isl_regionprops["Object"].tolist(), 'Object_Immune': imm_keys}
    inter_panel_neighbors = pd.DataFrame(d)
    inter_panel_neighbors["islet_roi_id"] = islet_roi_id
    inter_panel_neighbors["immune_roi_id"] = immune_roi_id

    ## Apply the transformation to the Islet Point Cloud -> These are the corrected Islet COORDINATES!
    src = np.array([a.T], copy=True).astype(np.float32)
    coord_icp = cv.transform(src, T_opt[0:2,:])
    
    ## Apply the transformation to the globally registered Islet image -> This is the corrected Islet IMAGE!
    isl_img_reg_icp = cv.warpAffine(isl_img_reg_global, T_opt[0:2,:], (isl_img_reg_global.shape[1], isl_img_reg_global.shape[0]))

    isl_regionprops = isl_regionprops.assign(centroid_1_dx = coord_icp[0,:,0], centroid_0_dy = coord_icp[0, :, 1])
    isl_regionprops["Panel"] = "Islet"

    # These are the corrected coordinates for both Panels.
    df_post_icp = pd.concat([isl_regionprops, imm_regionprops], axis=0)
    df_post_icp = df_post_icp.assign(panel_islet_parent=df_post_icp["Panel"] + "_" + df_post_icp["islet_parent"].astype(str))
    
    # could also return df_padded, T_opt, 
    return df_original, df_padded, df_pre_icp, df_post_icp, inter_panel_neighbors, indices_opt, isl_img, isl_img_reg_global, isl_img_reg_icp, imm_img, error_max


# ## Perform Registration and save output.
# 
# Iterate through all images perform the registration and write the outputs to csvs.
# 
# After Registration we obtain:
# - New set of cell-centroid coordiantes for the Islet and Immune panel (`df_post_icp`).
# - Inter panel neighbors (`inter_panel_neighbors`)
# We have to write them out to .csv.
# 
# ### Create Output Directories
# 
# First, create new output directories:
# - Registered cell-centroid coordinates: `data_cells/{panel}/registered_regionprops`
# - Inter-panel neighbors: `data_cells/{panel}/inter_panel_neighbors`.

# Cell data
measurements = ["inter_panel_neighbors", "registered_regionprops"]

for panel_name in panel_names:
    output_dir_cells = folders["data_cells"] / panel_name
    output_dir_cells.mkdir(exist_ok=True)
    
    for meas_type in measurements:
        meas_dir_cells = output_dir_cells / meas_type
        meas_dir_cells.mkdir(exist_ok=True)

## For QC: Store the residual errors and suspicious images for QC.
bin_arr = []; full_arr = []; idx_arr = []; sus = []

## Iterate through all images.
for idx in range(len(df)):
    idx_arr.append(idx)

    ## Register the binary segmentation masks.
    b_df_original, b_df_padded, b_df_pre_icp, b_df_post_icp, b_inter_panel_neighbors, b_indices_opt, b_isl_img, b_isl_img_reg_global, b_isl_img_reg_icp, b_imm_img, b_error_max = register_images(df, idx, sus = sus, upsample_factor=10, binary= True)
    bin_arr.append(b_error_max)

    ## Register the merged islet channel images (CHGA+SYP).
    df_original, df_padded, df_pre_icp, df_post_icp, inter_panel_neighbors, f_indices_opt, isl_img, isl_img_reg_global, isl_img_reg_icp, imm_img, error_max = register_images(df, idx, sus = sus, upsample_factor=10, binary= False)
    full_arr.append(error_max)
    
    # Choose the registration with lower residual error.
    if (b_error_max < error_max):
        lower_error = b_error_max
        indices_opt = b_indices_opt
        inter_panel_neighbors = b_inter_panel_neighbors
        df_original = b_df_original
        df_padded = b_df_padded
        df_pre_icp = b_df_pre_icp
        df_post_icp = b_df_post_icp
        isl_img_reg_global = b_isl_img_reg_global
        isl_img_reg_icp = b_isl_img_reg_icp
        imm_img = b_imm_img
        isl_img = b_isl_img
    else:
        lower_error = error_max  
    
    ## Very high error is often indicative of problematic alignment.
    ## Could think of normalizing this by the number of cells in the image.
    ## Just for QC.
    if lower_error > 1e6:
        sus.append(idx)
    
    #------------------------------------------------------------------------------------------------------------------
    ## Save the registered images.
    write_out_dfs(df, panel_names, idx, df_post_icp, inter_panel_neighbors)


## Make Dataframe:
df_error = pd.DataFrame({"idx": idx_arr, "binary": bin_arr, "full": full_arr})
df_error.astype({"idx": int, "binary": int, "full": int})

sus = np.unique(sus)
print(f"Number of suspect images: {len(sus)}")
sus

## Select images where sus matches the idx.
df_error_sus = df_error[df_error["idx"].isin(sus)]
print(len(df_error_sus))
print(df_error_sus)
print(df_error)


# # Plot images for QC:
# 
# Select a random image and plot its registraiton.

## Select random image.
# idx = 2423
idx = np.random.randint(0, len(df))
print(f"Select random image: {idx}")

b_df_original, b_df_padded, b_df_pre_icp, b_df_post_icp, b_inter_panel_neighbors, b_indices_opt, b_isl_img_reg_global, b_isl_img_reg_icp, b_imm_img, b_isl_img, b_error_max = register_images(df, idx, upsample_factor=10, binary= True)

df_original, df_padded, df_pre_icp, df_post_icp, inter_panel_neighbors, indices_opt, isl_img_reg_global, isl_img_reg_icp, imm_img, isl_img, error_max = register_images(df, idx, upsample_factor=10, binary= False)

# Choose the one with the lower error.
if (b_error_max < error_max):
    lower_error = b_error_max
    indices_opt = b_indices_opt
    inter_panel_neighbors = b_inter_panel_neighbors
    df_original = b_df_original
    df_padded = b_df_padded
    df_pre_icp = b_df_pre_icp
    df_post_icp = b_df_post_icp
    isl_img_reg_global = b_isl_img_reg_global
    isl_img_reg_icp = b_isl_img_reg_icp
    imm_img = b_imm_img
    isl_img = b_isl_img
else:
    lower_error = error_max  


# ### Check Global Registration Visually
# 
# Use the returned dataframes containing the coordinates to observe the registration.

# Subplot 1
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
fig, axes = plt.subplots(2, 2, figsize=(10, 10))

# Plot all point clouds.
plot_centroids(axes[0, 0], df_original, "Original")
plot_centroids(axes[0, 1], df_padded, title="After padding")
# Plot corrected_img to axes[0,2]
plot_centroids(axes[1, 0], df_pre_icp, "After Phase-correlation")
plot_centroids(axes[1, 1], df_post_icp, "After ICP")

plt.tight_layout()
plt.show()

## Plot the registered images.
fig, axes = plt.subplots(3, 2, figsize=(10, 10))

# First row: Original images
axes[0, 0].imshow(isl_img, vmin=0, vmax=150)
axes[0, 0].set_title("Islet Original")
axes[0, 1].imshow(imm_img, vmin=0, vmax=150)
axes[0, 1].set_title("Immune Original")

# Second row: Globally registered images
axes[1, 0].imshow(isl_img_reg_global, vmin=0, vmax=150)
axes[1, 0].set_title("Islet Global Registration")
axes[1, 1].imshow(imm_img, vmin=0, vmax=150)
axes[1, 1].set_title("Immune")

# Third row: ICP registered images
axes[2, 0].imshow(isl_img_reg_icp, vmin=0, vmax=150)
axes[2, 0].set_title("Islet ICP Registration")
axes[2, 1].imshow(imm_img, vmin=0, vmax=150)
axes[2, 1].set_title("Immune")

plt.tight_layout()
plt.show()

# Use idx = 5405 as image for edge case (negative coordinates after translation)

## Write success log-message to 04_Measurements.out
base_dir = Path("/home/processing/")
with open(base_dir / "txt_output" / "05_ImageRegistration.out", "w") as f:
    f.write("05_ImageRegistration.py completed successfully!")
