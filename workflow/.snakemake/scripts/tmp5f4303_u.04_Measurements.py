
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/mnt/snakemake', '/scratch/nsteen/T1D_preprocessing/processing']); import pickle; snakemake = pickle.loads(b'\x80\x04\x95\x01\x05\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94\x8c?/scratch/nsteen/T1D_preprocessing/processing/04_Measurements.py\x94a}\x94(\x8c\x06_names\x94}\x94\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x10\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x16)}\x94\x8c\x05_name\x94h\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94\x8c9/scratch/nsteen/processing/txt_output/04_Measurements.out\x94a}\x94(h\x0c}\x94h\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bub\x8c\x07threads\x94K\x08\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x08K\x01J\xc0\xd4\x01\x00M\xe8\x03\x8c\r/sctmp/nsteen\x94\x8c\x0810:00:00\x94e}\x94(h\x0c}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06mem_mb\x94K\x02N\x86\x94\x8c\x07disk_mb\x94K\x03N\x86\x94\x8c\x06tmpdir\x94K\x04N\x86\x94\x8c\x04time\x94K\x05N\x86\x94uh\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bhWK\x08hYK\x01h[J\xc0\xd4\x01\x00h]M\xe8\x03h_hS\x8c\x04time\x94hTub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bub\x8c\x06config\x94}\x94(\x8c\x08data_dir\x94\x8c\x10/scratch/nsteen/\x94\x8c\x08home_dir\x94\x8c"/home/nsteen/scripts/T1D_analysis/\x94\x8c\rcontainer_dir\x94\x8c&/home/nsteen/data/singularity/sandbox/\x94u\x8c\x04rule\x94\x8c\x12Measurements_04_py\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8c,/scratch/nsteen/T1D_preprocessing/processing\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = True; __real_file__ = __file__; __file__ = '/scratch/nsteen/T1D_preprocessing/processing/04_Measurements.py';
######## snakemake preamble end #########
#!/usr/bin/env python
# coding: utf-8

# # **Measurements**
# 
# ## **Introduction**
# 
# This is the fourth and last script in the processing pipeline for IMC data.
# 
# The goal is to extract measurements from the multichannel images generated in the first script using the islet and cell masks generated in the second and third scripts. This is performed using functions from the `steinbock` package. Documentation can be found here: https://bodenmillergroup.github.io/steinbock/latest/cli/measurement.
# 
# The following measurements are performed: 
# 
# **Measure intensities**  
# - Average marker intensities for islets.
# - Average marker intensities for single cells.
# 
# **Measure region properties**  
# - Islet-level spatial measurements, such as area and position.
# - Cell-level spatial measurements.
# - Distance of cells to segmented islets.
# 
# **Cell neighbors**  
# - Which cells are neighboring each other.
# 
# In case one of the three data files is not generated (for instance due to a missing mask), the corresponding files in other data folders are deleted at the end of this script.

# ## **Configuration**
# 
# ### **Import packages**

import logging
import numpy as np
import pandas as pd
from pathlib import Path
import pickle
import sys

from scipy import ndimage as ndi
from skimage import measure
from skimage.segmentation import expand_labels
from skimage.util import invert

from steinbock import io
from steinbock.measurement import intensities, neighbors, regionprops

logger = logging.getLogger(__name__)
print(sys.path)
print(sys.executable)
base_dir = Path("/home/processing/")

# ### **Load directories and panels**
# 
# Paths to input and output folders as well as antibody panels were exported by the first script (`01_Preprocessing.ipynb`). Here they are imported again.

with open(base_dir / "variables/folders.txt", "rb") as handle:
    folders = pickle.loads(handle.read())
folders

with open(base_dir / "variables/panels.txt", "rb") as handle:
    panels = pickle.loads(handle.read())

for panel_name, panel in panels.items():
    print(panel_name, "\n", panel.head())
panel_names = list(panels.keys())

# #### **Create output directories**

# Cell data
measurements = ["intensities", "regionprops", "neighbors"]

for panel_name in panel_names:
    output_dir_cells = folders["data_cells"] / panel_name
    output_dir_cells.mkdir(exist_ok=True)
    
    for meas_type in measurements:
        meas_dir_cells = output_dir_cells / meas_type
        meas_dir_cells.mkdir(exist_ok=True)
        
# Islet data
measurements_islets = ["intensities", "regionprops"]

for panel_name in panel_names:
    output_dir_islets = folders["data_islets"] / panel_name
    output_dir_islets.mkdir(exist_ok=True)
    
    for meas_type in measurements_islets:
        meas_dir_islets = output_dir_islets / meas_type
        meas_dir_islets.mkdir(exist_ok=True)        

# #### **Select the cell segmentation type to use**

segmentation_type = "whole-cell"

# ## **Measure intensities**
# 
# Here, the mean marker expression over the cell area, respectively the islet area, is measured.
# 
# Full documentation: https://bodenmillergroup.github.io/steinbock/latest/cli/measurement/#object-intensities.
# 
# ### **Islet intensities**

for panel_name in panel_names:
    img_subdir = folders["img"] / panel_name
    masks_subdir = folders["masks_islets"] / panel_name
    intensities_dir = folders["data_islets"] / panel_name / "intensities"
    
    for img_path, mask_path, intens in intensities.try_measure_intensities_from_disk(
        img_files = io.list_image_files(img_subdir),
        mask_files = io.list_image_files(masks_subdir),
        channel_names = panels[panel_name]["name"],
        intensity_aggregation = intensities.IntensityAggregation.MEAN
    ):
        intensities_file = img_path.name.replace('.tiff', '.csv')
        pd.DataFrame.to_csv(intens, Path(intensities_dir) / intensities_file)

# ### **Single cell intensities**

for panel_name in panel_names:
    img_subdir = folders["img"] / panel_name
    masks_subdir = folders["masks_cells"] / panel_name / segmentation_type
    intensities_dir = folders["data_cells"] / panel_name / "intensities"
    
    for img_path, mask_path, intens in intensities.try_measure_intensities_from_disk(
        img_files = io.list_image_files(img_subdir),
        mask_files = io.list_image_files(masks_subdir),
        channel_names = panels[panel_name]["name"],
        intensity_aggregation = intensities.IntensityAggregation.MEAN
    ):
        intensities_file = img_path.name.replace('.tiff', '.csv')
        pd.DataFrame.to_csv(intens, Path(intensities_dir) / intensities_file)

# ## **Measure region properties**
# 
# Documentation for region properties measurements: https://bodenmillergroup.github.io/steinbock/latest/cli/measurement/#region-properties.
# 
# ### **Properties to measure**
# 
# For a full list of measurable properties, refer to https://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.regionprops.

skimage_regionprops = [
        "area",
        "centroid",
        "major_axis_length",
        "minor_axis_length",
        "eccentricity",
    ]

# ### **Islet-level spatial measurements**
# 
# Measurement of spatial properties for islets.  
# 

for panel_name in panel_names:
    img_subdir = folders["img"] / panel_name
    islet_masks_subdir = folders["masks_islets"] / panel_name
    regions_dir = folders["data_islets"] / panel_name / "regionprops"
    
    # Measure cell regions props
    for img_path, mask_path, region_props in regionprops.try_measure_regionprops_from_disk(
        img_files = io.list_image_files(img_subdir),
        mask_files = io.list_image_files(islet_masks_subdir),
        skimage_regionprops = skimage_regionprops
    ):
        regprop_file = img_path.name.replace('.tiff', '.csv')
        pd.DataFrame.to_csv(region_props, Path(regions_dir) / regprop_file)

# ### **Cell-level spatial measurements**
# 
# Measurement of spatial properties for single cells.  
# 
# In addition, a binary transform function is applied to islet masks so that the transformed masks contain the distance of each pixel to the islet edge. The distance of a single cell to the islet is defined as the median distance of all the pixels that compose the cell. Cells located within islets are attributed positive distances while those located outside are attributed negative distances. The resulting distances are stored in the `distance_to_islet` column.
# 
# The `islet_closet` value defines which islet is located the closest to the considered cell, whereas the `islet_parent` value defines in which islet the considered cell is located (`0` if the cell is not located in any islet.

def transform_islet_masks(islet_mask):
    
     # Invert the islet mask (to measure cells located outside islets)
    inverted_mask = invert(islet_mask)
    inverted_mask[inverted_mask < 0] = 0

    # Apply binary transformation to the mask
    positive_mask = ndi.distance_transform_edt(islet_mask)
    negative_mask = ndi.distance_transform_edt(inverted_mask)
    
    # Store the masks in an array
    transformed_masks = np.zeros((2, islet_mask.shape[1], islet_mask.shape[2]))
    transformed_masks[0,...] = positive_mask
    transformed_masks[1,...] = negative_mask
    
    return transformed_masks

for panel_name in panel_names:
    img_subdir = folders["img"] / panel_name
    cell_masks_subdir = folders["masks_cells"] / panel_name / segmentation_type
    islet_masks_subdir = folders["masks_islets"] / panel_name
    regions_dir = folders["data_cells"] / panel_name / "regionprops"
    
    img_files = io.list_image_files(img_subdir)
    
    for img_file in img_files:
        try:
            # Load images and masks
            cell_mask_file = cell_masks_subdir / img_file.name
            islet_mask_file = islet_masks_subdir / img_file.name
            
            if (cell_mask_file.exists() and islet_mask_file.exists()):
                img = io.read_image(img_file)
                cell_mask = io.read_mask(cell_mask_file)
                islet_mask = io.read_image(islet_mask_file) # read as an image, not a mask

                # Measure cell-level region props
                region_props = regionprops.measure_regionprops(img, cell_mask, skimage_regionprops)

                # Apply binary transformation to the islet mask
                transformed_mask = transform_islet_masks(islet_mask)
                
                # Measure cell-to-islet distances
                islet_distances = intensities.measure_intensites(
                    transformed_mask, cell_mask,
                    channel_names = (["DistPOS", "DistNEG"]),
                    intensity_aggregation = intensities.IntensityAggregation.MEDIAN
                )
                islet_distances["distance_to_islet"] = np.where(islet_distances["DistPOS"] > 0,
                                                                islet_distances["DistPOS"],
                                                                -islet_distances["DistNEG"])
                islet_distances.drop(columns=["DistPOS", "DistNEG"], inplace=True)
                
                # Find closest islet
                expanded_mask = expand_labels(islet_mask, np.amax(islet_mask.shape))
                
                islet_closest = intensities.measure_intensites(
                    expanded_mask, cell_mask,
                    channel_names = (["islet_closest"]),
                    intensity_aggregation = intensities.IntensityAggregation.MIN
                )
                islet_closest["islet_closest"] = islet_closest["islet_closest"].astype("int")

                # Concatenate region props and islet distances
                region_props = pd.concat([region_props, islet_distances, islet_closest], axis=1)
                
                # Find parent islet
                region_props["islet_parent"] = np.where(region_props["distance_to_islet"] > 0,
                                                         region_props["islet_closest"], 0)

                # Save measurements as CSV files
                regprop_file = img_file.name.replace(".tiff", ".csv")
                pd.DataFrame.to_csv(region_props, Path(regions_dir) / regprop_file)
            
        except:
            logger.exception(f"Error measuring regionprops in {img_file}")    

# ## **Measure cell neighbors**
# 
# Cell located next to each other are recorded here. Several options are available (see below), but here we use the eculidean distance between cell borders to define neighbors.
# 
# Documentation: https://bodenmillergroup.github.io/steinbock/latest/cli/measurement/#object-neighbors.
# 
# ### **Settings**
# 
# **Neighborhood types:**
# + `NeighborhoodType.CENTROID_DISTANCE`
# + `NeighborhoodType.EUCLIDEAN_BORDER_DISTANCE`
# + `NeighborhoodType.EUCLIDEAN_PIXEL_EXPANSION`
# 
# **Thresholding:**
# + `dmax` (max distance between centroids)
# + `kmax` (k-nearest neighbors)

neighborhood_type = neighbors.NeighborhoodType.EUCLIDEAN_BORDER_DISTANCE
dmax = 10
kmax = 5

# ### **Measure neighbors**

for panel_name in panel_names:
    img_subdir = folders["img"] / panel_name
    masks_subdir = folders["masks_cells"] / panel_name / segmentation_type
    neighbors_dir = folders["data_cells"] / panel_name / "neighbors"

    for mask_path, neighb in neighbors.try_measure_neighbors_from_disk(
        mask_files = io.list_image_files(masks_subdir),
        neighborhood_type = neighborhood_type,
        metric = "euclidean",
        dmax = dmax,
        kmax = kmax
    ):
        neighb_file = mask_path.name.replace(".tiff", ".csv")
        neighb_file = f"{neighb_file}"
        pd.DataFrame.to_csv(neighb, Path(neighbors_dir) / neighb_file,
                            index=False)

# ## **Catch unmatched data files**
# 
# ### **Flag and delete unmatched data files**
# For each image, three data files should be generated corresponding to intensities, region props, and neighbors. If for one image one of these files is missing, the other ones are removed in order to avoid conflicts when importing data into R. Data files that do not have a matching file in every measurement folders are deleted.

delete_unmatched_files = True

for panel_name in panel_names:
    missing = set()
    
    # List intensity files
    intensity_dir = folders["data_cells"] / panel_name / "intensities"
    intensity_files = Path(intensity_dir).rglob("[!.]*.csv")
    intensity_files = frozenset(file.name for file in intensity_files)

    # Find matched data files in the other cell data folders
    for meas_type in ["regionprops", "neighbors"]:
        cur_dir = folders["data_cells"] / panel_name / meas_type
        cur_files = set([file.name for file in Path.iterdir(cur_dir)])

        missing.add(frozenset(intensity_files.difference(cur_files)))
        missing.add(frozenset(cur_files.difference(intensity_files)))
        
    # Find matched data files in the islet data folders
    for meas_type in ["intensities", "regionprops"]:
        cur_dir = folders["data_islets"] / panel_name / meas_type
        cur_files = set([file.name for file in Path.iterdir(cur_dir)])

        missing.add(frozenset(intensity_files.difference(cur_files)))
        missing.add(frozenset(cur_files.difference(intensity_files)))
        
    # Print out all missing images
    missing = [list(x) for x in missing]
    missing = [x for xs in missing for x in xs]
    print("Images with missing corresponding files:\n", missing)
    
    # Delete unmatched data files
    if delete_unmatched_files:
        unmatched_files = []
        for meas_type in ["intensities", "regionprops", "neighbors"]:
            cur_dir = folders["data_cells"] / panel_name / meas_type
            unmatched_files.extend([cur_dir / file for file in missing])

        print("\nDeleted files:")
        for file in unmatched_files:
            if file.is_file():
                print(file)
                Path.unlink(file, missing_ok=True)

## Write success log-message to 04_Measurements.out
with open(base_dir / "txt_output" / "04_Measurements.out", "w") as f:
    f.write("04_Measurements.py completed successfully!")

# ## **Next step**
# 
# This notebook is the last one in this processing pipeline. The next step is to load the measurements extracted here in R for data analysis. All data analysis for the current project is performed using the [T1D_analysis](https://github.com/BodenmillerGroup/T1D_analysis) repository.