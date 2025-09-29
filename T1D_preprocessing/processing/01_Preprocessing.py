#!/usr/bin/env python
# coding: utf-8

# # **IMC data processing pipeline**
# 
# This pipeline has been specifically developed for the Imaging Mass Cytometry (IMC) - Type 1 Diabetes (T1D) project.
# 
# ## **Introduction**
# 
# This pipeline extracts image data from Imaging Mass Cytometry aquisitions, performs islet- and cell-level image segmentation and extracts measurements from the segmented objects.  
# This pipeline is designed to work with two antibody panels applied to two consecutive tissue sections.
# 
# As input, the user should provide zipped folders containing IMC acquisition (one `.mcd` file with the associated `.txt` files), and a panel file (`panel.csv`) for each antibody panel that indicates the channels that were measured and the channels that should be used for segmentation. Detailed information about zipped folders and panel files can be found below.
# 
# This pipeline is based on functions from the [steinbock package](https://github.com/BodenmillerGroup/steinbock), full steinbock documentation can be found here: https://bodenmillergroup.github.io/steinbock.
# 
# 
# ### **Steps**
# 
# The full pipeline contains four notebooks that should be run sequentially:
#   
# **1. Preprocessing** *(current notebook)*
# - Process zipped folders.
# - Extract images from IMC acquisitions.
# - Catch unmatched images.
#   
# **2. Islet segmentation**
# - Prepare islet segmentation.
# - Segment islets.
# - Merge islets and match islets from consecutive sections.
# - Perform visual quality checks.
#   
# **3. Cell segmentation**
# - Prepare cell segmentation.
# - Segment cells.
#   
# **4. Measurements**
# - Measure cell intensities.
# - Measure region properties and cell distances to islets.
# - Measure cell neighbors.
# - Catch unmatched data files.

# # **Preprocessing**
# 
# ## **Configuration**
# 
# ### **Import packages**

import pandas as pd
import pickle
import re
import sys
from pathlib import Path
from steinbock import io
from steinbock.preprocessing import imc


print(sys.path)
print(sys.executable)


# ### Helper functions

# Helper functions to get unique elements from 2 lists with elements of list 1 having .tiff suffixes.
def get_unique_elements(list1, list2):
    # Include some " - split here"
    return [element for element in list1 if not any(element.replace(" - split.tiff", "").replace(".tiff", "") in item for item in list2)]

# Helper function to get all duplicates in 2 lists.
def get_duplicates(lst):
    unique_elements = set()
    duplicates = []
    for item in lst:
        if item in unique_elements:
            duplicates.append(item)
        else:
            unique_elements.add(item)
    
    return duplicates


# ### **Define input and output directories**
# 
# *Manual step:* enter the path to the directory where the data will be saved (named `folder_data` from here on).

# Data folder
folder_data = Path("/home/processing/")
Path(folder_data).mkdir(parents=True, exist_ok=True)
assert Path.exists(folder_data), f"{folder_data} does not exist"
print("Data folder:", folder_data)

# Git folder (folder containing the current notebook)
folder_git = Path.cwd()
assert Path.exists(folder_git), f"{folder_git} does not exist"
print("Git folder:", folder_git)


# #### **Create folders for intermediate processing steps**
# - `raw`: should contain user-provided zipped `.mcd` and `.txt` acquisitions.
# - `img`: store extracted images in `.tiff` format.
# - `seg_cells`: store image stacks for cell segmentation.
# - `seg_islets`: store image stacks for islet segmentation.
# - `masks_cells`: store cell segmentation masks.
# - `masks_islets`: store islet segmentation masks.
# - `data_cells`: store generated single cell-level data.
# - `data_islets`: store generated islet-level data.

folders = {
    "raw": folder_data / "raw",
    "img": folder_data / "img",
    "seg_cells": folder_data / "seg_cells",
    "seg_islets": folder_data / "seg_islets",
    "masks_cells": folder_data / "masks_cells",
    "masks_islets": folder_data / "masks_islets",
    "data_cells": folder_data / "data_cells",
    "data_islets": folder_data / "data_islets",
    "variables": folder_data / "variables"
}
# Make directories (if they do not exist)
for folder in folders.values():
    folder.mkdir(exist_ok=True)
    
# Add base previously defined data and git folders
folders["data"] = folder_data
folders["git"] = folder_git

# Export folder names for use in downstream notebooks
with open(folder_data / "variables" / "folders.txt", "wb") as handle:
    pickle.dump(folders, handle)


# ### **Antibody panels**
# 
# Panel files are user-provided `.csv` files that should be located in `folder_data` and contain the following columns:
# - `channel`: unique channel ID.
# - `metal`: metal isotope to which the antibody is conjugated (format: `Nd144`, `Ir191`).
# - `name`: name of the target marker.
# - `keep` should the channel be retained for processing and analysis? (1 = yes, 0 = no).
# - `deepcell` should the channel be used for cell segmentation and in which compartment the marker is expressed? (1 = nucleus, 2 = membrane, empty if the channel should not be used).
# - `isletseg` should the channel be used for islet segmentation? (1 = yes, empty = channel not used).

# Columns required in the panel file(s)
panel_cols = {
    "col_channel": "channel",
    "col_metal": "metal",
    "col_name": "name",
    "col_keep": "keep",
    "col_deeepcell": "deepcell",
    "col_islet_seg": "isletseg"
}


# *Manual step:* adapt the panel names and panel file names if needed. 

# List panel files
panels = {
    "Islet": (folder_data / 'panel_Islet.csv'),
    "Immune": (folder_data / 'panel_Immune.csv')
}


# #### **Load and display the panels**
# Loop through the panels
for panel_name, panel_path in panels.items():
    print("Panel:", panel_name)
    
    # Load the panel file
    assert Path.exists(panel_path), f"{panel_path} does not exist"
    cur_panel = pd.read_csv(panel_path, sep = ',', index_col = False)

    # Make sure that the required columns exist
    for col in panel_cols.values():
        assert(col in cur_panel.columns), f"Column {col} missing from panel"
    
    # Subset the panel
    cur_panel = cur_panel[cur_panel[panel_cols["col_keep"]]==1]
    panels[panel_name] = cur_panel
    
    # Display the panel
    print(panels[panel_name].head())
    
# Export the panels for use in downstream scripts
with open(folder_data / "variables" / "panels.txt", "wb") as handle:
     pickle.dump(panels, handle)


# ## **Process zipped folders**
# 
# IMC acquisitions generate `.mcd` and `.txt` files. Each acquisition session (corresponding ot one `.mcd` file) should be zipped in a folder containing:
# - The `.mcd` file.
# - All the associated `.txt` files generated during the acquisition (do not change any of the file names).
# The `.txt` files are used as a backup in case the data cannot be extracted from the `.mcd` file.
# 
# All the zipped folders should be stored in subfolders of the `raw` folder (in the `folder_data` directory). The subfolders should be named exactly like the panels in `panels` (see "List panel files" above).  
# 
# For the current dataset, the folder structure is the following, with zipped MCD and TXT files stored in `raw/Immune` and `raw/Islet`:
#folder_data
#    |_ data_cells
#|_ data_islets
#|_ img
#|_ masks_cells
#|_ masks_islets
#|_ raw
#    |_ Immune <- ZIP files from the Immune panel stored here
#    |_ Islet  <- ZIP files from the Islet panel stored here
#|_ seg_cells
#|
#|_ panel_Immune.csv <- Panel file (Immune panel)
#|_ panel_Islet.csv  <- Panel file (Islet panel)
####|_ seg_islets
# ### **List `.zip` folders**
# *Manual step:* define a regular expression to identify the naming scheme of `.zip` files.

# Part that all zipped files need to have in common
file_regex = '(?P<caseid>[0-9]{3,4})_(?P<panel>[a-zA-Z0-9]{5,6})*'
# List all zip folders that match the regular expression.

# List zip folders. This should 190 folders.
re_fn = re.compile(file_regex)
zip_folders = [f for f in folders["raw"].rglob("*") if 
               re.match(file_regex, f.name)]
print("\nTotal number of '.zip' folders:", len(zip_folders))
print("Zipped folders:")
for file in zip_folders:
    print(file.name)

# List all case IDs and panels
case_list = []
panel_list = []

for file in zip_folders:
    case_list.append(re_fn.search(file.name).group("caseid"))
    panel_list.append(re_fn.search(file.name).group("panel"))

# to remove pot. duplicates.
case_list = list(set(case_list))
panel_list = list(set(panel_list))

# Generate a table with case IDs as indexes and panels as columns
zip_table = pd.DataFrame(dtype=str, columns=panel_list, index=case_list)

for file in zip_folders:
    cur_case = re_fn.search(file.name).group("caseid")
    cur_panel = re_fn.search(file.name).group("panel")
    zip_table.loc[cur_case, cur_panel] = file.name
zip_table


# ## **Extract images from IMC acquisitions**
# 
# Here, images are extracted from raw IMC files and saved in the `img` folder. Each image corresponds to one acquisition in one file, with the image channels filtered (`keep` column in antibody panel) and sorted according to the the panel file.  
# 
# In case an `.mcd` file is corrupted, the steinbock function tries to extract missing acquisitions from matching `.txt` files. In a second step, images from unmatched `.txt` files are extracted as well.  
# 
# See the full documentation here: https://bodenmillergroup.github.io/steinbock/latest/cli/preprocessing/#image-conversion
# 
# ### **Settings**
# After image extraction, hot pixel filtering is performed using the threshold defined by the `hot_pixel_filtering` variable.
hot_pixel_filtering = 50


# ### **Fix mismatched regions**
# Due to an issue with region selection, the regions of interest (ROI) acquired on consecutive sections were mismatched (i.e., had different acquisition numbers) for a three of the 95 measured samples.  
# Here, we use manually-generated dictionaries to re-map the ROI numbers of the *Immune* panel, so that they correspond to the numbers of the *Islet* panel.
mismatched_cases = ("6043", "6048", "6429")
map_dict = [
    {"name": "6043",
     "001": "074", "002": "055", "003": "059", "004": "035", "005": "003", "006": "023", "007": "010", "008": "008", "009": "021", "010": "052",
     "011": "049", "012": "012", "013": "053", "014": "028", "015": "018", "016": "048", "017": "017", "018": "026", "019": "065", "020": "046", 
     "021": "068", "022": "016", "023": "054", "024": "063", "025": "025", "026": "073", "027": "042", "028": "007", "029": "044", "030": "027",
     "031": "005", "032": "031", "033": "014", "034": "020", "035": "064", "036": "036", "037": "039", "038": "076", "039": "061", "040": "050", 
     "041": "032", "042": "066", "043": "029", "044": "037", "045": "070", "046": "019", "047": "047", "048": "033", "049": "009", "050": "067", 
     "051": "051", "052": "011", "053": "006", "054": "043", "055": "001", "056": "022", "057": "034", "058": "062", "059": "045", "060": "058", 
     "061": "060", "062": "002", "063": "071", "064": "038", "065": "075", "066": "041", "067": "056", "068": "072", "069": "069", "070": "040", 
     "071": "013", "072": "057", "073": "004", "074": "030", "075": "015", "076": "024"},
    {"name": "6048",
     "001": "001", "002": "067", "003": "034", "004": "052", "005": "005", "006": "030", "007": "055", "008": "007", "009": "063", "010": "010",
     "011": "011", "012": "025", "013": "013", "014": "014", "015": "015", "016": "019", "017": "021", "018": "058", "019": "020", "020": "032", 
     "021": "045", "022": "002", "023": "023", "024": "024", "025": "047", "026": "018", "027": "027", "028": "028", "029": "056", "030": "008", 
     "031": "078", "032": "057", "033": "068", "034": "022", "035": "035", "036": "036", "037": "037", "038": "051", "039": "004", "040": "040", 
     "041": "033", "042": "042", "043": "043", "044": "050", "045": "054", "046": "046", "047": "075", "048": "048", "049": "049", "050": "016", 
     "051": "060", "052": "029", "053": "053", "054": "039", "055": "009", "056": "006", "057": "003", "058": "038", "059": "059", "060": "017", 
     "061": "061", "062": "062", "063": "065", "064": "064", "065": "071", "066": "066", "067": "069", "068": "073", "069": "012", "070": "070", 
     "071": "026", "072": "072", "073": "044", "074": "074", "075": "031", "076": "076", "077": "077", "078": "041"},
    {"name": "6429",
     "001": "023", "002": "056", "003": "017", "004": "040", "005": "013", "006": "053", "007": "079", "008": "073", "009": "036", "010": "051", 
     "011": "021", "012": "068", "013": "070", "014": "081", "015": "077", "016": "066", "017": "054", "018": "072", "019": "041", "020": "047", 
     "021": "033", "022": "008", "023": "063", "024": "042", "025": "014", "026": "022", "027": "045", "028": "061", "029": "009", "030": "004", 
     "031": "059", "032": "052", "033": "067", "034": "035", "035": "028", "036": "024", "037": "032", "038": "071", "039": "062", "040": "055",
     "041": "030", "042": "038", "043": "050", "044": "043", "045": "007", "046": "074", "047": "019", "048": "039", "049": "020", "050": "064", 
     "051": "058", "052": "076", "053": "002", "054": "075", "055": "011", "056": "029", "057": "044", "058": "001", "059": "012", "060": "069", 
     "061": "026", "062": "078", "063": "057", "064": "031", "065": "049", "066": "010", "067": "005", "068": "018", "069": "046", "070": "003",
     "071": "037", "072": "080", "073": "027", "074": "065", "075": "048", "076": "015", "077": "016", "078": "025", "079": "006", "080": "060", 
     "081": "034"}
]


# ### **Image conversion**
# Extract image stacks from IMC acquisitions. 
# Image and acquisition metadata are exported to `folder_data` as `images.csv`.
panels["Islet"]

for panel_name, panel in panels.items():
    print("Processing", panel_name, "panel")
    
    # Input and output folders
    image_info = []
    raw_subdir = folders["raw"] / panel_name
    img_subdir = folders["img"] / panel_name
    img_subdir.mkdir(exist_ok = True)  
    
    # List zipped files
    cur_mcd_files = imc.list_mcd_files(raw_subdir, unzip=True)
    cur_txt_files = imc.list_txt_files(raw_subdir, unzip=True)
    
    # Process files
    for (mcd_file, acquisition, img, matched_txt, recovered) in \
    imc.try_preprocess_images_from_disk(
        cur_mcd_files, cur_txt_files,
        hpf = hot_pixel_filtering,
        channel_names = panels[panel_name]["metal"],
        unzip = True
    ):
        cur_desc = acquisition.description
        cur_case = re_fn.search(mcd_file.name).group("caseid")
        
        # Renumber mismatched images
        if (panel_name == "Immune" and cur_case in mismatched_cases):
            cur_dict = next(item for item in map_dict if item["name"] == cur_case)
            cur_desc = cur_desc[:-3] + cur_dict.get(cur_desc[-3:])
            
        img_file = f"{mcd_file.stem}_{cur_desc}.tiff"
        io.write_image(img, img_subdir / img_file)

        # Save acquisition metadata
        image_info_row = imc.create_image_info(
            mcd_file, acquisition, img, matched_txt, recovered, img_file
        )
        
        # Renumber mismatched images in image metadta file
        if (panel_name == "Immune" and cur_case in mismatched_cases):
            image_info_row["acquisition_description"] = cur_desc
        
        image_info_row["panel"] = panel_name
        image_info.append(image_info_row)

    image_info = pd.DataFrame(image_info)
    image_meta_file = f"images_{panel_name}.csv"
    image_info.to_csv(folders["data"] / image_meta_file, index = False)


# ## **Catch unmatched images**
# 
# ### **Flag missing images**
# The ablated regions should be the same on all consecutive sections. Here, we attempt to match images from different panels, based on the ROI number. Images from one panel that do not have a corresponding image in the other panel(s) are flagged.

panel_names = list(panels.keys())
missing = set()

# List files for the first panel
images_panel0 = sorted([img.name.replace(panel_names[0], "") \
                        for img in Path.iterdir(folders["img"] / panel_names[0])])
images_panel0 = frozenset(images_panel0)

# Find matched images in the other panels
for panel_name in panel_names[1:]:
    cur_images = [img.name for img in Path.iterdir(folders["img"] / panel_name)]
    cur_list = set([img.replace(panel_name, "") for \
                    img in cur_images])
    
    missing.add(frozenset(images_panel0.difference(cur_list)))
    missing.add(frozenset(cur_list.difference(images_panel0)))

# Print out all missing images
missing = [list(x) for x in missing]
missing = sorted([x for xs in missing for x in xs])
print("Images with missing corresponding images (", len(missing),
      "missing images ):\n", missing)


# ### **Rescue split images**
# 
# Sometimes the laser stopped to ablate during acquisition. After restart, two split .txt-files were outputted. Here, split acquisitions which fully contain the Islet are retained (`keep_files`). 
# 
# Acquisitions, which splitted the islet are deleted below.

# These acquisitions were split but will be retained as Islet is visible in acquisitions.
keep_files = ['6547_Immune_ROI_046','6520_Immune_ROI_032', '6289_Immune_ROI_034', \
                '6147_Immune_ROI_011', '6538_Islet_ROI_055']

it = []

# Iterate through panels.
for panel_name in panel_names:
    cur_images = [img for img in Path.iterdir(folders["img"] / panel_name)]

    # Paths of split top and split bottom .tiffs
    # Example: split bottom: 6520_Immune_ROI_032 - split.tiff
    # Example: split top: 6520_Immune_ROI_032.tiff
    spl_top = [img for img in cur_images for file in keep_files if file in str(img.with_suffix("").name)]
    spl_bottom = [Path(re.sub(r" - split.tiff", ".tiff", string = str(x))) for x in spl_top]
    dict_paths = dict(zip(spl_top, spl_bottom))
    
    # used to update "missing" list 
    not_missing = [im.with_suffix("").name.replace(("_" + panel_name + "_"), "__") for im in spl_bottom]
    it.extend(get_unique_elements(missing, not_missing))
    
    # Delete split .tiff not containing Islet.
    # Replace file name if "- split.tiff" contains file name.
    for key, value in dict_paths.items():
        # Keep bottom split file.
        key.replace(value)
            
# Updated missing list
missing = get_duplicates(it)
print("Images with missing corresponding images (", len(missing),
      "updated missing images ):\n", missing)

# List files for the Islet panel
images_panel0 = sorted([img.name.replace(panel_names[0], "") \
                        for img in Path.iterdir(folders["img"] / panel_names[0])])
images_panel0 = frozenset(images_panel0)
images_panel0

# List files for the Immune panel
for panel_name in panel_names[1:]:
    cur_images = [img.name for img in Path.iterdir(folders["img"] / panel_name)]
    cur_list = set([img.replace(panel_name, "") for img in cur_images])
    print(cur_list)


# ### **Delete unmatched images**
# 
# Images that do not have a matching image in all the other panels are deleted.
delete_unmatched_images = True

if missing and delete_unmatched_images:
    for panel_name in panel_names:
        cur_dir = folders["img"] / panel_name
        unmatched_images = [
            cur_dir / im.replace("__", ("_" + panel_name + "_")) \
            for im in missing]
        
        for image in unmatched_images:
            print(f"Deleting {image}")
            Path.unlink(image, missing_ok=True)


# ## **Next steps**
# 
# The next step in this pipeline is islet segmentation, which is performed with the `02_IsletSegmentation.ipynb` notebook.
# get_ipython().system('conda list')
## Write success log-message to 01_Preprocessing.out
base_dir = Path("/home/processing/")
with open(base_dir / "txt_output" / "01_Preprocessing.out", "w") as f:
    f.write("01_Preprocessing.py completed successfully!")