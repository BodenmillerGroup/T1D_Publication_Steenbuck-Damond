import os
import shutil
from ij import IJ
from register_virtual_stack import Transform_Virtual_Stack_MT


# Parameters

## Change the input directory ("BASE" folder), and the case and panel lists if needed
input_dir = "R:\\Data3\\acquisition\\Batch4\\"
caseList = ("6036", "6055", "6090", "6147", "6225", "6228",
			"6303", "6321", "6388", "6396", "6414", "6421", "6428",
			"6437", "6458", "6510", "6519", "6521", "6522", "6526",
			"6532", "6547", "6550", "6553", "6558", "6562", "6563") # Use this format to process multiple cases: ("6289", "6310", "6328") # Use this format to process a single case: (["6450"])
panels = (["Immune"])  # Use this format to process multiple panels (do not include the base panel): ("Islet2", "Lympho", "Myelo")  # Use this format to process a single panel: (["Immune"])
panel_base = "Islet" # Base panel (used to define the ROIs)
img_ext = ".tif" # Image extension

###############

# MAIN

for case in caseList:
	print(case)
	dir_case = os.path.join(input_dir, case)
	dir_coord = os.path.join(dir_case, "coord" + os.sep)
	dir_transf = os.path.join(dir_case, "original" + os.sep)

	for panel in panels:
		subdir_coord = os.path.join(dir_coord, panel + os.sep)
		subdir_transf = os.path.join(dir_transf, panel + os.sep)
		if not os.path.exists(subdir_coord):
			os.mkdir(subdir_coord)
		
		# Transform the color mask
		print ("Transforming: " + dir_case + "..." + panel)
		Transform_Virtual_Stack_MT.exec(subdir_coord, subdir_coord, subdir_transf, True)
		IJ.run("Close")

		# Rename the transformed mask
		coord_base = os.path.join(subdir_coord, case + "_Merge_" + panel_base + img_ext)
		coord_panel = os.path.join(subdir_coord, case + "_Merge_" + panel + img_ext)
		os.remove(coord_panel)
		os.rename(coord_base, coord_panel)
		
print("done")