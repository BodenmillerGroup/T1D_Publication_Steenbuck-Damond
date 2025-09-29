import os
import re
from ij import IJ
from register_virtual_stack import Register_Virtual_Stack_MT


# Parameters

## Change the input directory ("BASE" folder), and the case and panel lists if needed
input_dir = "/Users/nathan/Downloads/region_selection"
p = re.compile("^\d{4}$")
## List directory
cases = os.listdir(input_dir)
caseList = [x for x in cases if p.match(x)]

# caseList = ("6034", "6044", "6048")
			
#"6055", "6063", "6090", "6123", "6126", "6135",
#"6145", "6147", "6151", "6171", "6178", "6181", "6197", "6208", "6209",
#"6224", "6228", "6229", "6234", "6247", "6264", "6285", "6289", "6301",
#"6303", "6310", "6314", "6321", "6324", "6328", "6347", "6362", "6367",
#"6380", "6388", "6396", "6397", "6400", "6405", "6414", "6418", "6420",
#"6421", "6422", "6424", "6428", "6429", "6433", "6437", "6447", "6449",
#"6450", "6456", "6458", "6469", "6483", "6488", "6494", "6496", "6505",
#"6526", "6531", "6532", "6533", "6534", "6538", "6547", "6549", "6550",
#"6509", "6510", "6512", "6514", "6517", "6518", "6519", "6520", "6521",
#"6551", "6553", "6558", "6562", "6563", "8002", "8011")
# Pancreatitis: "6036","6043","6061","6150","6225","6505","6522"

#caseList = ("6036", "6055", "6090", "6147", "6225", "6228", "6303",
#			"6321", "6388", "6396", "6414", "6421", "6428", "6437",
#			"6458", "6510", "6519", "6521", "6522", "6526", "6532",
#			"6547", "6550", "6553", "6558", "6562", "6563", "TMA") 

# # Use this format to process multiple cases: ("6289", "6310", "6328") # Use this format to process a single case: (["6450"])
panels = (["Immune"])  # Use this format to process multiple panels (do not include the base panel): ("Islet2", "Lympho", "Myelo")  # Use this format to process a single panel: (["Immune"])
panel_base = "Islet" # Base panel (used to define the ROIs)
img_ext = ".tif" # Image extension

###############

# Registration settings

## Shrinkage option (False = 0)
use_shrinking_constraint = 0

## Save transforms
save_transforms = 1
		 
p = Register_Virtual_Stack_MT.Param()

## SIFT parameters:
p.sift.maxOctaveSize = 1024
p.sift.fdSize = 8
p.sift.initialSigma = 1.6
p.sift.steps = 3	
p.maxEpsilon = 25

## 0=TRANSLATION, 1=RIGID, 2=SIMILARITY, 3=AFFINE
p.featuresModelIndex = 3
## 0=TRANSLATION, 1=RIGID, 2=SIMILARITY, 3=AFFINE, 4=ELASTIC, 5=MOVING_LEAST_SQUARES
p.registrationModelIndex = 3

## The "inlier ratio":
p.minInlierRatio = 0.05

###############

# MAIN

for case in caseList:
	dir_case = os.path.join(input_dir, case + os.sep)
	dir_input = os.path.join(dir_case, "original")
	dir_output = os.path.join(dir_case, "registered")
	
	for panel in panels:
		if (str(panel) != str(panel_base)):	
			subdir_input = os.path.join(dir_input, panel + os.sep)
			subdir_output = os.path.join(dir_output, panel + os.sep)
			dir_transf = subdir_input
			if not os.path.exists(subdir_output):
				os.mkdir(subdir_output)
				
			print("Processing case: " + str(case) + " - Panel " + str(panel))
			input_files = [fl for fl in os.listdir(subdir_input) if os.path.isfile(os.path.join(subdir_input, fl))]

			# Select the reference image (not the image corresponding to the base panel)
			for cur_file in input_files:
				if (((panel_base + "-DAPI.") not in cur_file) and (cur_file.endswith(img_ext))):
					reference_image = cur_file
					print(reference_image)

			# Execute alignment.
			print ("Aligning..." + reference_image)
			Register_Virtual_Stack_MT.exec(subdir_input, subdir_output, dir_transf, reference_image, p, use_shrinking_constraint)

			# Close the alignment window (comment out to check the alignment)
			IJ.run("Close")

print("done")
