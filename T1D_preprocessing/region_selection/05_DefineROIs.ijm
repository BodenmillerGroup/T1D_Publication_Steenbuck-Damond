/*
 * DEFINE ROIs
 * Nicolas Damond - 10 November 2021
 */
 

// *** GOALS ***
// Interactive script to select region of interests (ROIs) that will be measured by imaging mass cytometry.
// Merged images are opened and ROIs that were automatically selected at the previous step are shown.
// It is possible to interactively delete, add, and modify ROIs using ImageJ rectangle selection and ROI manager.
// Once all desired regions have been selected in the ROI manager, press "OK" to process the next image.
// In the end, the ROIs are saved as "CaseID_iBoxes_Original.zip". This file can be opened using imageJ ROI manager.


// *** INPUT SETTINGS ***

// Parameters
panel_to_crop = "Islet"; // Base panel (on which the ROIs are selected)

// Prompt the user to provide an input folder 
input_dir = getDirectory("Input Directory");  // The input folder should be the "BASE" directory

// *** MAIN ***

// Get the list of "caseID" folders in the input directory
folderlist = getFileList(input_dir);
processFolderNames(input_dir, folderlist);

// Remove tagged (NaN) files from the list (folders that don't correspond to a caseID)
folderlist = Array.delete(folderlist, NaN);

// Loop through files, load the segmented islets as ROI and save the ROI set.
for(i=0; i < folderlist.length; i++) {
	// Defining directories
	cur_case = "" + folderlist[i];
	dir_case = input_dir + File.separator + cur_case;
	dir_merged = dir_case + File.separator + "merged";
	filelist2 = getFileList(dir_merged);
	filelist3 = getFileList(dir_case);

	// Call function to select ROIs
	for(j=0; j < filelist2.length; j++) {
		selectROIS(dir_merged, dir_case, filelist2[j], cur_case, filelist3, panel_to_crop);
	}
	close("*");
	run("Collect Garbage");	
}

print("done");


// *** FUNCTIONS ***

// Function to open the merged image and interactively select ROIs
function selectROIS(dir_merged, dir_case, cur_img, cur_case, filelist3, panel){
	if(cur_img.contains(cur_case + "_Merge_" + panel + ".tif")) {
		// Open the image
		cur_img_path = dir_merged + File.separator + cur_img;
		print("Processing: " + cur_img);
		open(cur_img_path);

		// Open the ROI set
		for(k=0; k < filelist3.length; k++) {
			cur_rois = filelist3[k];

			if(cur_rois.contains(cur_case + "_iBoxes_Original")) {
				cur_rois_path = dir_case + File.separator + cur_rois;
				roiManager("Open", cur_rois_path);
				roicount = roiManager("Count");
				
				// Show all the ROIs and leave time to adapt the ROIs if needed
				roiManager("Show All");
				waitForUser("Modify the ROIs and press" + "'u'" + "to update, then click OK");

				run("Select All");
				roiManager("Save", dir_case + File.separator + cur_case + "_iBoxes_Original" + ".zip");
				roiManager("reset");
			}
		}
	}
}

// List all files in the input folder.
//     -> if the file is not a folder, modify to NaN (tag for deletion).
//     ->  if the folder name doesn't correspond to a caseID, modify to NaN (tag for deletion).
function processFolderNames(input_dir, folderlist) {
	for (i=0; i < folderlist.length; i++){
		if(File.isDirectory(input_dir + folderlist[i])){
			// test if file name == 5 (e.g. 6362/)
			// if it is == 5, convert numbers 1 to 4 to an integer (returns NaN if it doesn't work)
			if(lengthOf(folderlist[i]) == 5)
				folderlist[i] = parseInt(substring(folderlist[i], 0, 4));
			else if(lengthOf(folderlist[i]) == 4)
				folderlist[i] = parseInt(substring(folderlist[i], 0, 3));
			else
				folderlist[i] = NaN;
		}
		else
			folderlist[i] = NaN;

		//if (folderlist[i] != 6160)
			//folderlist[i] = NaN;
	}
}