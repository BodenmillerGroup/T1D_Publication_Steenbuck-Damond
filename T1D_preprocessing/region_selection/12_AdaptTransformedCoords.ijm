/*
 * TRANSFORM ROI COORDINATES
 * Nicolas Damond - 14 Februrary 2020
 */

// *** GOALS ***
// Interactive script to adapt the ROI from non-base panels after automated registration.
// 1. Load the registered coordinate masks for non-base panels. 
// 2. Threshold the mask and identify the coordinates.
// 4. Interactively adapt the coordinates.
// 5. Save the adapted coordinates.


// *** INPUT SETTINGS ***

// Parameters (adapt if needed)
panel_base = "Islet" // panel used to define the ROIs
panels = newArray("Islet", "Immune"); // enter all panels, including the base panel
file_ext = ".tif";

// Prompt the user to provide an input folder
input = getDirectory("Input Directory");


// *** MAIN ***

// Get the list of "caseID" folders in the input directory
folderlist = getFileList(input);
processFolderNames(input, folderlist);
// Remove tagged (NaN) files from the list (folders that don't correspond to a caseID)
folderlist = Array.delete(folderlist, NaN);


// Loop through files, define input/output directories, and call functions
for(i=0; i < folderlist.length; i++) {
	// Define directories
	cur_case = "" + folderlist[i];
	dir_case = input + File.separator + cur_case;
	dir_coord = dir_case + File.separator + "coord";
	dir_registered = dir_case + File.separator + "registered";

	for(j=0; j < panels.length; j++) {
		cur_panel = panels[j];
		cur_mask_fn = cur_case + "_Merge_" + cur_panel + file_ext;
		cur_mask_path = dir_coord + File.separator + cur_panel + File.separator + cur_mask_fn;

		if(!(cur_panel.matches(panel_base))) {
			print("Processing: " + cur_mask_fn);
			convertCoordinates(cur_mask_fn, cur_mask_path, dir_case, cur_case, cur_panel, panel_base, file_ext);
		}
	}
}

print("done");


// *** FUNCTIONS ***

function convertCoordinates(cur_mask_fn, cur_mask_path, dir_case, cur_case, cur_panel, panel_base, file_ext) {

	// Open the mask and threshold it
	open(cur_mask_path);
	setAutoThreshold("Default dark no-reset");
	run("Threshold...");
	setThreshold(1, 255);
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Close");
	
	// Identify the ROIs
	run("Analyze Particles...", "exclude include add");

	// Close the mask and reopen the transformed mask (not thresholded)
	close("*");
	run("Collect Garbage");
	open(cur_mask_path);
	roicount = roiManager("Count");
	print("Number of coordinates identified: " + roicount);

	// Measure all ROIs
	rois = newArray;
	for (k=0 ; k < roicount ; k++){
		rois = Array.concat(rois, k);
	}
	roiManager("select", rois);
	roiManager("Measure");

	// Store the mean intensity in a new array
	intensities = newArray(roicount);
	for (k=0 ; k < roicount ; k++) {
		intensities[k] = getResult("Mean", k);
	}
	roiManager("deselect");


	// Open the merged image
	merge_img_fn = cur_case + "_Merge_" + cur_panel + file_ext;
	merge_img_path = dir_case + File.separator + "registered" + File.separator + merge_img_fn;
	open(merge_img_path);
	
	roicount2 = roiManager("Count");
	print("Processing " + merge_img_fn + ": " + roicount2 + " ROIs");

	// Select ROIs one by one and get their coordinates
	xs = newArray(roicount2);
	ys = newArray(roicount2);

	// Get the ROIs coordinates
	for (k=0 ; k < roicount2 ; k++) {
		roiManager("select", k);
		Roi.getBounds(x,y,w,h);
		xs[k] = x;
		ys[k] = y;
	}

	// Convert the ROIs to points
	for (k=0 ; k < roicount2 ; k++) {
		roiManager("select", k);
		Roi.getBounds(x,y,w,h);
		xs[k] = x + (w/2);
		ys[k] = y + (h/2);
		makePoint(xs[k], ys[k]);
		roiManager("Add");
	}

	// Delete the rectangular ROIs
	for (k=0 ; k < roicount2 ; k++) {
		roiManager("select", newArray(roicount2));
		roiManager("delete");
	}

	// Interactive correction
	roiManager("Show All");
	waitForUser("Modify the coordinates and press " + "'u'" + " to update, then click OK");

	// Select all ROIs
	roiManager("deselect");
	roicount3 = roiManager("Count");
	rois3 = newArray;
	for (k=0 ; k < roicount3 ; k++){
		rois3 = Array.concat(rois3, k);
	}
	roiManager("select", rois3);
	roiManager("Measure");
	
	xt = newArray(roicount3);
	yt = newArray(roicount3);
	
	// Select ROIs one by one and get their coordinates
	for (k=0 ; k < roicount3 ; k++) {
		roiManager("select", k);
		Roi.getBounds(x, y, w, h);
		xt[k] = x;
		yt[k] = y;
	}

	// Save the ROIs
	fn_rois = dir_case + File.separator + cur_case + "_Coordinates_" + cur_panel + ".zip";
	roiManager("select", rois3);
	roiManager("Save", fn_rois);

	roiManager("deselect");
	roiManager("delete");
	close("*");
	run("Collect Garbage");
}

// List all files in the input folder.
//     -> if the file is not a folder, modify to NaN (tag for deletion).
//     ->  if the folder name doesn't correspond to a caseID, modify to NaN (tag for deletion).
function processFolderNames(input, folderlist) {
	for (i=0; i < folderlist.length; i++){
		if(File.isDirectory(input + folderlist[i])){
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