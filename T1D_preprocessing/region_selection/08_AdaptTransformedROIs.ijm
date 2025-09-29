/*
 * TRANSFORM ROI COORDINATES
 * Nicolas Damond - 14 Februrary 2020
 */
 
// *** GOALS ***
// Interactive script to adapt the ROI from non-base panels after automated registration.
// 1. Load the registered masks for non-base panels. 
// 2. Threshold the mask and identify the ROIs using "Analyze particles". Measure the intensity of each ROI (corresponds to colors, used to match ROIs from different panels).
// 3. Enlarge each ROI by a defined value.
// 4. Interactively adpat the ROIs.
// 5. Export all ROI coordinates to CSV files.

// *** INPUT SETTINGS ***

// Parameters (adapt if needed)
panel_base = "Islet"; // panel used to define the ROIs
panels = newArray("Islet", "Immune"); // enter all panels, including the base panel
panels_islet = newArray(); // enter other "Islet" type panels here in case there are any (otherwise leave the array empty). DO NOT include the "base" panel here!

enlarge_Islet = 50;  // Enlargement factor for the islet panel(s) (this enlargement factor will be applied to any panel in the "panels_islet" array above.
enlarge_Immune = 80;  // Enlargement factor for the other panel(s) (this enlargement factor will be applied to all the other panels).
file_ext = ".tif";


// *** MAIN ***

// Prompt the user to provide an input folder
input = getDirectory("Input Directory");

// Get the list of "caseID" folders in the input directory
folderlist = getFileList(input);
processFolderNames(input, folderlist);
// Remove tagged (NaN) files from the list (folders that don't correspond to a caseID)
folderlist = Array.delete(folderlist, NaN);

// Loop through files, define input/output directories, and call functions
for(i=0; i < folderlist.length; i++) {
	// Defining directories
	cur_case = "" + folderlist[i];
	dir_base = input + File.separator + cur_case;
	dir_mask = dir_base + File.separator + "mask";
	dir_registered = dir_base + File.separator + "registered";

	for(j=0; j < panels.length; j++) {
		cur_panel = panels[j];
		cur_mask_fn = cur_case + "_Merge_" + cur_panel + file_ext;
		cur_mask_path = dir_mask + File.separator + cur_panel + File.separator + cur_mask_fn;
		
		print("Processing: " + cur_mask_fn);
		convertCoordinates(cur_mask_fn, cur_mask_path, dir_base, cur_case, cur_panel, panel_base, file_ext);
	}
}

print("done");


// *** FUNCTIONS ***

function convertCoordinates(cur_mask_fn, cur_mask_path, dir_base, cur_case, cur_panel, panel_base, file_ext) {

	// Open the registered mask and retrieve the ROIs
	if(!(cur_panel.matches(panel_base))) {

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

		// Select the bounding box of each ROI, then enlarge it by a user-defined value
		for (k=0 ; k < roicount ; k++) {
			roiManager("select", k);
			run("Select Bounding Box (guess background color)");
			roiManager("Update");
			roiManager("select", k);

			// Enlarge the ROIs by the enlarge_Immune factor, unless the current panel is in the panels_islet list.
			for(m=0; m < panels.length; m++) {
				enlargement = enlarge_Immune;
				for(n=0; n < panels_islet.length; n++) {
					if(cur_panel.matches(panels_islet[n])) {
						enlargement = enlarge_Islet;
					}
				}
			}
			run("Enlarge...", "enlarge=enlargement");
			roiManager("Update");
		}

		// Open the merged image
		merge_img_fn = cur_case + "_Merge_" + cur_panel + file_ext;
		if(!(cur_panel.matches(panel_base))) {
			merge_img_path = dir_base + File.separator + "registered" + File.separator + merge_img_fn;
		} else if (cur_panel.matches(panel_base)) {
			merge_img_path = dir_base + File.separator + "merged" + File.separator + merge_img_fn;
		}
		open(merge_img_path);

		// Adapt the ROIs (for non-base panels)
		if(!(cur_panel.matches(panel_base))) {
			// Show all the ROIs and leave time to adapt the ROIs if needed
			roiManager("Show All");
			waitForUser("Modify the ROIs and press" + "'u'" + "to update, then click OK");
		} else if(cur_panel.matches(panel_base)) {
			// Open the ROIs
			fn_rois_base = dir_base + File.separator + cur_case + "_Boxes_" + cur_panel + ".zip";
			roiManager("Open", fn_rois_base);	
		}
	
		// Close the merged image and re-open the transformed mask

		// Save the ROIs
		fn_rois = dir_base + File.separator + cur_case + "_Boxes_" + cur_panel + ".zip";
		rois = newArray;
		for (k=0 ; k < roicount ; k++){
			rois = Array.concat(rois, k);
		}
		roiManager("select", rois);
		roiManager("Save", fn_rois);
		
		close("*");
		run("Collect Garbage");
		open(cur_mask_path);

		roicount2 = roiManager("Count");
		print("Processing " + merge_img_fn + ": " + roicount2 + " ROIs");
	
		if (roicount != roicount2)
			print("The number of ROIs has been modified");
	
		// Store the ROI coordinates in new arrays
		xs = newArray(roicount2);
		ys = newArray(roicount2);
		ws = newArray(roicount2);
		hs = newArray(roicount2);
	
		for (k=0 ; k < roicount2 ; k++) {
			roiManager("select", k);
			Roi.getBounds(x,y,w,h);
			xs[k] = x;
			ys[k] = y;
			ws[k] = w;
			hs[k] = h;
		}
	
		// Measure all ROIs
		rois = newArray;
		for (k=0 ; k < roicount2 ; k++){
			rois = Array.concat(rois, k);
		}
		roiManager("select", rois);
		roiManager("Measure");
	
		// Store the ROI color intensity in a new array
		intensities = newArray(roicount2);
		for (k=0 ; k < roicount2 ; k++) {
			intensities[k] = getResult("Max", k);
		}

		roiManager("deselect");
		roiManager("delete");
	
		// Export to CSV
		fn_csv = dir_base + File.separator + "TEMPBoxes_" + cur_panel + ".csv";
		run("Clear Results");
		for (k=0 ; k < roicount2 ; k++) {
	    	setResult("X", k, xs[k]);
	    	setResult("Y", k, ys[k]);
	    	setResult("W", k, ws[k]);
	    	setResult("H", k, hs[k]);
    		setResult("Color", k, intensities[k]);
		}
		updateResults();
		saveAs("Results", fn_csv);
		close("*");
	}

	// Generate CSV files to record the coordinates
	run("Clear Results");
	fn_csv = dir_base + File.separator + cur_case + "_Coordinates_" + cur_panel + ".csv";
	setResult("CoordNb", 0, 0);
	setResult("SlideX", 0, 0);
    setResult("SlideY", 0, 0);
    setResult("IMCX", 0, 0);
    setResult("IMCY", 0, 0);
	saveAs("Results", fn_csv);
	run("Clear Results");
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