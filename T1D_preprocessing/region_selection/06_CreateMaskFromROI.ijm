/*
 * GENERATE MASK FROM ROIs
 * Nicolas Damond - 13 Februrary 2020
 */
 

// *** GOALS ***
// 1. Read-in the ROI set defined at the previous step.
// 2. Enlarge each ROI by a fixed value.
// 3. Create a color mask from the ROI set.
// 4. Measure all ROIs, save the enlarged ROI set, and export the results as CSV.
// 5. Create a dummy image that will be used for mask transformation at the next step.


// *** INPUT SETTINGS ***

// Parameters (adapt if needed)
enlargement = 50;  // Pixel expansion factor (40 pixels ~= 25 micrometers)
panel_base = "Islet"; // panel used to define the original ROIs
panels = newArray("Islet", "Immune"); // all panels

// Prompt the user to provide an input folder 
input_dir = getDirectory("Input Directory");  // The input folder should be the "BASE" directory

setBatchMode(true);


// *** MAIN ***

// Get the list of "caseID" folders in the input directory
folderlist = getFileList(input_dir);
processFolderNames(input_dir, folderlist);
// Remove tagged (NaN) files from the list (folders that don't correspond to a caseID)
folderlist = Array.delete(folderlist, NaN);

// Loop through files, define input/output directories and call functions
for(i=0; i < folderlist.length; i++) {
	cur_case = "" + folderlist[i];
	dir_base = input_dir + File.separator + cur_case;
	dir_in = dir_base + File.separator + "merged";
	dir_out = dir_base + File.separator + "mask";
	filelist = getFileList(dir_in);
	
	for(k=0; k < filelist.length; k++) {
		dotIndex = indexOf(filelist[k], "." );
		cur_file = filelist[k].substring(0, dotIndex);
		cur_panel = getPanel(filelist[k], panels);
		extension = filelist[k].substring(dotIndex, filelist[k].length);

		if(extension.contains(".tif")) {
			print("Processing: " + cur_file);
			base_file = cur_case + "_Merge_" + cur_panel + extension;
			base_file_path = dir_in + File.separator + base_file;
			open(base_file_path);
			
			if(cur_panel.matches(panel_base)) {		
				// Open the base panel image and the original ROIs (not expanded)
				fn_roi = dir_base + File.separator + cur_case + "_iBoxes_Original.zip";
				roiManager("Open", fn_roi);
				
				// Create and save the ROI mask
				createColorMask(dir_out, panels, cur_file, cur_panel, extension);
			
				// Enlarge the rois and save the ROI set
				roicount = roiManager("Count");
				enlargeROIs(dir_base, cur_case, cur_panel, roicount);
				
				// Measure all ROIs and save the ROI coordinates in new arrays + Export to CSV
				measureROIs(dir_base, cur_panel, roicount);
			}
			
			if(!(cur_panel.matches(panel_base))) {
				// Create a temporary image that will be use for mask transformation
				createDummyImage(dir_out, cur_panel, cur_file, panel_base, extension);
			}
			close("*");
			run("Collect Garbage");
		}
	}
}

print("done");


// *** FUNCTIONS ***

// Create dummy images for transformation
function createDummyImage(dir_out, cur_panel, cur_file, panel_base, extension) {
	newImage("Dummy", "8-bit black", getWidth(), getHeight(), 1);
		
	folder_dummy = dir_out + File.separator + cur_panel;
	if (!(File.exists(folder_dummy)))
		File.makeDirectory(folder_dummy);
	fn_dummy = folder_dummy +  File.separator + cur_file + extension;
	selectWindow("Dummy");
	saveAs("Tiff", fn_dummy);
}

// Measure and save ROIs, export measurements to CSV
function measureROIs(dir_base, cur_panel, roicount) {
	rois = newArray;
	for (k=0 ; k < roicount ; k++){
		rois = Array.concat(rois, k);
	}
	roiManager("select", rois);
	roiManager("Measure");

	// Store the ROI coordinates in new arrays
	xs = newArray(roicount);
	ys = newArray(roicount);
	ws = newArray(roicount);
	hs = newArray(roicount);
	for (k=0 ; k < roicount ; k++) {
		roiManager("select", k);
		Roi.getBounds(x,y,w,h);
		xs[k] = x;
		ys[k] = y;
		ws[k] = w;
		hs[k] = h;
	}

	// Store the ROI color intensity in a new array
	intensities = newArray(roicount);
	
	for (k=0 ; k < roicount ; k++) {
		intensities[k] = getResult("Max", k);
	}
	roiManager("deselect");
	roiManager("delete");
	close("*");
	
	// Export to CSV
	fn_csv = dir_base + File.separator + "TEMPBoxes_" + cur_panel + ".csv";
	run("Clear Results");
	for (k=0 ; k < roicount ; k++) {
    	setResult("X", k, xs[k]);
    	setResult("Y", k, ys[k]);
    	setResult("W", k, ws[k]);
    	setResult("H", k, hs[k]);
    	setResult("Color", k, intensities[k]);
	}
	updateResults();
	saveAs("Results", fn_csv);
	run("Clear Results");	
}

// Function to enlarge the ROIs and save the ROI set
function enlargeROIs(dir_base, cur_case, cur_panel, roicount) {
	for (k=0; k < roicount; k++){ 
		roiManager("select", k);
		run("Enlarge...", "enlarge=enlargement");
		roiManager("Update");
	}
	fn_rois = dir_base + File.separator + cur_case + "_Boxes_" + cur_panel + ".zip";
	roiManager("Save", fn_rois);
}

// Function to create and save a color ROI mask
function createColorMask(dir_out, panels, cur_file, cur_panel, extension){
	newImage("Mask", "8-bit black", getWidth(), getHeight(), 1);
	
	for (i=0; i < roiManager("count"); i++) {
		roiManager("select", i);
		setColor(i+1);
		fill();
	}
	resetMinAndMax();
	run("glasbey");
	selectWindow("Mask");
	
	for (k=0; k < panels.length; k++) {
		folder_mask = dir_out + File.separator + panels[k];
		if (!(File.exists(folder_mask)))
			File.makeDirectory(folder_mask);
		fn_mask = folder_mask +  File.separator + cur_file + extension;
		saveAs("Tiff", fn_mask);
	}

}

// Identify panel from filename
function getPanel(cur_file, panels){
	for (i=0; i < panels.length; i++){
		if(cur_file.contains(panels[i] + ".")){
			cur_panel = panels[i];
			return cur_panel;
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