/*
 * GENERATE MASK FROM ROIs
 * Nicolas Damond - 13 Februrary 2020
 */

// *** GOALS ***
// 1. Read in the coordinates that were defined at the previous step.
// 2. Create a color mask from the coordinates selected at the previous step.
// 3. Create a dummy image that will be used for mask transformation at the next step.


// *** INPUT SETTINGS ***

// Panels (adapt if needed)
panel_base = "Islet" // panel used to define the original ROIs
panels = newArray("Islet", "Immune"); // all panels

// Prompt the user to provide an input folder
input = getDirectory("Input Directory");

setBatchMode(true);


// *** MAIN ***

// Get the list of "CaseID" folders in the input directory
folderlist = getFileList(input);
processFolderNames(input, folderlist);
// Remove tagged (NaN) files from the list (folders that don't correspond to a Case ID)
folderlist = Array.delete(folderlist, NaN);


// Loop through files, define input/output directories and call functions
for(i=0; i < folderlist.length; i++) {
	// Define working directories
	cur_case = "" + folderlist[i];
	dir_base = input + File.separator + cur_case;
	dir_in = dir_base + File.separator + "merged";
	dir_out = dir_base + File.separator + "coord";
	if(!(File.exists(dir_out)))
		File.makeDirectory(dir_out);
		
	// Open the merged image corresponding to the base panel
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
	
			// Create the mask
			if(cur_panel.matches(panel_base)) {
				// Open the previously defined set of coordinates
				fn_roi = dir_base + File.separator + cur_case + "_Coordinates_Islet.zip";
				roiManager("Open", fn_roi);
				
				// Create and save the coordinate mask
				createColorMask(dir_out, panels, cur_file, cur_panel, extension);
				roiManager("deselect");
				roiManager("delete");
			}
	
			if(!(cur_panel.matches(panel_base))) {
				// Create a temporary image that will be use for mask transformation
				createDummyImage(dir_out, cur_panel, cur_file, panel_base, extension);
			}
		}

		// Close the opened images
		close("*");
		run("Collect Garbage");
	}
}

print("done");


// *** FUNCTIONS ***

// Function to create and save a color mask for coordinates
function createColorMask(dir_out, panels, cur_file, cur_panel, extension){
	newImage("Mask", "8-bit black", getWidth(), getHeight(), 1);
	
	for (i=0; i < roiManager("count"); i++) {
		roiManager("select", i);
		setColor(i*10+1);
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