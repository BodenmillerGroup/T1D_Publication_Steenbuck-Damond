/*
 * LOAD SEGMENTED ISLETS
 * Nicolas Damond - 10 November 2021
 */
 
 
// *** GOALS ***
// 1. Identify the islets from the islet masks (exported from CellProfiler) and load them as a set of ROIs.
// 2. Get islets bounding boxes to define regions of interest.
// 3. Randomly subset the ROIs


// *** INPUT SETTINGS ***

// Parameters
panel_base = "Islet"; // Panel on which the islets are selected
nb_rois = 125;  // Number of ROIs to return (randomly subsetted from all identified ROIs)
file_ext = ".tiff";  // extension of the single channel slide scanner images
delimiter = "_";  // file name delimiter

// Prompt the user to provide an input folder 
input_dir = getDirectory("Input Directory");  // The input folder should be the "BASE" directory

setBatchMode(true);


// *** MAIN ***

// Get the list of folders in the input directory. Remove the folders that don't correspond to a case ID (tagged as NaN)
folderlist = getFileList(input_dir);
folderlist = processFolderNames(input_dir, folderlist);
folderlist = Array.delete(folderlist, NaN);

// Loop through files
for(i=0; i < folderlist.length; i++) {
	// Define directories
	cur_donor = "" + folderlist[i];
	dir_donor = input_dir + File.separator + cur_donor;
	dir_mask = dir_donor + File.separator + "mask";
	filelist2 = getFileList(dir_mask);
	
	// Process the binary masks exported from CellProfiler
	for(j=0; j < filelist2.length; j++) {
		processMasks(input_dir, filelist2[j], dir_mask, dir_donor, panel_base, delimiter, file_ext, nb_rois);
	}
}

print("done");


// *** FUNCTIONS ***

// Function to get ROIs from masks
function processMasks (input_dir, cur_file, dir_mask, dir_donor, panel_base, delim, file_ext, nb_rois) {
	cur_img = cur_file;
	parts = split(cur_img, delim);
	
	if(cur_img.endsWith("_binmask" + file_ext) && parts[1] == panel_base) {
		// Open the image
		cur_img_path = dir_mask + File.separator + cur_img;
		print("Processing: " + cur_img);
		open(cur_img_path);
		
		// Rotate and scale-up to revert the initial binning
		run("Rotate... ", "angle=180 grid=1 interpolation=Bilinear");
		run("Scale...", "x=4 y=4 interpolation=Bilinear create");
		
		// Identify the segmented objects
		run("Convert to Mask");
		run("Analyze Particles...", "pixel show=Outlines clear add");

		// Convert to bounding box
		roicount = roiManager("Count");
		for (n=0; n<roicount; n++){ 
			roiManager("select", n)
			run("To Bounding Box");
			roiManager("Update");
		}

		// Randomly subset ROIs
		roi_array = Array.getSequence(roicount);
		shuffle(roi_array);
		roi_array = Array.slice(roi_array, 0, roicount - nb_rois);
		Array.sort(roi_array);
		Array.reverse(roi_array);

		for (h=0; h<roi_array.length; h++) {
			roiManager("select", roi_array[h])
				roiManager("delete")
		}

		// Save the array
		run("Select All");
		roiManager("Save", dir_donor + File.separator + cur_donor + "_iBoxes_Original.zip");
		roiManager("reset");
		
		close("*");
		run("Collect Garbage");
	}
}

// Function to shuffle an array (from imagej macro examples)
function shuffle(array) {
   n = array.length;  // The number of items left to shuffle (loop invariant).
   while (n > 1) {
      k = randomInt(n);     // 0 <= k < n.
      n--;                  // n is now the last pertinent index;
      temp = array[n];  // swap array[n] with array[k] (does nothing if k==n).
      array[n] = array[k];
      array[k] = temp;
   }
}

// returns a random number, 0 <= k < n
function randomInt(n) {
   return n * random();
}

// List all files in the input folder.
//     -> if the file is not a folder, modify to NaN (tag for deletion).
//     ->  if the folder name doesn't correspond to a donorID, modify to NaN (tag for deletion).
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
	}
	return folderlist;
}