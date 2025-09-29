/*
 * IMPORT AND PROCESS SLIDE SCANNER FILES
 * Nicolas Damond - 06 Februrary 2020
 */

// *** GOALS ***
// 1. Generate DAPI images for registration
// 2. Generate Merged RGB images
// 3. Generate binned images for "mcd" files


// *** INPUT SETTINGS ***

// Images are expected to be named as "CaseID_Panel-Channel.tif" (e.g. "1234_Islet-DAPI.tif")
// Change the file extension and delimiters if needed
filext = ".tif";
delimiter1 = "_"; // Character that separates the CaseID from the Panel in the .tif image names.
delimiter2 = "-"; // Character that separates the Panel from the Channel in the .tif image names.
panel_base = "Islet"; // Panel used to define the ROIs

// Change the channel order and names if needed
redch = "Cy5";
greench = "Cy3";
bluech = "DAPI";
isletchannel = "Cy5"; // Channel containing the islet marker.

// Prompt the user to provide an input folder
input = getDirectory("Input Directory"); // "BASE" folder
input_dir = input + File.separator + "czi_if" + File.separator;

// Modify the string if needed (or use a prompt)
addstring = "_Merge_"; 
//Dialog.create("Suffix of merged image");
//Dialog.addString("Add suffix: ", "");
//Dialog.show();
//addstring = Dialog.getString();

setBatchMode(true);


// *** MAIN ***

// Get the list of files in the input folder. Remove the files that don't have the right extension (tagged as NaN)
filelist = getFileList(input_dir);
processFileNames(input_dir, filelist, filext);
filelist = Array.delete(filelist, NaN);

// Split the filenames into case IDs and channel names (stored in three new arrays: "cases", "panels" and "channels")
//     -> this considers that the file names are in the format: "CaseID + delimiter1 + Panel + delimiter2 + ChannelName + .fileExtension"
//     -> should be modified if the files are names differently
cases = newArray(filelist.length);
panels = newArray(filelist.length);
case_panel = newArray(filelist.length);
channels = newArray(filelist.length);
splitFileNames(filelist, cases, panels, case_panel, channels, delimiter1, delimiter2);

// Get unique case and channel IDs
caseID = ArrayUnique(cases);
panelID = ArrayUnique(panels);
case_panelID = ArrayUnique(case_panel);
channelID = ArrayUnique(channels);

// Process all files
for(i=0; i < case_panelID.length; i++) {
	// Parse file names
	cur_case_panel = split(case_panelID[i], delimiter1);
	cur_case = cur_case_panel[0];
	cur_panel = cur_case_panel[1];
	print("Processing case", cur_case, " - Panel", cur_panel);
	
	// Retrieve the names of individual channel images (imagenames[0] = red, imagenames[1] = green, imagenames[2] = blue)
	imagenames = parseFileNames(input_dir, filelist, filext, cur_case, cur_panel, delimiter1, delimiter2);
	if (imagenames[2] == "0")
		break;
	
	// Process the DAPI image (will be used later for registration)
	ProcessDAPI(input_dir, imagenames[2], cur_case, cur_panel, panelID, panel_base);

	// Process the RGB images
	ProcessRGB(input_dir, imagenames, cur_case, cur_panel, panel_base, addstring);
	
	// Generate binned images of the Islet channel (will be used for islet segmentation)
	createBinImages(input_dir, cur_case, cur_panel, panel_base, isletchannel, delimiter1, delimiter2);
}

print("done");



// *** FUNCTIONS ***

// Function to create a binned image of the islet channel
function createBinImages(input_dir, cur_case, cur_panel, panel_base, isletchannel, delimiter1, delimiter2) {
	
	// Open the image
	cur_imagename = cur_case + delimiter1 + cur_panel + delimiter2 + isletchannel + ".tif";
	cur_imagepath = input_dir + File.separator + cur_imagename;
	run("Bio-Formats Importer", "open=[cur_imagepath] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT stitch_tiles");
	
	// Create a binned image
	run("Rotate... ", "angle=180 grid=1 interpolation=Bilinear");
	run("Enhance Contrast", "saturated=0.35");
	//setMinAndMax(10, 50);
	run("Bin...", "x=4 y=4 bin=Average");
	
	// Save the binned image to the "mcd" subfolder
	casedir = File.getParent(input_dir) + File.separator + cur_case;
	mcd_dir = casedir + File.separator + "mcd";
	if(!File.exists(mcd_dir))
		File.makeDirectory(mcd_dir);
	out_imagename = cur_case + delimiter1 + cur_panel + delimiter2 + isletchannel + delimiter1 + "BIN" + ".tif";
	out_imagepath = mcd_dir + File.separator + out_imagename;
	saveAs("Tiff", out_imagepath);
	
	close("*");
	run("Collect Garbage");
}
			
// Function to open and merge individual channel images
function ProcessRGB(input_dir, imagenames, cur_case, cur_panel, panel_base, addstring) {
	// Open and rotate individual channels images
	for (k=0; k<imagenames.length; k++) {
		cur_imagename = imagenames[k];
		cur_imagepath = input_dir + File.separator + cur_imagename;
		run("Bio-Formats Importer", "open=[cur_imagepath] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT stitch_tiles");
		run("Rotate... ", "angle=180 grid=1 interpolation=Bilinear");
	}
	
	// Merge the channels
	merged_imagename = cur_case + addstring + cur_panel + ".tif";
	mergechannels = "c1=" + imagenames[0] + " " + "c2=" + imagenames[1] + " " + "c3=" + imagenames[2];
	run("Merge Channels...", mergechannels);

	// Adjust color balance for red (4), green (2) and blue (1) channels
	selectWindow("RGB");
	setMinAndMax(1, 40, 4);
	setMinAndMax(10, 50, 2);
	setMinAndMax(5, 25, 1);

	// Save the merged image in the "merged" folder
	casedir = File.getParent(input_dir) + File.separator + cur_case;
	outdir_merged = casedir + File.separator + "merged";
	if(!File.exists(outdir_merged))
		File.makeDirectory(outdir_merged);
	out_fn = outdir_merged + File.separator + merged_imagename;
	saveAs("Tiff", out_fn);
	
	close("*");
	run("Collect Garbage");
}

// Function to process DAPI images (used later for registration)
function ProcessDAPI(input_dir, dapi_imagename, cur_case, cur_panel, panelID, panel_base) {
	
	// Open the DAPI image, rotate it 180°, and adjust the contrast
	dapi_image = input_dir + dapi_imagename;
	run("Bio-Formats Importer", "open=[dapi_image] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT stitch_tiles");
	run("Rotate... ", "angle=180 grid=1 interpolation=Bilinear");
	run("Enhance Contrast", "saturated=0.35");
	
	// Create one directory the "original" folder that will be used to store output images
	casedir = File.getParent(input_dir) + File.separator + cur_case;
	if(!File.exists(casedir))
		File.makeDirectory(casedir);
			
	// Create "original" folder that will be used to store output images
	outdir_original = casedir + File.separator + "original";
	if(!File.exists(outdir_original))
		File.makeDirectory(outdir_original);
	
	// In the "original" folder, create one subfolder per panel:
	// 	- If the processed DAPI image corresponds to the "base panel" (used as a base for registration), save a copy of this image in each panel subfolder
	// 	- If the processed DAPI image corresponds to another panel, save the image in the panel subfolder corresponding to the panel of the current image.
	// Directory structure:
	// ¬ original
	// ¬¬ 1234_PanelBase_DAPI.tif
	// ¬¬ PanelX
	// ¬¬¬ 1234_PanelX_DAPI.tif
	// ¬¬¬ 1234_PanelBase_DAPI.tif
	// ¬¬ PanelY
	// ¬¬¬ 1234_PanelY_DAPI.tif
	// ¬¬¬ 1234_PanelBase_DAPI.tif
	
	if (!(cur_panel.matches(panel_base))) {
		outdir_panel = outdir_original + File.separator + cur_panel;
		if(!File.exists(outdir_panel))
			File.makeDirectory(outdir_panel);
		out_fn = outdir_panel + File.separator + dapi_imagename;
		saveAs("Tiff", out_fn);
	} else if (cur_panel.matches(panel_base)) {
		for (k=0 ; k < panelID.length ; k++) {
			if (!(panelID[k].matches(panel_base))) {
				outdir_panel = outdir_original + File.separator + panelID[k];
				if(!File.exists(outdir_panel))
					File.makeDirectory(outdir_panel);
				out_fn = outdir_panel + File.separator + dapi_imagename;
				saveAs("Tiff", out_fn);
			}
		}
		// Save in addition a copy of the base panel DAPI image in the "original" folder
		//out_fn = outdir_original + File.separator + dapi_imagename;
		//saveAs("Tiff", out_fn);
	}
	close("*");
	run("Collect Garbage");
}

// Go through the file list and parse file names
function parseFileNames(input_dir, filelist, filext, cur_case, cur_panel, delimiter1, delimiter2) {
	for (j=0; j < filelist.length; j++){
		// Create a folder for each case
		outdir = File.getParent(input_dir) + File.separator + cur_case;			
		if(!File.exists(outdir))
			File.makeDirectory(outdir);
		
		// Parse current file name 
		parts = split(filelist[j], delimiter2);

		// Retrieve image names corresponding to individual channels
		imagenames = newArray(3);
		if(parts[0] == cur_case + delimiter1 + cur_panel){
			if(parts[1] == redch)
				red = filelist[j]+filext;
			if(parts[1] == greench)
				green = filelist[j]+filext;
			if(parts[1] == bluech)
				blue = filelist[j]+filext;
		}
	}
	// Store image names in an array and return it
	redgreenblue = newArray(3);
	redgreenblue[0] = red;
	redgreenblue[1] = green;
	redgreenblue[2] = blue;
	return redgreenblue;
}

// List all files in the input folder.
//     -> if the extension is not right, modify to NaN (tag for deletion).
//     -> if the extension is right, remove it from the filename.
function processFileNames(input_dir, filelist, filext) {
	for (i=0; i < filelist.length; i++){
		if(!endsWith(filelist[i], filext))
			filelist[i] = NaN;
		else
			filelist[i] = replace(filelist[i], filext, "");
	}
}

// Split the file name based on delimiters, store each part in a separate array
function splitFileNames(filelist, cases, panels, case_panel, channels, delimiter1, delimiter2){
	for (i=0; i < filelist.length; i++){
		p1 = split(filelist[i], delimiter2);
		case_panel[i] = p1[0];
		channels[i] = p1[1];
		
		p2 = split(p1[0], delimiter1);
 		cases[i] = p2[0];
 		panels[i] = p2[1];
	}
}

function ArrayUnique(array) {
	array 	= Array.sort(array);
	array 	= Array.concat(array, 999999);
	uniqueA = newArray();
	i = 0;	
   	while (i<(array.length)-1) {
		if (array[i] == array[(i)+1]) {		
		} else {
			uniqueA = Array.concat(uniqueA, array[i]);
		}
   		i++;
   	}
	return uniqueA;
}