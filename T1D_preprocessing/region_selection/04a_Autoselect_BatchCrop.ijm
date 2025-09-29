/*
 * BATCH CROP SLIDE SCANNER IMAGES
 * Nicolas Damond - 09 November 2021
 */
 
// *** GOALS ***
// 1. Crop slide scanner images. The crops will be used for automated islet selection.


// *** INPUT SETTINGS ***
  
// Parameters (modify if needed)
islet_channel = "Cy5"; // Slide scanner channel containing the islet marker
panel_to_crop = "Islet"; // Crop only images from this panel
file_ext = ".tif";  // extension of the single channel slide scanner images
crop_nb = 2;  // number of crops to make per case
cropX = 1000;  // size of the crops (x-axis)
cropY = 1000;  // size of the crops (y-axis)
delimiter1 = "_";  // first delimiter in slide scanner file names
delimiter2 = "-";  // second delimiter in slide scanner file names

// Prompt the user to provide an input folder 
input_dir = getDirectory("Input Directory");  // The input folder should be the "BASE" directory

// Make folders
ssif_folder = input_dir + File.separator + "czi_if";
outdir1 = input_dir + File.separator + "islet_selection";
if(!File.exists(outdir1))
	File.makeDirectory(outdir1);

outdir = outdir1 + File.separator + "crops";
if(!File.exists(outdir))
	File.makeDirectory(outdir);

probabilities_dir = outdir1 + File.separator + "probabilities";
if(!File.exists(probabilities_dir))
	File.makeDirectory(probabilities_dir);
	
setBatchMode(true);


// *** MAIN ***

// Get the list of files in the input folder that end with the right extension
filelist = getFileList(ssif_folder);
filelist = Array.filter(filelist, file_ext);

// Split the filenames into case IDs and channel names (stored in three new arrays: "cases", "panels" and "channels")
//     -> this considers that the file names are in the format: "CaseID + delimiter1 + Panel + delimiter2 + ChannelName + .fileExtension"
//     -> should be modified if the files are names differently
cases = newArray(filelist.length);
panels = newArray(filelist.length);
case_panel = newArray(filelist.length);
channels = newArray(filelist.length);
splitFileNames(filelist, cases, panels, case_panel, channels, delimiter1, delimiter2, file_ext);

// Get unique case and channel IDs
caseID = ArrayUnique(cases);
panelID = ArrayUnique(panels);
case_panelID = ArrayUnique(case_panel);
channelID = ArrayUnique(channels);

// Process all case folders and call the crop function
for(i=0; i < case_panelID.length; i++){
	print("Processing ...", case_panelID[i]);
	
	// Parse file names
	for (j=0; j < filelist.length; j++){
		parts1 = split(filelist[j], delimiter2);
		parts2 = split(case_panelID[i], delimiter1);

		// If the image corresponds to the right case, channel and panel, open the file and bin it
		if((parts1[0] == case_panelID[i]) && (parts1[1] == islet_channel)  && (parts2[1] == panel_to_crop)) {
			imagename = ssif_folder + File.separator + filelist[j] + file_ext;
			cropImage(imagename, case_panelID[i], outdir, crop_nb, cropX, cropY);
		}
	}
}

print("done");


// *** FUNCTIONS ***

// Function to crop images
function cropImage(imagename, cur_case_panel, outdir, crop_nb, cropX, cropY) {
	// Open the image
	run("Bio-Formats Importer", "open=[imagename] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT stitch_tiles");
	run("Enhance Contrast", "saturated=0.35");
	run("Bin...", "x=4 y=4 bin=Average");
	height = getHeight;
	width = getWidth;
	
	// Generate random crops
	for(k=0; k < crop_nb; k++){
		// Pick random coordinates
		x_coord = getCropCoordinates(width, cropX);
		y_coord = getCropCoordinates(height, cropY);
		makeRectangle(x_coord, y_coord, cropX, cropY);
		roiManager("Add");
	}
	crop_name = outdir + File.separator + cur_case_panel + "_";
	RoiManager.multiCrop(crop_name, "save tif");
	
	// Delete ROIs from ROI manager
	run("Select All");
	roiManager("Delete");
	
	close("*");
	run("Collect Garbage");
}



// Function to return randomly-chosen crop coordinates
function getCropCoordinates(dim, cropDim){
	dimcr = dim - 1000;  // 1000 pixels padding to avoid image borders
	coord = 1000 + random() * (dimcr - cropDim);
	coord = floor(coord);
	return(coord);
}


// Split the file name based on delimiters, store each part in a separate array
function splitFileNames(filelist, cases, panels, case_panel, channels, delim1, delim2, file_ext){
	for (i=0; i < filelist.length; i++){
		filelist[i] = substring(filelist[i], 0, filelist[i].length - file_ext.length);
		
		p1 = split(filelist[i], delim2);
		case_panel[i] = p1[0];
		channels[i] = p1[1];
		
		p2 = split(p1[0], delim1);
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