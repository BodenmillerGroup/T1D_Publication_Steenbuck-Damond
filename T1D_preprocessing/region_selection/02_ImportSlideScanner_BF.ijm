/*
 * PROCESS BRIGHTFIELD SLIDE SCANNER IMAGES
 * Nicolas Damond - 13 Februrary 2020
 */

// *** GOALS ***
// 1. Process brightfield tissue images so that they can be loaded in .MCD files.
// 2. Process low resolution slide images so that they can be loaded in .MCD files.
// 3. Generates "mask" and "registered" subfolders that will be used in the next steps.


// *** INPUT SETTINGS ***

// Brightfield images are expected to be named as "case_BF_panel.tif" (e.g. "1234_BF_Islet.tif")
// Slide images are expected to be named as "case_BF_panel_2_pt2.tif" (e.g. "1234_BF_Islet_2_pt2.tif")
// If the name format is different, modify the file extension, delimiter and adapt the script
filext = ".tif";
delimiter = "_";

// Prompt the user to provide an input folder
input = getDirectory("Input Directory"); // "BASE" folder
input_dir = input + File.separator + "czi_bf" + File.separator;

setBatchMode(true);


// *** MAIN ***

// Get the list of files in the input folder that end with the right extension
filelist = getFileList(input_dir);
processFileNames(input_dir, filelist);
// Remove tagged (NaN) files from the list (files which don't have the right extension)
filelist = Array.delete(filelist, NaN);

// Create additional folders
createFolders(filelist, input_dir, delimiter);

// Process the files.
for (i=0; i < filelist.length; i++){
	parts = split(filelist[i], delimiter);
	// Define if the image is tissue or slide
	if(parts.length == 3)
		processTissueImages(filelist[i], filext, delimiter);
	if(parts.length == 5)
		if(parts[4] == "pt2")
			processSlideImages(filelist[i], filext, delimiter);
}


print("done");


// *** FUNCTIONS ***

// Function to process tissue images (rotate, convert to 8-bit and bin)
function processTissueImages(file, filext, delim) {
	// Open and process the file
	imagename = input_dir + file + filext;
	print(imagename);
	run("Bio-Formats Importer", "open=[imagename] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT stitch_tiles");
	run("RGB Color");
	run("8-bit");
	run("Rotate... ", "angle=180 grid=1 interpolation=Bilinear");
	
	// Bin the image
	run("Bin...", "x=2 y=2 bin=Average");

	// Save the processed file
	p1 = split(file, delim);
	outdir = File.getParent(input_dir) + File.separator + p1[0] + File.separator + "mcd";
	if(!File.exists(outdir))
		File.makeDirectory(outdir);
	outname = outdir + File.separator + file + filext;
	saveAs("Tiff", outname);
	close("*");
	run("Collect Garbage");
}

// Function to process slide images (rotate, convert to 8-bit and crop)
function processSlideImages(file, filext, delim) {
	// Open and process the file
	imagename = input_dir + file + filext;
	print(imagename);
	run("Bio-Formats Importer", "open=[imagename] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT stitch_tiles");
	run("RGB Color");
	run("8-bit");
	run("Rotate... ", "angle=180 grid=1 interpolation=Bilinear");
	setMinAndMax(170, 200);
	run("Apply LUT");
	
	// Crop the image
	slide_w = getWidth();
	slide_h = getHeight();
	makeRectangle(slide_w-2100, (slide_h-700)/2, 2100, 700);
	run("Crop");
	
	// Save the processed file
	p1 = split(file, delim);
	outdir = File.getParent(input_dir) + File.separator + p1[0] + File.separator + "mcd";
	if(!File.exists(outdir))
		File.makeDirectory(outdir);
	outname = outdir + File.separator + file + filext;
	saveAs("Tiff", outname);
	close("*");
	run("Collect Garbage");
}

// Function to create additional folders
function createFolders(file, input_dir, delim){
	for (i=0; i < file.length; i++){
		part = split(file[i], delim);
		dirbase = File.getParent(input_dir) + File.separator + part[0];
		dirmask = dirbase + File.separator + "mask";
		if(!File.exists(dirmask))
			File.makeDirectory(dirmask);
		dirreg = dirbase + File.separator + "registered";
		if(!File.exists(dirreg))
			File.makeDirectory(dirreg);
	}	
}
// List all files in the input folder.
//     -> if the extension is not right, modify to NaN (tag for deletion).
//     -> if the extension is right, remove it from the filename.
function processFileNames(input_dir, filelist) {
	for (i=0; i < filelist.length; i++){
		if(!endsWith(filelist[i], filext))
			filelist[i] = NaN;
		else
			filelist[i] = replace(filelist[i], filext, "");
	}
}