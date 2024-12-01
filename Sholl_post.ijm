/*
 * Macro template to process multiple images in a folder
 * 
 * Run Sholl analysis on cropped thresholded microglia images in Batch (all .tifs in a folder)
 */

#@ File (label = "Thresholded Microglia", style = "directory") input
#@ File (label = "Spreadsheets", style = "directory") output
#@ String (label = "File suffix", value = ".tif") suffix
#@ Boolean (label = "Soma analysis?", value=true, persist=false) option_soma

// See also Process_Folder.py for a version of this code
// in the Python scripting language.


processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
}

function processFile(input, output, file) {
	// Do the processing here by adding your own code.
	// print("Processing: " + file);
	setBatchMode("hide");
	open(input + File.separator + file);
setAutoThreshold("Default dark");
setThreshold(255, 255, "raw");
setOption("BlackBackground", true);
run("Convert to Mask");

	// Each image file (thresholded individual microglia) has an overlay saved with it.
	
	// IF: not doing soma analysis, overlay is a point selected at center of microglia. Use for Sholl
	if (option_soma == 0){
	roiManager("reset");
	nROIs = Overlay.size;
	// print("Number ROIs: " + d2s(nROIs,0));
	Overlay.activateSelection(nROIs-1); // Adds overlay to current selection
	roiManager("add"); // Adds current selection to ROI manager
	run("Select None"); // Unselect overlay object
	roiManager("select", 0); // Select the object in the ROI manager
	starting_rad = 0; // Start sholl from center
	}
	// IF: Soma analysis is selected: use polygon selection from previous macro to find center and starting radius
	if (option_soma == 1){
		roiManager("reset");
		nROIs = Overlay.size;
		Overlay.activateSelection(nROIs-1); // Adds overlay to current selection
		roiManager("add"); // Adds current selection to ROI manager
		run("Select None"); // Unselect overlay object
		roiManager("select", 0); // Select the object in the ROI manager
		// Get centroid and size of elipse axes of polygon
		List.setMeasurements;
		soma_x = List.getValue("X"); // centroid, x-coordinate
		soma_y = List.getValue("Y"); // centroid, y-coordinate
		soma_maj = List.getValue("Major"); // length of major axis (full length)
		soma_min = List.getValue("Minor"); // length of minor axis (full length)
		soma_area = List.getValue("Area"); // area of soma
		starting_rad = ((soma_min)/2); // Start Sholl from approximate edge of soma: radius is average of the short axis 
		roiManager("reset");
		toUnscaled(soma_x,soma_y);
		makePoint(soma_x,soma_y); // Make point at centroid of soma
		roiManager("add"); // Add center point to ROI manager
		roiManager("select", 0); // Select center point
		}
	
run("Legacy: Sholl Analysis (From Image)...", "starting="+d2s(starting_rad,2)+" ending=100 radius_step=2 #_samples=1 " + 
		"integration=Mean enclosing=1 #_primary=2 linear polynomial=[Best fitting degree] " + 
		"normalizer=Area save directory=["+output+"] do");
	
	close("*");
}




