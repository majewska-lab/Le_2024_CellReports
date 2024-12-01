// Brendan's Routine for Sholl Analysis  -Brendan Whitelaw, 2020
// Editted Linh Le - 2023 
// Part 1: Pre-processing Z-stacks to isolate, theshold, and smooth microglia before Sholl


/*
 * Macro template to process multiple images in a folder
 * 
 * Run Sholl analysis on cropped thresholded microglia images in Batch (all .tifs in a folder)
 */

// Dialog box made using Script Parameters: https://imagej.net/Script_Parameters
// These are 'global' parameters, so you can use them inside of functions 
// below without explicitly mentioning them as inputs

#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Zmax Output directory", style = "directory") dir_zmax
#@ File (label = "Cropped Microglia directory", style = "directory") dir_microglia
#@ String (label = "File suffix", value = ".tif") suffix

#@ Boolean (label = "Run Despeckle?", value=true, persist=false) option_despeckle
#@ Boolean (label = "Run Rolling-ball subtraction?", value=true, persist=false) option_subtraction
#@ Float (label = "Radius", value=50, persist=false) subtraction_radius
#@ Integer (label = "Small object exclusion size?", value=10, persist=false) option_exclusionsize
#@ Boolean (label = "Individual microglia thresholding?", value=true, persist=false) option_individualMG
#@ Boolean (label = "Soma analysis?", value=false, persist=false) option_soma

processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, dir_zmax, dir_microglia, suffix, list[i]);
	
}


// Runs analysis on each file in input folder with correct extension
function processFile(input, dir_zmax, dir_microglia,suffix, file) {
	// Do the processing here by adding your own code.
	// print("Processing: " + file);
	path = input + File.separator + file;
	print(file);
	// Run the isolate_microglia function
	isolate_microglia(path,dir_zmax,dir_microglia);
	
}

// Isolate microglia function to prepare images for Sholl analysis 
function isolate_microglia(path,dir_zmax,dir_microglia) {
	// Takes a z-stack (.nd2 as written) and performs semi-automated pre processing and cropping
	// Inputs: path (file-path of .nd2 z-stack); dir_zmax output folder for Maxz;
	// dir_microglia output for cropped microglia

close("*");
open(path); // open the file
stack_name = File.getName(path);
selectWindow(stack_name);

//// Choose channel for analysis and isolate that channel
//getDimensions(width, height, channels, slices, frames);
//if (channels > 1){
//Dialog.createNonBlocking("Choose the channel...");
//Dialog.addString("Channel: ", "");
//Dialog.show();
//new_channel = Dialog.getString();
//run("Split Channels");
//stack_name = "C"+new_channel+"-"+stack_name;
//selectWindow(stack_name);
//close("\\Others");
//}

// Despeckle and Smooth image if previously selected
if (option_despeckle == true){
	run("Despeckle", "stack");
}
if (option_subtraction == true){
	run("Subtract Background...", "rolling="+d2s(subtraction_radius,0) +" stack");
}
print("Done preprocessing");



// Generate a dialog box to determine slices for MaxZ projection and performe projections
Dialog.createNonBlocking("Choose the top and bottom Z slices...");
Dialog.addString("Lowest: ", "");
Dialog.addString("Highest: ", "");
Dialog.show();
lowest = Dialog.getString();
highest = Dialog.getString();
selectWindow(stack_name);
run("Z Project...", "start=" + lowest + " stop=" + highest + " projection=[Max Intensity]");



// Save Zmax projection 
zmax_name = "MAX_"+stack_name;
selectWindow(zmax_name); // Make sure to select the Max projection window
zmax_name = replace(zmax_name,suffix,".tif"); // rename max file name for saving
saveAs("tiff", dir_zmax + File.separator + zmax_name); // save MAX projection as .tif
// Note: zmax_name is also the name of the max projection window. stack_name is that of the z-stack


// Select microglial cell bodies (can use Z-stack to help find ones in the middle)
setTool("multipoint");
waitForUser("Click on cell bodies in middle of stack");
selectWindow(zmax_name);
Overlay.addSelection;
run("Select None");


// Use ROI manager to draw outline each microglia with polygons
close(stack_name);
roiManager("reset");
setOption("Show All", true);
setTool("polygon");
waitForUser("Outline microglia and add to ROI manager (press 't')");
ROIset_name = replace(zmax_name,".tif","_ROIset.zip");
roiManager("Save", dir_zmax+ File.separator+ ROIset_name);
run("Select None");
Overlay.remove // Removes the cell body overlays


// IF: Applying global threshold to image, cropping microglia, selecting cell center
if (option_individualMG == false){
	print("Global threshold");
	thresh_name = "thresholded_image";
	run("Duplicate...", "title="+thresh_name);
	selectWindow(thresh_name);
	
	run("Threshold...");
	setAutoThreshold("Default dark");
	waitForUser("Apply threshold... then click OK");
	run("Convert to Mask");
	
	// Remove small particles (size determined in options)
	run("Make Binary");
	roiManager("reset");
	run("Analyze Particles...", "size=0-"+d2s(option_exclusionsize,0)+" pixel show=Nothing add");
	roiManager("deselect");
	setColor(0);
	roiManager("fill");
	roiManager("reset");

	// Re-load microglia ROIs and then crop them out
	roiManager("open", dir_zmax+ File.separator+ ROIset_name)
	num_mg = roiManager("count");
	num_mg_string = d2s(num_mg, 0);
	setBackgroundColor(255, 255, 255); // Important after thresholiding the backround is white...
	// so when you 'Clear Outside' the ROI below it needs to also be white
	microglia_file_name = replace(zmax_name,".tif","_mgROI_");
	
	for(i = 0; i < num_mg; i++) {
		t_frame = i+1;
		roiManager("select",i);
		run("Duplicate...","cropped_microglia_temp");
		run("Clear Outside");
		run("Select None");
		setTool("point");
		waitForUser("Click center of cell body... then click OK");
		Overlay.addSelection;
		num_ROI = d2s(i+1,0);
		new_file_name = microglia_file_name + num_ROI + "of" + num_mg_string + ".tif";
		saveAs("tiff",dir_microglia + File.separator + new_file_name);
		close(new_file_name);
		selectWindow(thresh_name);
	}
	} // End if global threshold

// IF: Thresholding individual microglia, then cropping and selecting cell center
else{
	print("Individual MG Threshold");
	// Re-load microglia ROIs and then crop them out
	roiManager("reset");
	roiManager("open", dir_zmax+ File.separator+ ROIset_name);
	num_mg = roiManager("count");
	num_mg_string = d2s(num_mg, 0);
	setBackgroundColor(0, 0, 0); // To make background after clear outside black
	// so when you 'Clear Outside' the ROI below it needs to also be white
	microglia_file_name = replace(zmax_name,".tif","_mgROI_");

	for(i = 0; i < num_mg; i++) {
		t_frame = i+1;
		roiManager("reset");
		roiManager("open", dir_zmax+ File.separator+ ROIset_name);
		roiManager("select",i);
		run("Duplicate...", "title=cropped_microglia_temp");
		run("Clear Outside");
		run("Select None");
		selectWindow("cropped_microglia_temp");
		threshMG_name = "thresholded_MG";
		run("Duplicate...", "title="+threshMG_name);
		// Threshold individual microglia
		selectWindow(threshMG_name);
		run("Threshold...");
		setAutoThreshold("Default dark");
		waitForUser("Apply threshold... then click OK");
		run("Convert to Mask");

		// Remove small particles (size determined in options)
		selectWindow(threshMG_name);
		run("Make Binary");
		roiManager("reset");
		wait(200);
		selectWindow(threshMG_name);
		run("Analyze Particles...", "size=0-"+d2s(option_exclusionsize,0)+" pixel show=Nothing add");
		roiManager("deselect");
		setColor(0,0,0);
		roiManager("fill");
		roiManager("reset");

		// Select center of microglial cell body
		setTool("point");
		waitForUser("Click center of cell body... then click OK");
		Overlay.addSelection;
		num_ROI = d2s(i+1,0);
		new_file_name = microglia_file_name + num_ROI + "of" + num_mg_string + ".tif";
		saveAs("tiff",dir_microglia + File.separator + new_file_name);
		close(new_file_name);
		close("cropped_microglia_temp");
		selectWindow(zmax_name);
	}


}

// IF: Running soma analysis. Also includes using approximate soma radius as starting Sholl radius
if (option_soma == true) {
	selectWindow(zmax_name);
	print("Soma analysis");
	// Re-load microglia ROIs and then crop them out
	roiManager("reset");
	roiManager("open", dir_zmax+ File.separator+ ROIset_name);
	num_mg = roiManager("count");
	num_mg_string = d2s(num_mg, 0);
	microglia_file_name = replace(zmax_name,".tif","_mgROI_"); // Define microglia file name (maybe redundant)
	setBackgroundColor(0, 0, 0); // To make background after clear outside black
	// so when you 'Clear Outside' the ROI below it needs to also be white
	
	for(i = 0; i < num_mg; i++) {
		t_frame = i+1;
		selectWindow(zmax_name);
		roiManager("reset");
		roiManager("open", dir_zmax+ File.separator+ ROIset_name);
		// Select and isolate microglial cell
		roiManager("select",i);
		run("Duplicate...", "title=cropped_microglia_temp");
		run("Clear Outside");
		run("Select None");
		selectWindow("cropped_microglia_temp");
		// Draw polygon around soma and store the vertices
		setTool("polygon");
		waitForUser("Outline Soma... then click OK");
		getSelectionCoordinates(xpoints, ypoints);
		// Close image (non-thresholded microglia)
		close("cropped_microglia_temp");
		// Re-open thresholded microglia image, remove point overlay, and add soma outline overlay, and overwrite file
		num_ROI = d2s(i+1,0);
		new_file_name = microglia_file_name + num_ROI + "of" + num_mg_string + ".tif";
		new_file_path = dir_microglia + File.separator + new_file_name;
		open(new_file_path);
		Overlay.remove
		makeSelection("polygon", xpoints, ypoints);
		Overlay.addSelection;
		saveAs("tiff",new_file_path);
		close(new_file_name);
		selectWindow(zmax_name);
	}

}
roiManager("reset");
close("*");
print("done");

}


