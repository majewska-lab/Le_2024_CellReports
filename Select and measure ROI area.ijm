// Load and select ROI of images, measure ROI area for microglial density analysis 

// Dialog box made using Script Parameters: https://imagej.net/Script_Parameters
// These are 'global' parameters, so you can use them inside of functions 
// below without explicitly mentioning them as inputs

#@ File (label = "Image input directory", style = "directory") input
#@ File (label = "ROI output directory", style = "directory") dir_ROI

#@ String (label = "File suffix", value=".tif") suffix


processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input,dir_ROI, suffix, list[i]);
	
}


// Runs analysis on each file in input folder with correct extension
function processFile(input, dir_ROI, suffix, file) {
	// Do the processing here by adding your own code.
	// print("Processing: " + file);
	path = input + File.separator + file;
	print(file);
	// Run the ROIset function
	ROIset(path,dir_ROI);
	
}

function ROIset(path,dir_ROI) {
close("*");
open(path); // open the file
image_name = File.getName(path);
image_name_nosuff = File.getNameWithoutExtension(path);

run("Set Measurements...", "area redirect=None decimal=3");

selectWindow(image_name);
roiManager("reset");
ROIset_name = image_name_nosuff+"_ROIset.zip";
setOption("Show All", true);
setTool("polygon");
waitForUser("Outline Cortex/CA1/Corpus callosum and add to ROI manager (press 't')");
waitForUser("Rename all ROIs to reflect the corresponding brain region");
roiManager("Save", dir_ROI + File.separator+ ROIset_name); 
roiManager("reset");
if (File.exists(dir_ROI + File.separator + ROIset_name)==1) {
roiManager("Open", dir_ROI + File.separator+ ROIset_name); //Change to "Save" if need to create ROI
setBackgroundColor(0, 0, 0); // To make background after clear outside black
// so when you 'Clear Outside' the ROI below it needs to also be white
count = roiManager("count");
for (a = 0; a < count; a++) {
resultcounter = getValue("results.count");
roiManager("select", a)
ROIname = Roi.getName();
run("Measure");
setResult("Image",resultcounter, image_name);
setResult("ROI",resultcounter,ROIname);
updateResults();
}
run("Select None");

run("Close All");
print("done");
} 
else {print("Bad perfusion, Images were not processed");}

}
