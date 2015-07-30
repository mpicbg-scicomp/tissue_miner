/*  Draw and get ROI coordinates
 *  Just drag and drop the macro to Fiji and Run it
 *  Select the image file on which to define the ROIs
 *  Define or update the ROIs by first clicking on Polygon Selection in Fiji
 */


function define_ROI(listROI,outDir,file,pxsize){
	open(file);
	imID=getImageID();

	run("Set Scale...", "distance=1 known="+pxsize+" pixel=1 unit=microns");
	
	run("ROI Manager...");
	run("Labels...", "color=white font=12 show use draw bold");
	run("Colors...", "foreground=white background=white selection=green");
	//run("Colors...", "foreground=white background=white selection=red");
	roiManager("UseNames", "true");
	setTool("polygon");
	for (i=0;i<lengthOf(listROI);i++){
		waitForUser("Please, select "+listROI[i]+" and click OK to validate");
		roiManager("Add");
		roiManager("Select", i);
		roiManager("Rename", listROI[i]);
		roiManager("Set Color", "red");wait(400);
		roiManager("Deselect");
		roiManager("Show All"); 
		//roiManager("Show None");
	}
	
	roiManager("Save", outDir+"RoiSet.zip")
	selectWindow("ROI Manager");
	run("Close");

	selectImage(imID);
	close();

	//f=File.open(outDir+"ROI_log.txt");
	//print(f, "Image used to define ROIs: "+File.getName(file));
	//File.close(f);
}

function update_ROI(filePath, ROIpath,pxsize){
	open(filePath);
	imID=getImageID();

	run("Set Scale...", "distance=1 known="+pxsize+" pixel=1 unit=microns");
	
	roiManager("Open", ROIpath);
	waitForUser("Please, update all ROI + click OK to validate");
	roiManager("Save", ROIpath)
	selectWindow("ROI Manager");
	run("Close");

	selectImage(imID);
	close();
}

function extract_ROIcoord(filePath, ROIpath){
	roiManager("Open", ROIpath);
	
	setBatchMode(true);
	open(filePath);
	imID=getImageID();
	parentFolder=File.getParent(ROIpath)+File.separator();
	outputFile=parentFolder+"LastFrameRoi.txt";
	
	imgFileName=File.getName(filePath);
	timefield=substring(imgFileName, lastIndexOf(imgFileName, "_")+1, lengthOf(imgFileName)-4);
	frame=parseInt(replace(timefield, "[a-z]+", ""));
	
	// initialize the outputFile
	f=File.open(outputFile);
	print(f, "FileName frame ROIname index x y");
	File.close(f);

	// iterate over all ROI
	nbroi=roiManager("count");
	for (i=0; i<nbroi; i++){
		roiManager("Select", i);
		roiName=getInfo("roi.name");
		getSelectionCoordinates(x, y);
		for (j=0; j<x.length; j++) {
	        	File.append(imgFileName + " " + frame + " " + roiName+" "+j+" "+toString(parseInt(x[j]))+" "+toString(parseInt(y[j])), outputFile);
	        	print(imgFileName + " " + frame + " " + roiName+" "+j+" "+toString(parseInt(x[j]))+" "+toString(parseInt(y[j])));
		}
		roiManager("Deselect");
	}
	selectImage(imID);
	close();
	selectWindow("ROI Manager");
	run("Close");
	setBatchMode(false);
}





/*
 * ***************MAIN PROG**********************
 */
// Set up paths to files and ROI
inputFile = File.openDialog("Choose an image to draw ROI");
folderPath=File.getParent(inputFile)+File.separator();

// GUI
Dialog.create("Tasks to be performed");
data=newArray("pupal wing", "wing disk");
Dialog.addChoice("Data type: ",data,"pupal wing");
Dialog.addNumber("Pixel size (optional): ",0.207);
Dialog.show();
data=Dialog.getChoice();
pxsize=Dialog.getNumber();

// TODO: define data type and rois in a separate txt file called roi.conf
// pupal wing: blade, interL1-L2, ...
// wing disk: DV, AP, ...
// parse this file and generate the appropriate GUI
if (data=="pupal wing") {
	listROI=newArray("blade", "interL1-L2", "L2", "interL2-L3", "L3", "proxInterL3-L4", "distInterL3-L4", "L4", "proxInterL4-L5", "postCV", "distInterL4-L5", "L5", "postL5", "HBinterface", "hinge");	
}
else {
	listROI=newArray("DV", "AP", "DorsInLeft", "DorsInRight", "VentrInRight", "VentrInLeft", "DorsOutLeft", "DorsOutRight", "VentrOutRight", "VentrOutLeft");	
}

rep=0;

if (File.exists(folderPath+"RoiSet.zip")) {
	rep=getBoolean("RoiSet.zip already exists, update it ?");
	if (rep) {
		update_ROI(inputFile, toString(folderPath)+"RoiSet.zip",pxsize);
		extract_ROIcoord(inputFile, toString(folderPath)+"RoiSet.zip");
	}
	else {}
}
else {
	define_ROI(listROI,folderPath,inputFile,pxsize);
	update_ROI(inputFile, toString(folderPath)+"RoiSet.zip",pxsize);
	extract_ROIcoord(inputFile, toString(folderPath)+"RoiSet.zip");
}

exit();




