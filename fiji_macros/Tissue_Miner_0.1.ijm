// Tissue Miner 0.1 is tool kit to visualize and quantify the dynamic cell behavior in a 2D tissue



smFile="/media/project-raphael@mack/scripts/tissue_miner/workflow/tm.snkmk";

//print(exec("bash "+wf+" "+movieDbDir+" 2>&1"));


// GUI

Dialog.create("Tissue Miner 0.1");
Dialog.addMessage("Welcome to Tissue Miner 0.1 !\n\n \nPlease, click OK to analyze a time-lapse processed in Tissue Analyzer");
Dialog.addCheckbox("Use default workflow engine setup", true);
Dialog.show();
isDefaultSM=Dialog.getCheckbox();

if (!isDefaultSM) {
	Dialog.create("Tissue Miner 0.1");
	Dialog.addString("Path to snake file:", "", 50);
	Dialog.show();
	smFile=Dialog.getString();
	}

print(smFile);

transfoType=newArray("Draw axis","Ellipse Major Axis");
data=newArray("pupal wing", "wing disk");

Dialog.create("Tissue Miner 0.1");
Dialog.addCheckbox("Rotate ?", false);
Dialog.addChoice("Method: ",transfoType,"rotate_ellipse");
Dialog.addCheckbox("Draw Rois ?", false);
Dialog.addChoice("Data type: ",data,"pupal wing");
Dialog.addNumber("Pixel size (optional): ",0.207);
Dialog.show();
isRotate=Dialog.getCheckbox();
transfo=Dialog.getChoice();
isDrawRoi=Dialog.getCheckbox();
roiDef=Dialog.getChoice();
pxsize=Dialog.getNumber();


//isRotate=false;

if (isRotate)  {
	// Define trafo methods
	//transfoType=newArray("Draw axis","Ellipse Major Axis");

	// Select trafo method and image file to rotate
	//transfo=select_trafo(transfoType);
	path=select_image_file();

	// Define image directory and file name separately
	dir=File.getParent(path)+File.separator();
	print("dir",dir);
	fileName=File.getName(path);
	
	// Apply trafo
	if (transfo=="Draw axis"){
		rotate_drawAxis(dir, fileName); // rotate the tissue by drawing an axis
	}
	else if (transfo=="Ellipse Major Axis") {
		rotate_ellipse(dir, fileName); // rotate the tissue using fitted ellise
	}		
}

//isDrawRoi=true;
if (isDrawRoi) {
	// Set up paths to files and ROI
	inputFile = File.openDialog("Choose an image to draw ROI");
	folderPath=File.getParent(inputFile)+File.separator();
	
	// GUI
	//Dialog.create("Tasks to be performed");
	//data=newArray("pupal wing", "wing disk");
	//Dialog.addChoice("Data type: ",data,"pupal wing");
	//Dialog.addNumber("Pixel size (optional): ",0.207);
	//Dialog.show();
	//data=Dialog.getChoice();
	//pxsize=Dialog.getNumber();
	
	// TODO: define data type and rois in a separate txt file called roi.conf
	// pupal wing: blade, interL1-L2, ...
	// wing disk: DV, AP, ...
	// parse this file and generate the appropriate GUI
	if (roiDef=="pupal wing") {
		listROI=newArray("blade", "interL1-L2", "L2", "interL2-L3", "L3", "proxInterL3-L4", "distInterL3-L4", "L4", "proxInterL4-L5", "postCV", "distInterL4-L5", "L5", "postL5", "HBinterface", "hinge");	
	}
	else {
		listROI=newArray("DV", "AP", "DorsInLeft", "DorsInRight", "VentrInRight", "VentrInLeft", "DorsOutLeft", "DorsOutRight", "VentrOutRight", "VentrOutLeft");	
	}
	
	isUpdate=0;
	
	if (File.exists(folderPath+"RoiSet.zip")) {
		isUpdate=getBoolean("RoiSet.zip already exists, update it ?");
		if (isUpdate) {
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
	
	

}


// Run the workflow
movieDbDir=getDirectory("Select a time-lapse directory");
//movieDbDir="/home/etournay/RawData/Test_snakemake/demo";
workflow="/media/project-raphael@mack/scripts/flywing_shear/utils/epc_shear_workflow.sh";

print(exec("bash "+workflow+" "+movieDbDir+" 2>&1"));


exit();
//print(exec("snakemake --snakefile /media/project-raphael@mack/scripts/flywing_shear/utils/epcSnakeFile.snkmk -j 2 all"));




/*
 * Functions for image rotation
 */



function select_image_file(){
	inputpath = File.openDialog("Please, select an image to process");
	return inputpath;
}

function select_trafo(transfoType){
	Dialog.create("Type of transformation to be performed");
	Dialog.addChoice("Type: ",transfoType,"rotate_ellipse");
	Dialog.show();
	transfo=Dialog.getChoice();
	return transfo;
}


function rotate_ellipse(dir, file){
	print("Order of points: SO1 in proximal hinge ; SO2 on L3; SO3 on L3; any pt in post. compartment");
	open(dir + file); // for png files, better use open();
	//run("Bio-Formats Importer", "open="+ dir + file+" autoscale color_mode=Default view=[Standard ImageJ] stack_order=Default");
	imID=getImageID();
	
	//gets image dimensions of original
	width = getWidth();
	height = getHeight(); 

	// Perform point selection + transformation
	do {
		selectImage(imID);
		setTool("line");
		//gets points on wings
		leftButton=16;
		rightButton = 4;
		insideROI = 32;
		ctrl = 2;
		x2=-1; y2=-1; z2=-1; flags2=-1;
		printed = false;
		currentpoint=0;
		maxpoints=4;
		xCoords = newArray(maxpoints);
		yCoords = newArray(maxpoints);
		while (currentpoint<maxpoints){
			if (!printed){
				print("Please select point " + (currentpoint+1));
				printed = true;	
			}
			getCursorLoc(x,y,z,flags);
			if (x!=x2 || y!=y2 || z!=z2 || flags!=flags2) {
				if (flags&leftButton!=0){
					print(x+" "+y+" "+z+" "+flags);
					xCoords[currentpoint] = x;
					yCoords[currentpoint] = y;	
					printed = false;
					currentpoint++;
				}	
			}
			x2=x; y2=y; z2=z; flags2=flags;
			//print(currentpoint);
		}
		
		wait(200);//necessary, otherwise picture does not get rotated
		
		// Compute linear regression on the 3 selected SO: estimate a and b such that y=ax+b using the mean square method
		//a=cov(x,y)/var(x)
		//b=y_mean - x_mean*a
		x_mean=(1/3)*(xCoords[0]+xCoords[1]+xCoords[2]);
		y_mean=(1/3)*(yCoords[0]+yCoords[1]+yCoords[2]);
		cov_xy=(1/3)*((xCoords[0]-x_mean)*(yCoords[0]-y_mean)+(xCoords[1]-x_mean)*(yCoords[1]-y_mean)+(xCoords[2]-x_mean)*(yCoords[2]-y_mean));
		var_x =(1/3)*((xCoords[0]-x_mean)*(xCoords[0]-x_mean) + (xCoords[1]-x_mean)*(xCoords[1]-x_mean) + (xCoords[2]-x_mean)*(xCoords[2]-x_mean));
		a=cov_xy/var_x;
		b=y_mean - x_mean*a;
		//print("a",a);
		
		// Calculate directed angle between L3 and horizontal
		L3angleradians=atan2(cov_xy,var_x); // [-PI, +PI]
		L3angledeg = L3angleradians * 180 / PI;
		print("L3angledeg",L3angledeg);
		
		// As Fiji performs a clockwise rotation for positive angles, invert angle sign
		L3angledeg = -L3angledeg;
		L3angleradians = -L3angleradians;
		
		// Get center positions before rotation
		rotationcenterx = width / 2;
		rotationcentery = height / 2;
		
		// Clockwise rotation based on the image center
		print("Rotation angle =" + L3angledeg);
		run("Rotate... ", "angle=&L3angledeg grid=1 interpolation=Bilinear fill enlarge");
		
		
		// Caculate point coordinates in the new system of coordinates (new orientation and new boundaries)
		widthnew = getWidth();
		heightnew = getHeight();
		rotationcenternewx = widthnew / 2;
		rotationcenternewy = heightnew / 2;
		for (i=0; i<maxpoints; i++){
			//Verschiebung des Rotationszentrums zum Ursprung
			x = (xCoords[i]-rotationcenterx);
			y = (yCoords[i]-rotationcentery);
			//Rotation der Koordinaten im alten Koordinatensystem
			xnew = x * cos(L3angleradians) - y * sin(L3angleradians);
			ynew = x * sin(L3angleradians) + y * cos(L3angleradians);
			//Rückverschiebung der Koordinaten
			xCoords[i] = xnew + rotationcenternewx;
			yCoords[i] = ynew + rotationcenternewy;	
			print("i: "+i + " x: " + xCoords[i] + " y: " + yCoords[i]);	
		}

		// Put the anterior compartment up if it isn't
		isFlip=false;
		if (yCoords[3]<yCoords[0]){
			print("Turning image around");
			run("Flip Vertically");
			isFlip=true;
			for (i=0; i<maxpoints; i++){
				yCoords[i] = abs(yCoords[i]-heightnew);
				print("i: "+i + " x: " + xCoords[i] + " y: " + yCoords[i]);		
			}		
		}
	success=getBoolean("Rotation OK?");
	if (!success) {run("Undo");}
	} while (success==0);
	
	f=File.open(dir+"transformation.txt");
	//HEADER
	// Angle (rad) postive sign is anticlockwise: (invert the sign)
	// "Rotation center coordinates in old system: rotationCenterOld_x and rotationCenterOld_y
	// "Rotation center coordinates in new system: rotationCenterNew_x and rotationCenterNew_y
	// "Put the anterior compartment up if it isn't: isVerticalFlip"
	header="Angle_rad rotationCenterOld_x rotationCenterOld_y rotationCenterNew_x rotationCenterNew_y IsVerticalFlip newImHeight method";
	print(f, header);
	sep=" ";
	values=toString(L3angleradians) + sep + rotationcenterx + sep + rotationcentery + sep + rotationcenternewx + sep + rotationcenternewy + sep + isFlip + sep + heightnew + sep + "rotate_ellipse";
	print(f, values);
	File.close(f);
	selectImage(imID);
	close();
}


function rotate_drawAxis(dir, file){
	
	print("Draw a line to be aligned on the horizontal");
	
	// Perform axis selection + transformation
	do {
		// open selected file
		open(dir + file); 
		imID=getImageID();
		
		// get image dimensions of original
		width = getWidth();
		height = getHeight(); 
		
		// manually draw points on tissue
		setTool("line");
		run("ROI Manager...");
		waitForUser("Please, select a line and click OK to validate");
		roiManager("Add");
		
		
		// get axis angle 
		roiManager("Measure");
		rotation=getResult("Angle", 0);
		print("angle with respect to horizontal:", rotation);
	
		// Get center positions before rotation
		rotationcenterx = width/2;
		rotationcentery = height/2;
		
		// Clockwise rotation based on the image center (y-axis pointing downward)
		print("Rotation angle =" + rotation);

		// Rotation can only be enlarged if no roi on top of image: workaround = duplicate the image and rotate the duplicated image
		selectImage(imID); run("Duplicate...", "rot template");rotImID=getImageID();
		selectImage(imID);close(); 
		
		run("Rotate... ", "angle=&rotation grid=1 interpolation=Bilinear fill enlarge");	
		
		// Update coordinates for the new system of coordinates (new orientation => new boundaries)
		widthnew = getWidth();
		heightnew = getHeight();
		rotationcenternewx = widthnew / 2;
		rotationcenternewy = heightnew / 2;
		
		success=getBoolean("Rotation OK?");
		if (!success) {selectImage(rotImID);close();}
		selectWindow("ROI Manager");run("Close");
		selectWindow("Results");run("Close");	
		
	} while (success==0); //END OF LOOP

	
	// convert rotation angle in radian	
	rotationRad=rotation*PI/180;
	
	// Save transformation in transformation.txt
	// HEADER
	// Angle (rad) postive sign is anticlockwise: (invert the sign)
	// "Rotation center coordinates in old system: rotationCenterOld_x and rotationCenterOld_y
	// "Rotation center coordinates in new system: rotationCenterNew_x and rotationCenterNew_y
	// "Put the anterior compartment up if it isn't: isVerticalFlip"
	f=File.open(dir+"transformation.txt");
	header="Angle_rad rotationCenterOld_x rotationCenterOld_y rotationCenterNew_x rotationCenterNew_y IsVerticalFlip newImHeight method";
	print(f, header);
	sep=" "; isFlip=0;
	values=toString(rotationRad) + sep + rotationcenterx + sep + rotationcentery + sep + rotationcenternewx + sep + rotationcenternewy + sep + isFlip + sep + heightnew + sep + "rotate_drawAxis";
	print(f, values);
	File.close(f);

	// Clean up memory
	selectImage(rotImID);
	close();
}


/*
 * Functions for drawing ROIs
 */

function define_ROI(listROI,outDir,file,pxsize){
	open(file);
	imID=getImageID();

	run("Set Scale...", "distance=1 known="+pxsize+" pixel=1 unit=microns");
	
	run("ROI Manager...");
	run("Labels...", "color=white font=12 show use draw bold");
	run("Colors...", "foreground=white background=white selection=green");
	
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
	}
	
	roiManager("Save", outDir+"RoiSet.zip")
	selectWindow("ROI Manager");
	run("Close");

	selectImage(imID);
	close();
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



