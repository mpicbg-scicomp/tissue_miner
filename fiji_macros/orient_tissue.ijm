
/*
 * GUI MACRO to select one image on which to define the transformation
 * It must be the last image of the movie if rotate_wing_L4 method is used
 * It must be the first image of the movie if rotate_wing_SO method is used
 * 
 */

// TODO: make a flip function with possibility to redo the flip if the user isn't happy with what it did

// Select trafo method and image file to rotate
path=select_image_file();
// Define image directory and file name separately
dir=File.getParent(path)+File.separator();
//print("dir",dir);
fileName=File.getName(path);
open(dir + fileName);
imID=getImageID();

//transfoType=newArray("Draw an axis","Ellipse Major Axis");
transfoType=newArray("Draw a new X-axis", "Draw a new Y-axis");

//Dialog.create("Transformation options");
Dialog.create("Rotation options");
//Dialog.addCheckbox("Rotate ?", true);
Dialog.addChoice("Method: ",transfoType);
//Dialog.addCheckbox("Enable flip options ?", false);
Dialog.addNumber("Pixel size in micron (optional): ",0.207);
Dialog.show();
//isRotate=Dialog.getCheckbox();
isRotate=true;
rotMethod=Dialog.getChoice();
//isFlip=Dialog.getCheckbox();
isFlip=false;
pxsize=Dialog.getNumber();

selectImage(imID);
close();
	
// Cases: rotate only and rotate&flip
if (isRotate)  {

	// Apply trafo
	trafo=rotate_drawAxis(dir, fileName, pxsize, rotMethod); // rotate the tissue by drawing an axis
	print("trafo:", trafo);

	if (isFlip) {

		Dialog.create("Select the flips to apply");
		Dialog.addChoice("Flip: ", newArray("Horizontal (left-rigth exchange)","Vertical (up-side-down)", "Horizontal and Vertical", "None"));
		Dialog.show();
		flipMethod=Dialog.getChoice();

		if (flipMethod=="Horizontal (left-rigth exchange)") {
			run("Flip Horizontally");
			IsVerticalFlip=0;
			IsHorizontalFlip=1;
			trafo=update_trafo(dir, trafo, IsVerticalFlip, IsHorizontalFlip); print(trafo);
		}
		else if (flipMethod=="Vertical (up-side-down)") {
			run("Flip Vertically");
			IsVerticalFlip=1;
			IsHorizontalFlip=0;
			trafo=update_trafo(dir, trafo, IsVerticalFlip, IsHorizontalFlip); print(trafo);
		}
		else if (flipMethod=="Horizontal and Vertical") {
			run("Flip Horizontally");run("Flip Vertically");
			IsVerticalFlip=1;
			IsHorizontalFlip=1;
			trafo=update_trafo(dir, trafo, IsVerticalFlip, IsHorizontalFlip); print(trafo);
		}
		else if (flipMethod=="None") {
			IsVerticalFlip=0;
			IsHorizontalFlip=0;
			trafo=update_trafo(dir, trafo, IsVerticalFlip, IsHorizontalFlip); print(trafo);
		}
	}
	wait(3000); close();
}

// Case flip only
if (!isRotate && isFlip) {
	print("flip only");	
	open(dir + fileName); 
	imID=getImageID();
		
	// get image dimensions of original
	width = getWidth();
	height = getHeight();
	rotationcenterx = width/2;
	rotationcentery = height/2;
	rotationcenternewx = rotationcenterx;
	rotationcenternewy = rotationcentery;

	sep=" ";
	trafo= toString(0) + sep + rotationcenterx + sep + rotationcentery + sep + rotationcenternewx + sep + rotationcenternewy + sep + toString(0) + sep + toString(0) + sep + width + sep + height + sep + "FlipOnly"; 
	print(trafo);

	Dialog.create("Select the flips to apply");
	Dialog.addChoice("Flip: ", newArray("Horizontal (left-rigth exchange)","Vertical (up-side-down)", "Horizontal and Vertical","None"));
	Dialog.show();
	flipMethod=Dialog.getChoice();
		if (flipMethod=="Horizontal (left-rigth exchange)") {
		run("Flip Horizontally");
		IsVerticalFlip=0;
		IsHorizontalFlip=1;
		trafo=update_trafo(dir, trafo, IsVerticalFlip, IsHorizontalFlip); print(trafo);
	}
	else if (flipMethod=="Vertical (up-side-down)") {
		run("Flip Vertically");
		IsVerticalFlip=1;
		IsHorizontalFlip=0;
		trafo=update_trafo(dir, trafo, IsVerticalFlip, IsHorizontalFlip); print(trafo);
	}
	else if (flipMethod=="Horizontal and Vertical") {
		run("Flip Horizontally");run("Flip Vertically");
		IsVerticalFlip=1;
		IsHorizontalFlip=1;
		trafo=update_trafo(dir, trafo, IsVerticalFlip, IsHorizontalFlip); print(trafo);
	}
	else if (flipMethod=="None") {
		IsVerticalFlip=0;
		IsHorizontalFlip=0;
		trafo=update_trafo(dir, trafo, IsVerticalFlip, IsHorizontalFlip); print(trafo);
	}
	f=File.open(dir+"transformation.txt");
	header="Angle_rad rotationCenterOld_x rotationCenterOld_y rotationCenterNew_x rotationCenterNew_y IsVerticalFlip IsHorizontalFlip newWidth newImHeight method";
	print(f, header);
	print(f, trafo);
	File.close(f);

	wait(3000); close();
}

exit();



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

function update_trafo(dir, trafo, IsVerticalFlip,IsHorizontalFlip){
	sep=" ";
	fields=split(trafo, sep);
	fields[5]=IsVerticalFlip;
	fields[6]=IsHorizontalFlip;
	newTrafo = toString(fields[0]) + sep + fields[1] + sep + fields[2] + sep + fields[3] + sep + fields[4] + sep + fields[5] + sep + fields[6] + sep + fields[7] + sep + fields[8] + sep + fields[9];

	f=File.open(dir+"transformation.txt");
	header="Angle_rad rotationCenterOld_x rotationCenterOld_y rotationCenterNew_x rotationCenterNew_y IsVerticalFlip IsHorizontalFlip newImWidth newImHeight method";
	print(f, header);
	print(f, newTrafo);
	File.close(f);
	
	return newTrafo;
}
	

function rotate_drawAxis(dir, file, pxsize,rotMethod){
	
	//print("Draw a line to be aligned on the horizontal");
	
	// Perform axis selection + transformation
	do {
		// open selected file
		open(dir + file); 
		imID=getImageID();
		
		// get image dimensions of original
		width = getWidth();
		height = getHeight(); 

		// set scale in microns based on user-defined pixel size
		run("Set Scale...", "distance=1 known="+pxsize+" unit=micron");
		
		// manually draw points on tissue
		setTool("line");
		run("ROI Manager...");
		waitForUser("Please, select a line and click OK to validate");
		roiManager("Add");
	
		
		// get axis angle 
		roiManager("Measure");
		rotation=getResult("Angle", 0);
		if (rotMethod=="Draw a new X-axis") {
			print("angle with respect to horizontal:", rotation);
		}
		else {
			rotation=rotation-90;
			print("angle with respect to vertical:", rotation);
		}
		
	
		// Get center positions before rotation
		rotationcenterx = width/2;
		rotationcentery = height/2;
		
		// Clockwise rotation based on the image center (y-axis pointing downward)
		print("Rotation angle =" + rotation);

		// Rotation can only be enlarged if no roi on top of image: workaround = duplicate the image and rotate the duplicated image
		selectImage(imID); run("Duplicate...", "rot template");rotImID=getImageID();
		selectImage(imID);close(); 

		setBackgroundColor(0,0,0);
		run("Rotate... ", "angle=&rotation grid=1 interpolation=Bilinear fill enlarge");	
		
		// Update coordinates for the new system of coordinates (new orientation => new boundaries)
		// Rotate the four images corners A, B, C, D by the rotation matrix R(theta)
		// convert rotation angle in radian	
		rotationRad=rotation*PI/180;
		A=newArray(-width/2, -height/2);
		B=newArray(-1*-width/2, -height/2);//"-1*-" is needed in first place because of bug in macro language
		C=newArray(-1*-width/2, height/2);
		D=newArray(-width/2, height/2);

		// Clockwise rotation but bounded box would be the same by symmetry if counter-clockwise rotation
		Arot=newArray(A[0]*cos(rotationRad)-A[1]*sin(rotationRad),A[0]*sin(rotationRad)+A[1]*cos(rotationRad));
		Brot=newArray(B[0]*cos(rotationRad)-B[1]*sin(rotationRad),B[0]*sin(rotationRad)+B[1]*cos(rotationRad));
		Crot=newArray(C[0]*cos(rotationRad)-C[1]*sin(rotationRad),C[0]*sin(rotationRad)+C[1]*cos(rotationRad));
		Drot=newArray(D[0]*cos(rotationRad)-D[1]*sin(rotationRad),D[0]*sin(rotationRad)+D[1]*cos(rotationRad));

		// counter-clockwise rotation
		//Arot=newArray(A[0]*cos(rotation)+A[1]*sin(rotation),-A[0]*sin(rotation)+A[1]*cos(rotation));
		//Brot=newArray(B[0]*cos(rotation)+B[1]*sin(rotation),-B[0]*sin(rotation)+B[1]*cos(rotation));
		//Crot=newArray(C[0]*cos(rotation)+C[1]*sin(rotation),-C[0]*sin(rotation)+C[1]*cos(rotation));
		//Drot=newArray(D[0]*cos(rotation)+D[1]*sin(rotation),-D[0]*sin(rotation)+D[1]*cos(rotation));

		xCoord=newArray(Arot[0],Brot[0],Crot[0],Drot[0]); Array.sort(xCoord);
		yCoord=newArray(Arot[1],Brot[1],Crot[1],Drot[1]); Array.sort(yCoord);


		boundedBoxDim=newArray(-1*-floor(abs(xCoord[0])+xCoord[3]+1),floor(abs(yCoord[0])+yCoord[3])+1);		
		
		widthnew = getWidth();
		heightnew = getHeight();
		//rotationcenternewx = widthnew / 2;
		//rotationcenternewy = heightnew / 2;

		rotationcenternewx = boundedBoxDim[0] / 2;
		rotationcenternewy = boundedBoxDim[1] / 2;

		print("Fiji bounded box dim:         ", widthnew, heightnew);
		print("Calculated bounded box dim:", boundedBoxDim[0],boundedBoxDim[1]);
		
		success=getBoolean("Rotation OK?");
		if (!success) {selectImage(rotImID);close();}
		selectWindow("ROI Manager");run("Close");
		selectWindow("Results");run("Close");	
		
	} while (success==0); //END OF LOOP

	
	
	
	// Save transformation in transformation.txt
	// HEADER
	// Angle (rad) postive sign is anticlockwise: (invert the sign)
	// "Rotation center coordinates in old system: rotationCenterOld_x and rotationCenterOld_y
	// "Rotation center coordinates in new system: rotationCenterNew_x and rotationCenterNew_y
	// "Put the anterior compartment up if it isn't: isVerticalFlip"
	f=File.open(dir+"transformation.txt");
	header="Angle_rad rotationCenterOld_x rotationCenterOld_y rotationCenterNew_x rotationCenterNew_y IsVerticalFlip IsHorizontalFlip newImWidth newImHeight method";
	print(f, header);
	sep=" "; isVerticalFlip=0; isHorizontalFlip=0;
	values=toString(rotationRad) + sep + rotationcenterx + sep + rotationcentery + sep + rotationcenternewx + sep + rotationcenternewy + sep + isVerticalFlip + sep + isHorizontalFlip + sep + boundedBoxDim[0] + sep + boundedBoxDim[1] + sep + "drawAxis";
	print(f, values);
	File.close(f);

	// Clean up memory
	//selectImage(rotImID);
	//close();

	return values; 
}



exit();


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



// DEPRECATED CODE
transfoType=newArray("rotate_wing_L4","rotate_wing_SO","rotate_wingDisk");

transfo=transfo_gui(transfoType);
path=select_file_gui();
//path="/media/project_raphael@fileserver/FromNatalie/WingDisc_20120805_4.tif";

dir=File.getParent(path)+File.separator();
fileName=File.getName(path);

// DEBUG
//dir ="/media/project_raphael@fileserver/FromNatalie/";
//fileName="WingDisc_20120805.tif";



if (transfo=="rotate_wing_L4"){
	rotate_wing_L4(dir, fileName); // rotate the wing by selecting points on vein L4 + end of L5 on the last image of the movie
}
else if (transfo=="rotate_wing_SO") {
	rotate_wing_SO(dir, fileName); // rotate the wing by selecting 3 Sensory Organs on L3 on the first image of the movie + one point in the posterior compartment	
}
else if (transfo=="rotate_wingDisk"){
	rotate_wing_disk(dir, fileName); // rotate the wing disk by selecting a line through the AP boundary and one point in the dorsal compartment
}




function select_file_gui(){
	inputpath = File.openDialog("Please, select an image to process");
	return inputpath;
}

function transfo_gui(transfoType){
	Dialog.create("Type of transformation to be performed");
	Dialog.addChoice("Type: ",transfoType,"rotate_wing_SO");
	Dialog.show();
	transfo=Dialog.getChoice();
	return transfo;
}


function rotate_wing_SO(dir, file){
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
	values=toString(L3angleradians) + sep + rotationcenterx + sep + rotationcentery + sep + rotationcenternewx + sep + rotationcenternewy + sep + isFlip + sep + heightnew + sep + "rotate_wing_SO";
	print(f, values);
	File.close(f);
	selectImage(imID);
	close();
}


function rotate_wing_L4(dir, file){
	
	print("Order of points: ACV L4; end L4; PCV L4; end L5");

	run("Bio-Formats Importer", "open="+ dir + file+" autoscale color_mode=Default view=[Standard ImageJ] stack_order=Default");
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

		
		// Get L4 vector coordinates
		L4vectorx = xCoords[1]-xCoords[0];
		L4vectory = yCoords[1]-yCoords[0];

		// build a unitary vector out of L4
		L4vectorlength = sqrt(L4vectorx*L4vectorx+L4vectory*L4vectory);
		L4uniformvectorx = L4vectorx / L4vectorlength;
		L4uniformvectory = L4vectory / L4vectorlength;

		// Define horizontal unitary vector of the image system of coordinates
		horizontalx = 1;
		horizontaly = 0;

		// Calculate undirected angle between L4 and horizontal (from the scalar product)
		L4angleradians = acos(L4uniformvectorx * horizontalx + L4uniformvectory * horizontaly);
		L4angledeg = L4angleradians * 180 / PI;

		// As Fiji performs a clockwise rotation for positive angles, invert angle sign if end of L4 in pointing down
		if (yCoords[0] < yCoords[1]) {
			L4angledeg = -L4angledeg;
			L4angleradians = -L4angleradians;
		}
		
		// Get center positions before rotation
		rotationcenterx = width / 2;
		rotationcentery = height / 2;
		
		// Clockwise rotation based on the image center
		print("Rotation angle =" + L4angledeg);
		run("Rotate... ", "angle=&L4angledeg grid=1 interpolation=Bilinear fill enlarge");
		
		
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
			xnew = x * cos(L4angleradians) - y * sin(L4angleradians);
			ynew = x * sin(L4angleradians) + y * cos(L4angleradians);
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
	values=toString(L4angleradians) + sep + rotationcenterx + sep + rotationcentery + sep + rotationcenternewx + sep + rotationcenternewy + sep + isFlip + sep + heightnew + sep +"rotate_wing_L4";
	print(f, values);
	File.close(f);
	selectImage(imID);
	close();

}



function rotate_wing_disk(dir, file){
	print("Draw a line on the AP boundary starting in the dorsal compartment");
	open(dir + file); // for png files, better use open();
	//run("Bio-Formats Importer", "open="+ dir + file+" autoscale color_mode=Default view=[Standard ImageJ] stack_order=Default");
	imID=getImageID();
	
	//gets image dimensions of original
	width = getWidth();
	height = getHeight(); 

	// Perform point selection + transformations	
	do {
		selectImage(imID);
		setTool("line");
		
		//gets points on wings
		run("ROI Manager...");
		waitForUser("Please, select an AP line and click OK to validate");
		roiManager("Add");
		roiManager("Measure");
		angle=getResult("Angle", 0);
		print("angle with respect to horizontal and first point of the line:", angle);
		getSelectionCoordinates(x, y);
		for (j=0; j<x.length; j++) {
	        	//File.append(file + " " + roiName+" "+j+" "+toString(parseInt(x[j]))+" "+toString(parseInt(y[j])), outputFile);
	        	print(file +" "+j+" "+toString(parseInt(x[j]))+" "+toString(parseInt(y[j])));
		}
		roiManager("Show All");roiManager("Show None");
		selectWindow("Results");
		run("Close");
		selectWindow("ROI Manager");
		run("Close");
		
		// Calculate correction angle to set the line at -90 degrees (AP is vertical, dorsal compartment up)
		//if (angle<-90 && angle >=-180) {rotation=angle+90;}
		//else if(angle<=0 && angle >-90) {rotation=angle+90;}
		//else if(angle>0 && angle <=90) {rotation=angle+90;}
		//else if(angle>90 && angle <=180) {rotation=angle+90;}
		rotation=angle+90; // always valid because 'angle' is oriented
		rotationRad=rotation*PI/180;
		print("rotation",rotation);
	
		// Get center positions before rotation
		rotationcenterx = width / 2;
		rotationcentery = height / 2;
		
		// Clockwise rotation based on the image center
		print("Rotation angle =" + rotation);
		
		run("Rotate... ", "angle=&rotation grid=1 interpolation=Bilinear fill enlarge");
		
		// Caculate point coordinates in the new system of coordinates (new orientation and new boundaries)
		widthnew = getWidth();
		heightnew = getHeight();
		rotationcenternewx = widthnew / 2;
		rotationcenternewy = heightnew / 2;
		
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
	sep=" "; isFlip=0;
	values=toString(rotationRad) + sep + rotationcenterx + sep + rotationcentery + sep + rotationcenternewx + sep + rotationcenternewy + sep + isFlip + sep + heightnew + sep + "rotate_wing_disk";
	print(f, values);
	File.close(f);
	selectImage(imID);
	close();
}













