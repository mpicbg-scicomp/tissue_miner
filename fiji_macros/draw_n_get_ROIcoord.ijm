

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
	
	roiManager("Save", outDir+"RoiSet.zip");
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

function setDefaultTemplatesFile(templatesPathName)
{
	f = File.open(templatesPathName);
	print(f, "pupal wing,blade,interL1-L2,L2,interL2-L3,L3,proxInterL3-L4,distInterL3-L4,L4,proxInterL4-L5,postCV,distInterL4-L5,L5,postL5,HBinterface,hinge");
	//print(f, "wing disk,DV,AP,DorsInLeft,DorsInRight,VentrInRight,VentrInLeft,DorsOutLeft,DorsOutRight,VentrOutRight,VentrOutLeft");
	File.close(f);
}


function updateTemplateFile(templatesPathName, templates_str_list)
{
	// save the updated template file to disk 
	f = File.open(templatesPathName);
	for(i=0; i<templates_str_list.length; i++){
		print(f, templates_str_list[i]);
	}
	File.close(f);
}


// get a list of templates name

function getTemplatesName(templates_str_list){
	n_temp = templates_str_list.length;
	templatesName=newArray(n_temp);
	for(i=0; i<n_temp; i++){
	    chuncks=split(templates_str_list[i],",");
	    templatesName[i]=chuncks[0];
	}
	return templatesName;
}

// get the id of the template with name template name
function getTemplateId(templates_str_list, templateName)
{	
	id=-1;
	for(i=0; i<templates_str_list.length; i++){
	    chuncks=split(templates_str_list[i],",");
		if( chuncks[0]==templateName ){
			id = i;
		}
	}
	return id;
}

// create a dialog to input a string that can be given a default value
function updateTemplate(title, default)
{
	Dialog.create(title);
	Dialog.addMessage(title+":");
	Dialog.addString("", default, 100);			  
	Dialog.addMessage("please enter the template name and the name of the different\n"+
					  "structures separateded by a comma. For instance:\n"+
					  "pupal wing,blade,interL1-L2,L2,interL2-L3,L3 ... ");
	Dialog.show();
	return Dialog.getString();
}




// add item at the end of array
function appendToArray(array, new_item){
	aux_array = array;
	n= array.length;
	array = newArray( n+1 );
	for(i=0;i<n;i++){ array[i]=aux_array[i]; }
	array[n]=new_item;
	return array;
}

// remove item k from array and return array
function removeFromArray(array, k)
{
	aux = array;
	array = newArray(aux.length-1);
	for(i=0;i<k;i++){ array[i]=aux[i]; }
	for(i=k+1;i<aux.length;i++){ array[i-1]=aux[i]; }
	return array;
}




/*
 * ***************MAIN PROG**********************
 */
 
// Set up paths to files and ROI
fsep = File.separator();
inputFile = File.openDialog("Choose an image to draw ROI");
folderPath=File.getParent(inputFile)+fsep;

templatesPath= getDirectory("plugins");
templatesPathName = templatesPath+"Draw_andGetRoiCoord_templates.txt";


/////////////////////////////////////////////////////
// initialize the list of template //////////////////
/////////////////////////////////////////////////////

// if the templates file does not exist create one
if (File.exists(templatesPathName)!=1){
	setDefaultTemplatesFile(templatesPathName);
}
templates_fileStr = File.openAsString(templatesPathName);
templates_str_list = split(templates_fileStr, "\n");
templatesName = getTemplatesName(templates_str_list);

// if the file is empty set the file to default
if (templatesName.length==0){ 
	setDefaultTemplatesFile(templatesPathName);
	templates_fileStr = File.openAsString(templatesPathName);
	templatesName = getTemplatesName(templates_fileStr);
}
templates_str_list = split(templates_fileStr, "\n");


////////////////////////////////////////////////////////
// create the gui proposing the actions ////////////////
////////////////////////////////////////////////////////

actions = newArray("Draw rois for the selected template",
                   "Edit selected template",
                   "Remove selected template",
                   "Input a new template");
action = "";
pxsize = 0.207;
// redisplay the dialog as long as the action selection is related to adding, erasing, editing templates
while( (action!=actions[0]) | (templatesName.length==0) )
{
	// GUI
	Dialog.create("Tasks to be performed");
	Dialog.addChoice("Action: ",actions, actions[0]);
	Dialog.addChoice("Template : ",templatesName,templatesName[0]);
	Dialog.addNumber("Pixel size (optional): ",pxsize);
	Dialog.show();
	action = Dialog.getChoice();
	templateName = Dialog.getChoice();
	pxsize = Dialog.getNumber();

	if(action==actions[1]) //Edit the selected template
	{
		id = getTemplateId(templates_str_list, templateName);
		templates_str = updateTemplate( "Update template "+templateName, templates_str_list[id] );
		templates_str_list[id] = templates_str;
		IJ.log(templateName+" template was edited");
	}
	else if(action == actions[2]) //Erase the selected template
	{
		id = getTemplateId(templates_str_list, templateName);
		templates_str_list = removeFromArray(templates_str_list, id);
		IJ.log(templateName+" template was erased");
	}
	else if(action == actions[3]) //Add a new template
	{
		template_str = updateTemplate("Create a new template","");
		templates_str_list = appendToArray(templates_str_list, template_str);
		chuncks=split(template_str,",");
		IJ.log(chuncks[0]+" was created");
	}
	for(i=0; i<templates_str_list.length; i++){IJ.log("" + templates_str_list[i]); }
	updateTemplateFile(templatesPathName, templates_str_list);
	templatesName = getTemplatesName(templates_str_list);
}

// get the template selected
id = getTemplateId(templates_str_list, templateName);
template_str = templates_str_list[id];

// create an array from the template string
chuncks = split(template_str,",");
listROI = newArray(chuncks.length-1);
for(i=0; i<chuncks.length-1; i++){
	listROI[i] = chuncks[i+1];	
}


isRoiSet=0;

if (File.exists(folderPath+"RoiSet.zip")) {
	isRoiSet=getBoolean("RoiSet.zip already exists, update it ?");
	if (isRoiSet) {
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




