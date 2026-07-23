//directories
day = "2026_7_6";
basedir = "C:/Users/f008p47/Documents/Addition_Assay_Actin_Dynamics/";
onedrive_path = "C:/Users/f008p47/OneDrive/dartmouth_work/"


//ROI arrays
roiNames = newArray("TopLeft","TopRight","BottomLeft","BottomRight");
xCoords = newArray(0, 700, 0, 700);
yCoords = newArray(0, 0, 700, 700);
//mkdir day directory in onedrive
File.makeDirectory(onedrive_path + day);
// 1.INPUT...
inputdir = basedir + day + "/";

fileList = getFileList(inputdir);

// 2. Loop trough files
for ( i =0; i < fileList.length; i++) {
    if (endsWith(fileList[i], ".nd2")){
	     filePath = inputdir + fileList[i];
         // Open the image
         run ("Bio-Formats Importer", "open=[" + filePath + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
         Stack.setChannel(1);
	     run("Delete Slice", "delete=channel");
	     // set image variables 
	     mainImageID = getImageID();
	     mainTittle = getTitle();
	     imageName = File.getName(mainTittle);
	     baseName = replace(imageName, ".nd2", "");
	     if ( baseName == "A_FROG97-1_rGDB-GFP750nm_Utr594100nm_REP1") {
	     	run("Canvas Size...", "width=1400 height=1400 position=Top-Right zero");
	     }
	     else {
	     	run("Canvas Size...", "width=1400 height=1400 position=Center zero");
	     }
	     Stack.getDimensions(width, height, channels, slices, frames);
         // enchance contrast for each frame individually and rejoin them
         
     setBatchMode(true);
     for (f = 1; f <= frames; f++) {
         selectWindow(mainTittle);
         Stack.setFrame(f);
         
         
         padF = "" + f;
         while (lengthOf(padF) < 4) padF = "0" + padF;
         
         run("Duplicate...", "title=norm_frame_" + padF);
         run("Enhance Contrast", "saturated=0.4 normalize");
     }
     
     // Close original stack
     selectWindow(mainTittle);
     close();
     
     // join all frames to a new stack
     run("Images to Stack", "name=["+mainTittle+"] title=norm_frame_ use");
     setBatchMode(false);
       Stack.setFrame(1);
	   for (r=0; r < 4; r++) {
	     	 selectWindow(mainTittle);
	         makeRectangle(xCoords[r], yCoords[r], 700, 700);
	         run("Duplicate...", "title = [" + roiNames[r] + "] duplicate");
	         duplicateTittle = getTitle();
	     	 run("OrientationJ Vector Field", "tensor=5.0 gradient=0 radian=on vectorgrid=5 vectorscale=80.0 vectortype=0 vectoroverlay=on vectortable=on ");
	     	 csv_path = onedrive_path + day + "/" + baseName;
			 selectWindow("OJ-Table-Vector-Field-");
			 File.makeDirectory(csv_path);
			 saveAs("Results", csv_path + "/" + roiNames[r] + ".csv");
			 selectWindow(roiNames[r] + ".csv");
			 run("Close");
			 selectWindow(duplicateTittle);
			 run("Close");
	   }
	  run("Close All");
    }
}

showMessage("All files proccesed and csvs saved");
