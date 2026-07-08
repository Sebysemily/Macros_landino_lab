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
	     run("Canvas Size...", "width=1400 height=1400 position=Center zero");
	     Stack.getDimensions(width, height, channels, slices, frames);
	     quarterframe = floor(frames/4) + 1;
	     Stack.setFrame(quarterframe);
	     run("Enhance Contrast...", "saturated=0.40 normalize process_all");
	     mainImageID = getImageID();
	     mainTittle = getTitle();
	     imageName = File.getName(mainTittle);
	     baseName = replace(imageName, ".nd2", "");
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