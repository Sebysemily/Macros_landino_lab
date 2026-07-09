//This macro creates difference movies from multi-channel movies

//This macro creates difference movies from multi-channel movies

// 1. Input dir to iterate in the folder
basedir = "/home/sebas/Documents/fiji/Onedrive/macros/nd2_trial"; //change this to the path of the folder you want to iterate through
outputDir  = "/home/sebas/Documents/fiji/output"; //change this to the path of the folder you want to save the resliced images to


differenceNumber = getNumber("how many frames do you want to subtract?", 10) ; 
//differenceNumber = 10;
//asks the user for the number of frames to subtract
//Create inputdir variable to use in the for loop below
inputdir = basedir + "/";
// 2.For loop to iterate trough a folder
fileList = getFileList(inputdir); //creates a list of all the files in the specified folder
for ( i =0; i < fileList.length; i++) {
  if (endsWith(fileList[i], ".nd2")){
  	filePath = inputdir + fileList[i];
      // Open the image
      run ("Bio-Formats Importer", "open=[" + filePath + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
      fileName = getInfo("image.title"); 
      //gets and saves the file name for later
      dotIndex = indexOf(fileName, "."); 
      //gets dotIndex
      getDimensions(width, height, channels, slices, frames) ;		
      //gets and saves the movie dimensions for later use

      if (dotIndex > 0) {
        //this and the following line get the file name without the extension
        fileNameWithoutExtension = substring(fileName, 0, dotIndex); 
        //this and the above line get the file name without the extension
        newFileName = fileNameWithoutExtension + "_Diff" + differenceNumber + ".tif" ;
      } else {
        newFileName = fileName + "_Diff" + differenceNumber + ".tif" ;
      }

      counter = 1 ;//creates a counter variable that starts as 1 and increases by 1 with every trip through the loop
      while (counter <= channels) {  //runs a loop as long as the there are still channels left to duplicate
        Stack.setChannel(counter); //moves to channel x (whatever the x number through the loop is)
        run("Duplicate...", "title=Channel_" + counter + " duplicate channels=" + counter);
        //duplicates the active channel and renames it "Channel_X"
        setSlice(1);
        run("Duplicate...", "title=firstDup duplicate");
        run("Duplicate...", "title=lastDup duplicate");
        selectWindow("firstDup") ;
        for (i=1; i<=differenceNumber; i++) { //runs a loop
          run("Delete Slice");    //changes the slice for each position in the loop
        }
        selectWindow("lastDup") ;
        run("Reverse");
        for (i=1; i<=differenceNumber; i++) { //runs a loop
          run("Delete Slice");    //changes the slice for each position in the loop
        }
        run("Reverse");
        imageCalculator("Subtract create 32-bit stack", "firstDup","lastDup");
        selectWindow("firstDup") ;
        close() ; //closes first intermediate
        selectWindow("lastDup") ;
        close() ; //closes last intermediate
        selectWindow("Channel_" + counter);
        close(); //closes the duplicated channel selection
        selectWindow("Result of firstDup") ;
        rename("C" + counter) ; //renames difference movie to the channel it was generated from
        getMinAndMax(min, max) ;
        setMinAndMax(0, max) ; //thresholds the video 
        run("16-bit");
        run("Enhance Contrast", "saturated=0.35");
        selectWindow(fileName);
        counter += 1; //loops through again
      }

      if (channels == 2) {
        run("Merge Channels...", "c1=C1 c2=C2 create");
      } //merges two channels together

      if (channels == 3) {
        run("Merge Channels...", "c1=C1 c2=C2 c3=C3 create");
      } //merges 3 channels together

      if (channels == 4) {
        run("Merge Channels...", "c1=C1 c2=C2 c3=C3 c4=C4 create");
      } //merges 4 channels together

      if (channels == 1) {
        selectWindow("C1");
      } 

      rename(newFileName) ;
  saveAs("Tiff", outputDir + "/" + newFileName);
  run("Close All");
  }
showMessage("All files proccesed and saved to " + outputDir);

