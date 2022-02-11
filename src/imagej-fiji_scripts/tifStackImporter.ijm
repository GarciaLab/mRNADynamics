macro "tiffStackImporter" {
	inputDir = getDirectory("Select the folder containing your TIFF stacks: ");
	fileList = getFileList(inputDir);

	//setBatchMode(true);	//hides images as they're opened
	
	IJ.log("Found " + toString(fileList.length) + " TIF stacks");
	
	for (i=0; i<fileList.length; i++) {
	    showProgress(i, fileList.length);
	    if (filter(i, fileList[i])) {
	        tiffStackPath = ""+inputDir+fileList[i];
	        run("TIFF Virtual Stack...", "open="+tiffStackPath);
	    }
	}
	
	//setBatchMode(false);	//Show all further images that are opened
	
	// Concatenate all open TIF stacks - produces a Hyperstack
	run("Concatenate...", "all_open open");
}



function filter(i, name) {
    // ignore directories
    if (endsWith(name,"/")) return false;

    // ignore files that aren't TIFFs
    if (!endsWith(name,".tif")) return false;

    // ignore text files
    if (endsWith(name,".txt")) return false;

    // ignore His stack
    if (indexOf(name, "His") >= 0) return false;

    // does name contain both "Series002" and "ch01"
    // if (indexOf(name,"Series002")==-1) return false;
    // if (indexOf(name,"ch01")==-1) return false;

    // open only first 10 images
    // if (i>=10) return false;

    return true;
}