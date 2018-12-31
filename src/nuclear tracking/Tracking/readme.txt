Here's a list of the functions to call to do the full tracking :


* START segmentation and tracking :

	[nuclei, centers, ~, dataStructure] = mainTracking(names,'indMitosis',indMit,'embryoMask', embryo_mask);

names is a cell array containing the names of all frames in the movie in order.
indMitosis is an nx2 array containing the first and last frame of mitosis in every row.
embryoMask is the possible mask of the embryo. If no embryo edge is visible, true(size(an_image_from_the_movie)) can be given as input. To specify the space resolution or the time resolution of the movie, add the corresponding name/value pair to the inputs, i.e. mainTracking(...,'space resolution', xxx, 'time resolution' yyy, ...)


* PUT circles on the nuclei

	[ellipses] = putCirclesOnNuclei(centers,names,indMit);

* Convert nuclei structure into schnitzcell structure

	[Schnitzcells] = convertNucleiToSchnitzcells(nuclei);


* Perform corrections


* Convert back to nuclei and centers

	[centers,mapping,approvedCenters,nuclei] = convertSchnitzcellsToCentersAndMapping(schnitzcells,numberOfFrames,approvedSchnitzcells);


*** OPTIONAL : if corrections were made on the segmentation that were written on the Ellipses structure instead of the Schnitzcells structure, the 'centers' variable should be updated the following way:
centers = updateCentersFromElliipses(Ellipses, old_centers);


* Rerun tracking

	[nuclei2] = mainTracking(names,'data structure',dataStructure,'centers',centers,'mapping',mapp
ing,'approved',approvedCenters);


(If segmentation was corrected, re-run the function putCirclesOnNuclei to get a new 'ellipses' structure.)

* Convert nuclei structure into schnitzcell with the approved structure and get the mapping between the previous schnitzcell and the new one

	[schnitzcells, approvedSchnitz, map ] = convertNucleiToSchnitzcells( nuclei, approvedCenters, previously_approvedSchnitz, previous_schnitzcells )

etc...



------------------------------------
New stuff:

So I added a function called 'updateCentersFromElliipses' that does what I told you:

1. If you call it this way:

centers = updateCentersFromElliipses(Ellipses);

it just copies the position information from 'Ellipses' (the two first columns).

2. If you call it this way:

centers = updateCentersFromElliipses(Ellipses, old_centers);

it maps old_centers to Ellipses (in case Ellipse fitting changes the center) by matching the Ellipses centers to the closest old_centers element. If the distance between an Ellipses element and the closest element in old_centers is above the constant THRESHOLD_DISTANCE, then this is assumed to be a nucleus that was added and the center from the corresponding Ellipses element is copied.


-------------------------------------------

Now, every time you call 'mainTracking', you can get up to four outputs. The three first are the same as before, that is 'nuclei', 'xy' and 'approvedNuclei'. The fourth one is 'dataStructure'. I.e :

[nuclei, xy, approvedNuclei, dataStructure] = mainTracking(names,....)


This dataStructure stores the inputs that were given (e.g. indMitosis, space resolution, time resolution, embryoMask) and some temporary data for the tracking, such as the correlation data for each frame pair.

The result is a simplification and a dramatic speed up of subsequent calls to mainTracking, during the manual correction phase. You can basically replace almost all inputs by the data structure in the following way:

[...] = mainTracking(names,'data structure',dataStructure, 'centers', xy, 'mapping', manualMapping,...)
