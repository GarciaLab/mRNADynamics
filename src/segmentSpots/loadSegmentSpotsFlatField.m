function [ffim, doFF] = loadSegmentSpotsFlatField(PreProcPath, Prefix, FrameInfo)
  %Load flat-field. We need to process this file differently the images come
  %from a laser scanning or spinning disk microscope.
  doFF = 1;
  ffim = [];
  
  try 
    ffim = imread([PreProcPath, filesep, Prefix, filesep, Prefix, '_FF.tif']);
    %If we have a spinning disk confocal

    if strcmpi(FrameInfo(1).FileMode, 'dspin')
      disp('Assuming a spinning disk confocal for flat-field correction')
      %Note that we brought this back to the same parameters as for a
      %LSC. We need to figure out what's going on with the flatfields on
      %the spinning disk. If not, we can always crop the image.
      ffim = CPsmooth(ffim, 'Gaussian Filter', 256, 0);
      %If not, we assume we have a laser scanning confocal
    else 
      ffim = CPsmooth(ffim, 'Gaussian Filter', 256, 0);
    end 

    %Normalize the image
    ffim = double(ffim) / double(max(max(ffim)));
  catch 
    warning('Will not apply flat field correction');
    doFF = 0;
  end 

end 
