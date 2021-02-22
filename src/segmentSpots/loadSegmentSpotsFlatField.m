function [ffim, doFF] = loadSegmentSpotsFlatField(PreProcPath, Prefix, channelIndex)
  %Load flat-field
  
  doFF = 1;
  ffim = [];
  
  try
    prefixPath = [PreProcPath, filesep, Prefix];
    D = dir([prefixPath, filesep, Prefix, '*FF*.tif']);
    if length(D) > 1
        nameSuffix = ['_ch', iIndex(channelIndex, 2)];
        dCh = dir([prefixPath,filesep, Prefix, '*FF*', nameSuffix, '.tif']);
        ffim = imread([prefixPath,Prefix, filesep,dCh{1},name]);
    elseif length(D) == 1
        ffim = imread([prefixPath, filesep, Prefix, '_FF.tif']);
    else 
        doFF = 0; 
    end
    ffim = CPsmooth(ffim, 'Gaussian Filter', 256, 0);

    %Normalize the image
    ffim = double(ffim) / double(max(max(ffim)));
  catch 
    warning('Will not apply flat field correction');
    doFF = 0;
  end 

end 
