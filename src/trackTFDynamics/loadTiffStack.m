function tiff3DArray = loadTiffStack(tiffStackPath)
%
% function tiff3DArray = loadTiffStack(tiffStackPath)
%
% DESCRIPTION
% This function reads in a multi-page (i.e. 3D) TIFF stack (which is the 
% current file format for the mRNADynamics pipeline) and outputs it 
% as a 3D uint16 array, which is an easier format to use for image
% processing and analysis. Note, this is 10x faster than the BioFormats
% bfOpen3DVolume function.
%
% INPUT ARGUMENTS
% tiffStackPath: full path to a single 3D TIFF stack
% 
% OPTIONS
% N/A
%
% OUTPUT
% tiff3DArray: 3D (xyz) array containing the image from the TIFF stack
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 08/12/2020
% Last Updated: N/A
%

if ~exist(tiffStackPath,'file')
    error(['Cannot find file: ' tiffStackPath])
elseif ~(contains(tiffStackPath,'tif','IgnoreCase',true) || ...
            contains(tiffStackPath,'tiff','IgnoreCase',true))
    error('Input argument tiffStackPath must point to a TIF or TIFF file')
else
    imInfo = imfinfo(tiffStackPath);
end

xDim = imInfo(1).Width;
yDim = imInfo(1).Height;
zDim = numel(imInfo);


im3DArray = NaN(xDim,yDim,zDim);

% Read in each z slice individually
for i = 1:zDim
    currZSlice = imread(tiffStackPath,i);
    im3DArray(:,:,i) = currZSlice;
end

% Cast back to uint16 (changed to double in the process of adding to array)
tiff3DArray = uint16(im3DArray);



