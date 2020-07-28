function [MidImage,SurfImage] ...
    = rotateFullEmbryoLSM(rawDataPath, midFile, surfFile, HisChannel)
%% 
% function [MidImage,SurfImage] ...
%     = rotateFullEmbryoLSM(rawDataPath, midFile, surfFile, HisChannel)
%
% DESCRIPTION
% Rotates the full embryo image to match the zoomed-in times series images
% for Zeiss confocal data ('LSM' file mode, either LSM or CZI files)
%
% PARAMETERS
% rawDataPath: path for the folder containing the raw zoomed-in time series
%              images
% midFile: path for the midsagittal plane image
% surfFile: path for the surface image
% HisChannel: channel containing the nuclear marker that will be used to
%             align the full embyro and zoomed-in images
%
% OPTIONS
% N/A
%
% OUTPUT
% MidImage: rotated midsagittal plane image
% SurfImage: rotated surface image
%
%
% Author (contact): unknown (hggarcia@berkeley.edu)
% Created: unknown
% Last Updated: 2020-07-27
%
% Functionalized from code originally in FindAPAxisFullEmbryo, by Meghan
% Turner (meghan_turner@berkeley) on 2020-07-27
%

LSMMid = bfopen(midFile);
LSMSurf = bfopen(surfFile);
LSMMeta = LSMMid{:,4};
LSMMeta2 = LSMMid{:,2};

%Look for the image with the largest size. In this way, we avoid
%loading individual tiles in the case of a tile scan.
for i = 1:size(LSMMid,1)
    SizesImages(i) = size(LSMMid{i,1}{1,1},1);
end
[~,ImageCellToUse] = max(SizesImages);

%individual tiles if we're dealing with tile scan. Also, in CZI files,
%this seems to ensure a high-contrast image as well.
MidImage = LSMMid{ImageCellToUse,1}{HisChannel,:};
SurfImage = LSMSurf{ImageCellToUse,1}{HisChannel,:};


%Rotates the full embryo image to match the rotation of the zoomed
%time series
zoom_angle = NaN;    %Set defaults to NaN so we can tell if the 
mid_angle = NaN;     %angles were successfully extracted

%Figure out the rotation of the full embryo image
%This works for LSM files
mid_angle = LSMMeta2.get('Recording Rotation #1');
%If full_embryo_angle is empty, chances are we have a CZI file
if isempty(mid_angle)
    mid_angle = str2num(LSMMeta2.get('Global HardwareSetting|ParameterCollection|RoiRotation #1'));
%Check if angle was successfully extracted
elseif isnan(mid_angle)
    error('Could not extract rotation of FullEmbryo images')
end

%Figure out the rotation of the zoomed-in image. We need to check for
%both LSM and CZI files.
DLSMZoom = dir([rawDataPath,'*.lsm']);
DLSMZoom = [DLSMZoom,dir([rawDataPath,'*.czi'])];
LSMZoom = bfopen([rawDataPath,DLSMZoom(1).name]);
LSMMetaZoom2 = LSMZoom{:,2};

%This works for LSM files
zoom_angle = LSMMetaZoom2.get('Recording Rotation #1');
%If full_embryo_angle is empty, chances are we have a CZI file
if isempty(zoom_angle)
    zoom_angle = str2num(LSMMetaZoom2.get('Global HardwareSetting|ParameterCollection|RoiRotation #1'));
%Check if angle was successfully extracted
elseif isnan(zoom_angle)
    error('Could not extract rotation of FullEmbryo images')
end

MidImage = imrotate(MidImage, -zoom_angle + mid_angle);
SurfImage = imrotate(SurfImage, -zoom_angle + mid_angle);