function [MidImage,SurfImage] ...
    = rotateFullEmbryoLIF(Prefix, rawDataPath, midFile, surfFile, HisChannel)
%% 
% function [MidImage,SurfImage] ...
%    = rotateFullEmbryoLIF(dirFullEmbryo, fullEmbryoPath, MidFileIndex, SurfFileIndex, HisChannel)
%
% DESCRIPTION
% Rotates the full embryo image to match the zoomed-in times series images
% for Leica confocal data ('LIFExport' file mode)
% (Might also do some other necessary full embryo image processing.)
%
% PARAMETERS
% fullEmbryoPath: path for the folder containing the full embyro images
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

fullEmbryoPath = [rawDataPath,'FullEmbryo',filesep];
LIFMid = bfopen(midFile);
LIFSurf = bfopen(surfFile);

%By looking at the last image we make sure we're avoiding the
%individual tiles if we're dealing with tile scan
MidImage = LIFMid{end,1}{HisChannel,1};
SurfImage = LIFSurf{end,1}{HisChannel,1};
if size(MidImage) ~= size(SurfImage)
    MidImage = imresize(MidImage,length(SurfImage)/length(MidImage));
end

zoom_angle = 0;     % Set the default rotation to 0
mid_angle = 0;      % Set the default rotation to 0

zoom_angle = getZoomAngle(Prefix, rawDataPath);
mid_angle = getFullEmbryoAngle(fullEmbryoPath, midFile, Prefix);

MidImage = imrotate(MidImage, -zoom_angle + mid_angle);