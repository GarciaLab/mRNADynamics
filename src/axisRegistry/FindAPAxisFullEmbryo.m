function [coordA,coordP,xShift,yShift]=FindAPAxisFullEmbryo(Prefix, varargin)

% function [coordA,coordP,xShift,yShift]=FindAPAxisFullEmbryo(Prefix, varargin)
%
% DESCRIPTION
% Takes an image containing the whole embryo, and correlates the zoomed-in
% imaging field with the full embryo image to determine the shift and 
% then find the AP axis.
%
% PARAMETERS
% Prefix: Prefix of the dataset to analyze
%
% OPTIONS
% 'CorrectAxis': Opens to a GUI to allow manual correction of the xy
%                positions and orientations of the Anterior (A) and
%                Posterior (P) poles of the embryo following automatic
%                detection
% 'optionalResults': ???? Not sure exactly what this does ...
% 'FlipAP': Allows you to change the default Anterior/Posterior pole 
%           orientation. This script will look to the name of your MidSag 
%           image for an indication of the 'PA' orientation.
%
%
% OUTPUT
% coordA: xy coordinates of the anterior pole of the full embryo
% coordP: xy coordinates of the posterior pole of the full embryo
% xShift:
% yShift:
% 
%
% Author (contact): uknown (hggarcia@berkeley.edu)
% Created: XXXX-XX-XX
% Last Updated: 2020-07-27
% Documented by: Meghan Turner (meghan_turner@berkeley.edu)
%
%%

cleanupObj = onCleanup(@myCleanupFun);

% Process user options 
FlipAP = 0;
CorrectAxis = true;
optionalResults = '';

for i=1:length(varargin)
    if strcmpi(varargin{i},'FlipAP')
        FlipAP = 1;
    elseif strcmpi(varargin{i},'CorrectAxis')
        CorrectAxis = true;
    elseif strcmpi(varargin{i},'optionalResults')
        optionalResults = varargin{i+1};
    end
end
   
[SourcePath, ~, DefaultDropboxFolder, DropboxFolder, ~, ~, ~, ~] ... 
    = DetermineAllLocalFolders(Prefix, optionalResults);

%Find out the date it was taken
Dashes = strfind(Prefix,'-');
Date = Prefix(1:Dashes(3)-1);
EmbryoName = Prefix(Dashes(3)+1:end);

%Figure out what channels we have
[~, ~, ~, ~, ~, ~,Channel1, Channel2, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~,...
 Channel3,~,~, ~, ~]...
    = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);

% Figure out what type of data we're dealing with
rawDataPath = [SourcePath,filesep,Date,filesep,EmbryoName,filesep];
fullEmbryoPath = [rawDataPath,'FullEmbryo',filesep];
[~, FileMode] = DetermineFileMode(rawDataPath);

switch FileMode
    case 'TIF'
        dirFullEmbryo = dir([fullEmbryoPath,'*.tif']);
    case 'LIFExport' 
        dirFullEmbryo = dir([fullEmbryoPath,'*.lif']);
    case 'LSM'
        dirLSM = dir([fullEmbryoPath,'*.lsm']);
        dirCZI=dir([fullEmbryoPath,'*.czi']);
        dirFullEmbryo = [dirLSM, dirCZI];
    case 'SPIN' 
        dirFullEmbryo = dir([fullEmbryoPath,'*.nd']); 
    case 'ND2'
        %dirFullEmbryo = dir([fullEmbryoFolder,'*.nd2']); 
        error('Nikon point scanning .nd2 files not supported yet')
    case 'OMETIFF'
        error('OME-TIFF files not supported yet')
    case 'LAT'
        error('LAT (lattice light sheet) files not supported yet')
end

% Identify the midsagittal image
MidFileIndex=find(~cellfun('isempty',strfind(lower({dirFullEmbryo.name}),'mid')));
SurfFileIndex=find(~cellfun('isempty',strfind(lower({dirFullEmbryo.name}),'surf')));

if length(MidFileIndex)>1 && strcmpi(FileMode, ~'DSPIN')
    error('Too many midsagittal files in FullEmbryo folder')
end



% See if we don't want the default AP orientation
% MT 2020-07-27: Doesn't manual correction mode make this option
%                unnecessary? Maybe it can be removed?
if FlipAP == 1
    if strcmp(FileMode, 'TIF') || strcmp(FileMode, 'LIFExport') || strcmp(FileMode, 'LSM')
        if strcmp(dirFullEmbryo(MidFileIndex).name,'PA')
            FlipAP = 1;
        else
            FlipAP = 0;
        end
    end
end


% Grab the appropriate channel
ChannelToLoadTemp = contains([Channel1,Channel2,Channel3],'nuclear','IgnoreCase',true);

if sum(ChannelToLoadTemp) && sum(ChannelToLoadTemp)==1
    HisChannel = find(ChannelToLoadTemp);
elseif sum(ChannelToLoadTemp) && length(ChannelToLoadTemp)>=2
    ChannelToLoad = find(ChannelToLoadTemp);
    HisChannel = ChannelToLoad(1);
end

if isempty(HisChannel)
    error('LIF Mode error: Channel name not recognized. Check MovieDatabase.csv')
end

%% Rotate full embryo image and/or zoomed-in time series to match each other
% This is done differently for each type of microscopy data

% MT 2020-07-27: Functionalized each data type for more readable code
if strcmp(FileMode,'TIF')
    MidImage = imread([fullEmbryoPath,dirFullEmbryo(MidFileIndex).name],2);
    
elseif strcmp(FileMode,'LIFExport')
    midFile = [fullEmbryoPath,dirFullEmbryo(MidFileIndex).name];
    surfFile = [fullEmbryoPath,dirFullEmbryo(SurfFileIndex).name];
    [MidImage,SurfImage] = rotateFullEmbryoLIF(Prefix, rawDataPath, midFile, surfFile, HisChannel);
    
elseif strcmp(FileMode,'LSM')    
    midFile = [fullEmbryoPath,dirFullEmbryo(MidFileIndex).name];
    surfFile = [fullEmbryoPath,dirFullEmbryo(SurfFileIndex).name];
    [MidImage,SurfImage] = rotateFullEmbryoLSM(rawDataPath, midFile, surfFile, HisChannel);

elseif strcmp(FileMode, 'DSPIN')        
    [MidImage,SurfImage,xShift,yShift] = rotateFullEmbryoSPIN(dirFullEmbryo, fullEmbryoPath, Channel1, Channel2);

else
    warning('Image rotation correction not supported for this FileMode')
end

%%
%Save it to the Dropbox folder
mkdir([DropboxFolder,filesep,Prefix,filesep,'APDetection'])
imwrite(uint16(MidImage),[DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo.tif'],'compression','none');
imwrite(uint16(SurfImage),[DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryoSurf.tif'],'compression','none');

%Now, use them to find the embryo mask
liveExperiment = LiveExperiment(Prefix);
embMask = getEmbryoMaskLive(MidImage, liveExperiment.pixelSize_um);


%This code came from Michael's code
diagFigure = figure;
CC=bwconncomp(embMask);
if CC.NumObjects~=1
    warning('Failed to calculate embryo mask. Found more than one object in mask. Assigning arbitrary A and P positions.');
    coordA=[1,1];
    coordP=[1,1];
else
    
    % Rotate the mask to determine the AP axis as the extremal points of the mask
    Props=regionprops(CC,'Orientation');
    angle=Props.Orientation; % Angle is in DEGREES!


    I_mask_rot=imrotate(embMask,-angle);
    rotMatrix = [cosd(angle) sind(angle)
                -sind(angle) cosd(angle)];


    CC=bwconncomp(I_mask_rot);
    Props=regionprops(CC,'Centroid','MajorAxisLength', 'MinorAxisLength','Extrema');
    % After rotation, the major axis is aligned with x axis



    % for future diagnostic figures
    majorAxisBegin = Props.Centroid + [Props.MajorAxisLength/2,0];
    majorAxisEnd = Props.Centroid - [Props.MajorAxisLength/2,0];
    minorAxisBegin = Props.Centroid + [0, Props.MinorAxisLength/2];
    minorAxisEnd = Props.Centroid - [0, Props.MinorAxisLength/2];

    ext=Props.Extrema;
    coordP_rot=(ext(3,:)+ext(4,:))/2;
    coordA_rot=(ext(7,:)+ext(8,:))/2;

    % Calculate distances from centroid. In general Aneteriod will fall
    % further from centroid
    distA = sqrt(sum((coordA_rot-Props.Centroid).^2));
    distP = sqrt(sum((coordP_rot-Props.Centroid).^2));
    if distP > distA
        coordP_rot_temp = coordP_rot;
        coordP_rot = coordA_rot;
        coordA_rot = coordP_rot_temp;
    end
    
    if FlipAP
        temp = coordA_rot;
        coordA_rot = coordP_rot;
        coordP_rot = temp;
    end


    % coordA and coordP are the coordinates on the rotated image
    % We should rotate them back to the coordinates of the original picture
    % Remember that rotation was performed about the center of the image

    %coordinates of the center of the rotated image
    center_rot = 1/2*[size(I_mask_rot,2) size(I_mask_rot,1)];
    %coordinates of the center of the original image
    center = 1/2*[size(embMask,2) size(embMask,1)];

    coordA = center + (rotMatrix * (coordA_rot-center_rot)')';
    coordP = center + (rotMatrix * (coordP_rot-center_rot)')';
    
    % Save diagnostic figures to check the quality of axis determination
    imagesc(I_mask_rot);
    colormap(gray);
    axis image
    title('Anterior (green), posterior (red); rotated')
    hold on
    plot(coordA_rot(1),coordA_rot(2),'g.','MarkerSize',20);
    plot(coordP_rot(1),coordP_rot(2),'r.','MarkerSize',20);
    plot([majorAxisBegin(1),majorAxisEnd(1)],[majorAxisBegin(2),majorAxisEnd(2)],'b-');
    plot([minorAxisBegin(1),minorAxisEnd(1)],[minorAxisBegin(2),minorAxisEnd(2)],'b-');
    hold off
    saveas(gcf, [DropboxFolder,filesep,Prefix,filesep,'APMask.tif']);

end
    
%Save the AP and shift information
save([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'],'coordA','coordP');


clf
imagesc(MidImage)
axis image
axis off
title('Anterior (green), posterior (red); original')
hold on
plot(coordA(1),coordA(2),'g.','MarkerSize',20);
plot(coordP(1),coordP(2),'r.','MarkerSize',20);
hold off
saveas(gcf, [DropboxFolder,filesep,Prefix,filesep,'APEmbryo.tif']);
close(diagFigure);

% save embryo mask
save([DropboxFolder,filesep,Prefix,filesep,'EmbryoMask.mat'],'embMask');

if CC.NumObjects~=1 || CorrectAxis
    CorrectAPAxis(Prefix);
end
