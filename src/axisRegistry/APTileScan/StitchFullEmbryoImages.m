function StitchFullEmbryoImages(Prefix, varargin)
% author: Gabriella Martini
% date created: 12/30/19
% date last modified: 7/26/20

%% Parse Inputs 
if ~exist('Prefix')
    FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
    Dashes=strfind(FolderTemp,filesep);
    Prefix=FolderTemp((Dashes(end)+1):end);
end 

FullyAutomate = false;
StitchManually = false;
keepExistingStitching = false;
manualStitchOrder = false;
selectRegions = false;
useSurfStitching = false;
manualSeeding = false;
x = 1;
while x <= length(varargin)
    switch varargin{x}
        case{'keepExistingStitching'}
            keepExistingStitching = true;
        case {'FullyAutomate'}
            FullyAutomate = true;
            fprintf('Stitching fully automated.\n')
        case {'StitchManually'}
            StitchManually = true;
            fprintf('Stitching to be performed manually.\n')
        case {'MaxDeltaR'}
            MaxDeltaR = varargin{x+1};
            x = x+1;
            fprintf('Max change in row overlap to be used in stitching loop: %d\n', MaxDeltaR)
        case{'MaxDeltaC'}
            MaxDeltaC = varargin{x+1};
            x = x+1;
            fprintf('Max change in column overlap to be used in stitching loop: %d\n', MaxDeltaC)
        case{'manualStitchOrder'}
            manualStitchOrder = true;
        case{'selectStitchingRegions'}
            selectRegions=true;
        case{'useSurfStitchingInfo'}
            useSurfStitching=true;
        case{'manualSeeding'}
            manualSeeding=true;
        otherwise
            error(['Flag "', varargin{1}{x},'" not valid'])
    end
    x = x +1;
end

%% Stitch Surface and MidSaggital plane full embryo images


if ~keepExistingStitching
    varargin2 = {};
    if FullyAutomate
        varargin2{length(varargin2) + 1} = 'FullyAutomate';
    end
    if StitchManually 
        varargin2{length(varargin2) + 1} = 'StitchManually';
    end
    if manualStitchOrder
        varargin2{length(varargin2) + 1} = 'manualStitchOrder';
    end
    if selectRegions
        varargin2{length(varargin2) + 1} = 'selectStitchingRegions';
    end
    if exist('MaxDeltaR', 'var')
        varargin2{length(varargin2) + 1} = 'MaxDeltaR';
        varargin2{length(varargin2) + 1} = MaxDeltaR;
    end
    if exist('MaxDeltaC', 'var')
        varargin2{length(varargin2) + 1} = 'MaxDeltaC';
        varargin2{length(varargin2) + 1} = MaxDeltaC;
    end
    if useSurfStitching
        varargin2{length(varargin2) + 1} = 'useSurfStitchingInfo';
    end
    if manualSeeding
        varargin2{length(varargin2) + 1} = 'manualSeeding';
    end
    
    if length(varargin2) > 0
        EmbryoTileStitch(Prefix, 'Surf', varargin2);
        EmbryoTileStitch(Prefix, 'Mid', varargin2);
    else
        EmbryoTileStitch(Prefix, 'Surf');
        EmbryoTileStitch(Prefix, 'Mid');
    end
    
end
%% 

[SourcePath, ~, DefaultDropboxFolder, DropboxFolder, ~, ~,...
    ~, ~] = DetermineAllLocalFolders(Prefix);
stitchingDataFolder = [DropboxFolder,filesep,Prefix,filesep,'FullEmbryoStitching'];

if exist([stitchingDataFolder, filesep,'SurfTileArray.mat'], 'file')
     SurfData = load([stitchingDataFolder, filesep, 'SurfTileArray.mat']);
     SurfTileArray = SurfData.tile_array;
else
     error('No surface TileArray data stored. Seed a new TileArray using "NewTileArrayFromMetadata".')  
end



if exist([stitchingDataFolder, filesep,'MidTileArray.mat'], 'file')
     MidData = load([stitchingDataFolder, filesep, 'MidTileArray.mat']);
     MidTileArray = MidData.tile_array;
else
     error('No midsaggital plane TileArray data stored. Seed a new TileArray using "NewTileArrayFromMetadata".')  
end


%% 



%If the AP axis hasn't been specified check if it was specificed in the
%image names. Otherwise assume AP orientation

%Find out the date it was taken
Dashes=findstr(Prefix,'-');
Date=Prefix(1:Dashes(3)-1);
EmbryoName=Prefix(Dashes(3)+1:end);

%Figure out what type of experiment we have. Note: we name the var "DateFromDateColumn" to avoid shadowing previously defined "Date" var.
[DateFromDateColumn, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF, Channel3,~,~, ~, ~]...
    = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);

%Datatype is hardcoded in, unlike in FindAPAxisFullEmbryo
DLIF=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*.lif']);
FileMode='LIFExport';

DSTITCH=dir([stitchingDataFolder, filesep, '*.tif']);

%% 


% Identify the midsagittal image
MidLifFileIndex=find(~cellfun('isempty',strfind(lower({DLIF.name}),'midtile')));
SurfLifFileIndex=find(~cellfun('isempty',strfind(lower({DLIF.name}),'surftile')));
MidFileIndex=find(~cellfun('isempty',strfind(lower({DSTITCH.name}),'midtilestitch.tif')));
SurfFileIndex=find(~cellfun('isempty',strfind(lower({DSTITCH.name}),'surftilestitch.tif')));




%% 
%if strcmp(FileMode,'TIF')

       


%Rotates the full embryo image to match the rotation of the zoomed
%time series
zoom_angle = 0;
full_embryo_angle = 0;

LIFMid=bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,DLIF(MidLifFileIndex).name]);
LIFSurf=bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,DLIF(SurfLifFileIndex).name]);

MIDMeta = LIFMid{:, 4};
SURFMeta = LIFSurf{:,4};
%% 
[SurfNTiles, SurfNFrames, SurfNSlices, SurfNPlanes, SurfNChannels,...
    SurfFrame_Times] = getFrames(SURFMeta);

[MidNTiles, MidNFrames, MidNSlices, MidNPlanes, MidNChannels,...
    MidFrame_Times] = getFrames(MIDMeta);



PixelSize = double(MIDMeta.getPixelsPhysicalSizeX(1).value);% units: microns
PixelSize_m = double(PixelSize)*10^(-6);
surf_ypos = [];
surf_xpos = [];
for i=0:(SurfNTiles-1)
    surf_xpos(length(surf_xpos)+1) = -1*double(SURFMeta.getPlanePositionX(i,0).value);
    surf_ypos(length(surf_ypos)+1) = double(SURFMeta.getPlanePositionY(i,0).value);
end
surfxbound = find(surf_xpos == min(surf_xpos));
surfybound = find(surf_ypos == min(surf_ypos));
surfreftile = intersect(surfxbound, surfybound);

mid_ypos = [];
mid_xpos = [];
for i=0:(MidNTiles-1)
    mid_xpos(length(mid_xpos)+1) = -1*double(MIDMeta.getPlanePositionX(i,0).value);
    mid_ypos(length(mid_ypos)+1) = double(MIDMeta.getPlanePositionY(i,0).value);
end
midxbound = find(mid_xpos == min(mid_xpos));
midybound = find(mid_ypos == min(mid_ypos));
midreftile = intersect(midxbound, midybound);
DeltaR = round((surf_xpos(surfreftile)-mid_xpos(midreftile))/PixelSize_m, 0);
DeltaC = round((surf_ypos(surfreftile)-mid_ypos(midreftile))/PixelSize_m, 0);
%% 
% MaxProjFileName = [ID, 'TileStitch_Max.tif'];
% GaussFileName = [ID, 'TileStitch_MaxFiltered.tif'];
% ZStacksFileName = [ID, 'TileStitch.tif'];
MidImageNoPadding=imread([stitchingDataFolder,filesep,DSTITCH(MidFileIndex).name]);
SurfImageNoPadding=imread([stitchingDataFolder,filesep,DSTITCH(SurfFileIndex).name]);
MidImageMaxNoPadding = imread([stitchingDataFolder,filesep,'MidTileStitch_Max.tif']);
SurfImageMaxNoPadding = imread([stitchingDataFolder,filesep,'SurfTileStitch_Max.tif']);
MidImageFilteredNoPadding = imread([stitchingDataFolder,filesep,'MidTileStitch_MaxFiltered.tif']);
SurfImageFilteredNoPadding = imread([stitchingDataFolder,filesep,'SurfTileStitch_MaxFiltered.tif']);


if ~all(size(MidImageNoPadding) == size(SurfImageNoPadding))
    nrows = max(size(MidImageNoPadding, 1), size(SurfImageNoPadding, 1));
    ncols = max(size(MidImageNoPadding, 2), size(SurfImageNoPadding, 2));
    nz = size(MidImageNoPadding, 3);
    MidImage = zeros(nrows, ncols, nz, 'uint16');
    SurfImage = zeros(nrows, ncols, nz, 'uint16');
    MidMaxImage = zeros(nrows, ncols, 'uint16');
    SurfMaxImage = zeros(nrows, ncols, 'uint16');
    MidFilteredImage = zeros(nrows, ncols, 'uint16');
    SurfFilteredImage = zeros(nrows, ncols, 'uint16');
    [hm, wm] = size(MidImageNoPadding, 1:2);
    [hs, ws] = size(SurfImageNoPadding, 1:2);
    if DeltaR > 0
        if DeltaC > 0
            MidImage = MidImageNoPadding;
            SurfImage(DeltaR + 1:DeltaR + hs, DeltaC + 1:DeltaC + ws,:) = SurfImageNoPadding;
            MidMaxImage = MidImageMaxNoPadding;
            SurfMaxImage(DeltaR + 1:DeltaR + hs, DeltaC + 1:DeltaC + ws) = SurfImageMaxNoPadding;
            MidFilteredImage = MidImageFilteredNoPadding;
            SurfFilteredImage(DeltaR + 1:DeltaR + hs, DeltaC + 1:DeltaC + ws) = SurfImageFilteredNoPadding;
        else
            MidImage(1:hm, DeltaC + 1:DeltaC+wm,:) = MidImageNoPadding;
            SurfImage(DeltaR + 1:DeltaR + hs, 1:ws,:) = SurfImageNoPadding;
            MidMaxImage(1:hm, DeltaC + 1:DeltaC+wm) = MidImageMaxNoPadding;
            SurfMaxImage(DeltaR + 1:DeltaR + hs, 1:ws) = SurfImageMaxNoPadding;
            MidFilteredImage(1:hm, DeltaC + 1:DeltaC+wm) = MidImageFilteredNoPadding;
            SurfFilteredImage(DeltaR + 1:DeltaR + hs, 1:ws) = SurfImageFilteredNoPadding;
            
        end
    else
        if DeltaC > 0
            MidImage(DeltaR+1:DeltaR+hm, 1:wm,:) = MidImageNoPadding;
            SurfImage(1:hs, DeltaC+1:DeltaC+ws,:) = SurfImageNoPadding;
            MidMaxImage(DeltaR+1:DeltaR+hm, 1:wm) = MidImageMaxNoPadding;
            SurfMaxImage(1:hs, DeltaC+1:DeltaC+ws) = SurfImageMaxNoPadding;
            MidFilteredImage(DeltaR+1:DeltaR+hm, 1: wm) = MidImageFilteredNoPadding;
            SurfFilteredImage(1:hs, DeltaC+1:DeltaC+ws) = SurfImageFilteredNoPadding;
        else
            MidImage(DeltaR + 1:DeltaR + hm, DeltaC + 1:DeltaC + wm,:) = MidImageNoPadding;
            SurfImage(1:hs, 1:ws,:) = SurfImageNoPadding;
            MidMaxImage(DeltaR + 1:DeltaR + hm, DeltaC + 1:DeltaC + wm) = MidImageMaxNoPadding;
            SurfMaxImage(1:hs, 1:ws) = SurfImageMaxNoPadding;
            MidFilteredImage(DeltaR + 1:DeltaR + hm, DeltaC + 1:DeltaC + wm) = MidImageFilteredNoPadding;
            SurfFilteredImage(1:hs, 1:ws) = SurfImageFilteredNoPadding;
        end
    end
    
  
    
    %=MidImage = imresize(MidImage,length(SurfImage)/length(MidImage));
else
    MidImage = MidImageNoPadding;
    SurfImage = SurfImageNoPadding;
    MidMaxImage = MidImageMaxNoPadding;
    SurfMaxImage = SurfImageMaxNoPadding;
    MidFilteredImage = MidImageFilteredNoPadding;
    SurfFilteredImage = SurfImageFilteredNoPadding;
end

%% Save stitched data in tif files

outputFolder = [DropboxFolder,filesep,Prefix,filesep,'FullEmbryoStitching'];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
APDetectionFolder = [DropboxFolder,filesep,Prefix,filesep,'APDetection'];
if ~exist(APDetectionFolder, 'dir')
    mkdir([DropboxFolder,filesep,Prefix,filesep,'APDetection']);
end

MidMaxProjFileName = 'MidTileStitch_MaxPadded.tif';
MidGaussFileName = 'MidTileStitch_MaxFilteredPadded.tif';
MidZStacksFileName = 'MidTileStitchPadded.tif';
SurfMaxProjFileName = 'SurfTileStitch_MaxPadded.tif';
SurfGaussFileName = 'SurfTileStitch_MaxFilteredPadded.tif';
SurfZStacksFileName = 'SurfTileStitchPadded.tif';
     
imwrite(MidMaxImage,[outputFolder, filesep,MidMaxProjFileName],'compression','none');
imwrite(MidFilteredImage,[outputFolder, filesep,MidGaussFileName],'compression','none');
imwrite(MidImage,[outputFolder, filesep,MidZStacksFileName],'compression','none');
imwrite(SurfMaxImage,[outputFolder, filesep,SurfMaxProjFileName],'compression','none'); 
imwrite(SurfFilteredImage,[outputFolder, filesep,SurfGaussFileName],'compression','none');
imwrite(SurfImage,[outputFolder, filesep,SurfZStacksFileName],'compression','none');

end