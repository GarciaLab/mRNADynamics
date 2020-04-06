function GenerateStitchedData(Prefix, ID)
% author: Gabriella Martini
% date started: 12/29/19
% last modified: 12/29/19
%% Parse User Input 
if ~exist('Prefix')
     FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
     Dashes=strfind(FolderTemp,filesep);
     Prefix=FolderTemp((Dashes(end)+1):end);
end
if ~exist('ID', 'var')
    prompt = ['Enter an "ID" for stitching.',...
    'The standard inputs "Mid" and "Surf" will stitch the',...
    'Midsagittal and Surface Full Embryo images using the',...
    '"MidTile.lif" and "SurfTile.lif" files respectively.'];
    ID = input(prompt,'s');
    ID = [upper(ID(1)), ID(2:end)];
end
%% Determine local folders and extract metadata information 

[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
    DetermineLocalFolders(Prefix);
%Find out the date it was taken
Dashes=findstr(Prefix,'-');
Date=Prefix(1:Dashes(3)-1);
EmbryoName=Prefix(Dashes(3)+1:end);
if ~isempty(strfind(lower(ID), 'mid'))
    rawdatfile = 'MidTile';
elseif ~isempty(strfind(lower(ID), 'surf'))
    rawdatfile = 'SurfTile';
else
    rawdatfile = ID;
end
LIFPath = [SourcePath, filesep, Date, filesep, EmbryoName,filesep,...
'FullEmbryo\', rawdatfile,'.lif'];

r = bfGetReader();
% Decorate the reader with the Memoizer wrapper
r = loci.formats.Memoizer(r);
r.setId(LIFPath);
LIFImages = bfopen(LIFPath);
LIFMeta = LIFImages{:, 4};
r.close();

[Prefix, SkipFrames, ProjectionType, PreferredFileNameForTest, keepTifs,...
generateTifStacks, nuclearGUI, skipExtraction, rootFolder, zslicesPadding,...
lowbit] = exportDataForLivemRNA_processInputParameters(Prefix);
[rawDataPath, ~, DropboxFolder, ~, PreProcPath, rawDataFolder, Prefix, ExperimentType, Channel1, Channel2, ~,...
    Channel3] = readMovieDatabase(Prefix,'rootFolder', rootFolder);
[NSeries, NFrames, NSlices, NPlanes, NChannels, Frame_Times] = getFrames(LIFMeta);
if sum(NFrames)~=0
    [Frame_Times, First_Time] = obtainFrameTimes(XMLFolder, seriesPropertiesXML, NSeries, NFrames, NSlices, NChannels);
    [InitialStackTime, zPosition] = getFirstSliceTimestamp(NSlices, NSeries, NPlanes, NChannels, Frame_Times, XMLFolder, seriesXML);
else
    InitialStackTime = [];
    zPosition = [];
end
framesIndex  = 1;
NTiles = NSeries;
FrameInfo = recordFrameInfo(NFrames, NSlices, InitialStackTime, LIFMeta, zPosition);

[coatChannel, histoneChannel, fiducialChannel, inputProteinChannel, FrameInfo] =...
LIFExportMode_interpretChannels(ExperimentType, Channel1, Channel2, Channel3, FrameInfo);

%% Load stitching information saved in "tile_array" variable
if ~isempty(strfind(lower(ID), 'mid'))
    stitchingDatafile = 'MidTileArray.mat';
elseif ~isempty(strfind(lower(ID), 'surf'))
    stitchingDatafile = 'SurfTileArray.mat';
else
    stitchingDatafile = [ID, 'TileArray.mat'];
end
load([DropboxFolder, filesep, Prefix, filesep,'FullEmbryoStitching',filesep, stitchingDatafile])

%% Extract tiled data in z stacks, max projections, and gaussian filtered images

ZStacks = {};
for n=1:NTiles
    ImageSlices = generateHisSlicesTileScan(LIFImages, NSlices(n),...
        NChannels, fiducialChannel, framesIndex, n);
    ZStacks{n} = ImageSlices;
end

%MaxProjImgs = tile_array.tiles;
%GaussFiltImgs = tile_array.imgs;
%% Generate stitched images
% Stitched Max-projected image
StitchedMaxImm = imstitchTile(tile_array, 'NoFilter', true);
PixelSize = double(LIFMeta.getPixelsPhysicalSizeX(1).value);% units: microns
sigma = .6/PixelSize;
StitchedGaussImm = imgaussfilt(StitchedMaxImm, sigma);
StitchedZStacks = imstitchTile(tile_array, 'ZStack', true);
% Stitched raw data z stacks 
hs = [];
ws = [];
for n=1:NTiles
    hs(n) = size(tile_array.imgs{n}, 1);
    ws(n) = size(tile_array.imgs{n}, 2);
end

rs = [tile_array.rows{:}];
cs = [tile_array.cols{:}];
rbounds = rs + hs-1;
cbounds = cs + ws-1;
numrows = max(rbounds);
numcols = max(cbounds);



%% Save stitched data in tif files

outputFolder = [DropboxFolder,filesep,Prefix,filesep,'FullEmbryoStitching'];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
APDetectionFolder = [DropboxFolder,filesep,Prefix,filesep,'APDetection'];
if ~exist(APDetectionFolder, 'dir')
    mkdir([DropboxFolder,filesep,Prefix,filesep,'APDetection']);
end

MaxProjFileName = [ID, 'TileStitch_Max.tif'];
GaussFileName = [ID, 'TileStitch_MaxFiltered.tif'];
ZStacksFileName = [ID, 'TileStitch.tif'];
     
imwrite(StitchedMaxImm,[outputFolder, filesep,MaxProjFileName],'compression','none');
imwrite(StitchedGaussImm,[outputFolder, filesep,GaussFileName],'compression','none');
imwrite(StitchedZStacks,[outputFolder, filesep,ZStacksFileName],'compression','none');
imwrite(StitchedZStacks,[APDetectionFolder, filesep,ZStacksFileName],'compression','none'); 
end