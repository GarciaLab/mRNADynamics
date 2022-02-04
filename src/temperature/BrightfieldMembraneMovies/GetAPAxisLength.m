function [APLength, DVLength]=GetAPAxisLength(Prefix)
%%

[SourcePath, ~, DefaultDropboxFolder, DropboxFolder, ~, ~, ~, ~] ...
    = DetermineAllLocalFolders(Prefix, '');
liveExperiment = LiveExperiment(Prefix);
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
MidFileIndex=find(~cellfun('isempty',strfind(lower({dirFullEmbryo.name}),'mid.lif')));

if length(MidFileIndex)>1
    error('Too many midsagittal files in FullEmbryo folder')
end





% Grab the appropriate channel
ChannelToLoadTemp = contains([Channel1,Channel2,Channel3],'brightfield','IgnoreCase',true);

if sum(ChannelToLoadTemp) && sum(ChannelToLoadTemp)==1
    MemChannel = find(ChannelToLoadTemp);
elseif sum(ChannelToLoadTemp) && length(ChannelToLoadTemp)>=2
    ChannelToLoad = find(ChannelToLoadTemp);
    MemChannel = ChannelToLoad(1);
end

if isempty(MemChannel)
    error('LIF Mode error: Channel name not recognized. Check MovieDatabase.XLSX')
end

%% Rotate full embryo image and/or zoomed-in time series to match each other
% This is done differently for each type of microscopy data

% MT 2020-07-27: Functionalized each data type for more readable code
if strcmp(FileMode,'LIFExport')
    midFile = [fullEmbryoPath,dirFullEmbryo(MidFileIndex).name];
    
    LIFMid = bfopen(midFile);
    LIFMeta = LIFMid{:, 4};
    NSlicesAllChannels = size(LIFMid{end,1},1);
    IncludedSlices = [];
    for i = 1:NSlicesAllChannels
        fileString = LIFMid{end,1}{i,2};
        if contains(fileString,['C=',num2str(MemChannel)]) | ~contains(fileString,'C=')
            IncludedSlices=[IncludedSlices,i];
        end
    end
    MidImageSize = [size(LIFMid{end,1}{1,1}), length(IncludedSlices)];
    MidImageStack = zeros(MidImageSize, 'uint8');
    MedianImageStack = zeros(MidImageSize, 'uint8');
    RescaledMidImageStack = zeros(MidImageSize, 'uint8');
    for i = 1:length(IncludedSlices)
        TempImage = (LIFMid{end,1}{IncludedSlices(i),1}+1)/256-1;
        MidImageStack(:,:,i) = TempImage;
        MaxTempImage = max(max(TempImage));
        MinTempImage = min(min(TempImage));
        RescaledMidImageStack(:,:,i) = (TempImage-MinTempImage)*255/(MaxTempImage-MinTempImage);
    end
    MidImage = max(RescaledMidImageStack,[],3);
    MidImage2 = sum(uint16(MidImageStack),3);
    MidImage2 =  uint8((MidImage2+1)/256-1);
    MidImage3 = median(MidImageStack, 3);
    MaxTempImage = max(max(MidImage2));
    MinTempImage = min(min(MidImage2));
    MidImage2 = (MidImage2-MinTempImage)*255/(MaxTempImage-MinTempImage);

else
    warning('Image rotation correction not supported for this FileMode')
end

%%
%Save it to the Dropbox folder
mkdir([DropboxFolder,filesep,Prefix,filesep,'APDetection'])
imwrite(uint16(MidImage3),[DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo.tif'],'compression','none');



%This code came from Michael's code

coordA=[1,1];
coordP=[1,1];
coordD=[1,1];
coordV=[1,1];

%Save the AP and shift information

save([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'],'coordA','coordP', 'coordD', 'coordV');


CorrectAPAxis(Prefix);
close all
%%
load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])
if ~isempty(str2double(LIFMeta.getPixelsPhysicalSizeX(0))) &...
        ~isnan(str2double(LIFMeta.getPixelsPhysicalSizeX(0)))
    PixelSize = str2double(LIFMeta.getPixelsPhysicalSizeX(0));
else
    try
        PixelSize = str2double(LIFMeta.getPixelsPhysicalSizeX(0).value);
    catch %no idea man
        PixelSize = str2double(LIFMeta.getPixelsPhysicalSizeX(1).value);
    end
end

APLength= sqrt((coordA(1)-coordP(1))^2+(coordA(2)-coordP(2))^2)*PixelSize;
APSlope = (coordA(2)-coordP(2))/(coordA(1)-coordP(1));
APIntercept = coordA(2)-APSlope*coordA(1);
APtheta = atan(APSlope);

DVSlope = (coordD(2)-coordV(2))/(coordD(1)-coordV(1));
DVIntercept = coordD(2)-DVSlope*coordD(1);
DVtheta = atan(DVSlope);

DeltaTheta = APtheta-DVtheta;

DVLength = sqrt((coordD(1)-coordV(1))^2+(coordD(2)-coordV(2))^2)*PixelSize*abs(sin(DeltaTheta));


save([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'],'coordA','coordP','coordD','coordV', 'APLength','DVLength');