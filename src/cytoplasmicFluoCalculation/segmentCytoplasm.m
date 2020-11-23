function segmentCytoplasm(Prefix, varargin)
cleanupObj = onCleanup(@myCleanupFun);
%this function uses persistent (static) variables to speed computation.
%if not cleared, this could lead to errors 
clear fitSingleGaussian
warning('off', 'MATLAB:MKDIR:DirectoryExists');

optionalResults = [];

liveExperiment = LiveExperiment(Prefix);
yDim = liveExperiment.yDim;
xDim = liveExperiment.xDim;
zDim = liveExperiment.zDim;


inputChannels = liveExperiment.inputChannels;

[~, ~, DropboxFolder, ~, ~] = DetermineLocalFolders(Prefix, optionalResults);




FrameInfo = getFrameInfo(liveExperiment);
NFrames = length(FrameInfo);

[coatChannel, histoneChannel, fiducialChannel, inputProteinChannel, FrameInfo] =...
    LIFExportMode_interpretChannels(liveExperiment.experimentType, {liveExperiment.Channel1},...
    {liveExperiment.Channel2}, {liveExperiment.Channel3}, FrameInfo);

%% 


%Get the information about the AP axis as well as the image shifts
%used for the stitching of the two halves of the embryo
load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])


if ~exist('coordPZoom', 'var')
    warning('AddParticlePosition should have been run first. Running it now.')
    AddParticlePosition(Prefix, 'ManualAlignment')
    load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])
end



%Angle between the x-axis and the AP-axis
APAngle=atan2((coordPZoom(2)-coordAZoom(2)),(coordPZoom(1)-coordAZoom(1)));
DVAngle = APAngle + pi/2;
APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);
APResolution = liveExperiment.APResolution;

APPosImage=zeros(yDim,xDim);

for i=1:yDim
    for j=1:xDim
        Angle = atan2((i-coordAZoom(2)),(j-coordAZoom(1)));
        APDistance=sqrt((coordAZoom(2)-i).^2+(coordAZoom(1)-j).^2);
        APPosition=APDistance.*cos(Angle-APAngle);
        APPosImage(i,j)=APPosition/APLength;
    end
end

APbins = 0:APResolution:1;
APbinIDs=1:length(APbins);
APPosBinImage=zeros(size(APPosImage));
for i=APbinIDs(1:end-1)
    if i < APbinIDs(end-1)
        FilteredMask=(APbins(i)<=APPosImage)&(APbins(i+1)>APPosImage);
        APPosBinImage=APPosBinImage+FilteredMask*i;
    else
        FilteredMask=(APbins(i)<=APPosImage);
        APPosBinImage=APPosBinImage+FilteredMask*i;
    end
        
end

APPosBinImage3D = repmat(APPosBinImage, 1, 1, zDim);

MinAPbin = min(min(APPosBinImage));
MaxAPbin = max(max(APPosBinImage));

%% 


load([DropboxFolder, filesep, Prefix, filesep, 'cytoplasmFrames.mat']);



DogOutputFolder=[liveExperiment.procFolder,filesep,'dogs',filesep];

dogs = [];

DogOutputFolder = [liveExperiment.procFolder, 'dogs', filesep];

dogDir = dir([DogOutputFolder, '*_ch0', num2str(histoneChannel), '.*']);
dogStr = 'dogStack_';

if length(dir(DogOutputFolder)) <= 2
    error('Filtered movie files not found. Did you run FilterMovie?')
end

histoneStack = double(getMovieFrame(liveExperiment, cytoplasmFrames(1), histoneChannel));
histoneMaxProj = max(histoneStack, [], 3);
CytoplasmRegionOutfile = [DropboxFolder, filesep, Prefix, filesep, 'CytoplasmRegionBoundaries.mat'];
if isfile(CytoplasmRegionOutfile)
    load(CytoplasmRegionOutfile);
else
    [MinRow, MaxRow, MinColumn, MaxColumn] = SelectCytoplasmRegion(histoneMaxProj);
    save(CytoplasmRegionOutfile, 'MinRow', 'MaxRow', 'MinColumn', 'MaxColumn');
end


% Load cytoplasmFrame info 
load([DropboxFolder, filesep, Prefix, filesep, 'cytoplasmFrames.mat']);

ncs = [liveExperiment.nc9, liveExperiment.nc10, liveExperiment.nc11,...
    liveExperiment.nc12, liveExperiment.nc13, liveExperiment.nc14, NFrames];


%% 

NuclearCount = zeros(NFrames, zDim);
UseImagesForCalculation = zeros(NFrames, zDim);
MeanCytoplasmFluo = zeros(NFrames, zDim);
StdCytoplasmFluo = zeros(NFrames, zDim);
PixelCountCytoplasmFluo = zeros(NFrames, zDim);
APCytoplasmFluo = zeros(NFrames, length(APbinIDs), zDim);
APStdCytoplasmFluo = zeros(NFrames, length(APbinIDs), zDim);
APPixelCountCytoplasmFluo = zeros(NFrames, length(APbinIDs), zDim);
MeanCytoplasm3DFluo = zeros(NFrames);
StdCytoplasm3DFluo = zeros(NFrames);
PixelCountCytoplasm3DFluo = zeros(NFrames);
APCytoplasm3DFluo = zeros(NFrames, length(APbinIDs));
APStdCytoplasm3DFluo = zeros(NFrames, length(APbinIDs));
APPixelCountCytoplasm3DFluo = zeros(NFrames, length(APbinIDs));
MeanCytoplasmFluoWithNuclearMin = zeros(NFrames);
StdCytoplasmFluoWithNuclearMin = zeros(NFrames);
PixelCountCytoplasmFluoWithNuclearMin = zeros(NFrames);
APCytoplasmFluoWithNuclearMin = zeros(NFrames, length(APbinIDs));
APStdCytoplasmFluoWithNuclearMin = zeros(NFrames, length(APbinIDs));
APPixelCountCytoplasmFluoWithNuclearMin = zeros(NFrames, length(APbinIDs));
Cytoplasm3DMasks = zeros(NFrames, yDim, xDim,zDim);
waitbarFigure = waitbar(0, 'Segmenting cytoplasm');

for i = 1:length(inputChannels)
    ChN = inputChannels(i);
    for idx=1:length(cytoplasmFrames)
        CurrentFrame=cytoplasmFrames(idx);
        try waitbar(idx/length(cytoplasmFrames), waitbarFigure); catch; end
        % GM: added the uint16 conversion 

        imStack = double(getMovieFrame(liveExperiment, CurrentFrame, ChN));
        probStackFile = [DogOutputFolder, 'prob', Prefix, '_',...
            iIndex(CurrentFrame, 3),...
            '_ch', iIndex(histoneChannel, 2)];
        
        if exist([probStackFile, '.mat'], 'file')
            dogStack = load([probStackFile,'.mat'], 'dogStack');
            dogStack = dogStack.dogStack;
        elseif exist([probStackFile, '.tif'], 'file')
            dogStack = imreadStack2([probStackFile, '.tif'], yDim,...
                xDim, zDim);
        else
            error('Cannot find any dogs in the ProcessedData folder. Are they missing or incorrectly named?')
        end
        
        probStack = dogStack/max(max(max(dogStack)));
  
        for z=2:(1+zDim)
            nucleiImage = probStack(:,:,z);
            nucleiMask = zeros(size(nucleiImage));
            nucleiMask(nucleiImage > 0.5) = 1;
            LabledNuclei = bwlabel(nucleiMask);
            NuclearCount(CurrentFrame, z-1) = max(max(LabledNuclei));
            CytoplasmMask = zeros(size(nucleiImage));
            CytoplasmMask(nucleiImage < 0.2) = 1;
            if MinRow > 1
                CytoplasmMask(1:MinRow-1, :) = 0;
            end
            if MaxRow < yDim
                CytoplasmMask(MaxRow+1:end, :) = 0;
            end
            if MinColumn > 1
                CytoplasmMask(:,1:MinColumn-1) = 0;
            end
            if MaxColumn < xDim
                CytoplasmMask(:,MaxColumn+1:end) = 0;
            end
            Cytoplasm3DMasks(CurrentFrame, :,:,z-1) = CytoplasmMask;
            CytoplasmImage = imStack(:,:,z);
            CytoplasmFluoVector = CytoplasmImage((CytoplasmMask > 0)).';
            MeanCytoplasmFluo(CurrentFrame, z-1) = mean(CytoplasmFluoVector);
            StdCytoplasmFluo(CurrentFrame, z-1) = std(CytoplasmFluoVector);
            PixelCountCytoplasmFluo(CurrentFrame, z-1) = length(CytoplasmFluoVector);
            
            
            for APbin=MinAPbin:MaxAPbin
                CytoplasmImage = imStack(:,:,z);
                CytoplasmFluoVector = CytoplasmImage(((CytoplasmMask > 0) & (APPosBinImage == APbin))).';
                APCytoplasmFluo(CurrentFrame, APbin, z-1) = mean(CytoplasmFluoVector);
                APStdCytoplasmFluo(CurrentFrame, APbin, z-1) = std(CytoplasmFluoVector);
                APPixelCountCytoplasmFluo(CurrentFrame, APbin, z-1) = length(CytoplasmFluoVector);
                
            end
        end
        CytoplasmStack = imStack(:,:,2:zDim+1);
        Cytoplasm3DFluoVector = CytoplasmStack(Cytoplasm3DMasks(CurrentFrame, :,:,:)  > 0).';
        MeanCytoplasm3DFluo(CurrentFrame) = mean(Cytoplasm3DFluoVector);
        StdCytoplasm3DFluo(CurrentFrame) = std(Cytoplasm3DFluoVector);
        PixelCountCytoplasm3DFluo(CurrentFrame) = length(Cytoplasm3DFluoVector);
        for APbin=MinAPbin:MaxAPbin
            CytoplasmStack = imStack(:,:,2:zDim+1);
            CytoplasmFluoVector = CytoplasmStack((squeeze(Cytoplasm3DMasks(CurrentFrame, :,:,:))  > 0) & (APPosBinImage3D == APbin)).';
            APCytoplasm3DFluo(CurrentFrame, APbin) = mean(CytoplasmFluoVector);
            APStdCytoplasm3DFluo(CurrentFrame, APbin) = std(CytoplasmFluoVector);
            APPixelCountCytoplasm3DFluo(CurrentFrame, APbin) = length(CytoplasmFluoVector);
        end
        
       
        
    end  
    for nc_idx = find(ncs(1:6) > 0)
        NC = nc_idx + 8;
        MaxCycleNuclearCount = max(max(NuclearCount(ncs(nc_idx):(ncs(nc_idx+1)-1),:)));

        for frame=ncs(nc_idx):(ncs(nc_idx+1)-1)
            UseImagesForCalculation(frame,:) = NuclearCount(frame,:) >= .5*MaxCycleNuclearCount;
        end
    end
    for idx=1:length(cytoplasmFrames)
        CurrentFrame=cytoplasmFrames(idx);
        try waitbar(idx/length(cytoplasmFrames), waitbarFigure); catch; end
        if sum(UseImagesForCalculation(CurrentFrame,:) == 0)
            continue 
        end
        
        % GM: added the uint16 conversion 

        imStack = double(getMovieFrame(liveExperiment, CurrentFrame, ChN));
        probStackFile = [DogOutputFolder, 'prob', Prefix, '_',...
            iIndex(CurrentFrame, 3),...
            '_ch', iIndex(histoneChannel, 2)];
        
        if exist([probStackFile, '.mat'], 'file')
            dogStack = load([probStackFile,'.mat'], 'dogStack');
            dogStack = dogStack.dogStack;
        elseif exist([probStackFile, '.tif'], 'file')
            dogStack = imreadStack2([probStackFile, '.tif'], yDim,...
                xDim, zDim);
        else
            error('Cannot find any dogs in the ProcessedData folder. Are they missing or incorrectly named?')
        end
        
        probStack = dogStack/max(max(max(dogStack)));
       
        CytoplasmStack = imStack(:,:,2:zDim+1);
        Cyto3DMaskV2 = zeros(size(squeeze(Cytoplasm3DMasks(CurrentFrame,:,:,:))));
        Cyto3DMaskV2(squeeze(Cytoplasm3DMasks(CurrentFrame, :,:,:)) > 0) = 1;
        Cyto3DMaskV2(:,:,UseImagesForCalculation(CurrentFrame,:)==0) = 0;
        Cytoplasm3DFluoVector = CytoplasmStack(Cyto3DMaskV2 > 0).';
        MeanCytoplasmFluoWithNuclearMin(CurrentFrame) = mean(Cytoplasm3DFluoVector);
        StdCytoplasmFluoWithNuclearMin(CurrentFrame) = std(Cytoplasm3DFluoVector);
        PixelCountCytoplasmFluoWithNuclearMin(CurrentFrame) = length(Cytoplasm3DFluoVector);
        for APbin=MinAPbin:MaxAPbin
            CytoplasmStack = imStack(:,:,2:zDim+1);
            Cyto3DMaskV2 = zeros(size(squeeze(Cytoplasm3DMasks(CurrentFrame,:,:,:))));
            Cyto3DMaskV2((squeeze(Cytoplasm3DMasks(CurrentFrame, :,:,:)) > 0) & (APPosBinImage3D ~= APbin)) = 1;
            Cyto3DMaskV2(:,:,UseImagesForCalculation(CurrentFrame,:)==0) = 0;
            CytoplasmFluoVector = CytoplasmStack((Cyto3DMaskV2 > 0)).';
            APCytoplasmFluoWithNuclearMin(CurrentFrame, APbin) = mean(CytoplasmFluoVector);
            APStdCytoplasmFluoWithNuclearMin(CurrentFrame, APbin) = std(CytoplasmFluoVector);
            APPixelCountCytoplasmFluoWithNuclearMin(CurrentFrame, APbin) = length(CytoplasmFluoVector);
        end
        
       
        
    end  
end





%% 

outfile = [DropboxFolder, filesep, Prefix, filesep, 'CytoplasmFluorescence.mat'];
save(outfile, 'NuclearCount', 'MeanCytoplasmFluo', 'StdCytoplasmFluo', 'PixelCountCytoplasmFluo',...
    'APCytoplasmFluo','APStdCytoplasmFluo','APPixelCountCytoplasmFluo',...
    'MeanCytoplasm3DFluo', 'StdCytoplasm3DFluo', 'PixelCountCytoplasm3DFluo',...
    'APCytoplasm3DFluo','APStdCytoplasm3DFluo','APPixelCountCytoplasm3DFluo',...
    'MeanCytoplasmFluoWithNuclearMin', 'StdCytoplasmFluoWithNuclearMin',...
    'PixelCountCytoplasmFluoWithNuclearMin', 'APCytoplasmFluoWithNuclearMin',...
    'APStdCytoplasmFluoWithNuclearMin','APPixelCountCytoplasmFluoWithNuclearMin' );


