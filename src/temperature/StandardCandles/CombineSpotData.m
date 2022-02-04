% CombineSpotData.m
% author: Gabriella Martini
% date created: 1/26/22
% date last modified: 1/26/22
%
% Combines Data from individually segmented and picked frames

function CombineSpotData(Prefixes,outfolder,varargin)


UseGoodSpots = true;


k=1;
while k<=length(varargin)
    if strcmpi(varargin{k},'useallspots')
        UseGoodSpots=false;
        k=k+1;
    end
end

TimeStamp = datestr(now);
outdir = ['S:/Gabriella/Dropbox\StandardCandles\' outfolder];
mkdir(outdir);
plotdir =  ['S:/Gabriella/Dropbox/FixedTissueExperiments/StandardCandles/' outfolder];
mkdir(plotdir);
%%

for j = 1:length(Prefixes)
    liveExperiments{j} = LiveExperiment(Prefixes{j});
    if UseGoodSpots
        outpath = [liveExperiments{j}.resultsFolder, filesep, 'StoreSpotInfo.mat'];
        load(outpath, 'GoodSpots');
        EmbryoID = num2cell(j*ones(1,length(GoodSpots.Fits)));
        [GoodSpots.Fits(:).EmbryoIDs] = deal(EmbryoID{:});
        if j == 1
            
            
            Spots = GoodSpots;
            
            
        else
            Spots.Fits = [Spots.Fits GoodSpots.Fits];
        end
    else
        SpotsTemp = getSpots(liveExperiments{j});
        EmbryoID = num2cell(j*ones(1,length(SpotsTemp.Fits)));
        [SpotsTemp.Fits(:).EmbryoIDs] = deal(EmbryoID{:});
        if j == 1
            
            
            Spots = SpotsTemp;
            
            
        else
            Spots.Fits = [Spots.Fits SpotsTemp.Fits];
        end
    end
end

clear GoodSpots outpath EmbryoID

%%
NumSpots = length(Spots.Fits);


xDoG = [];
yDoG = [];
% x3D= [];
% y3D= [];
% z3D= [];
bwIntensity = [];
bwArea = [];
GaussianIntensityCorrected = [];
GaussianZCorrected = [];
GaussianIntensityCropped = [];
GaussianZCropped = [];
GaussianIntensitySmallSnip = [];

GaussianIntensity = [];
GaussianZ = [];
GaussianIntensityArea = [Spots.Fits(:).GaussianIntensityArea];
GaussianIntensityBackground = [Spots.Fits(:).GaussianIntensityBackground];
GaussianIntensityZBackground = [];
GaussianIntensityCorrectedArea = [Spots.Fits(:).GaussianIntensityCorrectedArea];
GaussianIntensityCorrectedBackground = [Spots.Fits(:).GaussianIntensityCorrectedBackground];
GaussianIntensityCorrectedZBackground = [];
GaussianIntensityCroppedArea = [Spots.Fits(:).GaussianIntensityCroppedArea];
GaussianIntensityCroppedBackground =  [Spots.Fits(:).GaussianIntensityCroppedBackground];
GaussianIntensityCroppedZBackground = [];
GaussianZSmallSnip = [Spots.Fits(:).GaussianZSmallSnip];
GaussianPeakSmallSnip =[Spots.Fits(:).GaussianPeakSmallSnip];
GaussianWidthSmallSnip = [Spots.Fits(:).GaussianWidthSmallSnip];
GaussianWidthXSmallSnip =  [Spots.Fits(:).GaussianWidthXSmallSnip];
GaussianWidthYSmallSnip = [Spots.Fits(:).GaussianWidthYSmallSnip];
GaussianRhoSmallSnip = [Spots.Fits(:).GaussianRhoSmallSnip];


GaussianAreaSmallSnip = [Spots.Fits(:).GaussianAreaSmallSnip];
GaussianZBackgroundSmallSnip = [];
GaussianBackgroundSmallSnip = [Spots.Fits(:).GaussianBackgroundSmallSnip];

GaussianKernelIntensity = [];
GaussianKernelZ = [Spots.Fits(:).GaussianKernelZ];
GaussianKernelArea = [Spots.Fits(:).GaussianKernelArea];
GaussianKernelZBackground = [];
GaussianKernelBackground = [Spots.Fits(:).GaussianKernelBackground];


IntegratedGaussianKernelIntensity = [];
IntegratedGaussianKernelZ = [Spots.Fits(:).IntegratedGaussianKernelZ];
IntegratedGaussianKernelArea = [Spots.Fits(:).IntegratedGaussianKernelArea];
IntegratedGaussianKernelZBackground = [];
IntegratedGaussianKernelBackground = [Spots.Fits(:).IntegratedGaussianKernelBackground];



GaussianInfo_A  = [];
GaussianError_A  = [];
GaussianInfo_x0  = [];
GaussianError_x0  =[];
GaussianInfo_y0  = [];
GaussianError_y0  = [];
GaussianInfo_rho  = [];
GaussianError_rho  = [];
GaussianInfo_sigx  = [];
GaussianError_sigx  = [];
GaussianInfo_sigy  = [];
GaussianError_sigy  = [];
GaussianInfo_offset  = [];
GaussianError_offset  =[];
GaussianInfo_offx  = [];
GaussianError_offx  = [];
GaussianInfo_offy  = [];
GaussianError_offy  = [];
zCount = [];
zCenter = [];

EmbryoIDs = [Spots.Fits(:).EmbryoIDs];
TotalPixels = [];


bwDiameter = [];
bwCircularity = [];
bwEccentricity = [];
bwMajorAxisLength = [];
bwMinorAxisLength = [];


dogFixedAreaIntensity =  [];
DOGIntensity = [];
FixedAreaIntensity= [];

%%
for SpotIndex = 1:NumSpots
    xDoG(SpotIndex)  =  mean(Spots.Fits(SpotIndex).xDoG);
    yDoG(SpotIndex)  =  mean(Spots.Fits(SpotIndex).yDoG);
    zCount(SpotIndex)  =  length(Spots.Fits(SpotIndex).z);
    zCenter(SpotIndex) = mean(Spots.Fits(SpotIndex).z);
    bwIntensity(SpotIndex)  = max(Spots.Fits(SpotIndex).bwIntensity);
    dogFixedAreaIntensity(SpotIndex)  = max(Spots.Fits(SpotIndex).dogFixedAreaIntensity);
    FixedAreaIntensity(SpotIndex)  = max(Spots.Fits(SpotIndex).FixedAreaIntensity);
    DOGIntensity(SpotIndex)  = max(Spots.Fits(SpotIndex).DOGIntensity);
    bwArea(SpotIndex)  = max(Spots.Fits(SpotIndex).bwArea);
    MatchIndex = find(single(nanmax(Spots.Fits(SpotIndex).GaussianIntensity)) == single(Spots.Fits(SpotIndex).GaussianIntensity), 1);
    
    [maxI, maxZindex] = max(Spots.Fits(SpotIndex).GaussianIntensity);
    GaussianIntensity(SpotIndex)  = maxI;
    GaussianZ(SpotIndex) = Spots.Fits(SpotIndex).z(maxZindex);
    GaussianIntensityZBackground(SpotIndex) = Spots.Fits(SpotIndex).GaussianIntensityZBackground(maxZindex);
    
    [maxI, maxZindex] = max(Spots.Fits(SpotIndex).GaussianIntensityCorrected);
    GaussianIntensityCorrected(SpotIndex)  = maxI;
    GaussianZCorrected(SpotIndex) = Spots.Fits(SpotIndex).z(maxZindex);
    GaussianIntensityCorrectedZBackground(SpotIndex) = Spots.Fits(SpotIndex).GaussianIntensityCorrectedZBackground(maxZindex);
    
    [maxI, maxZindex] = max(Spots.Fits(SpotIndex).GaussianIntensityCropped);
    GaussianIntensityCropped(SpotIndex)  = maxI;
    GaussianZCropped(SpotIndex) = Spots.Fits(SpotIndex).z(maxZindex);
    GaussianIntensityCroppedZBackground(SpotIndex) = Spots.Fits(SpotIndex).GaussianIntensityZBackground(maxZindex);
    
    [maxI, maxZindex]=  max(Spots.Fits(SpotIndex).GaussianIntensitySmallSnip);
    GaussianIntensitySmallSnip(SpotIndex) = maxI;
    GaussianZBackgroundSmallSnip(SpotIndex) = Spots.Fits(SpotIndex).GaussianZBackgroundSmallSnip(maxZindex);
    
    [maxI, maxZindex]= max(Spots.Fits(SpotIndex).GaussianKernelIntensity);
    GaussianKernelIntensity(SpotIndex) = maxI;
    GaussianKernelZBackground(SpotIndex) = Spots.Fits(SpotIndex).GaussianKernelZBackground(maxZindex);
    
    
    [maxI, maxZindex]= max(Spots.Fits(SpotIndex).IntegratedGaussianKernelIntensity);
    IntegratedGaussianKernelIntensity(SpotIndex) = maxI;
    IntegratedGaussianKernelZBackground(SpotIndex) = Spots.Fits(SpotIndex).IntegratedGaussianKernelZBackground(maxZindex);
    


    bwDiameter(SpotIndex)  = Spots.Fits(SpotIndex).bwDiameter(MatchIndex);
    bwCircularity(SpotIndex)  = Spots.Fits(SpotIndex).bwCircularity(MatchIndex);
    bwEccentricity(SpotIndex)  = Spots.Fits(SpotIndex).bwEccentricity(MatchIndex);
    bwMajorAxisLength(SpotIndex)  = Spots.Fits(SpotIndex).bwMajorAxisLength(MatchIndex);
    bwMinorAxisLength(SpotIndex)  = Spots.Fits(SpotIndex).bwMinorAxisLength(MatchIndex);
    GaussianInfo_A(SpotIndex)= Spots.Fits(SpotIndex).GaussianInfo(MatchIndex).A;
    GaussianError_A(SpotIndex)= Spots.Fits(SpotIndex).GaussianError(MatchIndex).A;
    GaussianInfo_x0(SpotIndex)= Spots.Fits(SpotIndex).GaussianInfo(MatchIndex).x0;
    GaussianError_x0(SpotIndex)= Spots.Fits(SpotIndex).GaussianError(MatchIndex).x0;
    GaussianInfo_y0(SpotIndex) = Spots.Fits(SpotIndex).GaussianInfo(MatchIndex).y0;
    GaussianError_y0(SpotIndex) = Spots.Fits(SpotIndex).GaussianError(MatchIndex).y0;
    GaussianInfo_rho(SpotIndex) = Spots.Fits(SpotIndex).GaussianInfo(MatchIndex).rho;
    GaussianError_rho(SpotIndex) = Spots.Fits(SpotIndex).GaussianError(MatchIndex).rho;
    GaussianInfo_sigx(SpotIndex) = Spots.Fits(SpotIndex).GaussianInfo(MatchIndex).sigma_x;
    GaussianError_sigx(SpotIndex) = Spots.Fits(SpotIndex).GaussianError(MatchIndex).sigma_x;
    GaussianInfo_sigy(SpotIndex) = Spots.Fits(SpotIndex).GaussianInfo(MatchIndex).sigma_y;
    GaussianError_sigy(SpotIndex) = Spots.Fits(SpotIndex).GaussianError(MatchIndex).sigma_y;
    GaussianInfo_offset(SpotIndex) = Spots.Fits(SpotIndex).GaussianInfo(MatchIndex).offset;
    GaussianError_offset(SpotIndex) = Spots.Fits(SpotIndex).GaussianError(MatchIndex).offset;
    GaussianInfo_offx(SpotIndex) = Spots.Fits(SpotIndex).GaussianInfo(MatchIndex).offset_x;
    GaussianError_offx(SpotIndex) = Spots.Fits(SpotIndex).GaussianError(MatchIndex).offset_x;
    GaussianInfo_offy(SpotIndex)= Spots.Fits(SpotIndex).GaussianInfo(MatchIndex).offset_y;
    GaussianError_offy(SpotIndex) = Spots.Fits(SpotIndex).GaussianError(MatchIndex).offset_y;
    TotalPixels(SpotIndex) = sum(Spots.Fits(SpotIndex).bwArea);
    
    
    
    
end

%%
clear j k MatchIndex SpotIndex maxI maxZindex UseGoodSpots

try
    save([outdir '\SpotInfo.mat'])
catch
    save([outdir '\SpotInfo.mat'],...
        '-v7.3', '-nocompression');
end