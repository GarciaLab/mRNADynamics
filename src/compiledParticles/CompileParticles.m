function CompiledParticles = CompileParticles(varargin)
% CompileParticles(varargin)
%
% DESCRIPTION
% This function puts together all the information we have about particles.
%
% ARGUMENTS
% varargin: A cell in which the first element is the prefix string of the data set
%           to analyze. Subsequent elements can be the options below.
%
% OPTIONS
% 'ForceAP': Force AP detection even if it's there already.
%
% 'SkipTraces': Don't output the individual traces.
%
% 'SkipFluctuations': Don't generate the plots of the correlation of signal
%                   and offset.
%
% 'SkipFits': Don't do the fits
%
% 'SkipMovie': Don't do the movie
%
% 'SkipAll': Skip all that can be skipped
%
% 'ApproveAll': Approves all particles. This is useful if we want to do a
%             quick check of, for example, the AP profiles
%
% 'MinParticles', N: Set the threshold for the minimum number of particles (N) per
%               AP bin for compilation. Default is 4.
%
% 'MinTime', M: %Require particles to exist for time M or else discard.
%               Default is 1.
%
% 'ROI', ROI1, ROI2: For Region of Interest (ROI) data. Assume that the ROI is top half of the imaging window.
%           Note that the origin is the left top of the image.
%           ROI1 and ROI2 are the y-position of the ROI rectangle. ROI1 is the
%           lower boundary of ROI and ROI2 is the upper boundary of non-ROI
%           since there is almost always scattering in the middle
% 'noHist': Force the code to assume there's no nuclear channel.
% 'doSingleFits': Generate single trace fits. Added by EL&AR
% 'manualSingleFits' : Compile values from manually generated single trace fits to
%               CompiledParticles.mat, the field names are 'fittedSlope', 'fittedTon' for
%               the initial slope and T_on respectively.
% 'optionalResults' : if you want to use a different dropbox folder
% 'minBinSize': changes the minimum size of allowed AP bins
% 'edgeWidth': remove ellipses and particles close to the boundary of the
% field of view

% Author (contact): Hernan Garcia (hggarcia@berkeley.edu)
% Created:
% Last Updated: 6/17/17 (AR)
%
% Documented by: Hernan Garcia (hggarcia@berkeley.edu)

%%

cleanupObj = onCleanup(@myCleanupFun);
warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% INITIALIZE ALL SAVED VARIABLES
% Please initialize any new variables you have added and want to save!!!!

savedVariables = {};

APFilter  =  {};
APbinArea = [];
APbinID = [];
AllTracesAP  =  {};
AllTracesVector = {};
CompiledParticles = {};
DVFilter = {};
DVbinArea = [];
DVbinID = [];
ElapsedTime = [];
elapsedTime_min = [];
EllipsePos = {};
EllipsePos_DV = {};
EllipsesFilteredPos = [];
EllipsesOnAP = {};
EllipsesOnDV = {};
FilteredParticlesPos = [];
MaxAPIndex = [];
MaxCyto = [];
MaxDVIndex = [];
MaxFrame = {};
MeanCyto = [];
MeanOffsetVector = [];
MeanSlopeVectorAP = {};
MeanSlopeVectorDV = {};
MeanVectorAP = {};
MeanVectorAP_ROI = {};
MeanVectorAP_nonROI = {};
MeanVectorAll = {};
MeanVectorAllAP = {};
MeanVectorAllDV = {};
MeanVectorAnterior = {};
MeanVectorDV = {};
MeanVectorDV_ROI = {};
MeanVectorDV_nonROI = {};
MedianCyto = [];
MinAPIndex = [];
MinDVIndex = [];
NEllipsesAP = [];
NEllipsesDV = [];
NOffsetParticles = [];
NParticlesAP = {};
NParticlesAP_ROI = {};
NParticlesAP_nonROI = {};
NParticlesAll = {};
NParticlesDV = {};
NParticlesDV_ROI = {};
NParticlesDV_nonROI = {};
NSlopeAP = {};
NSlopeDV = {};
NewCyclePos = [];
OnRatioAP = {};
OnRatioDV = {};
ParticleCountAP = {};
ParticleCountDV = {};
ParticleCountProbAP = {};
ParticleCountProbDV = {};
SDCyto = [];
SDOffsetVector = [];
SDSlopeVectorAP = {};
SDSlopeVectorDV = {};
SDVectorAP = {};
SDVectorAP_ROI = {};
SDVectorAP_nonROI = {};
SDVectorAll = {};
SDVectorDV = {};
SDVectorDV_ROI = {};
SDVectorDV_nonROI = {};
SEVectorAllAP = {};
SEVectorAllDV = {};
StemLoopEnd = '';
TotalEllipsesAP = [];
TotalEllipsesDV = [];
fittedLineEquations = [];
rateOnAP = [];
timeOnOnAP = [];
rateOnAPCell = [];
timeOnOnAPCell = [];
rateOnAPManual = [];
timeOnOnAPManual = [];
rateOnAPCellManual = [];
timeOnOnAPCellManual = [];
fittedSlope = [];
fittedTon = [];
MeanVector3DAll = {};
MeanVector3DAP = {};

nc9 = [];
nc10 = [];
nc11 = [];
nc12 = [];
nc13 = [];
nc14 = [];
ncFilter = [];
ncFilterID = [];
%%


[Prefix, ForceAP, SkipTraces, SkipFluctuations, SkipFits, SkipMovie, ...
    shouldSkipAll, shouldApproveAll, MinParticles, minTime, ROI,  noHist, ...
    ROI1, ROI2, slimVersion, manualSingleFits,...
    optionalResults, yToManualAlignmentPrompt, minBinSize, edgeWidth] = determineCompileParticlesOptions(varargin);


liveExperiment = LiveExperiment(Prefix);
FilePrefix=[Prefix,'_'];

DropboxFolder = liveExperiment.userResultsFolder;
PreProcPath = liveExperiment.userPreFolder;
ExperimentType = liveExperiment.experimentType;
ExperimentAxis = liveExperiment.experimentAxis;
APResolution = liveExperiment.APResolution;
DVResolution = liveExperiment.DVResolution;
nc9=liveExperiment.nc9; nc10=liveExperiment.nc10; 
nc11=liveExperiment.nc11;nc12=liveExperiment.nc12;
nc13=liveExperiment.nc13;nc14=liveExperiment.nc14;

Channels = liveExperiment.Channels;
Channel1 = Channels{1};
Channel2 = Channels{2};
Channel3 = Channels{3};

ncFrames = [zeros(1,8), nc9, nc10, nc11, nc12, nc13, nc14];


FrameInfo = getFrameInfo(liveExperiment);
pixelSize = liveExperiment.pixelSize_nm;
numFrames = liveExperiment.nFrames;
coatChannels = getCoatChannel(Channel1, Channel2, Channel3);


elapsedTime_min = computeElapsedTime(FrameInfo);
%copying this for backwards compatibility
ElapsedTime = elapsedTime_min;


APExperiment = strcmpi(ExperimentAxis, 'AP');
DVExperiment = strcmpi(ExperimentAxis, 'DV');

%Load Spots and Particles
disp('Loading Particles.mat...');
[Particles, SpotFilter] = getParticles(liveExperiment);
disp('Particles loaded.');
disp('Loading Spots.mat...');
Spots = getSpots(liveExperiment);
disp('Spots loaded.');
if isempty(Particles)
    SkipTraces=1;
    SkipFluctuations=1;
    SkipFits=1;
    SkipMovie=1;
end

%Delete the files in folder where we'll write again.
if ~SkipTraces
    delete([DropboxFolder,filesep,Prefix,filesep,'ParticleTraces',filesep,'*.*'])
end
if ~SkipFits
    % deleting the figures saved from fitSingleTraces.m
    delete([DropboxFolder,filesep,Prefix,filesep,'Fits',filesep,'*.*'])
end

%See if we had any lineage/nuclear information
if exist([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'], 'file') && ~noHist
    load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'], 'schnitzcells')
    haveHistoneChannel=true;
else
    disp('No lineage / nuclear information found. Proceeding without it.');
    haveHistoneChannel=false;
    % initialize the variables so that they can be plugged into the
    % CompileTraces below.
    schnitzcells={};
    Ellipses = {};
end


%Do we need to convert any NaN chars into doubles?
if strcmpi(nc14,'nan')
    nc14=nan;
end
if strcmpi(nc13,'nan')
    nc13=nan;
end
if strcmpi(nc12,'nan')
    nc12=nan;
end
if strcmpi(nc11,'nan')
    nc11=nan;
end
if strcmpi(nc10,'nan')
    nc10=nan;
end
if strcmpi(nc9,'nan')
    nc9=nan;
end


NewCyclePos=[nc9,nc10,nc11,nc12,nc13,nc14];
NewCyclePos=NewCyclePos(~(NewCyclePos==0));
NewCyclePos=NewCyclePos(~isnan(NewCyclePos));



%Create the particle array. This is done so that we can support multiple
%channels. Also figure out the number of channels
if iscell(Particles)
    nSpotChannels=length(Particles);
else
    Particles={Particles};
    nSpotChannels=1;
    Spots={Spots};
    SpotFilter={SpotFilter};
end

%Add the APPosition to Particles if they don't exist yet. Do this only if
%we took AP data. Otherwise just add x and y pixel coordinates

addParticleArgs = {Prefix};
if ~isempty(optionalResults)
    addParticleArgs = [addParticleArgs, optionalResults];
end
if yToManualAlignmentPrompt
    addParticleArgs = [addParticleArgs, 'yToManualAlignmentPrompt'];
end
if ~haveHistoneChannel
    addParticleArgs = [addParticleArgs, 'SkipAlignment'];
end

if APExperiment || DVExperiment
    if ~isfield(Particles{1},'APpos') || ForceAP
        try [Particles, SpotFilter] = AddParticlePosition(addParticleArgs{:});
        catch warning('Failed to add particle position. Is there no full embryo?'); end
    else disp('Using saved AP information (results from AddParticlePosition)'); end
end

if haveHistoneChannel, Ellipses = getEllipses(liveExperiment); end

%Folders for reports
if APExperiment
    mkdir([DropboxFolder,filesep,Prefix,filesep,'APMovie'])
end

if ~shouldSkipAll
    mkdir([DropboxFolder,filesep,Prefix,filesep,'ParticleTraces'])
    mkdir([DropboxFolder,filesep,Prefix,filesep,'TracesFluctuations'])
    mkdir([DropboxFolder,filesep,Prefix,filesep,'Offset'])
    mkdir([DropboxFolder,filesep,Prefix,filesep,'Fits'])
    mkdir([DropboxFolder,filesep,Prefix,filesep,'Probabilities'])
    mkdir([DropboxFolder,filesep,Prefix,filesep,'Various']);
end



%% Binning the axes into AP and DV
fullEmbryoExists = exist([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'], 'file');
shouldConvertToAP =  haveHistoneChannel...
    && (APExperiment || DVExperiment)...
    && fullEmbryoExists;

%Figure out the AP position of each of the nuclei.
if shouldConvertToAP
   [EllipsePos, APAngle, APLength]...
   = convertToFractionalEmbryoLength(Prefix);
end

%Divide the AP and DV axes into bins for generating means
if APExperiment || DVExperiment && fullEmbryoExists
    [APbinID, APbinArea] = binAxis(APResolution, FrameInfo, ...
        coordAZoom, APAngle, APLength, minBinSize, 'AP');
end

if DVExperiment && fullEmbryoExists
    if isempty(DVResolution)
        DVResolution = 100;
    end
    [DVbinID, DVbinArea] = binAxis(DVResolution, FrameInfo, ...
        coordAZoom, APAngle, APLength, minBinSize, 'DV');
end


%% Put together CompiledParticles

CompiledParticles = cell(nSpotChannels,1);

Particles = approveParticles(Particles,...
    shouldApproveAll, nSpotChannels, haveHistoneChannel);

[Particles, CompiledParticles, ncFilter, ncFilterID] =...
    ...
    compileTraces(...
    ...
    nSpotChannels, Particles, haveHistoneChannel, ...
    schnitzcells, minTime, ExperimentAxis, APbinID, APbinArea, CompiledParticles, ...
    Spots, SkipTraces, ncFilterID, ncFilter, ...
    elapsedTime_min, Ellipses, EllipsePos, PreProcPath, ...
    FilePrefix, Prefix, DropboxFolder, numFrames,...
    manualSingleFits, edgeWidth, pixelSize, coatChannels, fullEmbryoExists, liveExperiment);

%% ROI option
% This option is separating the CompiledParticles defined above into
% CompiledParticles_ROI and COmpiledParticles_nonROI
% written by YJK on 10/24/2017
CompiledParticles_ROI = cell(1,nSpotChannels); CompiledParticles_nonROI = cell(1,nSpotChannels);

if ROI
    for ChN=1:nSpotChannels
        % separate the CompileParticles into CompiledParticles_ROI and
        % CompiledParticles_nonROI using the Threshold (y position)
        t=1;
        s=1;
        
        % Use the ROI1 and ROI2 to split the Particles
        for ParticleIndex=1:length(CompiledParticles{ChN})
            if nanmean(CompiledParticles{ChN}(ParticleIndex).yPos) < ROI1
                CompiledParticles_ROI{ChN}(t)=CompiledParticles{ChN}(ParticleIndex);
                t=t+1;
            elseif nanmean(CompiledParticles{ChN}(ParticleIndex).yPos) > ROI2
                CompiledParticles_nonROI{ChN}(s)=CompiledParticles{ChN}(ParticleIndex);
                s=s+1;
            end
        end
    end
end
%% Bin particles into AP, DV and NC bins. Create filters that store binning information.

APFilter_ROI = []; APFilter_nonROI = []; DVFilter_ROI = []; DVFilter_nonROI = [];

if fullEmbryoExists
    [ncFilterID, ncFilter, APFilter, APFilter_ROI, APFilter_nonROI, ...
        DVFilter, DVFilter_ROI, DVFilter_nonROI] =...
        ...
        binParticles(...
        ...
        nc9, nc10, nc11, nc12,...
        nc13, nc14, nSpotChannels, CompiledParticles, ExperimentAxis, ROI,...
        APbinID, DVbinID, CompiledParticles_ROI, CompiledParticles_nonROI);
    
end


%% Averaging data
if fullEmbryoExists
    [AllTracesVector, AllTracesAP, AllTracesDV, MeanVectorAP_ROI, ...
        SDVectorAP_ROI, NParticlesAP_ROI, MeanVectorAP_nonROI, SDVectorAP_nonROI, ...
        NParticlesAP_nonROI, MeanVectorAP, SDVectorAP, NParticlesAP, MeanVectorDV_ROI, ...
        SDVectorDV_ROI, NParticlesDV_ROI, MeanVectorDV_nonROI, SDVectorDV_nonROI, ...
        NParticlesDV_nonROI, MeanVectorDV, SDVectorDV, NParticlesDV, ...
        MeanVectorAnterior, MeanVectorAll, SDVectorAll, NParticlesAll, MaxFrame, MeanVector3DAll, MeanVector3DAP] =...
        ...
        getAxisStatistics(...
        ...
        nSpotChannels, CompiledParticles, FrameInfo, ExperimentAxis, ...
        APFilter, ROI, CompiledParticles_ROI, CompiledParticles_nonROI, ...
        APFilter_ROI, APFilter_nonROI, NewCyclePos, DVFilter_ROI, ...
        DVFilter_nonROI, DVFilter);
else
    AllTracesVector = {};
    AllTracesVector{1} =...
        createAllTracesVector(FrameInfo,CompiledParticles{1},'NoAP');
end

if ~slimVersion && fullEmbryoExists
    %% Instantaneous rate of change
    
    [CompiledParticles, MeanSlopeVectorAP, SDSlopeVectorAP, NSlopeAP]...
        = instantRateOfChange(nSpotChannels, CompiledParticles, elapsedTime_min, ...
        ExperimentAxis, APFilter, MeanVectorAP, MeanSlopeVectorAP, ...
        SDSlopeVectorAP, NSlopeAP);
    
    %% Integrating each particle
    
    CompiledParticles = integrateParticles(nSpotChannels, elapsedTime_min, CompiledParticles);
    
    %% Information about the cytoplasm
    %If the nuclear masks are present then use them. Otherwise just calculate
    %the median of the images as a function of time
    %     if ~SkipAll
    %         [MeanCyto, SDCyto, MaxCyto, MedianCyto] =...
    %             ...
    %             getCytoplasmStatistics(...
    %             ...
    %             APExperiment, HistoneChannel, Prefix, numFrames, PreProcPath,...
    %             FrameInfo, NChannels);
    %     end
    
    %% Offset and fluctuations
    
    if nSpotChannels == 1
        
        [MeanOffsetVector, SDOffsetVector, NOffsetParticles] =...
            ...
            offsetAndFlux(...
            ...
            SkipFluctuations, ncFilter, elapsedTime_min, CompiledParticles, DropboxFolder, ...
            Prefix, ExperimentAxis, pixelSize, MeanVectorAll, SDVectorAll, MaxFrame, numFrames, shouldSkipAll);
        
    end
    
    %% Rate of mRNA production
    
    CompiledParticles = mRNAProdRate(nSpotChannels, CompiledParticles, ...
        ncFilter, elapsedTime_min, SkipFits, DropboxFolder, Prefix);
    
    
    %% First frames
    if ~shouldSkipAll
        plotFirstFrames(nSpotChannels, haveHistoneChannel, nc13, nc14, ...
            CompiledParticles, DropboxFolder, Prefix, elapsedTime_min, ExperimentAxis);
    end
    
    %% AP position of particle vs nucleus
    
    if haveHistoneChannel && APExperiment
        CompiledParticles = APPosParticleVsNucleus(nSpotChannels, ...
            CompiledParticles, schnitzcells, EllipsePos, DropboxFolder, Prefix, shouldSkipAll);
    end
    
    
    %% Fitting shapes to single traces (includes time on and initial rate of loading)
    % This section of code is will fit line segments piece-wise to the
    % single traces. fittedLineEquations correspond to the stored fitted lines
    % of the particles, where the indexing is as follows:
    % fittedLineEquations(particleNumber) with fields: Coefficients,
    % ErrorEstimation, and frameIndex, which are described below. This currently
    % does not support more than one channel. Please contact Emma to work on
    % implementing it for two channels.
    if ~SkipFits
        [CompiledParticles, fittedLineEquations] = fitShapesToTraces(Prefix, ...
            Particles, schnitzcells, FrameInfo, elapsedTime_min, CompiledParticles,Spots);
    end
    
    %% Probability of being on
    
    if haveHistoneChannel && (APExperiment || DVExperiment) && fullEmbryoExists
        
        [NEllipsesAP, MeanVectorAllAP, SEVectorAllAP, EllipsesFilteredPos, ...
            FilteredParticlesPos, OnRatioAP, ParticleCountAP, ParticleCountProbAP, ...
            EllipsesOnAP, rateOnAP, rateOnAPCell, timeOnOnAP, timeOnOnAPCell,...
            TotalEllipsesAP, rateOnAPManual, rateOnAPCellManual, timeOnOnAPManual, timeOnOnAPCellManual...
            ]...
            ...
            = computeAPFractionOn(...
            ...
            nSpotChannels, Particles, schnitzcells,...
            CompiledParticles, Ellipses, APbinID, FrameInfo, elapsedTime_min, DropboxFolder, ...
            Prefix, EllipsePos, nc12, nc13, nc14, numFrames, SkipFits, shouldSkipAll, ...
            APbinArea,pixelSize, manualSingleFits, edgeWidth,  DVbinArea, DVbinID, EllipsePos_DV);
        
    end
    
    % DV version. This should instead be smoothly integrated with the AP
    % version since there's a lot of duplicate code here.
    if haveHistoneChannel && DVExperiment && fullEmbryoExists %JAKE: Need to change this later
        
        [NEllipsesDV, MeanVectorAllDV, SEVectorAllDV, OnRatioDV, ParticleCountDV, ...
            ParticleCountProbDV, TotalEllipsesDV, EllipsesOnDV, EllipsesFilteredPos, ...
            FilteredParticlesPos] =...
            ...
            DVProbOn(nSpotChannels, Particles, schnitzcells, ...
            CompiledParticles, Ellipses, FrameInfo, DropboxFolder, Prefix, ...
            elapsedTime_min, DVbinID, EllipsePos_DV, nc12, nc13, nc14, numFrames, ...
            DVbinArea,edgeWidth, pixelSize);
        
    end
    
    %% Calculation of particle speed
    try
        calcParticleSpeeds(nSpotChannels, Particles, ...
            Spots, elapsedTime_min, schnitzcells, Ellipses);
    end
    
    %% Movie of AP profile
    
    %I want to make a movie of the average fluorescence as a function of AP as
    %a function of time. In order to make life easier I'll just export to a
    %folder. I can then load everything in ImageJ.
    
    if ~SkipMovie && APExperiment && fullEmbryoExists
        APProfileMovie(MeanVectorAP, NParticlesAP, MinParticles, ...
            APbinID, SDVectorAP, FrameInfo, elapsedTime_min, DropboxFolder, Prefix, nc9, nc10, nc11, nc12, nc13, nc14)
    end
    
end


CompiledParticles = addCycle(CompiledParticles, ncFrames);

%% Input-output

%Compile the nuclear fluorescence information if we have the appropriate
%experiment type and axis
try
    if strcmpi(ExperimentType,'inputoutput')
        if ~ROI, CompileNuclearProtein(Prefix);
        else, CompileNuclearProtein(Prefix,'ROI',ROI1,ROI2); end
    end
catch, warning('Couldn''t run CompileNuclearProtein.'); end



%% Save everything

%Now save all the information

savedVariables = [savedVariables,'APFilter', 'APbinArea', 'APbinID', 'AllTracesAP',...
    'AllTracesVector', 'CompiledParticles', 'DVFilter', 'DVbinArea',...
    'DVbinID', 'ElapsedTime', 'EllipsePos','EllipsePos_DV', 'EllipsesFilteredPos',...
    'EllipsesOnAP', 'EllipsesOnDV', 'FilteredParticlesPos', 'MaxAPIndex',...
    'MaxCyto', 'MaxDVIndex', 'MaxFrame', 'MeanCyto',...
    'MeanOffsetVector', 'MeanSlopeVectorAP', 'MeanSlopeVectorDV', 'MeanVectorAP',...
    'MeanVectorAP_ROI', 'MeanVectorAP_nonROI', 'MeanVectorAll', 'MeanVectorAllAP',...
    'MeanVectorAllDV', 'MeanVectorAnterior', 'MeanVectorDV', 'MeanVectorDV_ROI',...
    'MeanVectorDV_nonROI', 'MedianCyto', 'MinAPIndex', 'MinDVIndex',...
    'NEllipsesAP', 'NEllipsesDV', 'NOffsetParticles', 'NParticlesAP',...
    'NParticlesAP_ROI', 'NParticlesAP_nonROI', 'NParticlesAll', 'NParticlesDV',...
    'NParticlesDV_ROI', 'NParticlesDV_nonROI', 'NSlopeAP', 'NSlopeDV',...
    'NewCyclePos', 'OnRatioAP', 'OnRatioDV', 'ParticleCountAP',...
    'ParticleCountDV', 'ParticleCountProbAP', 'ParticleCountProbDV', 'Prefix',...
    'SDCyto', 'SDOffsetVector', 'SDSlopeVectorAP', 'SDSlopeVectorDV',...
    'SDVectorAP', 'SDVectorAP_ROI', 'SDVectorAP_nonROI', 'SDVectorAll',...
    'SDVectorDV', 'SDVectorDV_ROI', 'SDVectorDV_nonROI', 'SEVectorAllAP',...
    'SEVectorAllDV', 'StemLoopEnd', 'TotalEllipsesAP', 'TotalEllipsesDV',...
    'fittedLineEquations', 'rateOnAP', 'timeOnOnAP','rateOnAPCell', 'timeOnOnAPCell','nc10', 'nc11', 'nc12',...
    'nc13', 'nc14', 'nc9', 'ncFilter',...
    'ncFilterID', 'rateOnAPManual', 'rateOnAPCellManual', 'timeOnOnAPManual', 'timeOnOnAPCellManual'...
    'MeanVector3DAP', 'MeanVector3DAll'];

CompiledParticlesFile = [DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'];

try
    save(CompiledParticlesFile, savedVariables{:},'-v6');
catch
    %save as 7.3 only if we really need to
    save(CompiledParticlesFile, savedVariables{:},'-v7.3', '-nocompression');
end

CompiledParticlesToken = now;
save([DropboxFolder,filesep,Prefix,filesep,'CompiledParticlesToken.mat'],'CompiledParticlesToken', '-v6')

disp('CompiledParticles.mat saved.');

end