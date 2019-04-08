function CompileParticles(varargin)
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
% 'intArea': Change the area (in pixels) of integration used in offset calculations
% 'noHist': Force the code to assume there's no nuclear channel.
% 'doSingleFits': Generate single trace fits. Added by EL&AR 
% 'manualSingleFits' : Compile values from manually generated single trace fits to
%               CompiledParticles.mat, the field names are 'fittedSlope', 'fittedTon' for
%               the initial slope and T_on respectively.
% 'optionalResults' : if you want to use a different dropbox folder

% Author (contact): Hernan Garcia (hggarcia@berkeley.edu)
% Created:
% Last Updated: 6/17/17 (AR)
%
% Documented by: Hernan Garcia (hggarcia@berkeley.edu)

close all;

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
EllipsePos = {};
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

nc9 = [];
nc10 = [];
nc11 = [];
nc12 = [];
nc13 = [];
nc14 = [];
ncFilter = [];
ncFilterID = [];
%%

%Information about about folders
[~,~,DefaultDropboxFolder,~,~]=...
    DetermineLocalFolders;

[Prefix, ForceAP, SkipTraces, SkipFluctuations, SkipFits, SkipMovie, ...
    SkipAll, ApproveAll, MinParticles, minTime, ROI, intArea, noHist, ...
    ROI1, ROI2, slimVersion, manualSingleFits, optionalResults, yToManualAlignmentPrompt] = determineCompileParticlesOptions(varargin);

FilePrefix=[Prefix,'_'];

%What type of experiment are we dealing with? Get this out of MovieDatabase
[~,~,DropboxFolder,~, PreProcPath,...
    ~, ~, ~, ~, ~,~] = readMovieDatabase(Prefix, optionalResults);

% refactor in progress, we should replace readMovieDatabase with getExperimentDataFromMovieDatabase
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoopEnd, APResolution,...
    Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
    nc9, nc10, nc11, nc12, nc13, nc14, CF,Channel3,prophase,metaphase, anaphase] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);


%Load all the information
disp('Loading Particles.mat...');
load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat']);
disp('Particles loaded.');
disp('Loading Spots.mat...');
load([DropboxFolder,filesep,Prefix,filesep,'Spots.mat']);
disp('Spots loaded.');
if isempty(Particles)
    SkipTraces=1;
    SkipFluctuations=1;
    SkipFits=1;
    SkipMovie=1;
end

%Check that FrameInfo exists
if exist([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'file')
    load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo');
    pixelSize = FrameInfo(1).PixelSize;
else
    warning('No FrameInfo.mat found. Trying to continue')
    %Adding frame information
    DHis=dir([PreProcPath,filesep,FilePrefix(1:end-1),filesep,'*His*.tif']);
    if ~isempty(DHis)
        FrameInfo(length(DHis)).nc=[];
    end
    
    %Adding information
    
    Dz=dir([PreProcPath,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'*001*.tif']);
    NumberSlices=length(Dz)-1;
    
    for i=1:numFrames
        FrameInfo(i).NumberSlices=NumberSlices;
    end
end

numFrames = length(FrameInfo);



%See how  many frames we have and adjust the index size of the files to
%load accordingly
if numFrames<1E3
    NDigits=3;
elseif numFrames<1E4
    NDigits=4;
else
    error('No more than 10,000 frames currently supported.')
end

%Create the particle array. This is done so that we can support multiple
%channels. Also figure out the number of channels
if iscell(Particles)
    NChannels=length(Particles);
else
    Particles={Particles};
    Spots={Spots};
    SpotFilter={SpotFilter};
    NChannels=1;
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
    load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])
    HistoneChannel=1;
else
    disp('No lineage / nuclear information found. Proceeding without it.');
    HistoneChannel=0;
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


%Add the APPosition to Particles if they don't exist yet. Do this only if
%we took AP data. Otherwise just add x and y pixel coordinates

addParticleArgs = {Prefix};
if ~isempty(optionalResults)
    addParticleArgs = [addParticleArgs, optionalResults];
end
if yToManualAlignmentPrompt
    addParticleArgs = [addParticleArgs, 'yToManualAlignmentPrompt'];
end

if strcmpi(ExperimentAxis,'AP')
    if (~isfield(Particles{1},'APpos')) || ForceAP
        if ~HistoneChannel
            addParticleArgs = [addParticleArgs, 'SkipAlignment'];
        end
        AddParticlePosition(addParticleArgs{:});
    else
        disp('Using saved AP information (results from AddParticlePosition)')
    end
elseif strcmpi(ExperimentAxis,'dv') && exist([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'],'file')
    AddParticlePosition(Prefix);
elseif strcmpi(ExperimentAxis,'dv') || strcmpi(ExperimentAxis,'NoAP')
    AddParticlePosition(Prefix,'NoAP');
else
    error('Experiment axis not recognized in MovieDatabase')
end

%Reload Particles.mat
load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
%Create the particle array. This is done so that we can support multiple
%channels. Also figure out the number of channels
if iscell(Particles)
    NChannels=length(Particles);
else
    Particles={Particles};
    NChannels=1;
end


if HistoneChannel
    load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'])
end


%Folders for reports
warning('off', 'MATLAB:MKDIR:DirectoryExists');
if strcmpi(ExperimentAxis,'AP')
    mkdir([DropboxFolder,filesep,Prefix,filesep,'APMovie'])
end
mkdir([DropboxFolder,filesep,Prefix,filesep,'ParticleTraces'])
mkdir([DropboxFolder,filesep,Prefix,filesep,'TracesFluctuations'])
mkdir([DropboxFolder,filesep,Prefix,filesep,'Offset'])
mkdir([DropboxFolder,filesep,Prefix,filesep,'Fits'])
mkdir([DropboxFolder,filesep,Prefix,filesep,'Probabilities'])
mkdir([DropboxFolder,filesep,Prefix,filesep,'Various']);


%% Put together CompiledParticles

CompiledParticles = cell(NChannels,1);
%Approve all particles if the mode has been selected
if ApproveAll
    for ChN=1:NChannels
        %Check that the approved field is present. If not include
        %it. This can occur if CheckParticleTracking is not run
        %first.
        if ~isfield(Particles{ChN},'Approved')
            for i=1:length(Particles{ChN})
                Particles{ChN}(i).Approved=0;
            end
        end
        
        for i=1:length(Particles{ChN})
            %Make sure the particle has an associated nucleus if we are in
            %HistoneChannel mode
            if HistoneChannel
                if ~isempty(Particles{ChN}(i).Nucleus)
                    %If a particle has been explicitly rejected then don't
                    %approve it!
                    if Particles{ChN}(i).Approved~=-1
                        Particles{ChN}(i).Approved=1;
                    end
                end
            else
                Particles{ChN}(i).Approved=1;
            end
        end
    end
end

%First, figure out the AP position of each of the nuclei.
if strcmpi(ExperimentAxis, 'AP') || strcmpi(ExperimentAxis, 'DV')
    %Load the AP detection information
    load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])
    %Angle between the x-axis and the AP-axis
    if exist('coordPZoom', 'var')
        APAngle=atan2((coordPZoom(2)-coordAZoom(2)),(coordPZoom(1)-coordAZoom(1)));
    else
        error('coordPZoom not defined. Was AddParticlePosition.m run?')
    end
    APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);
end

if HistoneChannel && (strcmpi(ExperimentAxis,'AP') || strcmpi(ExperimentAxis,'DV'))
    %The information in Ellipses is
    %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
    for i=1:length(Ellipses)
        for j=1:size(Ellipses{i},1)
            
            %Angle between the x-axis and the particle using the A position as a
            %zero

            Angles=atan2((Ellipses{i}(j,2)-coordAZoom(2)),(Ellipses{i}(j,1)-coordAZoom(1)));
            
            %Distance between the points and the A point
            Distances=sqrt((coordAZoom(2)-Ellipses{i}(j,2)).^2+(coordAZoom(1)-Ellipses{i}(j,1)).^2);
            APPositions=Distances.*cos(Angles-APAngle);
            EllipsePos{i}(j)=APPositions/APLength;
            
            DVPositions=Distances.*sin(Angles-APAngle);
            EllipsePos_DV{i}(j)=DVPositions;
        end
    end
end


%Get the actual time corresponding to each frame
if isfield(FrameInfo,'FileMode')
    if strcmp(FrameInfo(end).FileMode,'TIF')
        for j=1:numFrames
            ElapsedTime(j)=etime(datevec(FrameInfo(j).TimeString),datevec(FrameInfo(1).TimeString));
        end
    elseif strcmp(FrameInfo(end).FileMode,'LSM')||strcmp(FrameInfo(end).FileMode,'LSMExport')||...
            strcmp(FrameInfo(end).FileMode,'LIFExport')||strcmp(FrameInfo(end).FileMode,'LAT')
        for j=1:numFrames
            ElapsedTime(j)=FrameInfo(j).Time-FrameInfo(1).Time;
        end
    else
        error('File mode not supported. Cannot extract time information. Include format in ExportDataForLivemRNA.m')
    end
else
    warning('No FileMode information found. Assuming that this is TIF from the 2-photon.')
    for j=1:numFrames
        ElapsedTime(j)=etime(datevec(FrameInfo(j).TimeString),datevec(FrameInfo(1).TimeString));
    end
end

ElapsedTime=ElapsedTime/60;     %Time is in minutes

%Divide the AP and DV axes into bins for generating means, etc. 
if strcmpi(ExperimentAxis,'AP') || strcmpi(ExperimentAxis,'DV')  
    [APbinID, APbinArea] = binAPAxis(APResolution, FrameInfo, ...
        coordAZoom, APAngle, APLength);
end
if strcmpi(ExperimentAxis,'DV')
    [DVbinID, DVbinArea] = binDVAxis(FrameInfo, coordAZoom, APAngle);
end    
    

%Now get the particle information for those that were approved
[Particles, CompiledParticles, ncFilter, ncFilterID] =...
    compileTraces(NChannels, Particles, HistoneChannel, ...
    schnitzcells, minTime, ExperimentAxis, APbinID, APbinArea, CompiledParticles, ...
    Spots, SkipTraces, nc9, nc10, nc11, nc12, nc13, nc14, ncFilterID, ncFilter, ...
    ElapsedTime, intArea, Ellipses, EllipsePos, PreProcPath, ...
    FilePrefix, Prefix, DropboxFolder, NDigits, manualSingleFits);

%% ROI option
% This option is separating the CompiledParticles defined above into
% CompiledParticles_ROI and COmpiledParticles_nonROI
% written by YJK on 10/24/2017
CompiledParticles_ROI = {}; CompiledParticles_nonROI = {};
for ChN=1:NChannels
    if ROI
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
%% Create filters

% APFilter needs to be changed for ROI option, so that we can have two
% filters for each CompiledParticles (ROI and nonROI) (YJK on 10/24/2017)

%nc filters:
try
    
APFilter_ROI = []; APFilter_nonROI = []; DVFilter_ROI = []; DVFilter_nonROI = [];
[ncFilterID, ncFilter, APFilter, APFilter_ROI, APFilter_nonROI, ...
    DVFilter, DVFilter_ROI, DVFilter_nonROI] = createFilters(nc9, nc10, nc11, nc12, ...
    nc13, nc14, NChannels, CompiledParticles, ExperimentAxis, ROI, APFilter,...
    APFilter_ROI, APFilter_nonROI, APbinID, DVbinID, DVFilter, DVFilter_ROI, ...
    DVFilter_nonROI, CompiledParticles_ROI, CompiledParticles_nonROI);
end

%% Binning and averaging data

% Here I need to think about how to bin ROI and non-ROI particles.
% One way is having a function that splits the CompiledParticles into CompiledParticles_ROI
% and CompiledParticles_nonROI according to its y-position (Threshold). Then, pass them through
% the AverageTraces (YJK on 10/22/2017)
try
[AllTracesVector, AllTracesAP, AllTracesDV, MeanVectorAP_ROI, ...
    SDVectorAP_ROI, NParticlesAP_ROI, MeanVectorAP_nonROI, SDVectorAP_nonROI, ...
    NParticlesAP_nonROI, MeanVectorAP, SDVectorAP, NParticlesAP, MeanVectorDV_ROI, ...
    SDVectorDV_ROI, NParticlesDV_ROI, MeanVectorDV_nonROI, SDVectorDV_nonROI, ...
    NParticlesDV_nonROI, MeanVectorDV, SDVectorDV, NParticlesDV, ...
    MeanVectorAnterior, MeanVectorAll, SDVectorAll, NParticlesAll, MaxFrame] =...
    computeAPandDVStatistics(NChannels, CompiledParticles, FrameInfo, ExperimentAxis, ...
    APFilter, ROI, CompiledParticles_ROI, CompiledParticles_nonROI, ...
    APFilter_ROI, APFilter_nonROI, NewCyclePos, MaxFrame, ...
    MeanVectorAP_ROI, SDVectorAP_ROI, NParticlesAP_ROI, MeanVectorAP_nonROI, ...
    SDVectorAP_nonROI, NParticlesAP_nonROI, MeanVectorAP, SDVectorAP, ...
    NParticlesAP, MeanVectorDV_ROI, SDVectorDV_ROI, NParticlesDV_ROI, ...
    MeanVectorDV_nonROI, SDVectorDV_nonROI, NParticlesDV_nonROI, ...
    MeanVectorDV, SDVectorDV, NParticlesDV, MeanVectorAnterior, DVFilter_ROI, ...
    DVFilter_nonROI, DVFilter);
end
if ~slimVersion
    %% Instantaneous rate of change
    try
    [CompiledParticles, MeanSlopeVectorAP, SDSlopeVectorAP, NSlopeAP]...
        = instantRateOfChange(NChannels, CompiledParticles, ElapsedTime, ...
        ExperimentAxis, APFilter, MeanVectorAP, MeanSlopeVectorAP, ...
        SDSlopeVectorAP, NSlopeAP);
    end
    %% Integrating each particle
try
    CompiledParticles = integrateParticles(NChannels, ElapsedTime, CompiledParticles);
end

    %% Information about the cytoplasm
    %If the nuclear masks are present then use them. Otherwise just calculate
    %the median of the images as a function of time

    %HG on 8/6/16: Why was this commented out? Did I do this?


    % if NChannels==1
    %
    %     if HistoneChannel&strcmpi(ExperimentAxis,'AP')
    %         [MeanCyto,SDCyto,MedianCyto,MaxCyto]=GetCytoMCP(Prefix);
    %     else
    %         MeanCyto=[];
    %         SDCyto=[];
    %         MaxCyto=[];
    %
    %         h=waitbar(0,'Calculating the median cyto intentisy');
    %         for i=1:numFrames
    %             waitbar(i/numFrames,h)
    %             for j=1:FrameInfo(1).NumberSlices
    %                 Image(:,:,j)=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(i,3),'_z',iIndex(j,2),'.tif']);
    %             end
    %             ImageMax=max(Image,[],3);
    %             MedianCyto(i)=median(double(ImageMax(:)));
    %         end
    %         close(h)
    %     end
    % else
    MeanCyto=[];
    SDCyto=[];
    MaxCyto=[];
    MedianCyto=[];
    % end




    %% Offset and fluctuations
try
    [MeanOffsetVector, SDOffsetVector, NOffsetParticles] = offsetAndFlux(NChannels, ...
        SkipFluctuations, ncFilter, ElapsedTime, CompiledParticles, DropboxFolder, ...
        Prefix, ExperimentAxis, intArea, MeanVectorAll, SDVectorAll, MaxFrame, numFrames);

    %% Rate of mRNA production

    CompiledParticles = mRNAProdRate(NChannels, CompiledParticles, ...
        ncFilter, ElapsedTime, SkipFits, DropboxFolder, Prefix);


    %% First frames

    firstFramesThing(NChannels, HistoneChannel, SkipAll, nc13, nc14, ...
        CompiledParticles, DropboxFolder, Prefix, ElapsedTime, ExperimentAxis);

    %% AP position of particle vs nucleus

    if HistoneChannel&&strcmpi(ExperimentAxis,'AP')
        CompiledParticles = APPosParticleVsNucleus(NChannels, ...
            CompiledParticles, schnitzcells, EllipsePos, DropboxFolder, Prefix);
    end


    %% Fitting shapes` to single traces (includes time on and initial rate of loading)
    % This section of code is will fit line segments piece-wise to the
    % single traces. fittedLineEquations correspond to the stored fitted lines
    % of the particles, where the indexing is as follows:
    % fittedLineEquations(particleNumber) with fields: Coefficients,
    % ErrorEstimation, and frameIndex, which are described below. This currently
    % does not support more than one channel. Please contact Emma to work on
    % implementing it for two channels.
    if ~SkipFits
        [CompiledParticles, fittedLineEquations] = fitShapesToTraces(Prefix, ...
            Particles, schnitzcells, FrameInfo, ElapsedTime, CompiledParticles,Spots);
    end

    %% Probability of being on

    if HistoneChannel&&strcmpi(ExperimentAxis,'AP') || strcmpi(ExperimentAxis,'DV')
        [NEllipsesAP, MeanVectorAllAP, SEVectorAllAP, EllipsesFilteredPos, ...
            FilteredParticlesPos, OnRatioAP, ParticleCountAP, ParticleCountProbAP, ...
            EllipsesOnAP, rateOnAP, rateOnAPCell, timeOnOnAP, timeOnOnAPCell, TotalEllipsesAP, rateOnAPManual, rateOnAPCellManual, timeOnOnAPManual, timeOnOnAPCellManual]...
            = computeAPFractionOn(NChannels, Particles, schnitzcells, ...
            CompiledParticles, Ellipses, APbinID, FrameInfo, ElapsedTime, DropboxFolder, ...
            Prefix, EllipsePos, nc12, nc13, nc14, numFrames, SkipFits, SkipAll, ...
            APbinArea,pixelSize, manualSingleFits);
    end

    % DV version. This should instead be smoothly integrated with the AP
    % version since there's a lot of duplicate code here. 
    if HistoneChannel&&strcmpi(ExperimentAxis,'DV') %JAKE: Need to change this later
        [NEllipsesDV, MeanVectorAllDV, SEVectorAllDV, OnRatioDV, ParticleCountDV, ...
            ParticleCountProbDV, TotalEllipsesDV, EllipsesOnDV, EllipsesFilteredPos, ...
            FilteredParticlesPos] = DVProbOn(NChannels, Particles, schnitzcells, ...
            CompiledParticles, Ellipses, FrameInfo, DropboxFolder, Prefix, ...
            ElapsedTime, DVbinID, EllipsePos_DV, nc12, nc13, nc14, numFrames, ...
            DVbinArea);
    end
end
    %% Calculation of particle speed
    try
        calcParticleSpeeds(NChannels, Particles, ...
            Spots, ElapsedTime, schnitzcells, Ellipses);
    catch
    end
    %% Movie of AP profile

    %I want to make a movie of the average fluorescence as a function of AP as
    %a function of time. In order to make life easier I'll just export to a
    %folder. I can then load everything in ImageJ.

    if ~SkipMovie&&strcmpi(ExperimentAxis,'AP')
        APProfileMovie(MeanVectorAP, NParticlesAP, MinParticles, ...
            APbinID, SDVectorAP, FrameInfo, ElapsedTime, DropboxFolder, Prefix, nc9, nc10, nc11, nc12, nc13, nc14)
    end

    %% Checking correlations with position and expression
    %
    % 1+1;
    %
    % i=100
    %
    % for j=1:length(CompiledParticles{1}(i).Frame)
    %
    %     CurrentFrame=CompiledParticles{1}(i).Frame(j);
    %     EllipsePos=find((schnitzcells(CompiledParticles{1}(i).Nucleus).frames)==CurrentFrame);
    %     CurrentEllipse=Ellipses{CompiledParticles{1}(i).Frame};
    %     CurrentEllipse=CurrentEllipse(EllipsePos,:);
    %
    %     Position(j)=(CompiledParticles{1}(i).xPos(j)-CurrentEllipse(1))^2+...
    %         (CompiledParticles{1}(i).yPos(j)-CurrentEllipse(2))^2;
    % end
    %
    % figure(9)
    % plot(  Position ,CompiledParticles{1}(i).Fluo,'.k')
    %
    %
end

%% Input-output

%Compile the nuclear fluorescence information if we have the appropriate
%experiment type and axis
%Simon: This script is not required for Dl-Venus experiments in DV so I
%added the ExperimentAxis condition.
% ROI option added (YJK on 10/27/2017) :
% When the data is acquired in ROI mode, I added the ROI option,
% two thresholds to be incorporated into CompileNuclearProtein
if strcmpi(ExperimentType,'inputoutput') 
    %&& (strcmpi(ExperimentAxis,'AP')||strcmpi(ExperimentAxis,'DV'))
    if ~ROI
        CompileNuclearProtein(Prefix)
    else
        CompileNuclearProtein(Prefix,'ROI',ROI1,ROI2)
    end
end



%% Save everything

%Now save all the information

savedVariables = [savedVariables,'APFilter', 'APbinArea', 'APbinID', 'AllTracesAP',...
    'AllTracesVector', 'CompiledParticles', 'DVFilter', 'DVbinArea',...
    'DVbinID', 'ElapsedTime', 'EllipsePos', 'EllipsesFilteredPos',...
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
    'ncFilterID', 'rateOnAPManual', 'rateOnAPCellManual', 'timeOnOnAPManual', 'timeOnOnAPCellManual'];

save([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'],...
    savedVariables{:},'-v7.3');
 
end