function [Particles, Spots, SpotFilter, schnitzcells] = CheckParticleTracking(varargin)
% function [Particles, Spots, SpotFilter, schnitzcells] = CheckParticleTracking(varargin)
%
% DESCRIPTION
% The point of this function is to check the tracking of particles. The
% logic of this function should be similar to schnitzcells: We want to be
% able to correct both the segmentation and tracking.
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
%
% OPTIONS
% NoSort : Flag to sort or not particles according to their starting frame
% sortByLength : Flag to have the particles with the most points (and thus
%       potentially most useful) tracked show up first
% ForCompileAll : Flag to just save the data. This is good for CompileAll
% speedmode : Flag to plot only ellipses for current particle & save time
% sistermode : Decide whether you want to do sister chromatid analysis
% nc, NC : Only look at particles that show up in nc13
%    % Currently this only starts at nc13...not restrict you to nc13 Added by Emma
%    % Also, this option shows you the max projection.
%
% CONTROLS
% Frame specific:
% . , Move a frame forward/backward
% > < Move five frames forward/backward
% ; ' Move to the next empty frame within a particle
% a z Move up/down in Z
% j Jump to a specified frame
% g b Increase/decrease histone channel contrast
% ! @ Change the contrast in transcription channel (! increases, @ resets
%       it back to the lowest)
%
% Particle specific:
% m Move to the next particle
% n Move to the previous particle
% k Jump to a specified particle by inputting particle number
% \ Jump to a specified particle by clicking
% c Connect two existing particle traces. This will join the current
%  particle's trace to the clicked particle's trace.
% d Separate traces forward. A new particle is created at the current frame
%  and this particle is disconnected from the current nucleus. If this is
%  done on a particle with only one frame then
%  it disconnects it from its nucleus.
%%
% q Cycle between approved status: green - approved; yellow - approved but
%  with conditions (drift of nucleus, for example)
% w Disapprove a trace
% p Identify a particle. It will also tell you the particle associated with
%  the clicked nucleus.
% e Approve/Disapprove a frame within a trace
% u Move a particle detected with Threshold2 into the our structure.
% i Move a particle detected with Threshold2 into the our structure and
%  connect it to the current particle. This is a combination of "u" and
%  "c". %AR 1/15/18: This is currently deprecated. Talk to HG if you need
%  this function.
% [ Add a spot that was not recognized originally by segmentSpots, creating
%    a new particle if you've used '+' to zoom in on the particle of
%    interest or adding to the current particle if you are zoomed out in
%    'o'. Note that the command forces ZoomMode. To toggle, use 'o' or '+'
%    depending on whether you're adding to an existing trace or creating
%    a new trace, respectively.
% { Same as [ but uses the exact pixel and z-plane that you click on.
%    Useful if the algorithms get the centroid positioning wrong (i.e. your
%    particle is put in the wrong place by '[').
% # remove a spot from Spots and erase its frame in Particles
%
%
%
% Nuclear tracking specific:
% l Modify a nuclear lineage and associate a particle with a nucleus.
%       Usage:
%       Click on one new nucleus + ENTER: Continue the schnitz with that nucleus.
%       Click on the current nucleus + ENTER: Split the schnitz. This time
%           point will be the first frame of the new schnitz.
%       Click on two nuclei: Split the current nucleus into two daughter
%       nuclei.
%       Click on the same nucleus twice: Split the current nucleus, but
%       with only one daughter nucleus.
% 2 set parent of current nucleus
% p Find the particle associated with the clicked nucleus. It will also tell
%  you the closest particle associated you clicked on.
%
%
% General:
% 8 Change channels
% t Show/hide particles from the second threshold
% s Save the current Particles structure
% x Save and exit
% h Show non-approved particles yellow or dissapproved particlesz
% y Input the frame/nc information again. This only works in the absence of
%  the histone channel
% r Reorder the particles according to initial frame
% f Redo tracking. It only gets done on the non-approved particles.
% o Zoom in/out around the particle's first frame.
% + Zoom anywhere button. Click with the mouse to specify the position to
%     to zoom in on after hitting this.
% -/= Change the zoom factor when in zoom mode.
% 0 Enter debug mode to fix things manually
% ~ Switch figure 1 from a single plane image to a z or time projection.
% 3 Fits a line to the polymerase loading regime of the trace.
% F Start a fitting mode for the single trace (in the current particle).
% (in progress)
%
% OUTPUT
% Particles: A modified Particles
% Spots: A modified Spots
% SpotFilter: A modified SpotFilter
% schnitzcells: A modified schnitzcells
%
% Author (contact): Hernan Garcia (hgarcia@berkeley.edu)
% Created:
% Last Updated: 1/13/2018

close all


warning('off','MATLAB:nargchk:deprecated')
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary')

%%Initialization
schnitzcells = [];
Ellipses = [];
correspondingNCInfo = [];
IntegrationArea=[]; %Initialized here to avoid dynamic assignment later in function

xForZoom = 0;
yForZoom = 0;

% Parameters for fitting
lineFit = 0;
fitApproved = 0;
FramesToFit = [];
FrameIndicesToFit = [];

%% Information about about folders

%Get the folders
[~,~,DefaultDropboxFolder,~,PreProcPath]=...
    DetermineLocalFolders;

if isempty(varargin)
    DataFolder=uigetdir(DefaultDropboxFolder,'Select data set to analyze');
else
    [~,~,DropboxFolder,~,~]=...
        DetermineLocalFolders(varargin{1});
    DataFolder=[DropboxFolder,filesep,varargin{1}];
end

[Prefix, Sort, sortByLength, ForCompileAll, SpeedMode, SisterMode, ...
    ncRange, projectionMode, plot3DGauss, intScale, NC, ...
    startNC, endNC] = determineCheckParticleTrackingOptions(varargin);

%%

FilePrefix=[DataFolder(length(DropboxFolder)+2:end),'_'];

%Now get the actual folders
[~,~,DropboxFolder,~,PreProcPath]=...
    DetermineLocalFolders(FilePrefix(1:end-1));

disp('Loading Particles.mat...')
load([DataFolder,filesep,'Particles.mat'], 'Particles', 'SpotFilter')
disp('Particles.mat loaded')
disp('Loading Spots.mat...')
load([DataFolder,filesep,'Spots.mat'], 'Spots')
disp('Spots.mat loaded')

%Check that FrameInfo exists
if exist([DataFolder,filesep,'FrameInfo.mat'], 'file')
    load([DataFolder,filesep,'FrameInfo.mat'], 'FrameInfo')
else
    warning('No FrameInfo.mat found. Trying to continue')
    %Adding frame information
    DHis=dir([PreProcPath,filesep,FilePrefix(1:end-1),filesep,'*His*.tif']);
    FrameInfo(length(DHis)).nc=[];
    %Adding information
    Dz=dir([PreProcPath,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'*001*.tif']);
    NumberSlices=length(Dz)-1;
    for i=1:numFrames
        FrameInfo(i).NumberSlices=NumberSlices;
    end
end
xSize = FrameInfo(1).PixelsPerLine;
ySize = FrameInfo(1).LinesPerFrame;
pixelSize = FrameInfo(1).PixelSize*1000; %nm
zStep = FrameInfo(1).ZStep;
snippet_size = 2*(floor(1300/(2*pixelSize))) + 1; % nm. note that this is forced to be odd
LinesPerFrame = FrameInfo(1).LinesPerFrame;
PixelsPerLine = FrameInfo(1).PixelsPerLine;
numFrames =length(FrameInfo);


%See how  many frames we have and adjust the index size of the files to
%load accordingly
if numFrames<1E3
    NDigits=3;
elseif numFrames<1E4
    NDigits=4;
else
    error('No more than 10,000 frames supported.')
end



%Some parameters:
SnippetEdge=13;     %Size of the snippets generated by Michael's code in pixels.
storedTimeProjection = []; % Don't need to wait for timeProjection to finish each time its called

%Create the particle array. This is done so that we can support multiple
%channels. Also figure out the number of channels
if iscell(Particles)
    NChannels=length(Particles);
else
    Particles={Particles};
    if ~iscell(Spots)
        Spots={Spots};
    end
    SpotFilter={SpotFilter};
    NChannels=1;
end


%Add FramesApproved where necessary
for NCh=1:NChannels
    if ~isfield(Particles{NCh},'FrameApproved')
        for i=1:length(Particles{NCh})
            Particles{NCh}(i).FrameApproved=true(size(Particles{NCh}(i).Frame));
        end
    else
        for i=1:length(Particles{NCh})
            if isempty(Particles{NCh}(i).FrameApproved)
                Particles{NCh}(i).FrameApproved=true(size(Particles{NCh}(i).Frame));
            end
        end
    end
end



%Check if we have the histone channel and we have done the nuclear
%segmentation.
if exist([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
        FilePrefix(1:end-1),'-His_',iIndex(1,NDigits),'.tif'], 'file')||...
        exist([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
        FilePrefix(1:end-1),'_His_',iIndex(1,NDigits),'.tif'], 'file')
    %(MT, 2018-02-11) Added support for lattice imaging with bad histone
    %channel, maybe temporary - FIX LATER
    if exist([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,'Ellipses.mat'], 'file')
        load([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,'Ellipses.mat'], 'Ellipses')
        UseHistoneOverlay=1;
    else
        warning('Ellipses.mat does not exist. Proceeding as though there is no Histone channel. If you expect a Histone channel, there is something wrong.')
        UseHistoneOverlay=0;
    end
else
    UseHistoneOverlay=0;
end

%Check that we have the nuclear tracking done using schnitzcells
if exist([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'], 'file')
    UseSchnitz=1;
else
    UseSchnitz=0;
end

% we name the variable DataFolderColumnValue to avoid shadowing previously defined DataFolder var, which is actually a subfolder inside dropbox
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
    Channel1, Channel2, Objective, Power, DataFolderColumnValue, DropboxFolderName, Comments,...
    nc9, nc10, nc11, nc12, nc13, nc14, CF,Channel3,prophase,metaphase] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);


if exist([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'], 'file')
    
    for i=1:numFrames
        if i<nc9
            FrameInfo(i).nc=8;
        elseif (i>=nc9)&(i<nc10)
            FrameInfo(i).nc=9;
        elseif (i>=nc10)&(i<nc11)
            FrameInfo(i).nc=10;
        elseif (i>=nc11)&(i<=nc12)
            FrameInfo(i).nc=11;
        elseif (i>=nc12)&(i<=nc13)
            FrameInfo(i).nc=12;
        elseif (i>=nc13)&(i<=nc14)
            FrameInfo(i).nc=13;
        elseif i>=nc14
            FrameInfo(i).nc=14;
        end
    end
else
    %     warning('No nuclear marker channel may result in strange behavior.');
    
    for i=1:numFrames
        if i<nc9
            FrameInfo(i).nc=8;
        elseif (i>=nc9)&(i<nc10)
            FrameInfo(i).nc=9;
        elseif (i>=nc10)&(i<nc11)
            FrameInfo(i).nc=10;
        elseif (i>=nc11)&(i<=nc12)
            FrameInfo(i).nc=11;
        elseif (i>=nc12)&(i<=nc13)
            FrameInfo(i).nc=12;
        elseif (i>=nc13)&(i<=nc14)
            FrameInfo(i).nc=13;
        elseif i>=nc14
            FrameInfo(i).nc=14;
        end
    end
end

% %Get the actual time corresponding to each frame
try
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
end
ElapsedTime=ElapsedTime/60;     %Time is in minutes
anaphase = [nc9,nc10,nc11,nc12,nc13,nc14];
anaphaseInMins = anaphase;
for i = 1:length(anaphase)
    if anaphase(i) > 0
        anaphaseInMins(i) = ElapsedTime(anaphase(i)); % in units of minutes
    end
end

%prophase and metaphase 
prophaseInMins=[]; metaphaseInMins=[];

try
    prophaseInMins = prophase;
    for i = 1:length(prophase)
        if prophase(i) > 0
            prophaseInMins(i) = ElapsedTime(prophase(i)); %mins
        end
    end
    metaphaseInMins = metaphase;
    for i = 1:length(metaphase)
        if metaphase(i) > 0
            metaphaseInMins(i) = ElapsedTime(metaphase(i)); %mins
        end
    end
end


try
    correspondingNCInfo = [FrameInfo.nc]; % the assigned nc of the frames
catch
end

save([DataFolder,filesep,'FrameInfo.mat'],'FrameInfo') %this is here so that a user will still get an updated
%frameinfo.mat even if they abort checkparticletracking without saving (to
%prevent issues with compileparticles)

%Check if we have already determined nc
if (~isfield(FrameInfo,'nc'))&&(~UseHistoneOverlay)
    %FrameInfo=DetermineNC(fad,Particles,FrameInfo);  AR 3/14/16: This
    %script seems to have disappeared.
    
elseif UseSchnitz
    load([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'])
    
    %Remove the schnitz fields that can give us problems potentially if
    %present. I don't know how this came to be, but it's for fields that
    %are not all that relevant. The fields are: approved, ang
    if isfield(schnitzcells,'approved')
        schnitzcells=rmfield(schnitzcells,'approved');
    end
    if isfield(schnitzcells,'ang')
        schnitzcells=rmfield(schnitzcells,'ang');
    end
end

% %Load the DoG images. Necessary for particle addition/subtraction
% dog = [];
% num_frames = numFrames;
% zSize = FrameInfo(1).NumberSlices + 2;
% OutputFolder1=[FISHPath,filesep,Prefix,'_',filesep,'dogs'];
% for current_frame = 1:num_frames
%     for i = 1:zSize
%         dog(:,:,i,current_frame) = double(imread([OutputFolder1, filesep,'DOG_',Prefix,'_',iIndex(current_frame,3),'_z',iIndex(i,2),'.tif']));
%     end
% end

%Order particles by the earliest frame they appear at. This makes the
%tracking a lot easier! Can also track by the number of spots in a trace
if Sort
    for ChN=1:NChannels
        nParticles = length(Particles{ChN});
        sortIndex = zeros(1, nParticles);
        for i=1:length(Particles{ChN})
            if sortByLength %sort by most points in particle
                sortIndex(i)=length(Particles{ChN}(i).Frame);
                direction = 'descend';
            else %Otherwise, sort by first frame as normal
                sortIndex(i)=Particles{ChN}(i).Frame(1);
                direction = 'ascend';
            end
        end
        [~,Permutations]=sort(sortIndex,direction);
        Particles{ChN}=Particles{ChN}(Permutations);
    end
end

%Some flags and initial parameters
ShowThreshold2=1;                    %Whether to show particles below the threshold
HideApprovedFlag=0;
ParticleToFollow=[];
ZSlices=FrameInfo(1).NumberSlices+2; %Note that the blank slices are included
CurrentZ=round(ZSlices/2);
ManualZFlag=0;
CurrentParticle=1;
PreviousParticle=1;
lastParticle = 0; %this gets flagged if there's a drop to one particle within the Particles structure.
CurrentFrameWithinParticle=1;
CurrentChannel=1;
PreviousChannel=CurrentChannel;
CurrentFrame=Particles{1}(1).Frame(1);
DisplayRange=[];
DisplayRangeSpot=[];
ZoomMode=0;
GlobalZoomMode=0;
ZoomRange=50;
nameSuffix='';

%Set up the default contrast settings for the MCP channel depending on the
%microscope that was used used
if strcmpi(FrameInfo(1).FileMode,'dspin')
    %For spinning disk, we set the contrast to the maximum and minimum
    minContrast=[];
    maxContrast=[];
else
    %For all other microscopes, we have a default. HG is not sure this will
    %actually work well beyond Leica SP8.
    minContrast = 0; % Default contrast settings for gfp channel
    maxContrast = 80;
end

% Changing the intial frames and particle if justNC13
if ncRange
    
    if strcmpi('nc15',endNC)
        lastNCFrame = numFrames;
    else
        lastNCFrame = eval(endNC)-1; % This will not include the 1st frame of the next NC
    end
    firstNCFrame = eval(startNC);
    particlesInRange = particlesWithinFrames(Prefix,firstNCFrame,lastNCFrame);
    CurrentParticle = particlesInRange(1);
    CurrentFrame = Particles{1}(CurrentParticle).Frame(1);
    %     ncRangeFigure = figure();
    %     set(gcf,'units', 'normalized', 'position',[0.35, 0.55, .2, .33])
    %     uicontrol('Parent', ncRangeFigure, 'Style', 'text','String','Implementing justNC13','Units','normalized', 'Position', [0.25 0.5 0.5 0.35])
    %     uicontrol('Parent', ncRangeFigure, 'Style', 'text','String',['nc13 : ' num2str(firstNCFrame)],'Units','normalized', 'Position', [0.25 0.40 0.5 0.35])
    %     uicontrol('Parent', ncRangeFigure, 'Style', 'text','String',['Number of Particles: ' num2str(length(particlesInRange))],'Units','normalized','Position', [0.25 0.30 0.5 0.35])
    %     uicontrol('Parent', ncRangeFigure, 'Style', 'text','String',['Particles in range: ' num2str(particlesInRange)],'Units','normalized','Position', [0.25 0.20 0.5 0.35])
    disp(['nc range: ' num2str(NC)])
    disp(['start frame: ' num2str(firstNCFrame)])
    disp(['end frame: ' num2str(lastNCFrame)])
    disp(['Particles in range: ' num2str(particlesInRange)])
    disp(['Number of Particles: ' num2str(length(particlesInRange))])
end

%Define the windows
Overlay=figure;

if UseHistoneOverlay
    HisOverlayFig=figure;
    HisOverlayFigAxes = axes(HisOverlayFig);
end

overlayAxes =subplot(1, 2, 1, 'Parent', Overlay);
traceFigAxes = subplot(1, 2, 2, 'Parent', Overlay);

zFig = figure;
zProfileFigAxes =subplot(1, 2, 1, 'Parent', zFig);
zTraceAxes  = subplot(1, 2, 2, 'Parent', zFig);

snipFig = figure();
snippetFigAxes =subplot(1, 3, 1, 'Parent', snipFig);
rawDataAxes =subplot(1, 3, 2, 'Parent', snipFig);
gaussianAxes =subplot(1, 3, 3, 'Parent', snipFig);

%Define the windows

set(Overlay,'units', 'normalized', 'position', [0.01, .45, .82, .33]);

if UseHistoneOverlay
    set(HisOverlayFig,'units', 'normalized', 'position',[0.01, 0.1, .33, .33]);
end

set(overlayAxes,'units', 'normalized', 'position', [-.25 .06 .9 .9])
set(traceFigAxes,'units', 'normalized', 'position', [.48 .17 .48 .63])
set(snipFig,'units', 'normalized', 'position',[0.355, 0.15, 3*(.2/2), .33/2]);
set(zFig,'units', 'normalized', 'position',[0.67, 0.15, .2, .33/2]);

%Define user interface
[controls, frame_num, z_num, particle_num, ...
    add_spot, smart_add_spot, delete_spot, ...
    fit_spot, averagingLength, approve_fit] = setupControls(Overlay);
set(0, 'CurrentFigure', Overlay);
import java.awt.Robot;
import java.awt.event.KeyEvent;
robot = Robot;
fake_event = KeyEvent.VK_T;
no_clicking = false;

frame_num.ValueChangedFcn = @frame_num_changed;
    function frame_num_changed(~, ~)

        figure(Overlay);
        numValidFrames = length({Spots{1}.Fits});
        [CurrentFrame, ManualZFlag] = changeFrame(str2double(frame_num.Value), numValidFrames);
        robot.keyPress(fake_event);
        robot.keyRelease(fake_event);
    end

z_num.ValueChangedFcn = @z_num_changed;
    function z_num_changed(~,~)

        figure(Overlay);
        [CurrentZ,ManualZFlag] = changeZSlice(str2double(z_num.Value), ZSlices);
        robot.keyPress(fake_event);
        robot.keyRelease(fake_event);

    end

particle_num.ValueChangedFcn = @particle_num_changed;
    function particle_num_changed(~, ~)

        figure(Overlay);
        [CurrentParticle,CurrentFrame, ManualZFlag] = changeParticle(...
            str2double(particle_num.Value), Particles, numParticles, CurrentChannel);
        robot.keyPress(fake_event);
        robot.keyRelease(fake_event);

    end

add_spot.ButtonPushedFcn = @add_spot_pushed;
    function add_spot_pushed(~,~)
        figure(Overlay);
        smart_add = '{';
        if smart_add_spot.Value
            smart_add = '[';
        end
        PathPart1 = [PreProcPath,filesep,FilePrefix(1:end-1),filesep,FilePrefix];
        PathPart2 = [nameSuffix,'.tif'];
        Path3 = [PreProcPath,filesep,Prefix,filesep,Prefix];
        no_clicking = true;
        [numParticles, SpotFilter, Particles, Spots, PreviousParticle] =...
            addSpot(ZoomMode, GlobalZoomMode, Particles, CurrentChannel, ...
            CurrentParticle, CurrentFrame, CurrentZ, Overlay, snippet_size, PixelsPerLine, ...
            LinesPerFrame, Spots, ZSlices, PathPart1, PathPart2, Path3, FrameInfo, pixelSize, ...
            SpotFilter, numParticles, smart_add, xSize, ySize, NDigits, intScale);

        robot.keyPress(fake_event);
        robot.keyRelease(fake_event);
        no_clicking = false;
    end

delete_spot.ButtonPushedFcn = @delete_spot_pushed;
    function delete_spot_pushed(~,~)
        no_clicking = true;
        figure(Overlay);
        [Spots, SpotFilter, ZoomMode, GlobalZoomMode, CurrentFrame, ...
            CurrentParticle, Particles, ManualZFlag, DisplayRange, lastParticle, PreviousParticle] =...
            removeSpot(ZoomMode, GlobalZoomMode, Frames, CurrentFrame, ...
            CurrentChannel, CurrentParticle, CurrentParticleIndex, Particles, Spots, SpotFilter, ...
            numParticles, ManualZFlag, DisplayRange, lastParticle, PreviousParticle);

        robot.keyPress(fake_event);
        robot.keyRelease(fake_event);
        no_clicking = false;
    end

% The part below is added by Yang Joon Kim, for single MS2 trace linear
% fitting (for the inital slope). Contact yjkim90@berkeley.edu for further
% discussion or improvement.
% Define the averaging window
averagingLength.ValueChangedFcn = @averagingLength_changed;
    function averagingLength_changed(~,~)
        averagingLength = str2double(averagingLength.Value);
    end

% Fit the initial slope, by clicking two points, you can define the window
% for fitting.
fit_spot.ButtonPushedFcn = @fit_spot_pushed;
    function fit_spot_pushed(~,~)
        %lineFit = 0;
        clear fit1E;
        figure(Overlay);
        
    [lineFit, Coefficients, fit1E, Particles, FramesToFit, FrameIndicesToFit] =...
        fitInitialSlope(CurrentParticle, Particles, Spots, CurrentChannel, schnitzcells, ...
        ElapsedTime, anaphaseInMins, correspondingNCInfo, traceFigAxes, Frames, anaphase, ...
        averagingLength, FramesToFit, FrameIndicesToFit, lineFit)
    end

approve_fit.ButtonPushedFcn = @fit_approve;
    function fit_approve(~,~)
        % Define the fitApproved as true
        fitApproved=1
        % save the fitted values (Slope and Time on) in Particles.mat
        if ~isempty(Coefficients)
            singleTraceLoadingRate = Coefficients(1,1); %au/min
            if singleTraceLoadingRate >= 0 %some easy quality control
                singleTraceTimeOn = roots(Coefficients(1,:));
                Particles{CurrentChannel}(CurrentParticle).fittedSlope =  singleTraceLoadingRate;
                Particles{CurrentChannel}(CurrentParticle).fittedTON =  singleTraceTimeOn;
                Particles{CurrentChannel}(CurrentParticle).fitApproved = 1;
            else     
                Particles{CurrentChannel}(CurrentParticle).fittedSlope =  NaN;
                Particles{CurrentChannel}(CurrentParticle).fittedTON =  NaN; 
                Particles{CurrentChannel}(CurrentParticle).fitApproved = 0;
            end
        else
            Particles{CurrentChannel}(CurrentParticle).fittedSlope =  NaN;
            Particles{CurrentChannel}(CurrentParticle).fittedTON =  NaN;   
            Particles{CurrentChannel}(CurrentParticle).fitApproved = 0;
        end
    end

cc=1;

% Create the approved field if it does not exist
for NCh=1:NChannels
    if ~isfield(Particles{NCh},'Approved')
        for i=1:length(Particles{NCh})
            Particles{NCh}(i).Approved=0;
        end
    end
end

%See if we just want to save the data
if ForCompileAll
    cc='x';
end


%Figure out channel-specific information
if NChannels==1
    if contains(Channel1{1}, 'MCP') || contains(Channel1{1}, 'PCP')
        nameSuffix=['_ch',iIndex(1,2)];
        coatChannel = 1;
    elseif contains(Channel2{1}, 'MCP') || contains(Channel2{1}, 'PCP')
        nameSuffix=['_ch',iIndex(2,2)];
        coatChannel = 2;
    end
elseif strcmpi(ExperimentType,'2spot2color')
    %We are assuming that channels 1 and 2 are assigned to coat
    %proteins. We should do a better job with this.
    coatChannels = [1,2];
    coatChannel=coatChannels(1);
else
    error('Experiment type not recognized')
end



%This flag allows the code to directly pass a command without waiting for
%the user to press a key or click on the figure
SkipWaitForButtonPress=[];
lineFit = 0; % the initial rise of the trace was not fitted

% these variables are only use when lineFit = 1
fit1E = []; 
Coefficients = []; 

while (cc~='x')
    
    %Update the name suffix
    if strcmpi(ExperimentType,'2spot2color')
        nameSuffix=['_ch',iIndex(coatChannel,2)];
    end
    
    numParticles = length(Particles{CurrentChannel});
    
    %Get the coordinates of all the spots in this frame
    [x,y,z]=SpotsXYZ(Spots{CurrentChannel}(CurrentFrame));
    
    %If the approved field does not exist create it
    if ~isfield(Particles{CurrentChannel},'Approved')
        for i=1:numParticles
            Particles{CurrentChannel}(i).Approved=0;
        end
    end
    
    %Pull out the right particle if it exists in this frame
    CurrentParticleIndex=...
        Particles{CurrentChannel}(CurrentParticle).Index(Particles{CurrentChannel}(CurrentParticle).Frame==...
        CurrentFrame);
    %This is the position of the current particle
    xTrace=x(CurrentParticleIndex);
    yTrace=y(CurrentParticleIndex);
    
    if (~isempty(xTrace))&&(~ManualZFlag)
        CurrentZ=z(CurrentParticleIndex);
        CurrentZIndex=find(...
            Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).z==...
            CurrentZ);
        ManualZFlag=0;
    end
    
    if NChannels==1 % inputoutput mode can also be in this case, changed CurrentChannel to the coatChannel (YJK : 1/15/2018)
        if strcmpi(projectionMode,'None (Default)')
            Image=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
                FilePrefix,iIndex(CurrentFrame,NDigits),'_z',iIndex(CurrentZ,2),nameSuffix,'.tif']);
        elseif strcmpi(projectionMode,'Max Z')
            [Image,~] = zProjections(Prefix, coatChannel, CurrentFrame, ZSlices, NDigits,DropboxFolder,PreProcPath);
        elseif strcmpi(projectionMode,'Median Z')
            [~,Image] = zProjections(Prefix, coatChannel, CurrentFrame, ZSlices, NDigits,DropboxFolder,PreProcPath);
        elseif strcmpi(projectionMode,'Max Z and Time')
            if isempty(storedTimeProjection)
                if ncRange
                    Image = timeProjection(Prefix, coatChannel,'nc',NC);
                    storedTimeProjection = Image;
                else
                    Image = timeProjection(Prefix, CurrentChannel);
                    storedTimeProjection = Image;
                end
            else
                Image = storedTimeProjection;
            end
        end
        
    elseif NChannels>1
        Image=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
            FilePrefix,iIndex(CurrentFrame,NDigits),'_z',iIndex(CurrentZ,2),...
            nameSuffix,'.tif']);
    else
        error('ExperimentType and/or channel not supported.')
    end
    
    set(0, 'CurrentFigure', Overlay);
    imshow(Image,DisplayRangeSpot,'Border','Tight','Parent',overlayAxes, ...
        'InitialMagnification', 'fit')
    hold(overlayAxes,'on')
    
    if UseHistoneOverlay
        HisPath1 = [PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
            FilePrefix(1:end-1),'-His_',iIndex(CurrentFrame,NDigits),'.tif'];
        HisPath2 = [PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
            FilePrefix(1:end-1),'_His_',iIndex(CurrentFrame,NDigits),'.tif'];
        [ImageHis, xForZoom, yForZoom] = plotFrame(overlayAxes, Image, SpeedMode, FrameInfo, Particles, ...
            Spots, CurrentFrame, ShowThreshold2, ...
            Overlay, CurrentChannel, CurrentParticle, ZSlices, CurrentZ, numFrames, ...
            schnitzcells, UseSchnitz, DisplayRange, Ellipses, SpotFilter, ZoomMode, GlobalZoomMode, ...
            ZoomRange,xForZoom,yForZoom,UseHistoneOverlay, HisOverlayFigAxes, HisPath1, HisPath2);

    else
        plotFrame(overlayAxes, Image, SpeedMode, ...
            FrameInfo, Particles, Spots, CurrentFrame, ShowThreshold2, ...
            Overlay, CurrentChannel, CurrentParticle, ZSlices, CurrentZ, numFrames, ...
            schnitzcells, UseSchnitz, DisplayRange, Ellipses, SpotFilter, ...
            ZoomMode, GlobalZoomMode, ZoomRange,xForZoom,yForZoom,UseHistoneOverlay);
    end
    
    if ~isempty(xTrace)
        MaxZIndex=find(...
            Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).z==...
            Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).brightestZ);
        CurrentZIndex=find(...
            Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).z==...
            CurrentZ);
        if isempty(CurrentZIndex)
            %             warning('This particle has a gap in its z-profile. This is
            %             highly suspect.'); %this if statement should only happen
            %             between two spots, not past the PSF boundaries
        end
    end
    
    %Check to see if spots structure contains multi-slice fields
    multi_slice_flag = isfield(Spots{CurrentChannel}(CurrentFrame).Fits...
        (CurrentParticleIndex),'IntegralZ');
    
    % PLOTS SNIPPET
    
    FullSlicePath = [PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(CurrentFrame,3)...
        ,'_z' iIndex(CurrentZ,2) '_ch' iIndex(coatChannel,2) '.tif'];
    if exist('CurrentSnippet', 'var')
        CurrentSnippet = plotSnippet(snippetFigAxes, rawDataAxes, gaussianAxes, xTrace, ...
            CurrentZIndex, FullSlicePath, Spots, CurrentChannel, CurrentFrame, ...
            CurrentParticleIndex, ExperimentType, intScale, snippet_size, xSize, ...
            ySize, SnippetEdge, CurrentSnippet);
    else
        CurrentSnippet = plotSnippet(snippetFigAxes, rawDataAxes, gaussianAxes, xTrace, ...
            CurrentZIndex, FullSlicePath, Spots, CurrentChannel, CurrentFrame, ...
            CurrentParticleIndex, ExperimentType, intScale, snippet_size, xSize, ...
            ySize, SnippetEdge);
    end
    
    % PLOTS TRACE OF CURRENT PARTICLE
    if exist('AmpIntegral', 'var')
        [Frames,AmpIntegral,GaussIntegral,AmpIntegral3,AmpIntegral5, ...
            ErrorIntegral, ErrorIntegral3, ErrorIntegral5,backGround3, ...
            AmpIntegralGauss3D, ErrorIntegralGauss3D, PreviousParticle] = plotTrace(traceFigAxes, ...
            FrameInfo, CurrentChannel, PreviousChannel, ...
        CurrentParticle, PreviousParticle, lastParticle, HideApprovedFlag, lineFit, anaphaseInMins, ...
        ElapsedTime, schnitzcells, Particles, plot3DGauss, anaphase, prophase, metaphase,prophaseInMins, metaphaseInMins,Prefix, ...
        DefaultDropboxFolder, numFrames, CurrentFrame, ZSlices, CurrentZ, Spots, ...
        correspondingNCInfo, fit1E, Coefficients, ExperimentType, Frames, AmpIntegral, ...
        GaussIntegral, AmpIntegral3, AmpIntegral5, ErrorIntegral, ErrorIntegral3, ...
        ErrorIntegral5, backGround3, AmpIntegralGauss3D, ErrorIntegralGauss3D);
    else
        [Frames,AmpIntegral,GaussIntegral,AmpIntegral3,AmpIntegral5, ...
            ErrorIntegral, ErrorIntegral3, ErrorIntegral5,backGround3, ...
            AmpIntegralGauss3D, ErrorIntegralGauss3D, PreviousParticle] = plotTrace(traceFigAxes, ...
            FrameInfo, CurrentChannel, PreviousChannel, ...
        CurrentParticle, PreviousParticle, lastParticle, HideApprovedFlag, lineFit, anaphaseInMins, ...
        ElapsedTime, schnitzcells, Particles, plot3DGauss, anaphase,prophase, metaphase, prophaseInMins, metaphaseInMins,Prefix, ...
        DefaultDropboxFolder, numFrames, CurrentFrame, ZSlices, CurrentZ, Spots, ...
        correspondingNCInfo, fit1E, Coefficients, ExperimentType);
    end 

    

    % PLOT Z SLICE RELATED FIGURES
    if exist('MaxZProfile', 'var')
        [MaxZProfile, Frames] = plotZFigures(zProfileFigAxes, zTraceAxes, ExperimentType, ...
            xTrace, Spots, CurrentFrame, CurrentChannel, CurrentParticleIndex, ZSlices, ...
            CurrentZ, CurrentZIndex, PreviousParticle, CurrentParticle, ...
            PreviousChannel, Particles, Frames, MaxZProfile);
    else
        [MaxZProfile, Frames] = plotZFigures(zProfileFigAxes, zTraceAxes, ExperimentType, ...
            xTrace, Spots, CurrentFrame, CurrentChannel, CurrentParticleIndex, ZSlices, ...
            CurrentZ, CurrentZIndex, PreviousParticle, CurrentParticle, ...
            PreviousChannel, Particles, Frames);
    end
    
    % UPDATE UICONTROLS
    updateControls(frame_num, z_num, particle_num, CurrentFrame, CurrentZ, CurrentParticle);
    
    set(0, 'CurrentFigure', Overlay);
    if isempty(SkipWaitForButtonPress)
        ct = waitforbuttonpress; % ct=0 for click and ct=1 for keypress
        cc = get(Overlay, 'CurrentCharacter');
        cm2 = get(overlayAxes, 'CurrentPoint');
           
        current_axes = get(Overlay, 'CurrentAxes');
        if strcmpi(cc, '') || ct == 0
            cc = 'donothing';
        end
        is_control = isa(get(Overlay, 'CurrentObject'), 'matlab.ui.control.UIControl');
        if ct == 0 && cm2(1,1) < xSize && current_axes == overlayAxes...
                && ~no_clicking && ~is_control
            
            [CurrentParticle, CurrentFrame, ManualZFlag] = toNearestParticle(Spots, ...
            Particles, CurrentFrame, CurrentChannel, UseHistoneOverlay, ...
            schnitzcells, [cm2(1,1), cm2(2,2)]);
        end
    else
        cc=SkipWaitForButtonPress;
        SkipWaitForButtonPress=[];
    end
    
    numValidFrames = length({Spots{1}.Fits});
    
    if strcmpi(cc, 'donothing')
        %do nothing
    elseif cc=='.' %Move forward one frame
        [CurrentFrame, ManualZFlag] = changeFrame(CurrentFrame + 1, numValidFrames);
    elseif (cc==',') %Move backward one frame
        [CurrentFrame, ManualZFlag] = changeFrame(CurrentFrame - 1, numValidFrames);
    elseif (cc=='>') %Move forward five frames
        [CurrentFrame, ManualZFlag] = changeFrame(CurrentFrame + 5, numValidFrames);
    elseif (cc=='<') %#ok<*AND2> %Move backward five frames
        [CurrentFrame, ManualZFlag] = changeFrame(CurrentFrame - 5, numValidFrames);
    elseif (cc=='''')&(CurrentFrame<numValidFrames) %Move to the next skipped frame
        %within the particle
        CurrentFrame = nextSkippedFrame(Particles, CurrentChannel, ... 
            CurrentParticle, CurrentFrame);        
    elseif (cc==';')&(CurrentFrame>1) %Move to the previous skipped frame
        %within the particle
        CurrentFrame = previousSkippedFrame(Particles, CurrentChannel, ... 
            CurrentParticle, CurrentFrame);
    elseif (cc=='a')&(CurrentZ<ZSlices) %Move up in Z
        [CurrentZ,ManualZFlag] = changeZSlice(CurrentZ + 1, ZSlices);
    elseif (cc=='z')&(CurrentZ>1) %Move down in Z
        [CurrentZ,ManualZFlag] = changeZSlice(CurrentZ - 1, ZSlices);
    elseif cc=='j'
        try
            iJump= inputdlg('Frame to jump to:',...
                'Move to frame');
            iJump=str2double(iJump{1});
        catch
            iJump=CurrentFrame;
        end
        [CurrentFrame, ManualZFlag] = changeFrame(iJump, numValidFrames);
        DisplayRange=[];
    elseif cc=='k'
        try
            ParticleJump = inputdlg('Particle to jump to:',...
                'Move to particle');
            ParticleJump=str2double(ParticleJump{1});
        catch
            ParticleJump=CurrentParticle;
        end
        [CurrentParticle,CurrentFrame, ManualZFlag] = ...
            changeParticle(ParticleJump, Particles, numParticles, CurrentChannel);  
        DisplayRange=[];
    elseif cc=='g' & UseHistoneOverlay     %Increase histone channel contrast
        if isempty(DisplayRange)
            DisplayRange=[min(min(ImageHis)),max(max(ImageHis))/1.5];
        else
            DisplayRange=[DisplayRange(1),DisplayRange(2)/1.5];
        end
        
    elseif cc=='b' & UseHistoneOverlay     %Decrease histone channel contrast
        DisplayRange=[min(min(ImageHis)),max(max(ImageHis))*1.5];
        
    elseif cc=='#' %remove a spot from Spots and erase its frame in Particles
        [Spots, SpotFilter, ZoomMode, GlobalZoomMode, CurrentFrame, ...
            CurrentParticle, Particles, ManualZFlag, DisplayRange, lastParticle, PreviousParticle] =...
            removeSpot(ZoomMode, GlobalZoomMode, Frames, CurrentFrame, ...
            CurrentChannel, CurrentParticle, CurrentParticleIndex, Particles, Spots, SpotFilter, ...
            numParticles, ManualZFlag, DisplayRange, lastParticle, PreviousParticle);
    elseif cc=='[' | cc=='{' %#ok<*OR2> %Add particle and all of its shadows to Spots.
        PathPart1 = [PreProcPath,filesep,FilePrefix(1:end-1),filesep,FilePrefix];
        PathPart2 = [nameSuffix,'.tif'];
        Path3 = [PreProcPath,filesep,Prefix,filesep,Prefix];
        [numParticles, SpotFilter, Particles, Spots,PreviousParticle] =...
            addSpot(ZoomMode, GlobalZoomMode, Particles, CurrentChannel, ...
            CurrentParticle, CurrentFrame, CurrentZ, Overlay, snippet_size, PixelsPerLine, ...
            LinesPerFrame, Spots, ZSlices, PathPart1, PathPart2, Path3, FrameInfo, pixelSize, ...
            SpotFilter, numParticles, cc, xSize, ySize, NDigits, intScale);
    elseif cc=='r'
        Particles = orderParticles(numParticles, CurrentChannel, Particles);       
    elseif cc=='f'
        [Particles, schnitzcells] = redoTracking(DataFolder, ...
            UseHistoneOverlay, FrameInfo, DropboxFolder, FilePrefix, schnitzcells, ...
            Particles, NChannels, CurrentChannel, numParticles);
    elseif cc=='c'
        [PreviousParticle, Particles] = combineTraces(Spots, ...
            CurrentChannel, CurrentFrame, Particles, CurrentParticle);
    elseif cc=='p' %Identify a particle. It will also tell you the particle associated with
        %  the clicked nucleus.
        identifyParticle(Spots, Particles, CurrentFrame, ...
            CurrentChannel, UseHistoneOverlay, schnitzcells);       
    elseif cc=='\' %Moves to clicked particle.
        [CurrentParticle, CurrentFrame, ManualZFlag] = toNearestParticle(Spots, ...
            Particles, CurrentFrame, CurrentChannel, UseHistoneOverlay, schnitzcells);
    elseif cc=='u'
        [x2,y2,z2]=SpotsXYZ(Spots{CurrentChannel}(CurrentFrame));
        if ~isempty(x2)
            ClickedSpot=ginput(1);
            
            UnfilterSpot(Spots{CurrentChannel},SpotFilter{CurrentChannel},...
                ClickedSpot,Particles{CurrentChannel},CurrentFrame)
        end
    elseif cc=='i'
        warning(' AR 1/15/18: This is currently deprecated. Talk to HG if you need this function.')
    elseif cc=='d' || cc=='v'  %d Separate traces forward at the current frame.
        [Particles, PreviousParticle] = separateTraces(Particles, ...
            CurrentChannel, CurrentFrame, CurrentParticle);
    elseif cc=='q'      %Approve a trace
        if Particles{CurrentChannel}(CurrentParticle).Approved==1
            Particles{CurrentChannel}(CurrentParticle).Approved=2;
        elseif Particles{CurrentChannel}(CurrentParticle).Approved==0
            Particles{CurrentChannel}(CurrentParticle).Approved=1;
        elseif Particles{CurrentChannel}(CurrentParticle).Approved==2
            Particles{CurrentChannel}(CurrentParticle).Approved=0;
        end
    elseif cc=='w'      %Disapprove a trace
        if Particles{CurrentChannel}(CurrentParticle).Approved==-1
            Particles{CurrentChannel}(CurrentParticle).Approved=0;
        else
            Particles{CurrentChannel}(CurrentParticle).Approved=-1;
        end
        
    elseif cc=='s'  
        saveChanges(NChannels, Particles, Spots, SpotFilter, DataFolder, ...
            FrameInfo, UseHistoneOverlay, FilePrefix, ...
            schnitzcells, DropboxFolder);      
    elseif cc=='t'
        ShowThreshold2=~ShowThreshold2;
        
        %     elseif (cc=='y')&(~UseHistoneOverlay)
        %             FrameInfo=DetermineNC(fad,Particles{CurrentChannel},FrameInfo);
        %
        %AR 9/5/18- this button is deprecated. leaving this comment in
        %case we want to replace the functionality.
        %
    elseif cc=='h'
        if HideApprovedFlag==0
            HideApprovedFlag=1;         %Show only non-approved traces
        elseif HideApprovedFlag==1
            HideApprovedFlag=2;         %Show only yellow and red traces
        elseif HideApprovedFlag==2
            HideApprovedFlag=0;
        end
        
        %HideApprovedFlag=~HideApprovedFlag;
    elseif cc=='o'
        if ~GlobalZoomMode
            ZoomMode=~ZoomMode;
        else
            disp('Try again after exiting global zoom mode by hitting ''+''')
        end
        
    elseif cc=='+'
        if ~ZoomMode
            if ~GlobalZoomMode
                %AR 12/14/2018 ginputc doesn't work at the moment. not sure
                %why. maybe related to the fact that this is an image
                %within a figure that has multiple axes. 
%                 [ConnectPositionx,ConnectPositiony]=ginputc(1,'color', 'r', 'linewidth',1);
                [ConnectPositionx,ConnectPositiony]=ginput(1);
                xForZoom = round(ConnectPositionx);
                yForZoom = round(ConnectPositiony);
            else
            end
            GlobalZoomMode =~ GlobalZoomMode;
        else
            disp('Try again after exiting zoom mode by hitting ''o''')
        end
        
    elseif (cc=='m')&(CurrentParticle<numParticles)
        [lineFit, CurrentParticle, CurrentFrame, ManualZFlag, DisplayRange] =...
            goNextParticle(CurrentParticle, CurrentChannel, HideApprovedFlag, Particles);
    elseif (cc=='n')&(CurrentParticle>1)
        [lineFit, CurrentParticle, CurrentFrame, ManualZFlag, DisplayRange] =...
            goPreviousParticle(CurrentParticle, CurrentChannel, HideApprovedFlag, Particles);    
    elseif cc=='e'
        Particles{CurrentChannel}(CurrentParticle).FrameApproved(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame)=...
            ~Particles{CurrentChannel}(CurrentParticle).FrameApproved(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame);

        %Schnitzcells specific
        
    elseif cc=='l' %Split a nucleus and select one or two daughter cells or stop the lineage
        [Particles, PreviousParticle, schnitzcells] = splitNuclei(schnitzcells, ...
            CurrentFrame, CurrentChannel, CurrentParticle, Particles);    
    elseif cc=='2' %2 set parent of current nucleus
        schnitzcells = setParentNucleus(schnitzcells, ...
            CurrentFrame, CurrentChannel, CurrentParticle, Particles);
    elseif cc=='8' && NChannels > 1      %Switch channels
        [CurrentChannel, PreviousChannel, coatChannel, CurrentParticle] =...
            switchChannels(CurrentChannel, CurrentParticle, Particles, ...
            UseHistoneOverlay, coatChannels, NChannels);       
    elseif cc=='~'      %Switch projection mode
        projectionMode = chooseProjection;
        disp(['projectionMode : ' projectionMode])
        
    elseif cc=='!' %Increase contrast in the Overlay figure
        if isempty(DisplayRangeSpot)
            DisplayRangeSpot=[min(min(Image)),max(max(Image))/1.5];
        else
            DisplayRangeSpot=[DisplayRangeSpot(1),DisplayRangeSpot(2)/1.5];
        end
    elseif cc=='@'      %Decrease spot channel contrast
        DisplayRangeSpot=[min(min(Image)),max(max(Image))*1.5];
    elseif cc=='0'      %Debugging mode
        keyboard;
        
    elseif cc== '3'
        [lineFit, Coefficients, fit1E, Particles] =...
            fitLine(CurrentParticle, Particles, Spots, CurrentChannel, schnitzcells, ...
            ElapsedTime, anaphaseInMins, correspondingNCInfo, traceFigAxes, Frames, anaphase);
    end
end


save([DataFolder,filesep,'FrameInfo.mat'],'FrameInfo')

%If we only have one channel bring Particles back to the legacy
%format without any cells
if NChannels==1
    Particles=Particles{1};
    Spots=Spots{1};
    SpotFilter=SpotFilter{1};
end


if UseHistoneOverlay
    save([DataFolder,filesep,'Particles.mat'],'Particles','SpotFilter', '-v7.3')
    save([DataFolder,filesep,'Spots.mat'],'Spots', '-v7.3')
    save([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'],'schnitzcells')
else
    save([DataFolder,filesep,'Particles.mat'],'Particles','SpotFilter', '-v7.3')
    save([DataFolder,filesep,'Spots.mat'],'Spots','-v7.3')
end

close all
if ishandle(controls)
    close(controls)
end

disp('Particles saved.')
disp(['(Left off at Particle #', num2str(CurrentParticle), ')'])

%% Extra stuff that is useful in debug mode

%Reset approve status of all approved particles in a certain nc
% nc=13;
%
% %Determine the start and end frame of the nc
% if nc==14
%     disp('Do this')
% else
%     eval(['ncStart=nc',num2str(nc),';']);
%     eval(['ncEnd=nc',num2str(nc+1),';']);
% end
%
% for i=1:length(Particles)
%     if Particles(i).Approved==1
%
%         if (min(Particles(i).Frame(Particles(i).FrameApproved))>=ncStart)&...
%                 (min(Particles(i).Frame(Particles(i).FrameApproved))<ncEnd)
%             Particles(i).Approved=0;
%         end
%     end
% end


%This is if the reference from schnitz to nucleus is screwed up. It
%manifests itself in that cenx and ceny of the schnitz will not coincide
%with the actual ellipse. In this case we delete the schnitz (emptying it)
%and the nuclear reference and use '1' to recreate the schnitz.

% ParticlesToDecouple=[314:316,318,320:325,328,329,331:333,335,337,339,341,342,344:346]
%
% for ParticleToDecouple=ParticlesToDecouple
%
%     SchnitzToDecouple=Particles(ParticleToDecouple).Nucleus;
%
%     %Delete the schnitz by emptying it
%     schnitzcells(SchnitzToDecouple).P=[];
%     schnitzcells(SchnitzToDecouple).E=[];
%     schnitzcells(SchnitzToDecouple).D=[];
%     schnitzcells(SchnitzToDecouple).frames=[];
%     schnitzcells(SchnitzToDecouple).cenx=[];
%     schnitzcells(SchnitzToDecouple).ceny=[];
%     schnitzcells(SchnitzToDecouple).len=[];
%     schnitzcells(SchnitzToDecouple).cellno=[];
%
%     %Decouple particle and schnitz
%     Particles(ParticleToDecouple).Nucleus=[];
% end


end
%
%     function plotButtonPushed(src,event)
%         bar(randn(1,5));
%     end
