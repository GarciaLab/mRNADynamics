function [Particles, Spots, SpotFilter, schnitzcells] = CheckParticleTracking(Prefix, varargin)
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
% plot3DGauss: plot 3D gaussian fitting intensities in tracefig
% sortByLength: sort particles by duration instead of time observed
%
% CONTROLS
% Frame specific:
% ./, :Move a frame forward/backward
% >/< :Move five frames forward/backward
% ;/' :Move to the next/previous empty frame within a particle
% a/z :Move up/down in z-slice
%  t  :Jump to a specific z-slice
%  j  :Jump to a specified frame
% g/b :Increase/decrease histone channel contrast
% !/@ :Change the contrast in transcription channel (! increases/@ resets
%       it back to the lowest)
%
% Particle specific:
% m :Move to the next particle
% n :Move to the previous particle
% k :Jump to a specified particle by inputting particle number
% \ :Jump to a specified particle by clicking
% c :Connect two existing particle traces. This will join the current
%     particle's trace to the clicked particle's trace.
% d :Separate traces forward. A new particle is created at the current frame
%     and this particle is disconnected from the current nucleus. If this is
%     done on a particle with only one frame then it disconnects
%     it from its nucleus.
% q :Cycle between approved status: green - approved; yellow - approved but
%     with conditions (drift of nucleus, for example)
% w :Disapprove a trace
% p :Identify a particle. It will also tell you the particle associated with
%     the clicked nucleus.
% e :Approve/Disapprove a frame within a trace
% u :Move a particle detected with Threshold2 into the our structure.
% i :Move a particle detected with Threshold2 into the our structure and
%     connect it to the current particle. This is a combination of "u" and
%     "c". %AR 1/15/18: This is currently deprecated. Talk to HG if you need
%     this function.
% [ :Add a spot that was not recognized originally by segmentSpots, creating
%     a new particle if you've used '+' to zoom in on the particle of
%     interest or adding to the current particle if you are zoomed out in
%     'o'. Note that the command forces ZoomMode. To toggle, use 'o' or '+'
%     depending on whether you're adding to an existing trace or creating
%     a new trace, respectively.
% { :Same as "[" but uses the exact pixel and z-plane that you click on.
%     Useful if the algorithms get the centroid positioning wrong (i.e. your
%     particle is put in the wrong place by '[').
% # :Remove a spot from Spots and erase its frame in Particles
%
%
% Nuclear tracking specific:
% C :Add a nucleus to the specified location. Re-run nuclear tracking
%       after if this is done.
% V :Remove a nucleus from the specified location. Re-run nuclear tracking
%       after if this is done.
% l :Modify a nuclear lineage and associate a particle with a nucleus.
%       Usage:
%       Click on one new nucleus + ENTER: Continue the schnitz with that nucleus.
%       Click on the current nucleus + ENTER: Split the schnitz. This time
%           point will be the first frame of the new schnitz.
%       Click on two nuclei: Split the current nucleus into two daughter
%       nuclei.
%       Click on the same nucleus twice: Split the current nucleus, but
%       with only one daughter nucleus.
% 2 :Set parent of current nucleus
% p :Find the particle associated with the clicked nucleus. It will also tell
%     you the closest particle associated you clicked on.
% $ :Add particle to the nearest nucleus
%
%
% General:
% 8 :Change channels
% t :Show/hide particles from the second threshold
% s :Save the current Particles structure
% x :Save and exit
% h :Show non-approved particles yellow or dissapproved particles
% y :Input the frame/nc information again. This only works in the absence of
%     the histone channel
% r :Reorder the particles according to initial frame
% f :Redo tracking. It only gets done on the non-approved particles.
% o :Zoom in/out around the particle's first frame.
% + :Zoom anywhere button. Click with the mouse to specify the position to
%     zoom in on after hitting this.
% -/= :Change the zoom factor when in zoom mode.
% 0 :Enter debug mode to fix things manually
% ~ :Switch figure 1 from a single plane image to a z or time projection.

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

warning('off', 'MATLAB:nargchk:deprecated')
warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary')

%% Initialization
schnitzcells = [];
Ellipses = [];
correspondingNCInfo = [];
IntegrationArea = []; %Initialized here to avoid dynamic assignment later in function
fish = false;
nucleiModified = false;

xForZoom = 0;
yForZoom = 0;

% Parameters for fitting
lineFitted = 0; % equals 1 if a line has been fitted
fitApproved = 0;
FramesToFit = []; % actual frames of the movie that were used for fitting
FrameIndicesToFit = 0; % index of the current particle that were used for fitting
Coefficients = []; % coefficients of the fitted line

SnippetEdge = 13; %Size of the snippets generated by Michael's code in pixels.
storedTimeProjection = []; % Don't need to wait for timeProjection to finish each time its called

[sortByFrames, sortByLength, ForCompileAll, SpeedMode, ~, ...
    ncRange, projectionMode, plot3DGauss, NC, ...
    startNC, endNC, optionalResults, nWorkers, fish, noHisOverlay] = determineCheckParticleTrackingOptions(varargin{:});

if fish
    noHisOverlay = true;
    projectionMode = 'Max Z';
end

%% Information about about folders

% Get the folders
[~,ProcPath,DropboxFolder,~, PreProcPath,~, Prefix, ~,~,~,~,~, ~, ~, movieDatabase]...
    = readMovieDatabase(Prefix, optionalResults);

DataFolder = [DropboxFolder, filesep, Prefix];
FilePrefix = [Prefix, '_'];


[Particles, SpotFilter, Spots, FrameInfo, Spots3D] = loadCheckParticleTrackingMats(DataFolder, PreProcPath);
[xSize, ySize, pixelSize, zStep, snippet_size, LinesPerFrame, PixelsPerLine,...
    numFrames] = getFrameInfoParams(FrameInfo);

  %See how  many frames we have and adjust the index size of the files to load accordingly
  if numFrames < 1E3
    NDigits = 3;
  elseif numFrames < 1E4
    NDigits = 4;
  else
    error('No more than 10,000 frames supported.')
  end
  %Create the particle array. This is done so that we can support multiple
  %channels. Also figure out the number of channels
  if iscell(Particles)
      NChannels = length(Particles);
  else
      Particles = {Particles};
      if ~iscell(Spots)
          Spots = {Spots};
      end
      SpotFilter = {SpotFilter};
      NChannels = 1;
  end

%Add FramesApproved where necessary
Particles = addFrameApproved(NChannels, Particles);

[Ellipses, UseHistoneOverlay, UseSchnitz] = checkHistoneAndNuclearSegmentation(PreProcPath, FilePrefix, NDigits, DropboxFolder, noHisOverlay);

if fish
    UseSchnitz = false;
end

[~, ExperimentType, ~, ~, ~, ~, Channel1, Channel2, ~, ~, ~, ~, ~, ...
    nc9, nc10, nc11, nc12, nc13, nc14, ~, Channel3, prophase, metaphase] =...
    getExperimentDataFromMovieDatabase(Prefix, movieDatabase);

Channels = {Channel1{1}, Channel2{1}, Channel3{1}};

for i = 1:numFrames
    if i < nc9
        FrameInfo(i).nc = 8;
    elseif (i >= nc9) & (i < nc10)
        FrameInfo(i).nc = 9;
    elseif (i >= nc10) & (i < nc11)
        FrameInfo(i).nc = 10;
    elseif (i >= nc11) & (i <= nc12)
        FrameInfo(i).nc = 11;
    elseif (i >= nc12) & (i <= nc13)
        FrameInfo(i).nc = 12;
    elseif (i >= nc13) & (i <= nc14)
        FrameInfo(i).nc = 13;
    elseif i >= nc14
        FrameInfo(i).nc = 14;
    end
end

%Get the actual time corresponding to each frame, in minutes
ElapsedTime = getFrameElapsedTime(FrameInfo, numFrames);

anaphase = [nc9, nc10, nc11, nc12, nc13, nc14];
[anaphaseInMins, prophaseInMins, metaphaseInMins] = getPhasesDurationsInMins(anaphase, prophase, metaphase, ElapsedTime);

try
    correspondingNCInfo = [FrameInfo.nc]; % the assigned nc of the frames
end

% this save() is here so that a user will still get an updated
% frameinfo.mat even if they abort checkparticletracking without saving
% (to prevent issues with compileparticles)
save([DataFolder, filesep, 'FrameInfo.mat'], 'FrameInfo');

schnitzPath = [DropboxFolder, filesep, FilePrefix(1:end - 1), filesep, FilePrefix(1:end - 1), '_lin.mat'];

if  exist(schnitzPath, 'file')
    
    disp('Loading schnitzcells...')
    load(schnitzPath, 'schnitzcells');
    disp('schnitzcells loaded.')
    
    %Remove the schnitz fields that can give us problems potentially if
    %present. I don't know how this came to be, but it's for fields that
    %are not all that relevant. The fields are: approved, ang
    if isfield(schnitzcells, 'approved')
        schnitzcells = rmfield(schnitzcells, 'approved');
    end
    
    if isfield(schnitzcells, 'ang')
        schnitzcells = rmfield(schnitzcells, 'ang');
    end
    
end

[Particles] = sortParticles(sortByFrames, sortByLength, NChannels, Particles);
        
%Some flags and initial parameters
ShowThreshold2 = 1; %Whether to show particles below the threshold
HideApprovedFlag = 0;
ParticleToFollow = [];
PreviousParticle = 1;
lastParticle = 0; %this gets flagged if there's a drop to one particle within the Particles structure.
CurrentFrameWithinParticle = 1;

cptState = CPTState(Spots, Particles, 0, 0, FrameInfo(1).NumberSlices, 1, 1); %work in progress, 2019-12, JP.

if ~isempty(cptState.Particles{cptState.CurrentChannel})
    cptState.CurrentFrame = cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Frame(CurrentFrameWithinParticle);
else
    error('Looks like the Particles structure is empty. There''s nothing to check.');
end
DisplayRange = [];
DisplayRangeSpot = [];
ZoomMode = 0;
GlobalZoomMode = 0;
ZoomRange = 50;
nameSuffix = '';
hImage = [];
CurrentSnippet = [];
oim = [];
hIm = [];
Image = [];
ellipseHandles = {};
spotHandles = {};
ellipseHisHandles = {};


Frames = [];

% Changing the intial frames and particle if justNC13
if ncRange
    
    if strcmpi('nc15', endNC)
        lastNCFrame = numFrames;
    else
        lastNCFrame = eval(endNC) - 1; % This will not include the 1st frame of the next NC
    end
    
    firstNCFrame = eval(startNC);
    particlesInRange = particlesWithinFrames(Prefix, firstNCFrame, lastNCFrame);
    cptState.CurrentParticle = particlesInRange(1);
    cptState.CurrentFrame = cptState.Particles{1}(cptState.CurrentParticle).Frame(1);
    disp(['nc range: ' num2str(NC)])
    disp(['start frame: ' num2str(firstNCFrame)])
    disp(['end frame: ' num2str(lastNCFrame)])
    disp(['Particles in range: ' num2str(particlesInRange)])
    disp(['Number of Particles: ' num2str(length(particlesInRange))])
end

%Define user interface
[Overlay, overlayAxes, snippetFigAxes, rawDataAxes, gaussianAxes, traceFigAxes, zProfileFigAxes,...
    zTraceAxes,HisOverlayFig,HisOverlayFigAxes] = checkParticleTracking_drawGUI(UseHistoneOverlay, fish);

[controls, frame_num, z_num, particle_num, ...
    add_spot, smart_add_spot, delete_spot, ...
    fit_spot, averagingLength, approve_fit] = setupControls(Overlay);

set(0, 'CurrentFigure', Overlay);
import java.awt.Robot;
import java.awt.event.KeyEvent;
robot = Robot;
fake_event = KeyEvent.VK_T;
no_clicking = false;

coatChannels = [1, 2]; % JP temporary, will be used only if 2spot2color, could be refactored into cptState

[frameChangeTextInput, frameChangeKeyInput] = FrameChangeEventHandler(cptState, robot, fake_event);
frame_num.ValueChangedFcn = frameChangeTextInput;

[zSliceChangeTextInput, zSliceChangeKeyInput] = ZSliceChangeEventHandler(cptState, robot, fake_event);
z_num.ValueChangedFcn = zSliceChangeTextInput;

[particleChangeTextInput, particleChangeKeyInput] = ParticleChangeEventHandler(cptState, robot, fake_event);
particle_num.ValueChangedFcn = particleChangeTextInput;

channelSwitchKeyInput = ChannelSwitchEventHandler(cptState, NChannels, coatChannels, UseHistoneOverlay);

add_spot.ButtonPushedFcn = @add_spot_pushed;

    function add_spot_pushed(~, ~)
        figure(Overlay);
        smart_add = '{';
        
        if smart_add_spot.Value
            smart_add = '[';
        end
        PathPart1 = [PreProcPath, filesep, FilePrefix(1:end - 1), filesep, FilePrefix];
        PathPart2 = [nameSuffix, '.tif'];
        Path3 = [PreProcPath, filesep, Prefix, filesep, Prefix];
        no_clicking = true;
        [SpotFilter, cptState.Particles, cptState.Spots, PreviousParticle, cptState.CurrentParticle] = ...
            addSpot(ZoomMode, GlobalZoomMode, cptState.Particles, cptState.CurrentChannel, ...
            cptState.CurrentParticle, cptState.CurrentFrame, cptState.CurrentZ, Overlay, snippet_size, PixelsPerLine, ...
            LinesPerFrame, cptState.Spots, cptState.ZSlices, PathPart1, PathPart2, Path3, FrameInfo, pixelSize, ...
            SpotFilter,smart_add, xSize, ySize, NDigits, ...
            cptState.coatChannel, UseHistoneOverlay, schnitzcells, nWorkers, plot3DGauss);
        robot.keyPress(fake_event);
        robot.keyRelease(fake_event);
        no_clicking = false;
    end

delete_spot.ButtonPushedFcn = @delete_spot_pushed;

    function delete_spot_pushed(~, ~)
        no_clicking = true;
        figure(Overlay);
        [cptState.Spots, SpotFilter, ZoomMode, GlobalZoomMode, cptState.CurrentFrame, ...
            cptState.CurrentParticle, cptState.Particles, cptState.ManualZFlag, DisplayRange, lastParticle, PreviousParticle] = ...
            removeSpot(ZoomMode, GlobalZoomMode, Frames, cptState.CurrentFrame, ...
            cptState.CurrentChannel, cptState.CurrentParticle, CurrentParticleIndex, cptState.Particles, cptState.Spots, SpotFilter, ...
             cptState.ManualZFlag, DisplayRange, lastParticle, PreviousParticle);
        
        robot.keyPress(fake_event);
        robot.keyRelease(fake_event);
        no_clicking = false;
    end

% The part below is added by Yang Joon Kim, for single MS2 trace linear
% fitting (for the inital slope). Contact yjkim90@berkeley.edu for further
% discussion or improvement.
% Define the averaging window
AveragingLength = 1; % Default
averagingLength.ValueChangedFcn = @averagingLength_changed;

    function averagingLength_changed(~, ~)
        AveragingLength = str2double(averagingLength.Value);
    end

% Fit the initial slope, by clicking two points, you can define the window
% for fitting.
fit_spot.ButtonPushedFcn = @fit_spot_pushed;

    function fit_spot_pushed(~, ~)
        %lineFit = 0;
        %clear lineFitHandle;
        figure(Overlay);
        
        [lineFitted, Coefficients, FramesToFit, FrameIndicesToFit] = ...
            fitInitialSlope(cptState.CurrentParticle, cptState.Particles, cptState.Spots, cptState.CurrentChannel, schnitzcells, ...
            ElapsedTime, anaphaseInMins, correspondingNCInfo, traceFigAxes, Frames, anaphase, ...
            AveragingLength, FramesToFit, FrameIndicesToFit, lineFitted);
    end

approve_fit.ButtonPushedFcn = @fit_approve;

    function fit_approve(~, ~)
        % Define the fitApproved as true
        fitApproved = 1;
        % save the fitted values (Coefficients and lineFitHandle) in Particles.mat
        % For now, I will save lineFitHandle (which is a plot for the initial
        % slope), but it might take up too much space, then I need a better
        % way to regenerate the plot, using the Coefficients and fittedFrames
        % Quality control : Check whether the slope is positive
        if Coefficients(1, 1) > 0
            cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).fitApproved = 1;
            cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Coefficients = Coefficients;
            %cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).lineFitHandle =  lineFitHandle;
            cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).fittedFrames = FrameIndicesToFit; % use the index of particle trace for convenience
            %cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).fittedYSegment = currentYSegment;
        else
            cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).fitApproved = 0;
            cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Coefficients = [];
            %cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).lineFitHandle =  [];
            cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).fittedFrames = [];
            %cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).fittedYSegment = [];
        end
        
    end

% Create the approved field if it does not exist
for NCh = 1:NChannels
    
    if ~isfield(cptState.Particles{NCh}, 'Approved')
        
        for i = 1:length(cptState.Particles{NCh})
            cptState.Particles{NCh}(i).Approved = 0;
        end
        
    end
    
end

cc = 1;

if ForCompileAll
    %we just want to save the data, we set 'x' as user command
    cc = 'x';
end

%This flag allows the code to directly pass a command without waiting for
%the user to press a key or click on the figure
SkipWaitForButtonPress = [];

while (cc ~= 'x')
    
    %% Main loop - start
    %%
    %Figure out channel-specific information
    if NChannels == 1
        
        for ch = 1:length(Channels)
            if contains(Channels{ch}, 'MCP') || contains(Channels{ch}, 'PCP') || contains(Channels{ch}, 'spot') 
                nameSuffix = ['_ch', iIndex(ch, 2)];
                cptState.coatChannel = ch;
            end
        end

    elseif strcmpi(ExperimentType, '2spot2color')
        %We are assuming that channels 1 and 2 are assigned to coat
        %proteins. We should do a better job with this.
        cptState.coatChannel = coatChannels(cptState.CurrentChannel);
    else
        error('Experiment type not recognized')
    end

    %Update the name suffix
    if strcmpi(ExperimentType, '2spot2color')
        nameSuffix = ['_ch', iIndex(cptState.coatChannel, 2)];
    end
    
    %Get the coordinates of all the spots in this frame
    [x, y, z] = SpotsXYZ(cptState.Spots{cptState.CurrentChannel}(cptState.CurrentFrame));
    
    %If the approved field does not exist create it
    if ~isfield(cptState.Particles{cptState.CurrentChannel}, 'Approved')
        
        for i = 1:cptState.numParticles()
            cptState.Particles{cptState.CurrentChannel}(i).Approved = 0;
        end
        
    end
    
    %Pull out the right particle if it exists in this frame
    CurrentParticleIndex = ...
        cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Index(cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Frame == ...
        cptState.CurrentFrame);
    %This is the position of the current particle
    xTrace = x(CurrentParticleIndex);
    yTrace = y(CurrentParticleIndex);
    
    if (~isempty(xTrace)) && (~cptState.ManualZFlag)
        cptState.CurrentZ = z(CurrentParticleIndex);
        CurrentZIndex = find(...
            cptState.Spots{cptState.CurrentChannel}(cptState.CurrentFrame).Fits(CurrentParticleIndex).z == ...
            cptState.CurrentZ);
        cptState.ManualZFlag = 0;
    end
    
    
    if strcmpi(projectionMode, 'None')
        Image = imread([PreProcPath, filesep, FilePrefix(1:end - 1), filesep, ...
            FilePrefix, iIndex(cptState.CurrentFrame, NDigits), '_z', iIndex(cptState.CurrentZ, 2), nameSuffix, '.tif']);
    elseif strcmpi(projectionMode, 'Max Z')
        Image = zProjections(Prefix, cptState.coatChannel, cptState.CurrentFrame, cptState.ZSlices, NDigits, DropboxFolder, PreProcPath, FrameInfo, 'max', nWorkers);
    elseif strcmpi(projectionMode, 'Median Z')
        Image = zProjections(Prefix, cptState.coatChannel, cptState.CurrentFrame, cptState.ZSlices, NDigits, DropboxFolder, PreProcPath, FrameInfo, 'median', nWorkers);
    elseif strcmpi(projectionMode, 'Max Z and Time')
        
        if isempty(storedTimeProjection)
            
            if ncRange
                Image = timeProjection(Prefix, cptState.coatChannel,FrameInfo, DropboxFolder,PreProcPath, 'nc', NC);
                storedTimeProjection = Image;
            else
                Image = timeProjection(Prefix, cptState.CurrentChannel,FrameInfo, DropboxFolder,PreProcPath);
                storedTimeProjection = Image;
            end
            
        else
            Image = storedTimeProjection;
        end
        
    end
    
    
    
    set(0, 'CurrentFigure', Overlay);
%     if isempty(hIm)
        hIm = imshow(Image, DisplayRangeSpot, 'Border', 'Tight', 'Parent', overlayAxes, ...
            'InitialMagnification', 'fit');
%     else
%         hIm.CData = Image;
%     end
%     
    hold(overlayAxes, 'on')
    
    if UseHistoneOverlay
        HisPath1 = [PreProcPath, filesep, FilePrefix(1:end - 1), filesep, ...
            FilePrefix(1:end - 1), '-His_', iIndex(cptState.CurrentFrame, NDigits), '.tif'];
        HisPath2 = [PreProcPath, filesep, FilePrefix(1:end - 1), filesep, ...
            FilePrefix(1:end - 1), '_His_', iIndex(cptState.CurrentFrame, NDigits), '.tif'];
        
        [ImageHis, xForZoom, yForZoom, oim,ellipseHandles]= displayOverlays(overlayAxes, Image, SpeedMode, FrameInfo, cptState.Particles, ...
            cptState.Spots, cptState.CurrentFrame, ShowThreshold2, ...
            Overlay, cptState.CurrentChannel, cptState.CurrentParticle, cptState.ZSlices, cptState.CurrentZ, numFrames, ...
            schnitzcells, UseSchnitz, DisplayRange, Ellipses, SpotFilter, ZoomMode, GlobalZoomMode, ...
            ZoomRange, xForZoom, yForZoom, fish, UseHistoneOverlay, HisOverlayFigAxes, HisPath1, HisPath2, oim, ellipseHandles);
        
    else
        displayOverlays(overlayAxes, Image, SpeedMode, ...
            FrameInfo, cptState.Particles, cptState.Spots, cptState.CurrentFrame, ShowThreshold2, ...
            Overlay, cptState.CurrentChannel, cptState.CurrentParticle, cptState.ZSlices, cptState.CurrentZ, numFrames, ...
            schnitzcells, UseSchnitz, DisplayRange, Ellipses, SpotFilter, ...
            ZoomMode, GlobalZoomMode, ZoomRange, xForZoom, yForZoom, fish, UseHistoneOverlay);
    end
    
    if ~isempty(xTrace)
        MaxZIndex = find(...
            cptState.Spots{cptState.CurrentChannel}(cptState.CurrentFrame).Fits(CurrentParticleIndex).z == ...
            cptState.Spots{cptState.CurrentChannel}(cptState.CurrentFrame).Fits(CurrentParticleIndex).brightestZ);
        CurrentZIndex = find(...
            cptState.Spots{cptState.CurrentChannel}(cptState.CurrentFrame).Fits(CurrentParticleIndex).z == ...
            cptState.CurrentZ);
        
        if isempty(CurrentZIndex)
            %             warning('This particle has a gap in its z-profile. This is
            %             highly suspect.'); %this if statement should only happen
            %             between two spots, not past the PSF boundaries
        end
        
    end
    
    %Check to see if spots structure contains multi-slice fields
    multi_slice_flag = isfield(cptState.Spots{cptState.CurrentChannel}(cptState.CurrentFrame).Fits ...
        (CurrentParticleIndex), 'IntegralZ');
    
    % PLOT SNIPPET
    
    FullSlicePath = [PreProcPath, filesep, Prefix, filesep, Prefix, '_', iIndex(cptState.CurrentFrame, 3) ...
        , '_z' iIndex(cptState.CurrentZ, 2) '_ch' iIndex(cptState.coatChannel, 2) '.tif'];
    
        [CurrentSnippet, hImage] = plotSnippet(snippetFigAxes, rawDataAxes, gaussianAxes, xTrace, ...
            CurrentZIndex, FullSlicePath, cptState.Spots, cptState.CurrentChannel, cptState.CurrentFrame, ...
            CurrentParticleIndex, ExperimentType, snippet_size, xSize, ...
            ySize, SnippetEdge, FrameInfo, CurrentSnippet, hImage, pixelSize);

    
    % PLOTS TRACE OF CURRENT PARTICLE
    if ~fish
        plottrace_argin = {};
        if exist('AmpIntegral', 'var')
            plottrace_argin = [plottrace_argin, AmpIntegral, GaussIntegral, AmpIntegral3,  ...
                ErrorIntegral, ErrorIntegral3, backGround3, ...
                AmpIntegralGauss3D, ErrorIntegralGauss3D, FrameIndicesToFit];
        end
        if ~isempty(Spots3D)
            plottrace_argin = [plottrace_argin, Spots3D];
        end
        [Frames, AmpIntegral, GaussIntegral, AmpIntegral3, ...
            ErrorIntegral, ErrorIntegral3,  backGround3, ...
            AmpIntegralGauss3D, ErrorIntegralGauss3D, PreviousParticle] =...
            ...
            plotTrace(traceFigAxes, ...
            ...
            FrameInfo, cptState.CurrentChannel, cptState.PreviousChannel, ...
            cptState.CurrentParticle, PreviousParticle, lastParticle, HideApprovedFlag, lineFitted, anaphaseInMins, ...
            ElapsedTime, schnitzcells, cptState.Particles, plot3DGauss, anaphase, prophase, metaphase, prophaseInMins, metaphaseInMins, Prefix, ...
            numFrames, cptState.CurrentFrame, cptState.ZSlices, cptState.CurrentZ, cptState.Spots, ...
            correspondingNCInfo, Coefficients, ExperimentType,cptState.PreviousFrame, Frames,...
            Channels, PreProcPath, DropboxFolder, plottrace_argin{:});
    end

    
    % PLOT Z SLICE RELATED FIGURES
    plotzvars = {zProfileFigAxes, zTraceAxes, ExperimentType, ...
    xTrace, cptState.Spots, cptState.CurrentFrame, cptState.CurrentChannel, CurrentParticleIndex, cptState.ZSlices, ...
    cptState.CurrentZ, CurrentZIndex, PreviousParticle, cptState.CurrentParticle, ...
    cptState.PreviousChannel, cptState.Particles, Frames, fish};
        if exist('MaxZProfile', 'var')
            plotzvars = [plotzvars, MaxZProfile];
        end
        [MaxZProfile, Frames] = plotZFigures(plotzvars{:});

    
    % UPDATE UICONTROLS
    updateControls(frame_num, z_num, particle_num, cptState.CurrentFrame, cptState.CurrentZ, cptState.CurrentParticle);
    
    set(0, 'CurrentFigure', Overlay);
    
    if isempty(SkipWaitForButtonPress)
        % Wait for user input to select command to execute
        ct = waitforbuttonpress; % ct=0 for click and ct=1 for keypress
        cc = get(Overlay, 'CurrentCharacter');
        cm2 = get(overlayAxes, 'CurrentPoint');
        
        current_axes = get(Overlay, 'CurrentAxes');
        
        if strcmpi(cc, '') || ct == 0
            cc = 'donothing';
        end
        
        is_control = isa(get(Overlay, 'CurrentObject'), 'matlab.ui.control.UIControl');
        
%         if ct == 0 && cm2(1, 1) < xSize && current_axes == overlayAxes ...
%                 &&~no_clicking &&~is_control
%             
%             [cptState.CurrentParticle, cptState.CurrentFrame, cptState.ManualZFlag] = toNearestParticle(cptState.Spots, ...
%                 cptState.Particles, cptState.CurrentFrame, cptState.CurrentChannel, UseHistoneOverlay, ...
%                 schnitzcells, [cm2(1, 1), cm2(2, 2)]);
%             
%         end
    else
        cc = SkipWaitForButtonPress;
        SkipWaitForButtonPress = [];
    end
    
    numValidFrames = length({cptState.Spots{1}.Fits});
    
    frameChangeKeyInput(cc);
    zSliceChangeKeyInput(cc);
    particleChangeKeyInput(cc);
    channelSwitchKeyInput(cc);
    
    if strcmpi(cc, 'donothing')
        %do nothing
        
    elseif cc == 'j'
        DisplayRange = []; %Temporary, sould be moved to FrameChangeHandler
        
    elseif (cc == '''') & (cptState.CurrentFrame < numValidFrames)%Move to the next skipped frame
        %within the particle
        cptState.PreviousFrame = cptState.CurrentFrame;
        cptState.CurrentFrame = nextSkippedFrame(cptState.Particles, cptState.CurrentChannel, ...
            cptState.CurrentParticle, cptState.CurrentFrame);
    elseif (cc == ';') & (cptState.CurrentFrame > 1)%Move to the previous skipped frame
        %within the particle
        cptState.PreviousFrame = cptState.CurrentFrame;
        cptState.CurrentFrame = previousSkippedFrame(cptState.Particles, cptState.CurrentChannel, ...
            cptState.CurrentParticle, cptState.CurrentFrame);
        
    elseif cc == 't'
        DisplayRange = []; %Temporary, sould be moved to ZSliceChangeHandler
    
    elseif cc == 'k'
        DisplayRange = []; %Temporary, sould be moved to ParticleChangeHandler

    elseif cc == 'g' & UseHistoneOverlay %Increase histone channel contrast
        
        if isempty(DisplayRange)'
            DisplayRange = [min(min(ImageHis)), max(max(ImageHis)) / 1.5];
        else
            DisplayRange = [DisplayRange(1), DisplayRange(2) / 1.5];
        end
            disp('increased nuclear contrast');

        
    elseif cc == 'b' & UseHistoneOverlay %Decrease histone channel contrast
        DisplayRange = [min(min(ImageHis)), max(max(ImageHis)) * 1.5];
        disp('decreased nuclear contrast');
        
    elseif cc == '#' %remove a spot from cptState.Spots and erase its frame in Particles
        [cptState.Spots, SpotFilter,cptState.CurrentFrame, ...
            cptState.CurrentParticle, cptState.Particles, cptState.ManualZFlag, lastParticle, PreviousParticle] = ...
            removeSpot(Frames, cptState.CurrentFrame, ...
            cptState.CurrentChannel, cptState.CurrentParticle, CurrentParticleIndex, cptState.Particles, cptState.Spots, SpotFilter ...
            );
        
    elseif cc == '^' %remove a whole trace from cptState.Spots and Particles. AR 7/9/2019 a work in progress
        for f = 1:length(Frames)
            [cptState.Spots, SpotFilter,cptState.CurrentFrame, ...
                cptState.CurrentParticle, cptState.Particles, cptState.ManualZFlag, lastParticle, PreviousParticle] = ...
                removeSpot(Frames, f, ...
                cptState.CurrentChannel, cptState.CurrentParticle, CurrentParticleIndex, cptState.Particles, cptState.Spots, SpotFilter, ...
                cptState.numParticles());
            
        end
    elseif cc == '[' | cc == '{' %#ok<*OR2> %Add particle and all of its shadows to cptState.Spots.
        PathPart1 = [PreProcPath, filesep, FilePrefix(1:end - 1), filesep, FilePrefix];
        PathPart2 = [nameSuffix, '.tif'];
        Path3 = [PreProcPath, filesep, Prefix, filesep, Prefix];
        [SpotFilter, cptState.Particles, cptState.Spots,...
            PreviousParticle, cptState.CurrentParticle, ZoomMode, GlobalZoomMode] = ...
            addSpot(ZoomMode, GlobalZoomMode, cptState.Particles, cptState.CurrentChannel, ...
            cptState.CurrentParticle, cptState.CurrentFrame, cptState.CurrentZ, Overlay, snippet_size, PixelsPerLine, ...
            LinesPerFrame, cptState.Spots, cptState.ZSlices, PathPart1, PathPart2, Path3, FrameInfo, pixelSize, ...
            SpotFilter, cc, xSize, ySize, NDigits,...
           Prefix, PreProcPath, ProcPath, cptState.coatChannel, UseHistoneOverlay, schnitzcells, nWorkers, plot3DGauss);
    elseif cc == 'r'
        cptState.Particles = orderParticles(cptState.numParticles(), cptState.CurrentChannel, cptState.Particles);
    elseif cc == 'f'
        [cptState.Particles, schnitzcells] = redoTracking(DataFolder, ...
            UseHistoneOverlay, FrameInfo, DropboxFolder, FilePrefix, schnitzcells, ...
            cptState.Particles, NChannels, cptState.CurrentChannel, cptState.numParticles());
    elseif cc == 'c'
        [PreviousParticle, cptState.Particles] = combineTraces(cptState.Spots, ...
            cptState.CurrentChannel, cptState.CurrentFrame, cptState.Particles, cptState.CurrentParticle);
     elseif cc == 'C'  %add ellipse
           
        [ConnectPositionx,ConnectPositiony] = ginput(1);
        
        cm = [round(ConnectPositionx), round(ConnectPositiony)];
        [Rows, Cols] = size(ImageHis);
         if (cm(1,2)>0)&(cm(1,1)>0)&(cm(1,2)<=Rows)&(cm(1,1)<=Cols)
            
            %Add a circle to this location with the mean radius of the
            %ellipses found in this frame
            
            %(x, y, a, b, theta, maxcontourvalue, time,
            %particle_id)
            if ~isempty(Ellipses{CurrentFrame})
                MeanRadius=mean((Ellipses{CurrentFrame}(:,3)+Ellipses{CurrentFrame}(:,4))/2);
            elseif ~isempty(Ellipses{CurrentFrame+1})
                MeanRadius=mean((Ellipses{CurrentFrame+1}(:,3)+Ellipses{CurrentFrame+1}(:,4))/2);
            elseif ~isempty(Ellipses{CurrentFrame-1})
                MeanRadius=mean((Ellipses{CurrentFrame-1}(:,3)+Ellipses{CurrentFrame-1}(:,4))/2);
            end
            
            try
                Ellipses{CurrentFrame}(end+1,:)=...
                    [cm(1,1),cm(1,2),MeanRadius,MeanRadius,0,0,0,0,0];
            catch
                Ellipses{CurrentFrame}(end+1,:)=...
                    [cm(1,1),cm(1,2),MeanRadius,MeanRadius,0,0,0,0];
            end
        end
    
         nucleiModified = true;
         
    elseif cc == 'V'
        %remove ellipse
           [ConnectPositionx,ConnectPositiony] = ginput(1);
        
        cm = [round(ConnectPositionx), round(ConnectPositiony)];
        [Rows, Cols] = size(ImageHis);
        
          if (cm(1,2)>0)&(cm(1,1)>0)&(cm(1,2)<=Rows)&(cm(1,1)<=Cols)
            %Find out which ellipses we clicked on so we can delete it
            
            %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
            Distances=sqrt((Ellipses{CurrentFrame}(:,1)-cm(1,1)).^2+...
                (Ellipses{CurrentFrame}(:,2)-cm(1,2)).^2);
            [~,MinIndex]=min(Distances);
            
            Ellipses{CurrentFrame}=[Ellipses{CurrentFrame}(1:MinIndex-1,:);...
                Ellipses{CurrentFrame}(MinIndex+1:end,:)];
          end
          
        nucleiModified = true;
        
    elseif cc == 'p' % Identify a particle. It will also tell you the particle associated with
        %  the clicked nucleus.
        identifyParticle(cptState.Spots, cptState.Particles, cptState.CurrentFrame, ...
            cptState.CurrentChannel, UseHistoneOverlay, schnitzcells);
    elseif cc == '\' %Moves to clicked particle.
        [cptState.CurrentParticle, cptState.CurrentFrame, cptState.ManualZFlag] = toNearestParticle(cptState.Spots, ...
            cptState.Particles, cptState.CurrentFrame, cptState.CurrentChannel, UseHistoneOverlay, schnitzcells);
        
    elseif cc == 'u'
        [x2, y2, z2] = SpotsXYZ(cptState.Spots{cptState.CurrentChannel}(cptState.CurrentFrame));
        
        if ~isempty(x2)
            ClickedSpot = ginput(1);
            
            UnfilterSpot(cptState.Spots{cptState.CurrentChannel}, SpotFilter{cptState.CurrentChannel}, ...
                ClickedSpot, cptState.Particles{cptState.CurrentChannel}, cptState.CurrentFrame)
        end
        
    elseif cc == 'i'
        warning(' AR 1/15/18: This is currently deprecated. Talk to HG if you need this function.')
    elseif cc == 'd' || cc == 'v' %d Separate traces forward at the current frame.
        [cptState.Particles, PreviousParticle] = separateTraces(cptState.Particles, ...
            cptState.CurrentChannel, cptState.CurrentFrame, cptState.CurrentParticle);
    elseif cc == 'q' %Approve a trace
        
        if cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Approved == 1
            cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Approved = 2;
        elseif cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Approved == 0
            cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Approved = 1;
        elseif cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Approved == 2
            cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Approved = 0;
        end
        
    elseif cc == 'w' %Disapprove a trace
        
        if cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Approved == -1
            cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Approved = 0;
        else
            cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Approved = -1;
        end
        
    elseif cc == 's'
        saveChanges(NChannels, cptState.Particles, cptState.Spots, SpotFilter, DataFolder, ...
            FrameInfo, UseHistoneOverlay, FilePrefix, ...
            schnitzcells, DropboxFolder);

    elseif cc == 'h'
        
        if HideApprovedFlag == 0
            HideApprovedFlag = 1; %Show only non-approved traces
        elseif HideApprovedFlag == 1
            HideApprovedFlag = 2; %Show only yellow and red traces
        elseif HideApprovedFlag == 2
            HideApprovedFlag = 0;
        end
        
        %HideApprovedFlag=~HideApprovedFlag;
    elseif cc == 'o'
        
        if ~GlobalZoomMode
            ZoomMode = ~ZoomMode;
        elseif GlobalZoomMode & ~ZoomMode
            GlobalZoomMode = false;
        end
        
        
    elseif cc == '+'
        
        if ~ZoomMode & ~GlobalZoomMode
                [ConnectPositionx, ConnectPositiony] = ginput(1);
                xForZoom = round(ConnectPositionx);
                yForZoom = round(ConnectPositiony);
            GlobalZoomMode = true;
        elseif ZoomMode & ~GlobalZoomMode
            ZoomMode = false;
        elseif ~ZoomMode & GlobalZoomMode
            GlobalZoomMode = false;
        end
        
    elseif (cc == 'm') & (cptState.CurrentParticle < cptState.numParticles())
        [lineFitted, cptState.CurrentParticle, cptState.CurrentFrame, cptState.ManualZFlag, DisplayRange] = ...
            goNextParticle(cptState.CurrentParticle, cptState.CurrentChannel, HideApprovedFlag, cptState.Particles);
        
    elseif (cc == 'n') & (cptState.CurrentParticle > 1)
        [lineFitted, cptState.CurrentParticle, cptState.CurrentFrame, cptState.ManualZFlag, DisplayRange] = ...
            goPreviousParticle(cptState.CurrentParticle, cptState.CurrentChannel, HideApprovedFlag, cptState.Particles);
        
    elseif cc == 'e'
        cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).FrameApproved(cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Frame == cptState.CurrentFrame) = ...
            ~cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).FrameApproved(cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Frame == cptState.CurrentFrame);
        
        %Schnitzcells specific
        
    elseif cc == 'l' %Split a nucleus and select one or two daughter cells or stop the lineage
        [cptState.Particles, PreviousParticle, schnitzcells] = splitNuclei(schnitzcells, ...
            cptState.CurrentFrame, cptState.CurrentChannel, cptState.CurrentParticle, cptState.Particles);
    elseif cc == '2' %2 set parent of current nucleus
        schnitzcells = setParentNucleus(schnitzcells, ...
            cptState.CurrentFrame, cptState.CurrentChannel, cptState.CurrentParticle, cptState.Particles);
    elseif cc == '~' %Switch projection mode
        projectionMode = chooseProjection;
        disp(['projectionMode : ' projectionMode])
        
    elseif cc == '!' %Increase contrast in the Overlay figure
        
        if isempty(DisplayRangeSpot)
            DisplayRangeSpot = [min(min(Image)), max(max(Image)) / 1.5];
        else
            DisplayRangeSpot = [DisplayRangeSpot(1), DisplayRangeSpot(2) / 1.5];
        end
        
        disp('increased spot contrast');
        
    elseif cc == '@' %Decrease spot channel contrast
        DisplayRangeSpot = [min(Image(:)), max(Image(:)) * 1.5];
        
         disp('decreased spot contrast');
         
    elseif cc == '$' %add particle to nucleus
            cptState.Particles = addNucleusToParticle(cptState.Particles, cptState.CurrentFrame, ...
                cptState.CurrentChannel, UseHistoneOverlay, schnitzcells, cptState.CurrentParticle);
    elseif cc == '0' %Debugging mode
        keyboard;
        
    end
    
    %% Main loop - end
end

save([DataFolder, filesep, 'FrameInfo.mat'], 'FrameInfo')

%If we only have one channel bring Particles back to the legacy
%format without any cells
if NChannels == 1
    cptState.Particles = cptState.Particles{1};
    cptState.Spots = cptState.Spots{1};
    SpotFilter = SpotFilter{1};
end

Spots = cptState.Spots; % store Spots as a local variable so we can save it
Particles = cptState.Particles; % store Particles as a local variable so we can save it

if UseHistoneOverlay
    save([DataFolder, filesep, 'Particles.mat'], 'Particles', 'SpotFilter', '-v7.3')
    save([DataFolder, filesep, 'Spots.mat'], 'Spots', '-v7.3')
    save([DropboxFolder, filesep, FilePrefix(1:end - 1), filesep, FilePrefix(1:end - 1), '_lin.mat'], 'schnitzcells')
else
    save([DataFolder, filesep, 'Particles.mat'], 'Particles', 'SpotFilter', '-v7.3')
    save([DataFolder, filesep, 'Spots.mat'], 'Spots', '-v7.3')
end

close all

if ishandle(controls)
    close(controls)
end

disp('Particles saved.')
disp(['(Left off at Particle #', num2str(cptState.CurrentParticle), ')'])

if nucleiModified 
    
   save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses');

    %Decide whether we need to re-track
    userPrompt = 'Did you make changes to nuclei and thus require re-tracking? (y/n)';
    reTrackAnswer = inputdlg(userPrompt);
    if contains(reTrackAnswer,'n')
        disp('Ellipses saved. Per user input, not re-tracking. Exiting.')
    else
      opts = {};
       if fish
           opts = [opts, 'markandfind'];
       end
       disp('Ellipses saved. Running TrackNuclei to incorporate changes.')
       TrackNuclei(Prefix,'NoBulkShift','ExpandedSpaceTolerance', 1.5, 'retrack', 'nWorkers', 1, opts{:}); 
    end
end

end
