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
% . , Move a frame forward/backward
% > < Move five frames forward/backward
% ; ' Move to the next empty frame within a particle
% a z Move up/down in Z
% t jump to a specific z-slice
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
% $ add particle to the nearest nucleus
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

xForZoom = 0;
yForZoom = 0;

% Parameters for fitting
lineFitted = 0; % equals 1 if a line has been fitted
fitApproved = 0;
FramesToFit = []; % actual frames of the movie that were used for fitting
FrameIndicesToFit = []; % index of the current particle that were used for fitting
Coefficients = []; % coefficients of the fitted line

SnippetEdge = 13; %Size of the snippets generated by Michael's code in pixels.
storedTimeProjection = []; % Don't need to wait for timeProjection to finish each time its called

[sortByFrames, sortByLength, ForCompileAll, SpeedMode, ~, ...
    ncRange, projectionMode, plot3DGauss, intScale, NC, ...
    startNC, endNC, optionalResults] = determineCheckParticleTrackingOptions(varargin{:});


%% Information about about folders

% Get the folders

[rawDataPath,ProcPath,DropboxFolder,MS2CodePath, PreProcPath,...
    rawDataFolder, Prefix, ExperimentType,Channel1,Channel2,OutputFolder,...
    Channel3, spotChannels, MovieDataBaseFolder, movieDatabase]...
    = readMovieDatabase(Prefix, optionalResults);

DataFolder = [DropboxFolder, filesep, Prefix];
FilePrefix = [Prefix, '_'];


[Particles, SpotFilter, Spots, FrameInfo] = loadCheckParticleTrackingMats(DataFolder, PreProcPath);

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
[Particles] = addFrameApproved(NChannels, Particles);

[Ellipses, UseHistoneOverlay, UseSchnitz] = checkHistoneAndNuclearSegmentation(PreProcPath, FilePrefix, NDigits, DropboxFolder);

% we name the variable DataFolderColumnValue to avoid shadowing previously defined DataFolder var, which is actually a subfolder inside dropbox
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution, ...
    Channel1, Channel2, Objective, Power, DataFolderColumnValue, ~, Comments, ...
    nc9, nc10, nc11, nc12, nc13, nc14, CF, Channel3, prophase, metaphase] = getExperimentDataFromMovieDatabase(Prefix, movieDatabase);

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
catch
    %
end

% this save() is here so that a user will still get an updated
% frameinfo.mat even if they abort checkparticletracking without saving
% (to prevent issues with compileparticles)
save([DataFolder, filesep, 'FrameInfo.mat'], 'FrameInfo');

%Check if we have already determined nc
if (~isfield(FrameInfo, 'nc')) && (~UseHistoneOverlay)
    %do nothing
    
elseif UseSchnitz
    load([DropboxFolder, filesep, FilePrefix(1:end - 1), filesep, FilePrefix(1:end - 1), '_lin.mat']);
    
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
ZSlices = FrameInfo(1).NumberSlices + 2; %Note that the blank slices are included
CurrentZ = round(ZSlices / 2);
ManualZFlag = 0;
CurrentParticle = 1;
PreviousParticle = 1;
lastParticle = 0; %this gets flagged if there's a drop to one particle within the Particles structure.
CurrentFrameWithinParticle = 1;
CurrentChannel = 1;
PreviousChannel = CurrentChannel;
if ~isempty(Particles{CurrentChannel})
    CurrentFrame = Particles{CurrentChannel}(CurrentParticle).Frame(CurrentFrameWithinParticle);
else
    error('Looks like the Particles structure is empty. There''s nothing to check.');
end
DisplayRange = [];
DisplayRangeSpot = [];
ZoomMode = 0;
GlobalZoomMode = 0;
ZoomRange = 50;
nameSuffix = '';


% Changing the intial frames and particle if justNC13
if ncRange
    
    if strcmpi('nc15', endNC)
        lastNCFrame = numFrames;
    else
        lastNCFrame = eval(endNC) - 1; % This will not include the 1st frame of the next NC
    end
    
    firstNCFrame = eval(startNC);
    particlesInRange = particlesWithinFrames(Prefix, firstNCFrame, lastNCFrame);
    CurrentParticle = particlesInRange(1);
    CurrentFrame = Particles{1}(CurrentParticle).Frame(1);
    disp(['nc range: ' num2str(NC)])
    disp(['start frame: ' num2str(firstNCFrame)])
    disp(['end frame: ' num2str(lastNCFrame)])
    disp(['Particles in range: ' num2str(particlesInRange)])
    disp(['Number of Particles: ' num2str(length(particlesInRange))])
end

%Define user interface
[Overlay, overlayAxes, snippetFigAxes, rawDataAxes, gaussianAxes, traceFigAxes, zProfileFigAxes,...
    zTraceAxes,HisOverlayFig,HisOverlayFigAxes] = checkParticleTracking_drawGUI(UseHistoneOverlay);

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

    function z_num_changed(~, ~)
        
        figure(Overlay);
        [CurrentZ, ManualZFlag] = changeZSlice(str2double(z_num.Value), ZSlices);
        robot.keyPress(fake_event);
        robot.keyRelease(fake_event);
        
    end

particle_num.ValueChangedFcn = @particle_num_changed;

    function particle_num_changed(~, ~)
        
        figure(Overlay);
        [CurrentParticle, CurrentFrame, ManualZFlag] = changeParticle(...
            str2double(particle_num.Value), Particles, numParticles, CurrentChannel);
        robot.keyPress(fake_event);
        robot.keyRelease(fake_event);
        
    end

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
        [numParticles, SpotFilter, Particles, Spots, PreviousParticle, CurrentParticle] = ...
            addSpot(ZoomMode, GlobalZoomMode, Particles, CurrentChannel, ...
            CurrentParticle, CurrentFrame, CurrentZ, Overlay, snippet_size, PixelsPerLine, ...
            LinesPerFrame, Spots, ZSlices, PathPart1, PathPart2, Path3, FrameInfo, pixelSize, ...
            SpotFilter, numParticles, smart_add, xSize, ySize, NDigits, intScale,...
            coatChannel, UseHistoneOverlay, schnitzcells);
        
        robot.keyPress(fake_event);
        robot.keyRelease(fake_event);
        no_clicking = false;
    end

delete_spot.ButtonPushedFcn = @delete_spot_pushed;

    function delete_spot_pushed(~, ~)
        no_clicking = true;
        figure(Overlay);
        [Spots, SpotFilter, ZoomMode, GlobalZoomMode, CurrentFrame, ...
            CurrentParticle, Particles, ManualZFlag, DisplayRange, lastParticle, PreviousParticle] = ...
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
            fitInitialSlope(CurrentParticle, Particles, Spots, CurrentChannel, schnitzcells, ...
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
            Particles{CurrentChannel}(CurrentParticle).fitApproved = 1;
            Particles{CurrentChannel}(CurrentParticle).Coefficients = Coefficients;
            %Particles{CurrentChannel}(CurrentParticle).lineFitHandle =  lineFitHandle;
            Particles{CurrentChannel}(CurrentParticle).fittedFrames = FrameIndicesToFit; % use the index of particle trace for convenience
            %Particles{CurrentChannel}(CurrentParticle).fittedYSegment = currentYSegment;
        else
            Particles{CurrentChannel}(CurrentParticle).fitApproved = 0;
            Particles{CurrentChannel}(CurrentParticle).Coefficients = [];
            %Particles{CurrentChannel}(CurrentParticle).lineFitHandle =  [];
            Particles{CurrentChannel}(CurrentParticle).fittedFrames = [];
            %Particles{CurrentChannel}(CurrentParticle).fittedYSegment = [];
        end
        
    end

% Create the approved field if it does not exist
for NCh = 1:NChannels
    
    if ~isfield(Particles{NCh}, 'Approved')
        
        for i = 1:length(Particles{NCh})
            Particles{NCh}(i).Approved = 0;
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

        if contains(Channel1{1}, 'MCP') || contains(Channel1{1}, 'PCP')
            nameSuffix = ['_ch', iIndex(1, 2)];
            coatChannel = 1;
        elseif contains(Channel2{1}, 'MCP') || contains(Channel2{1}, 'PCP')
            nameSuffix = ['_ch', iIndex(2, 2)];
            coatChannel = 2;
        end

    elseif strcmpi(ExperimentType, '2spot2color')
        %We are assuming that channels 1 and 2 are assigned to coat
        %proteins. We should do a better job with this.
        coatChannels = [1, 2];
        coatChannel = coatChannels(CurrentChannel);
    else
        error('Experiment type not recognized')
    end

    %Update the name suffix
    if strcmpi(ExperimentType, '2spot2color')
        nameSuffix = ['_ch', iIndex(coatChannel, 2)];
    end
    %%
    numParticles = length(Particles{CurrentChannel});
    
    %Get the coordinates of all the spots in this frame
    [x, y, z] = SpotsXYZ(Spots{CurrentChannel}(CurrentFrame));
    
    %If the approved field does not exist create it
    if ~isfield(Particles{CurrentChannel}, 'Approved')
        
        for i = 1:numParticles
            Particles{CurrentChannel}(i).Approved = 0;
        end
        
    end
    
    %Pull out the right particle if it exists in this frame
    CurrentParticleIndex = ...
        Particles{CurrentChannel}(CurrentParticle).Index(Particles{CurrentChannel}(CurrentParticle).Frame == ...
        CurrentFrame);
    %This is the position of the current particle
    xTrace = x(CurrentParticleIndex);
    yTrace = y(CurrentParticleIndex);
    
    if (~isempty(xTrace)) && (~ManualZFlag)
        CurrentZ = z(CurrentParticleIndex);
        CurrentZIndex = find(...
            Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).z == ...
            CurrentZ);
        ManualZFlag = 0;
    end
    
    
    if strcmpi(projectionMode, 'None')
        Image = imread([PreProcPath, filesep, FilePrefix(1:end - 1), filesep, ...
            FilePrefix, iIndex(CurrentFrame, NDigits), '_z', iIndex(CurrentZ, 2), nameSuffix, '.tif']);
    elseif strcmpi(projectionMode, 'Max Z')
        Image = zProjections(Prefix, coatChannel, CurrentFrame, ZSlices, NDigits, DropboxFolder, PreProcPath, FrameInfo, 'max');
    elseif strcmpi(projectionMode, 'Median Z')
        Image = zProjections(Prefix, coatChannel, CurrentFrame, ZSlices, NDigits, DropboxFolder, PreProcPath, FrameInfo, 'median');
    elseif strcmpi(projectionMode, 'Max Z and Time')
        
        if isempty(storedTimeProjection)
            
            if ncRange
                Image = timeProjection(Prefix, coatChannel,FrameInfo, DropboxFolder,PreProcPath, 'nc', NC);
                storedTimeProjection = Image;
            else
                Image = timeProjection(Prefix, CurrentChannel,FrameInfo, DropboxFolder,PreProcPath);
                storedTimeProjection = Image;
            end
            
        else
            Image = storedTimeProjection;
        end
        
    end
    
    
    
    set(0, 'CurrentFigure', Overlay);
    imshow(Image, DisplayRangeSpot, 'Border', 'Tight', 'Parent', overlayAxes, ...
        'InitialMagnification', 'fit')
    hold(overlayAxes, 'on')
    
    if UseHistoneOverlay
        HisPath1 = [PreProcPath, filesep, FilePrefix(1:end - 1), filesep, ...
            FilePrefix(1:end - 1), '-His_', iIndex(CurrentFrame, NDigits), '.tif'];
        HisPath2 = [PreProcPath, filesep, FilePrefix(1:end - 1), filesep, ...
            FilePrefix(1:end - 1), '_His_', iIndex(CurrentFrame, NDigits), '.tif'];
        
        [ImageHis, xForZoom, yForZoom] = displayOverlays(overlayAxes, Image, SpeedMode, FrameInfo, Particles, ...
            Spots, CurrentFrame, ShowThreshold2, ...
            Overlay, CurrentChannel, CurrentParticle, ZSlices, CurrentZ, numFrames, ...
            schnitzcells, UseSchnitz, DisplayRange, Ellipses, SpotFilter, ZoomMode, GlobalZoomMode, ...
            ZoomRange, xForZoom, yForZoom, UseHistoneOverlay, HisOverlayFigAxes, HisPath1, HisPath2);
        
    else
        displayOverlays(overlayAxes, Image, SpeedMode, ...
            FrameInfo, Particles, Spots, CurrentFrame, ShowThreshold2, ...
            Overlay, CurrentChannel, CurrentParticle, ZSlices, CurrentZ, numFrames, ...
            schnitzcells, UseSchnitz, DisplayRange, Ellipses, SpotFilter, ...
            ZoomMode, GlobalZoomMode, ZoomRange, xForZoom, yForZoom, UseHistoneOverlay);
    end
    
    if ~isempty(xTrace)
        MaxZIndex = find(...
            Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).z == ...
            Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).brightestZ);
        CurrentZIndex = find(...
            Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).z == ...
            CurrentZ);
        
        if isempty(CurrentZIndex)
            %             warning('This particle has a gap in its z-profile. This is
            %             highly suspect.'); %this if statement should only happen
            %             between two spots, not past the PSF boundaries
        end
        
    end
    
    %Check to see if spots structure contains multi-slice fields
    multi_slice_flag = isfield(Spots{CurrentChannel}(CurrentFrame).Fits ...
        (CurrentParticleIndex), 'IntegralZ');
    
    % PLOTS SNIPPET
    
    FullSlicePath = [PreProcPath, filesep, Prefix, filesep, Prefix, '_', iIndex(CurrentFrame, 3) ...
        , '_z' iIndex(CurrentZ, 2) '_ch' iIndex(coatChannel, 2) '.tif'];
    
    if exist('CurrentSnippet', 'var')
        CurrentSnippet = plotSnippet(snippetFigAxes, rawDataAxes, gaussianAxes, xTrace, ...
            CurrentZIndex, FullSlicePath, Spots, CurrentChannel, CurrentFrame, ...
            CurrentParticleIndex, ExperimentType, intScale, snippet_size, xSize, ...
            ySize, SnippetEdge, FrameInfo, CurrentSnippet);
    else
        CurrentSnippet = plotSnippet(snippetFigAxes, rawDataAxes, gaussianAxes, xTrace, ...
            CurrentZIndex, FullSlicePath, Spots, CurrentChannel, CurrentFrame, ...
            CurrentParticleIndex, ExperimentType, intScale, snippet_size, xSize, ...
            ySize, SnippetEdge, FrameInfo);
    end
    
    % PLOTS TRACE OF CURRENT PARTICLE
    if exist('AmpIntegral', 'var')
        [Frames, AmpIntegral, GaussIntegral, AmpIntegral3, AmpIntegral5, ...
            ErrorIntegral, ErrorIntegral3, ErrorIntegral5, backGround3, ...
            AmpIntegralGauss3D, ErrorIntegralGauss3D, PreviousParticle] = plotTrace(traceFigAxes, ...
            FrameInfo, CurrentChannel, PreviousChannel, ...
            CurrentParticle, PreviousParticle, lastParticle, HideApprovedFlag, lineFitted, anaphaseInMins, ...
            ElapsedTime, schnitzcells, Particles, plot3DGauss, anaphase, prophase, metaphase, prophaseInMins, metaphaseInMins, Prefix, ...
            numFrames, CurrentFrame, ZSlices, CurrentZ, Spots, ...
            correspondingNCInfo, Coefficients, ExperimentType, ...
            Frames, AmpIntegral, GaussIntegral, AmpIntegral3, AmpIntegral5, ...
            ErrorIntegral, ErrorIntegral3, ErrorIntegral5, backGround3, ...
            AmpIntegralGauss3D, ErrorIntegralGauss3D, FrameIndicesToFit);
    else
        [Frames, AmpIntegral, GaussIntegral, AmpIntegral3, AmpIntegral5, ...
            ErrorIntegral, ErrorIntegral3, ErrorIntegral5, backGround3, ...
            AmpIntegralGauss3D, ErrorIntegralGauss3D, PreviousParticle] = plotTrace(traceFigAxes, ...
            FrameInfo, CurrentChannel, PreviousChannel, ...
            CurrentParticle, PreviousParticle, lastParticle, HideApprovedFlag, lineFitted, anaphaseInMins, ...
            ElapsedTime, schnitzcells, Particles, plot3DGauss, anaphase, prophase, metaphase, prophaseInMins, metaphaseInMins, Prefix, ...
            numFrames, CurrentFrame, ZSlices, CurrentZ, Spots, ...
            correspondingNCInfo, Coefficients, ExperimentType);
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
        % Wait for user input to select command to execute
        ct = waitforbuttonpress; % ct=0 for click and ct=1 for keypress
        cc = get(Overlay, 'CurrentCharacter');
        cm2 = get(overlayAxes, 'CurrentPoint');
        
        current_axes = get(Overlay, 'CurrentAxes');
        
        if strcmpi(cc, '') || ct == 0
            cc = 'donothing';
        end
        
        is_control = isa(get(Overlay, 'CurrentObject'), 'matlab.ui.control.UIControl');
        
        if ct == 0 && cm2(1, 1) < xSize && current_axes == overlayAxes ...
                &&~no_clicking &&~is_control
            
            [CurrentParticle, CurrentFrame, ManualZFlag] = toNearestParticle(Spots, ...
                Particles, CurrentFrame, CurrentChannel, UseHistoneOverlay, ...
                schnitzcells, [cm2(1, 1), cm2(2, 2)]);
        end
    else
        cc = SkipWaitForButtonPress;
        SkipWaitForButtonPress = [];
    end
    
    numValidFrames = length({Spots{1}.Fits});
    
    if strcmpi(cc, 'donothing')
        %do nothing
        
    elseif cc == '.'%Move forward one frame
        [CurrentFrame, ManualZFlag] = changeFrame(CurrentFrame + 1, numValidFrames);
    elseif (cc == ',')%Move backward one frame
        [CurrentFrame, ManualZFlag] = changeFrame(CurrentFrame - 1, numValidFrames);
    elseif (cc == '>')%Move forward five frames
        [CurrentFrame, ManualZFlag] = changeFrame(CurrentFrame + 5, numValidFrames);
    elseif (cc == '<')%#ok<*AND2>%Move backward five frames
        [CurrentFrame, ManualZFlag] = changeFrame(CurrentFrame - 5, numValidFrames);
    elseif cc == 'j'
        
        try
            iJump = inputdlg('Frame to jump to:', ...
                'Move to frame');
            iJump = str2double(iJump{1});
        catch
            iJump = CurrentFrame;
        end
        
        [CurrentFrame, ManualZFlag] = changeFrame(iJump, numValidFrames);
        DisplayRange = [];
        
    elseif (cc == '''') & (CurrentFrame < numValidFrames)%Move to the next skipped frame
        %within the particle
        CurrentFrame = nextSkippedFrame(Particles, CurrentChannel, ...
            CurrentParticle, CurrentFrame);
    elseif (cc == ';') & (CurrentFrame > 1)%Move to the previous skipped frame
        %within the particle
        CurrentFrame = previousSkippedFrame(Particles, CurrentChannel, ...
            CurrentParticle, CurrentFrame);
        
    elseif (cc == 'a') & (CurrentZ < ZSlices)%Move up in Z
        [CurrentZ, ManualZFlag] = changeZSlice(CurrentZ + 1, ZSlices);
    elseif (cc == 'z') & (CurrentZ > 1)%Move down in Z
        [CurrentZ, ManualZFlag] = changeZSlice(CurrentZ - 1, ZSlices);
    elseif cc == 't'
        
        try
            iJump = inputdlg('z-slice to jump to:', ...
                'Move to z-slice');
            iJump = str2double(iJump{1});
        catch
            iJump = CurrentFrame;
        end
        [CurrentZ, ManualZFlag] = changeZSlice(iJump, ZSlices);
        DisplayRange = [];
    
    elseif cc == 'k'
        
        try
            ParticleJump = inputdlg('Particle to jump to:', ...
                'Move to particle');
            ParticleJump = str2double(ParticleJump{1});
        catch
            ParticleJump = CurrentParticle;
        end
        
        [CurrentParticle, CurrentFrame, ManualZFlag] = ...
            changeParticle(ParticleJump, Particles, numParticles, CurrentChannel);
        DisplayRange = [];
    elseif cc == 'g' & UseHistoneOverlay%Increase histone channel contrast
        
        if isempty(DisplayRange)'
            DisplayRange = [min(min(ImageHis)), max(max(ImageHis)) / 1.5];
        else
            DisplayRange = [DisplayRange(1), DisplayRange(2) / 1.5];
        end
        
    elseif cc == 'b' & UseHistoneOverlay%Decrease histone channel contrast
        DisplayRange = [min(min(ImageHis)), max(max(ImageHis)) * 1.5];
        
    elseif cc == '#'%remove a spot from Spots and erase its frame in Particles
        [Spots, SpotFilter,CurrentFrame, ...
            CurrentParticle, Particles, ManualZFlag, lastParticle, PreviousParticle] = ...
            removeSpot(Frames, CurrentFrame, ...
            CurrentChannel, CurrentParticle, CurrentParticleIndex, Particles, Spots, SpotFilter, ...
            numParticles);
    elseif cc == '^'%remove a whole trace from Spots and Particles. AR 7/9/2019 a work in progress
        for f = 1:length(Frames)
            [Spots, SpotFilter,CurrentFrame, ...
                CurrentParticle, Particles, ManualZFlag, lastParticle, PreviousParticle] = ...
                removeSpot(Frames, f, ...
                CurrentChannel, CurrentParticle, CurrentParticleIndex, Particles, Spots, SpotFilter, ...
                numParticles);
        end
    elseif cc == '[' | cc == '{'%#ok<*OR2> %Add particle and all of its shadows to Spots.
        PathPart1 = [PreProcPath, filesep, FilePrefix(1:end - 1), filesep, FilePrefix];
        PathPart2 = [nameSuffix, '.tif'];
        Path3 = [PreProcPath, filesep, Prefix, filesep, Prefix];
        [numParticles, SpotFilter, Particles, Spots,...
            PreviousParticle, CurrentParticle, ZoomMode, GlobalZoomMode] = ...
            addSpot(ZoomMode, GlobalZoomMode, Particles, CurrentChannel, ...
            CurrentParticle, CurrentFrame, CurrentZ, Overlay, snippet_size, PixelsPerLine, ...
            LinesPerFrame, Spots, ZSlices, PathPart1, PathPart2, Path3, FrameInfo, pixelSize, ...
            SpotFilter, numParticles, cc, xSize, ySize, NDigits,...
            intScale, Prefix, PreProcPath, ProcPath, coatChannel, UseHistoneOverlay, schnitzcells);
    elseif cc == 'r'
        Particles = orderParticles(numParticles, CurrentChannel, Particles);
    elseif cc == 'f'
        [Particles, schnitzcells] = redoTracking(DataFolder, ...
            UseHistoneOverlay, FrameInfo, DropboxFolder, FilePrefix, schnitzcells, ...
            Particles, NChannels, CurrentChannel, numParticles);
    elseif cc == 'c'
        [PreviousParticle, Particles] = combineTraces(Spots, ...
            CurrentChannel, CurrentFrame, Particles, CurrentParticle);
    elseif cc == 'p'%Identify a particle. It will also tell you the particle associated with
        %  the clicked nucleus.
        identifyParticle(Spots, Particles, CurrentFrame, ...
            CurrentChannel, UseHistoneOverlay, schnitzcells);
    elseif cc == '\'%Moves to clicked particle.
        [CurrentParticle, CurrentFrame, ManualZFlag] = toNearestParticle(Spots, ...
            Particles, CurrentFrame, CurrentChannel, UseHistoneOverlay, schnitzcells);
    elseif cc == 'u'
        [x2, y2, z2] = SpotsXYZ(Spots{CurrentChannel}(CurrentFrame));
        
        if ~isempty(x2)
            ClickedSpot = ginput(1);
            
            UnfilterSpot(Spots{CurrentChannel}, SpotFilter{CurrentChannel}, ...
                ClickedSpot, Particles{CurrentChannel}, CurrentFrame)
        end
        
    elseif cc == 'i'
        warning(' AR 1/15/18: This is currently deprecated. Talk to HG if you need this function.')
    elseif cc == 'd' || cc == 'v'%d Separate traces forward at the current frame.
        [Particles, PreviousParticle] = separateTraces(Particles, ...
            CurrentChannel, CurrentFrame, CurrentParticle);
    elseif cc == 'q'%Approve a trace
        
        if Particles{CurrentChannel}(CurrentParticle).Approved == 1
            Particles{CurrentChannel}(CurrentParticle).Approved = 2;
        elseif Particles{CurrentChannel}(CurrentParticle).Approved == 0
            Particles{CurrentChannel}(CurrentParticle).Approved = 1;
        elseif Particles{CurrentChannel}(CurrentParticle).Approved == 2
            Particles{CurrentChannel}(CurrentParticle).Approved = 0;
        end
        
    elseif cc == 'w'%Disapprove a trace
        
        if Particles{CurrentChannel}(CurrentParticle).Approved == -1
            Particles{CurrentChannel}(CurrentParticle).Approved = 0;
        else
            Particles{CurrentChannel}(CurrentParticle).Approved = -1;
        end
        
    elseif cc == 's'
        saveChanges(NChannels, Particles, Spots, SpotFilter, DataFolder, ...
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
        
    elseif (cc == 'm') & (CurrentParticle < numParticles)
        [lineFitted, CurrentParticle, CurrentFrame, ManualZFlag, DisplayRange] = ...
            goNextParticle(CurrentParticle, CurrentChannel, HideApprovedFlag, Particles);
    elseif (cc == 'n') & (CurrentParticle > 1)
        [lineFitted, CurrentParticle, CurrentFrame, ManualZFlag, DisplayRange] = ...
            goPreviousParticle(CurrentParticle, CurrentChannel, HideApprovedFlag, Particles);
    elseif cc == 'e'
        Particles{CurrentChannel}(CurrentParticle).FrameApproved(Particles{CurrentChannel}(CurrentParticle).Frame == CurrentFrame) = ...
            ~Particles{CurrentChannel}(CurrentParticle).FrameApproved(Particles{CurrentChannel}(CurrentParticle).Frame == CurrentFrame);
        
        %Schnitzcells specific
        
    elseif cc == 'l'%Split a nucleus and select one or two daughter cells or stop the lineage
        [Particles, PreviousParticle, schnitzcells] = splitNuclei(schnitzcells, ...
            CurrentFrame, CurrentChannel, CurrentParticle, Particles);
    elseif cc == '2'%2 set parent of current nucleus
        schnitzcells = setParentNucleus(schnitzcells, ...
            CurrentFrame, CurrentChannel, CurrentParticle, Particles);
    elseif cc == '8' && NChannels > 1%Switch channels
        [CurrentChannel, PreviousChannel, coatChannel, CurrentParticle] = ...
            switchChannels(CurrentChannel, CurrentParticle, Particles, ...
            UseHistoneOverlay, coatChannels, NChannels);
    elseif cc == '~'%Switch projection mode
        projectionMode = chooseProjection;
        disp(['projectionMode : ' projectionMode])
        
    elseif cc == '!'%Increase contrast in the Overlay figure
        
        if isempty(DisplayRangeSpot)
            DisplayRangeSpot = [min(min(Image)), max(max(Image)) / 1.5];
        else
            DisplayRangeSpot = [DisplayRangeSpot(1), DisplayRangeSpot(2) / 1.5];
        end
        
    elseif cc == '@'%Decrease spot channel contrast
        DisplayRangeSpot = [min(Image(:)), max(Image(:)) * 1.5];
    elseif cc == '$' %add particle to nucleus
            Particles = addNucleusToParticle(Particles, CurrentFrame, ...
                CurrentChannel, UseHistoneOverlay, schnitzcells, CurrentParticle);
    elseif cc == '0'%Debugging mode
        keyboard;
        
    end
    
    %% Main loop - end
end

save([DataFolder, filesep, 'FrameInfo.mat'], 'FrameInfo')

%If we only have one channel bring Particles back to the legacy
%format without any cells
if NChannels == 1
    Particles = Particles{1};
    Spots = Spots{1};
    SpotFilter = SpotFilter{1};
end

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
disp(['(Left off at Particle #', num2str(CurrentParticle), ')'])

end