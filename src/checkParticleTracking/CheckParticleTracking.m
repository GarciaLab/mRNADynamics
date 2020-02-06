function [outs, movieMat, hisMat] = CheckParticleTracking(Prefix, varargin)
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
% outs.Particles: A modified Particles
% outs.Spots: A modified Spots
% outs.SpotFilter: A modified SpotFilter
% outs.schnitzcells: A modified schnitzcells
% outs.FrameInfo: A modified FrameInfo
% movieCell: TBD
% hisCell: TBD
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

% Parameters for fitting
FramesToFit = []; % actual frames of the movie that were used for fitting

SnippetEdge = 13; %Size of the snippets generated by Michael's code in pixels.
storedTimeProjection = []; % Don't need to wait for timeProjection to finish each time its called

[sortByFrames, sortByLength, ForCompileAll, SpeedMode, ~, ...
    ncRange, projectionMode, plot3DGauss, NC, ...
    startNC, endNC, optionalResults, nWorkers, fish,...
    noHisOverlay, multiView, preStructs, preMovie, movieMat, hisMat] = determineCheckParticleTrackingOptions(varargin{:});

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

if isempty(preStructs)
    [Particles, SpotFilter, Spots, FrameInfo, schnitzcells, Spots3D] = loadCheckParticleTrackingMats(DataFolder, PreProcPath, FilePrefix);
else
    Particles = preStructs{1};
    SpotFilter = preStructs{3};
    Spots = preStructs{2};
    schnitzcells = preStructs{4};
    FrameInfo = preStructs{5};
    Spots3D = {};
end


[xSize, ySize, pixelSize, zStep, snippet_size,...
    nFrames, nSlices, nDigits] = getFrameInfoParams(FrameInfo);

nSlices = nSlices + 2; %due to padding;
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

[Ellipses, UseHistoneOverlay, UseSchnitz] = checkHistoneAndNuclearSegmentation(PreProcPath, FilePrefix, nDigits, DropboxFolder, noHisOverlay, fish);

[~, ExperimentType, ~, ~, ~, ~, Channel1, Channel2, ~, ~, ~, ~, ~, ...
    nc9, nc10, nc11, nc12, nc13, nc14, ~, Channel3, prophase, metaphase] =...
    getExperimentDataFromMovieDatabase(Prefix, movieDatabase);

Channels = {Channel1{1}, Channel2{1}, Channel3{1}};
coats = getCoatChannel(Channel1, Channel2, Channel3);

ch = coats(1); %assumes the experiment is _not_ 2spot2color
maxTimeCell = [];
if preMovie
    if isempty(movieMat)
        [movieMat, hisMat, maxMat, ~, ~]...
            = makeMovieMats(Prefix, PreProcPath, nWorkers, FrameInfo, Channels);
        movieMat = squeeze(movieMat(ch,:,:,:,:));
        maxMat = squeeze(maxMat(ch,:,:,:));
    else
        maxMat = squeeze(max(movieMat(:,:,:,:), [],1));
    end
else
    movieMat = []; hisMat = []; maxMat = []; 
end



for i = 1:nFrames
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
ElapsedTime = getFrameElapsedTime(FrameInfo, nFrames);

ncFrames = [nc9, nc10, nc11, nc12, nc13, nc14]; %anaphases
ncFramesFull = [zeros(1,8), ncFrames, nFrames]; %more useful for some things
[anaphaseInMins, prophaseInMins, metaphaseInMins] = getPhasesDurationsInMins(ncFrames, prophase, metaphase, ElapsedTime);

try
    correspondingNCInfo = [FrameInfo.nc]; % the assigned nc of the frames
end

% this save() is here so that a user will still get an updated
% frameinfo.mat even if they abort checkparticletracking without saving
% (to prevent issues with compileparticles)
save([DataFolder, filesep, 'FrameInfo.mat'], 'FrameInfo');


[Particles] = sortParticles(sortByFrames, sortByLength, NChannels, Particles);

%Some flags and initial parameters
ShowThreshold2 = 1; %Whether to show particles below the threshold
ParticleToFollow = [];
CurrentFrameWithinParticle = 1;

cptState = CPTState(Spots, Particles, SpotFilter, schnitzcells, Ellipses, FrameInfo, UseHistoneOverlay, nWorkers, plot3DGauss, projectionMode); %work in progress, 2019-12, JP.

if ~isempty(cptState.Particles{cptState.CurrentChannel})
    cptState.CurrentFrame = cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Frame(CurrentFrameWithinParticle);
else
    error('Looks like the Particles structure is empty. There''s nothing to check.');
end

ZoomRange = 50;
snipImageHandle = [];
CurrentSnippet = [];
oim = [];
ImageHandle = [];
ellipseHandles = {};
spotHandles = {};
ellipseHisHandles = {};



% Changing the intial frames and particle if justNC13
if ncRange
    
    if strcmpi('nc15', endNC)
        lastNCFrame = nFrames;
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
    zTraceAxes,HisOverlayFig,HisOverlayFigAxes] = checkParticleTracking_drawGUI(cptState.UseHistoneOverlay, fish, cptState.plot3DGauss, ExperimentType);

if multiView
    multiFig = figure;
    blankImage = zeros(cptState.FrameInfo(1).LinesPerFrame, cptState.FrameInfo(1).PixelsPerLine);
end


set(0, 'CurrentFigure', Overlay);

frameChangeKeyInput = FrameChangeEventHandler(cptState);
zSliceChangeKeyInput = ZSliceChangeEventHandler(cptState);
particleChangeKeyInput = ParticleChangeEventHandler(cptState);
channelSwitchKeyInput = ChannelSwitchEventHandler(cptState, NChannels, cptState.UseHistoneOverlay);
zoomParticleToggleKeyInput = ZoomParticleToggleEventHandler(cptState);
zoomAnywhereKeyInput = ZoomAnywhereEventHandler(cptState);
histoneContrastKeyInput = HistoneContrastChangeEventHandler(cptState);
addSpotKeyInput = AddSpotEventHandler(cptState, PreProcPath, ProcPath, Prefix);
deleteSpotKeyInput = DeleteSpotEventHandler(cptState);
ellipsesKeyInput = EllipsesEventHandler(cptState);
tracesKeyInput = TracesEventHandler(cptState);
nuclearTrackingKeyInput = NuclearTrackingEventHandler(cptState);
generalKeyInput = GeneralEventHandler(cptState, DataFolder, DropboxFolder, FilePrefix, NChannels);

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

%% Main loop - start
while (cc ~= 'x')
    
    cptState.coatChannel = getCoatChannel(Channel1, Channel2, Channel3);
    if length(cptState.coatChannel) == 1
        cptState.nameSuffix = ['_ch', iIndex(cptState.coatChannel, 2)];
    else
        cptState.nameSuffix = ['_ch', iIndex(cptState.currentChannel, 2)];
    end
    
    inds = find(cptState.CurrentFrame > ncFramesFull);
    currentNC = inds(end);
    
    %Get the coordinates of all the spots in this frame
    [x, y, z] = SpotsXYZ(cptState.Spots{cptState.CurrentChannel}(cptState.CurrentFrame));
    
    %If the approved field does not exist create it
    if ~isfield(cptState.Particles{cptState.CurrentChannel}, 'Approved')
        
        for i = 1:cptState.numParticles()
            cptState.Particles{cptState.CurrentChannel}(i).Approved = 0;
        end
        
    end
    
    %Pull out the right particle if it exists in this frame
    cptState.CurrentParticleIndex = ...
        cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Index(cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Frame == ...
        cptState.CurrentFrame);
    %This is the position of the current particle
    xTrace = x(cptState.CurrentParticleIndex);
    yTrace = y(cptState.CurrentParticleIndex);
    
    if (~isempty(xTrace)) && (~cptState.ManualZFlag)
        cptState.CurrentZ = z(cptState.CurrentParticleIndex);
        CurrentZIndex = find(...
            cptState.Spots{cptState.CurrentChannel}(cptState.CurrentFrame).Fits(cptState.CurrentParticleIndex).z == ...
            cptState.CurrentZ);
        cptState.ManualZFlag = 0;
    end
    
    multiImage = {};
    if strcmpi(cptState.projectionMode, 'None')
        if ~preMovie
            cptState.ImageMat = imread([PreProcPath, filesep, FilePrefix(1:end - 1), filesep, ...
                FilePrefix, iIndex(cptState.CurrentFrame, nDigits), '_z', iIndex(cptState.CurrentZ, 2), cptState.nameSuffix, '.tif']);
        else
            cptState.ImageMat = squeeze(movieMat(cptState.CurrentZ, cptState.CurrentFrame, :, :));
        end
        if multiView
%             if cptState.CurrentParticle == 2
                1
%             end
            for z = 1:-1:-1
                for f = -1:1
                    if any( 1:nSlices == cptState.CurrentZ + z) && any( 1:nFrames == cptState.CurrentFrame + f)
                        multiImage{z+2, f+2} = squeeze(movieMat(cptState.CurrentZ+z, cptState.CurrentFrame+f,:,:));
                    else
                        multiImage{z+2, f+2} = blankImage;
                    end
                end
            end
        end
    elseif strcmpi(cptState.projectionMode, 'Max Z')
        if preMovie
            if nFrames > 1
                cptState.ImageMat = squeeze(maxMat(cptState.CurrentFrame,:,:));
            else
                cptState.ImageMat = maxMat;
            end
        else
            cptState.ImageMat = zProjections(Prefix, cptState.coatChannel, cptState.CurrentFrame, cptState.ZSlices, nDigits, DropboxFolder, PreProcPath, cptState.FrameInfo, 'max', cptState.nWorkers);
        end
    elseif strcmpi(cptState.projectionMode, 'Max Z and Time')
        if preMovie
            if isempty(maxTimeCell)
                cptState.ImageMat = squeeze(max(max(movieMat(:,ncFramesFull(currentNC):ncFramesFull(currentNC+1),:,:), [], 3), [], 2)); % ch z t x y
            end
        else
            if isempty(storedTimeProjection)
                
                if ncRange
                    cptState.ImageMat = timeProjection(Prefix, cptState.coatChannel, cptState.FrameInfo, DropboxFolder,PreProcPath, 'nc', NC);
                    storedTimeProjection = cptState.ImageMat;
                else
                    cptState.ImageMat = timeProjection(Prefix, cptState.CurrentChannel, cptState.FrameInfo, DropboxFolder,PreProcPath);
                    storedTimeProjection = cptState.ImageMat;
                end
                
            else
                cptState.ImageMat = storedTimeProjection;
            end
            
        end
    end
    
    if multiView && ~exist('subAx', 'var')
        tiles = tiledlayout(multiFig, 3, 3, 'TileSpacing', 'none', 'Padding', 'none');
        subAx = cell(3);
        n = 0;
        for z = 1:3
            for f = 1:3
                n = n + 1;
                subAx{z, f} = nexttile(tiles, n);
                title(subAx{z,f},['z: ', num2str(cptState.CurrentZ + z - 2), ' frame: ', num2str(cptState.CurrentFrame + f - 2)])
            end
        end
    end
    
    if ~multiView
        subAx = {};
    end
    
    
    set(0, 'CurrentFigure', Overlay);
    %     if isempty(hIm)
    ImageHandle = imshow(cptState.ImageMat, cptState.DisplayRangeSpot, 'Border', 'Tight', 'Parent', overlayAxes, ...
        'InitialMagnification', 'fit');
    if multiView
        for z = 1:size(multiImage, 1)
            for f= 1:size(multiImage, 2)
                if ~isempty(subAx{z,f}.Children)
                    subAx{z,f}.Children.CData = multiImage{z, f};
                else
                    imshow(multiImage{z, f}, cptState.DisplayRangeSpot, 'Border', 'Tight', 'Parent', subAx{z,f},...
                        'InitialMagnification', 'fit');
                end
                title(subAx{z,f},['z: ', num2str(cptState.CurrentZ + z - 2), ' frame: ', num2str(cptState.CurrentFrame + f - 2)])
            end
        end
    end
    %     else
    %         hIm.CData = Image;
    %     end
    %
    hold(overlayAxes, 'on')
    
    if cptState.UseHistoneOverlay
        if ~preMovie
            HisPath = [PreProcPath, filesep, FilePrefix(1:end - 1), filesep, ...
                FilePrefix(1:end - 1), '-His_', iIndex(cptState.CurrentFrame, nDigits), '.tif'];
            hisImage=imread(HisPath);
        else
            hisImage = squeeze(hisMat(cptState.CurrentFrame, :,:));
        end
        [cptState.ImageHis, cptState.xForZoom, cptState.yForZoom, oim,ellipseHandles]= displayOverlays(overlayAxes, cptState.ImageMat, SpeedMode, cptState.FrameInfo, cptState.Particles, ...
            cptState.Spots, cptState.CurrentFrame, ShowThreshold2, ...
            Overlay, cptState.CurrentChannel, cptState.CurrentParticle, cptState.ZSlices, cptState.CurrentZ, nFrames, ...
            cptState.schnitzcells, UseSchnitz, cptState.DisplayRange, cptState.Ellipses, cptState.SpotFilter, cptState.ZoomMode, cptState.GlobalZoomMode, ...
            ZoomRange, cptState.xForZoom, cptState.yForZoom, fish, cptState.UseHistoneOverlay, subAx, HisOverlayFigAxes,...
            oim, ellipseHandles, hisImage);
        
    else
        displayOverlays(overlayAxes, cptState.ImageMat, SpeedMode, ...
            cptState.FrameInfo, cptState.Particles, cptState.Spots, cptState.CurrentFrame, ShowThreshold2, ...
            Overlay, cptState.CurrentChannel, cptState.CurrentParticle, cptState.ZSlices, cptState.CurrentZ, nFrames, ...
            cptState.schnitzcells, UseSchnitz, cptState.DisplayRange, cptState.Ellipses, cptState.SpotFilter, ...
            cptState.ZoomMode, cptState.GlobalZoomMode, ZoomRange, cptState.xForZoom, cptState.yForZoom, fish, cptState.UseHistoneOverlay, subAx);
    end
    
    if ~isempty(xTrace)
        MaxZIndex = find(...
            cptState.Spots{cptState.CurrentChannel}(cptState.CurrentFrame).Fits(cptState.CurrentParticleIndex).z == ...
            cptState.Spots{cptState.CurrentChannel}(cptState.CurrentFrame).Fits(cptState.CurrentParticleIndex).brightestZ);
        CurrentZIndex = find(...
            cptState.Spots{cptState.CurrentChannel}(cptState.CurrentFrame).Fits(cptState.CurrentParticleIndex).z == ...
            cptState.CurrentZ);
        
        if isempty(CurrentZIndex)
            %             warning('This particle has a gap in its z-profile. This is
            %             highly suspect.'); %this if statement should only happen
            %             between two spots, not past the PSF boundaries
        end
        
    end
    
    %Check to see if spots structure contains multi-slice fields
    multi_slice_flag = isfield(cptState.Spots{cptState.CurrentChannel}(cptState.CurrentFrame).Fits ...
        (cptState.CurrentParticleIndex), 'IntegralZ');
    
    % PLOT SNIPPET
    
    [CurrentSnippet, snipImageHandle] = plotSnippet(snippetFigAxes, rawDataAxes, gaussianAxes, xTrace, ...
        CurrentZIndex, cptState.ImageMat, cptState.Spots, cptState.CurrentChannel, cptState.CurrentFrame, ...
        cptState.CurrentParticleIndex, ExperimentType, snippet_size, xSize, ...
        ySize, SnippetEdge, cptState.FrameInfo, CurrentSnippet, snipImageHandle, pixelSize);
    
    
    % PLOTS TRACE OF CURRENT PARTICLE
    if ~fish
        plottrace_argin = {};
        if exist('AmpIntegral', 'var')
            plottrace_argin = [plottrace_argin, AmpIntegral, GaussIntegral, AmpIntegral3,  ...
                ErrorIntegral, ErrorIntegral3, backGround3, ...
                AmpIntegralGauss3D, ErrorIntegralGauss3D, cptState.FrameIndicesToFit];
        end
        if ~isempty(Spots3D)
            plottrace_argin = [plottrace_argin, Spots3D];
        end
        [cptState.Frames, AmpIntegral, GaussIntegral, AmpIntegral3, ...
            ErrorIntegral, ErrorIntegral3,  backGround3, ...
            AmpIntegralGauss3D, ErrorIntegralGauss3D, cptState.PreviousParticle] =...
            ...
            plotTrace(traceFigAxes, ...
            ...
            cptState.FrameInfo, cptState.CurrentChannel, cptState.PreviousChannel, ...
            cptState.CurrentParticle, cptState.PreviousParticle, cptState.lastParticle, cptState.HideApprovedFlag, cptState.lineFitted, anaphaseInMins, ...
            ElapsedTime, cptState.schnitzcells, cptState.Particles, cptState.plot3DGauss, ncFrames, prophase, metaphase, prophaseInMins, metaphaseInMins, Prefix, ...
            nFrames, cptState.CurrentFrame, cptState.ZSlices, cptState.CurrentZ, cptState.Spots, ...
            correspondingNCInfo, cptState.Coefficients, ExperimentType,cptState.PreviousFrame, cptState.Frames,...
            Channels, PreProcPath, DropboxFolder, plottrace_argin{:});
    end
    
    
    % PLOT Z SLICE RELATED FIGURES
    plotzvars = {zProfileFigAxes, zTraceAxes, ExperimentType, ...
        xTrace, cptState.Spots, cptState.CurrentFrame, cptState.CurrentChannel, cptState.CurrentParticleIndex, cptState.ZSlices, ...
        cptState.CurrentZ, CurrentZIndex, cptState.PreviousParticle, cptState.CurrentParticle, ...
        cptState.PreviousChannel, cptState.Particles, cptState.Frames, fish};
    if exist('MaxZProfile', 'var')
        plotzvars = [plotzvars, MaxZProfile];
    end
    [MaxZProfile, cptState.Frames] = plotZFigures(plotzvars{:});

    set(0, 'CurrentFigure', Overlay);

    % Wait for user input to select command to execute
    ct = waitforbuttonpress; % ct=0 for click and ct=1 for keypress
    cc = get(Overlay, 'CurrentCharacter');
    
    if strcmpi(cc, '') || ct == 0
        cc = 'donothing';
    end
    
    frameChangeKeyInput(cc);
    zSliceChangeKeyInput(cc);
    particleChangeKeyInput(cc);
    channelSwitchKeyInput(cc);
    zoomParticleToggleKeyInput(cc);
    zoomAnywhereKeyInput(cc);
    histoneContrastKeyInput(cc);
    addSpotKeyInput(cc);
    deleteSpotKeyInput(cc);
    ellipsesKeyInput(cc);
    tracesKeyInput(cc);
    nuclearTrackingKeyInput(cc);
    generalKeyInput(cc);
    
end
%% Main loop - end

% save after exiting the main loop - the user pressed 'x'
saveChanges(NChannels, cptState.Particles, cptState.Spots, cptState.SpotFilter, DataFolder, ...
                cptState.FrameInfo, cptState.UseHistoneOverlay, FilePrefix, ...
                cptState.schnitzcells, DropboxFolder);

close all

disp(['(Left off at Particle #', num2str(cptState.CurrentParticle), ')'])

outs = {cptState.Particles, cptState.Spots, cptState.SpotFilter, cptState.schnitzcells, cptState.FrameInfo};

CheckNucleiModified(cptState, DropboxFolder, Prefix, fish)

end
