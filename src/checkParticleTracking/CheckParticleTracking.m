function [outs, movieCell, hisCell] = CheckParticleTracking(Prefix, varargin)
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
lineFitted = 0; % equals 1 if a line has been fitted
FramesToFit = []; % actual frames of the movie that were used for fitting

SnippetEdge = 13; %Size of the snippets generated by Michael's code in pixels.
storedTimeProjection = []; % Don't need to wait for timeProjection to finish each time its called

[sortByFrames, sortByLength, ForCompileAll, SpeedMode, ~, ...
    ncRange, projectionMode, plot3DGauss, NC, ...
    startNC, endNC, optionalResults, nWorkers, fish,...
    noHisOverlay, multiView, preStructs, preMovie, movieCell, hisCell] = determineCheckParticleTrackingOptions(varargin{:});

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


[xSize, ySize, pixelSize, zStep, snippet_size, LinesPerFrame, PixelsPerLine,...
    nFrames, nSlices] = getFrameInfoParams(FrameInfo);

%See how  many frames we have and adjust the index size of the files to load accordingly
if nFrames < 1E3
    NDigits = 3;
elseif nFrames < 1E4
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

[Ellipses, UseHistoneOverlay, UseSchnitz] = checkHistoneAndNuclearSegmentation(PreProcPath, FilePrefix, NDigits, DropboxFolder, noHisOverlay, fish);

[~, ExperimentType, ~, ~, ~, ~, Channel1, Channel2, ~, ~, ~, ~, ~, ...
    nc9, nc10, nc11, nc12, nc13, nc14, ~, Channel3, prophase, metaphase] =...
    getExperimentDataFromMovieDatabase(Prefix, movieDatabase);

Channels = {Channel1{1}, Channel2{1}, Channel3{1}};

nPadding = 2; %normally we pad a blank image above and below the stack.
%if a movie is unpadded this will require modification
maxCell = [];
if preMovie
    if isempty(movieCell)
        coats = getCoatChannel(Channel1, Channel2, Channel3);
        movieCell = zeros(length(coats), nSlices, nFrames, xSize, ySize, 'uint16'); % ch z t x y
        if UseSchnitz
            hisCell =  zeros(nFrames, xSize, ySize, 'uint16'); %t x y
        end
        pth = [PreProcPath, filesep, Prefix, filesep,Prefix];
        for ch = 1:length(coats)
            c = coats(ch);
            for f = 1:nFrames
                parfor z = 1:nSlices+nPadding
                    %reminder to squeeze the array before accessing an image
                    movieCell(ch, z, f, :, :) = imread([pth,'_',iIndex(f, NDigits), '_z', iIndex(z, 2), ['_ch', iIndex(c, 2)], '.tif']);
                end
                if UseSchnitz
                    hisCell(f, :, :) = imread([pth,'-His_', iIndex(f, NDigits), '.tif']);
                end
            end
        end
        maxCell = squeeze(max(movieCell(:,:,:,:, :), [], 2));
    end
else
    movieCell = [];
    hisCell = [];
    maxCell = [];
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

anaphase = [nc9, nc10, nc11, nc12, nc13, nc14];
[anaphaseInMins, prophaseInMins, metaphaseInMins] = getPhasesDurationsInMins(anaphase, prophase, metaphase, ElapsedTime);

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
HideApprovedFlag = 0;
ParticleToFollow = [];
CurrentFrameWithinParticle = 1;

cptState = CPTState(Spots, Particles, SpotFilter, schnitzcells, Ellipses, FrameInfo, UseHistoneOverlay, nWorkers, plot3DGauss); %work in progress, 2019-12, JP.

if ~isempty(cptState.Particles{cptState.CurrentChannel})
    cptState.CurrentFrame = cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Frame(CurrentFrameWithinParticle);
else
    error('Looks like the Particles structure is empty. There''s nothing to check.');
end
DisplayRangeSpot = [];
ZoomRange = 50;
snipImageHandle = [];
CurrentSnippet = [];
oim = [];
ImageHandle = [];
ImageMat = [];
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
end


set(0, 'CurrentFigure', Overlay);

coatChannels = [1, 2]; % JP temporary, will be used only if 2spot2color, could be refactored into cptState

frameChangeKeyInput = FrameChangeEventHandler(cptState);
zSliceChangeKeyInput = ZSliceChangeEventHandler(cptState);
particleChangeKeyInput = ParticleChangeEventHandler(cptState);
channelSwitchKeyInput = ChannelSwitchEventHandler(cptState, NChannels, coatChannels, cptState.UseHistoneOverlay);
zoomParticleToggleKeyInput = ZoomParticleToggleEventHandler(cptState);
zoomAnywhereKeyInput = ZoomAnywhereEventHandler(cptState);
histoneContrastKeyInput = HistoneContrastChangeEventHandler(cptState);
addSpotKeyInput = AddSpotEventHandler(cptState, PreProcPath, ProcPath, Prefix);
deleteSpotKeyInput = DeleteSpotEventHandler(cptState);
ellipsesKeyInput = EllipsesEventHandler(cptState);

AveragingLength = 1; % Default

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
    
    cptState.coatChannel = getCoatChannel(Channel1, Channel2, Channel3);
    if length(cptState.coatChannel) == 1
        cptState.nameSuffix = ['_ch', iIndex(cptState.coatChannel, 2)];
    else
        cptState.nameSuffix = ['_ch', iIndex(cptState.currentChannel, 2)];
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
    if strcmpi(projectionMode, 'None')
        if ~preMovie
            ImageMat = imread([PreProcPath, filesep, FilePrefix(1:end - 1), filesep, ...
                FilePrefix, iIndex(cptState.CurrentFrame, NDigits), '_z', iIndex(cptState.CurrentZ, 2), cptState.nameSuffix, '.tif']);
        else
            ImageMat = squeeze(movieCell(cptState.CurrentChannel,cptState.CurrentZ, cptState.CurrentFrame, :, :));
        end
        if multiView
            for i = 1:-1:-1
                for j = -1:1
                    if any( 1:nSlices == cptState.CurrentZ + i) && any( 1:nFrames == cptState.CurrentFrame + j)
                        multiImage{i+2, j+2} = imread([PreProcPath, filesep, FilePrefix(1:end - 1), filesep, ...
                            FilePrefix, iIndex(cptState.CurrentFrame+j, NDigits), '_z', iIndex(cptState.CurrentZ+i, 2), cptState.nameSuffix, '.tif']);
                    else
                        multiImage{i+2, j+2} = zeros(cptState.FrameInfo(1).LinesPerFrame, cptState.FrameInfo(1).PixelsPerLine);
                    end
                end
            end
        end
    elseif strcmpi(projectionMode, 'Max Z')
        if preMovie
            if isempty(maxCell)
                maxCell = squeeze(max(movieCell(:,:,:,:, :), [], 2));
            end
            if nFrames > 1
                ImageMat = squeeze(maxCell(cptState.CurrentFrame,:,:));
            else
                ImageMat = maxCell;
            end
        else
            ImageMat = zProjections(Prefix, cptState.coatChannel, cptState.CurrentFrame, cptState.ZSlices, NDigits, DropboxFolder, PreProcPath, cptState.FrameInfo, 'max', cptState.nWorkers);
        end
    elseif strcmpi(projectionMode, 'Median Z')
        ImageMat = zProjections(Prefix, cptState.coatChannel, cptState.CurrentFrame, cptState.ZSlices, NDigits, DropboxFolder, PreProcPath, cptState.FrameInfo, 'median', cptState.nWorkers);
    elseif strcmpi(projectionMode, 'Max Z and Time')
        
        if isempty(storedTimeProjection)
            
            if ncRange
                ImageMat = timeProjection(Prefix, cptState.coatChannel, cptState.FrameInfo, DropboxFolder,PreProcPath, 'nc', NC);
                storedTimeProjection = ImageMat;
            else
                ImageMat = timeProjection(Prefix, cptState.CurrentChannel, cptState.FrameInfo, DropboxFolder,PreProcPath);
                storedTimeProjection = ImageMat;
            end
            
        else
            ImageMat = storedTimeProjection;
        end
        
    end
    
    
    if multiView && ~exist('subAx', 'var')
        tiles = tiledlayout(multiFig, 3, 3, 'TileSpacing', 'none', 'Padding', 'none');
        n = 0;
        for i = 1:3
            for j = 1:3
                n = n + 1;
                subAx{i, j} = nexttile(tiles, n);
                title(subAx{i,j},['z: ', num2str(cptState.CurrentZ + i - 2), ' frame: ', num2str(cptState.CurrentFrame + j - 2)])
            end
        end
    end
    
    if ~multiView
        subAx = {};
    end
    
    
    set(0, 'CurrentFigure', Overlay);
    %     if isempty(hIm)
    ImageHandle = imshow(ImageMat, DisplayRangeSpot, 'Border', 'Tight', 'Parent', overlayAxes, ...
        'InitialMagnification', 'fit');
    if multiView
        for i = 1:size(multiImage, 1)
            for j = 1:size(multiImage, 2)
                if ~isempty(subAx{i,j}.Children)
                    subAx{i,j}.Children.CData = multiImage{i, j};
                else
                    imshow(multiImage{i, j}, DisplayRangeSpot, 'Border', 'Tight', 'Parent', subAx{i,j},...
                        'InitialMagnification', 'fit');
                end
                title(subAx{i,j},['z: ', num2str(cptState.CurrentZ + i - 2), ' frame: ', num2str(cptState.CurrentFrame + j - 2)])
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
                FilePrefix(1:end - 1), '-His_', iIndex(cptState.CurrentFrame, NDigits), '.tif'];
            ImageHisMat=imread(HisPath);
        else
            ImageHisMat = squeeze(hisCell(cptState.CurrentFrame, :,:));
        end
        [cptState.ImageHis, cptState.xForZoom, cptState.yForZoom, oim,ellipseHandles]= displayOverlays(overlayAxes, ImageMat, SpeedMode, cptState.FrameInfo, cptState.Particles, ...
            cptState.Spots, cptState.CurrentFrame, ShowThreshold2, ...
            Overlay, cptState.CurrentChannel, cptState.CurrentParticle, cptState.ZSlices, cptState.CurrentZ, nFrames, ...
            cptState.schnitzcells, UseSchnitz, cptState.DisplayRange, cptState.Ellipses, cptState.SpotFilter, cptState.ZoomMode, cptState.GlobalZoomMode, ...
            ZoomRange, cptState.xForZoom, cptState.yForZoom, fish, cptState.UseHistoneOverlay, subAx, HisOverlayFigAxes,...
            oim, ellipseHandles, ImageHisMat);
        
    else
        displayOverlays(overlayAxes, ImageMat, SpeedMode, ...
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
        CurrentZIndex, ImageMat, cptState.Spots, cptState.CurrentChannel, cptState.CurrentFrame, ...
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
            cptState.CurrentParticle, cptState.PreviousParticle, cptState.lastParticle, HideApprovedFlag, lineFitted, anaphaseInMins, ...
            ElapsedTime, cptState.schnitzcells, cptState.Particles, cptState.plot3DGauss, anaphase, prophase, metaphase, prophaseInMins, metaphaseInMins, Prefix, ...
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
        
    else
        cc = SkipWaitForButtonPress;
        SkipWaitForButtonPress = [];
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
    
    if strcmpi(cc, 'donothing')
        %do nothing

    elseif cc == 'f'
        [cptState.Particles, cptState.schnitzcells] = redoTracking(DataFolder, ...
            cptState.UseHistoneOverlay, cptState.FrameInfo, DropboxFolder, FilePrefix, cptState.schnitzcells, ...
            cptState.Particles, NChannels, cptState.CurrentChannel, cptState.numParticles());
    elseif cc == 'c'
        [cptState.PreviousParticle, cptState.Particles] = combineTraces(cptState.Spots, ...
        
    elseif cc == 'p' % Identify a particle. It will also tell you the particle associated with
        %  the clicked nucleus.
        identifyParticle(cptState.Spots, cptState.Particles, cptState.CurrentFrame, ...
            cptState.CurrentChannel, cptState.UseHistoneOverlay, cptState.schnitzcells);
    elseif cc == '\' %Moves to clicked particle.
        [cptState.CurrentParticle, cptState.CurrentFrame, cptState.ManualZFlag] = toNearestParticle(cptState.Spots, ...
            cptState.Particles, cptState.CurrentFrame, cptState.CurrentChannel, cptState.UseHistoneOverlay, cptState.schnitzcells);
        
    elseif cc == 'u'
        [x2, y2, z2] = SpotsXYZ(cptState.Spots{cptState.CurrentChannel}(cptState.CurrentFrame));
        
        if ~isempty(x2)
            ClickedSpot = ginput(1);
            
            UnfilterSpot(cptState.Spots{cptState.CurrentChannel}, cptState.SpotFilter{cptState.CurrentChannel}, ...
                ClickedSpot, cptState.Particles{cptState.CurrentChannel}, cptState.CurrentFrame)
        end
        
    elseif cc == 'i'
        warning(' AR 1/15/18: This is currently deprecated. Talk to HG if you need this function.')
    elseif cc == 'd' || cc == 'v' %d Separate traces forward at the current frame.
        [cptState.Particles, cptState.PreviousParticle] = separateTraces(cptState.Particles, ...
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
        saveChanges(NChannels, cptState.Particles, cptState.Spots, cptState.SpotFilter, DataFolder, ...
            cptState.FrameInfo, cptState.UseHistoneOverlay, FilePrefix, ...
            cptState.schnitzcells, DropboxFolder);
    elseif cc == 'h'
        
        if HideApprovedFlag == 0
            HideApprovedFlag = 1; %Show only non-approved traces
        elseif HideApprovedFlag == 1
            HideApprovedFlag = 2; %Show only yellow and red traces
        elseif HideApprovedFlag == 2
            HideApprovedFlag = 0;
        end
        
    elseif (cc == 'm') & (cptState.CurrentParticle < cptState.numParticles())
        [lineFitted, cptState.CurrentParticle, cptState.CurrentFrame, cptState.ManualZFlag, cptState.DisplayRange] = ...
            goNextParticle(cptState.CurrentParticle, cptState.CurrentChannel, HideApprovedFlag, cptState.Particles);
        
    elseif (cc == 'n') & (cptState.CurrentParticle > 1)
        [lineFitted, cptState.CurrentParticle, cptState.CurrentFrame, cptState.ManualZFlag, cptState.DisplayRange] = ...
            goPreviousParticle(cptState.CurrentParticle, cptState.CurrentChannel, HideApprovedFlag, cptState.Particles);
        
    elseif cc == 'e'
        cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).FrameApproved(cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Frame == cptState.CurrentFrame) = ...
            ~cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).FrameApproved(cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Frame == cptState.CurrentFrame);
        
        %Schnitzcells specific
        
    elseif cc == 'l' %Split a nucleus and select one or two daughter cells or stop the lineage
        [cptState.Particles, cptState.PreviousParticle, cptState.schnitzcells] = splitNuclei(cptState.schnitzcells, ...
            cptState.CurrentFrame, cptState.CurrentChannel, cptState.CurrentParticle, cptState.Particles);
    elseif cc == '2' %2 set parent of current nucleus
        cptState.schnitzcells = setParentNucleus(cptState.schnitzcells, ...
            cptState.CurrentFrame, cptState.CurrentChannel, cptState.CurrentParticle, cptState.Particles);
    elseif cc == '~' %Switch projection mode
        projectionMode = chooseProjection;
        disp(['projectionMode : ' projectionMode])
        
    elseif cc == '!' %Increase contrast in the Overlay figure
        
        if isempty(DisplayRangeSpot)
            DisplayRangeSpot = [min(min(ImageMat)), max(max(ImageMat)) / 1.5];
        else
            DisplayRangeSpot = [DisplayRangeSpot(1), DisplayRangeSpot(2) / 1.5];
        end
        
        disp('increased spot contrast');
        
    elseif cc == '@' %Decrease spot channel contrast
        DisplayRangeSpot = [min(ImageMat(:)), max(ImageMat(:)) * 1.5];
        
        disp('decreased spot contrast');
        
    elseif cc == '$' %add particle to nucleus
            cptState.Particles = addNucleusToParticle(cptState.Particles, cptState.CurrentFrame, ...
                cptState.CurrentChannel, cptState.UseHistoneOverlay, cptState.schnitzcells, cptState.CurrentParticle);
    elseif cc == '0' %Debugging mode
        keyboard;
        
    end
    
    %% Main loop - end
end

% JP: we should reuse saveChanges function to do all the saving below

FrameInfo = cptState.FrameInfo;
save([DataFolder, filesep, 'FrameInfo.mat'], 'FrameInfo')

%If we only have one channel bring Particles back to the legacy
%format without any cells
if NChannels == 1
    cptState.Particles = cptState.Particles{1};
    cptState.Spots = cptState.Spots{1};
    cptState.SpotFilter = SpotFilter{1};
end

% store Spots, Paticles and SpotFilter as local variables so we can both save and return them
Spots = cptState.Spots;
Particles = cptState.Particles;
SpotFilter = cptState.SpotFilter;
schnitzcells = cptState.schnitzcells;

if cptState.UseHistoneOverlay
    save([DataFolder, filesep, 'Particles.mat'], 'Particles', 'SpotFilter', '-v7.3')
    save([DataFolder, filesep, 'Spots.mat'], 'Spots', '-v7.3')
    save([DropboxFolder, filesep, FilePrefix(1:end - 1), filesep, FilePrefix(1:end - 1), '_lin.mat'], 'schnitzcells')
else
    save([DataFolder, filesep, 'Particles.mat'], 'Particles', 'SpotFilter', '-v7.3')
    save([DataFolder, filesep, 'Spots.mat'], 'Spots', '-v7.3')
end

close all

try
    if ishandle(controls)
        close(controls)
    end
end

disp('Particles saved.')
disp(['(Left off at Particle #', num2str(cptState.CurrentParticle), ')'])

if cptState.nucleiModified
    Ellipses = cptState.Ellipses; % assign it to a local variable so we can save it to its .mat below
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

outs = {Particles, Spots, SpotFilter, schnitzcells, FrameInfo};

end
