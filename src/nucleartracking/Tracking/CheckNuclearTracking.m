function CheckNuclearTracking(Prefix, varargin)
% Adapted directly from CheckParticleTracking
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
% nc, NC : Only look at particles that show up in nc13
%    % Currently this only starts at nc13...not restrict you to nc13 Added by Emma
%    % Also, this option shows you the max projection.
% nWorkers: number of parallel pools

%
% CONTROLS
% Frame specific:
% ./, :Move a frame forward/backward
% >/< :Move five frames forward/backward
% ;/' :Move to the next/previous empty frame within a nucleus
% a/z :Move up/down in z-slice
%  t  :Jump to a specific z-slice
%  j  :Jump to a specified frame
% g/b :Increase/decrease histone channel contrast
% !/@ :Change the contrast in protein channel (! increases/@ resets
%       it back to the lowest)
%
% Particle specific:
% m :Move to the next nucleus
% n :Move to the previous nucleus
% k :Jump to a specified nucleus by inputting a nucleus number
% \ :Jump to a specified nucleus by clicking
% c :Connect two existing nuclear traces. This will join the current
%     nucleus's trace to the clicked nucleus's trace.
% d :Separate traces forward. A new nucleus is created at the current frame
%     and this nucleus is disconnected from the current nucleus.
% q :Cycle between approved status: green - approved; yellow - approved but
%     with conditions (drift of nucleus, for example)
% w :Disapprove a trace
% p :Identify a particle. It will also tell you the particle associated with
%     the clicked nucleus.
% e :Approve/Disapprove a frame within a trace
% [ :Add a spot that was not recognized originally by segmentSpots, creating
%     a new particle if you've used '+' to zoom in on the particle of
%     interest or adding to the current particle if you are zoomed out in
%     'o'. Note that the command forces ZoomMode. To toggle, use 'o' or '+'
%     depending on whether you're adding to an existing trace or creating
%     a new trace, respectively.
% { :Same as "[" but uses the exact pixel and z-plane that you click on.
%     Useful if the algorithms get the centroid positioning wrong (i.e. your
%     particle is put in the wrong place by '[').
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
%none
% Author (contact): Gabriella Martini (martini@berkeley.edu)
% Created: 9/7/2020
% Last Updated: 9/7/2020

cleanupObj = onCleanup(@myCleanupFun);
warning('off', 'MATLAB:nargchk:deprecated')
warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary')

%% Initialization

liveExperiment = LiveExperiment(Prefix);

schnitzcells = [];
Ellipses = [];
correspondingNCInfo = [];
IntegrationArea = []; %Initialized here to avoid dynamic assignment later in function
fish = false;
% Parameters for fitting
FramesToFit = []; % actual frames of the movie that were used for fitting
SnippetEdge = 13; %Size of the snippets generated by Michael Tikhonov's code in pixels.


[sortByFrames, sortByLength,...
    ncRange, projectionMode, NC, ...
    startNC, endNC, optionalResults, nWorkers,...
    noHisOverlay,preStructs] = determineCheckNuclearTrackingOptions(varargin{:});


%% Information about about folders

% Get the folders
[~,~,DropboxFolder,~, ~,~, ~, ~,~,~,~,~, ~, ~, movieDatabase]...
    = readMovieDatabase(Prefix, optionalResults);

PreProcPath = liveExperiment.userPreFolder;

DataFolder = [DropboxFolder, filesep, Prefix];
FilePrefix = [Prefix, '_'];


if isempty(preStructs)
    [FrameInfo, schnitzcells] =...
        loadCheckNuclearTrackingMats(DataFolder, PreProcPath, FilePrefix);
else
    FrameInfo = preStructs{1};
    schnitzcells = preStructs{2};
end

xSize = liveExperiment.xDim;
ySize = liveExperiment.yDim;
pixelSize_nm = liveExperiment.pixelSize_nm;
snippetSize_px = liveExperiment.snippetSize_px;
nFrames = liveExperiment.nFrames;
nSlices = liveExperiment.zDim;
nDigits = liveExperiment.nDigits;

nSlices = nSlices + 2; %due to padding;
%Create the particle array. This is done so that we can support multiple
%channels. Also figure out the number of channels

numInputChannels = 1;

[Ellipses, UseHistoneOverlay, UseSchnitz] =...
    checkHistoneAndNuclearSegmentation(...
    PreProcPath, FilePrefix, nDigits, DropboxFolder, noHisOverlay, false);

[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ...
    ~, ~, ~, ~, ~, ~, ~, ~, prophase, metaphase] =...
    getExperimentDataFromMovieDatabase(Prefix, movieDatabase);

ExperimentType = liveExperiment.experimentType;
Channels = liveExperiment.Channels;
Channel1 = liveExperiment.Channel1;
Channel2 = liveExperiment.Channel2;
Channel3 = liveExperiment.Channel3;
nc9 = liveExperiment.nc9;
nc10 = liveExperiment.nc10;
nc11 = liveExperiment.nc11;
nc12 = liveExperiment.nc12;
nc13 = liveExperiment.nc13;
nc14 = liveExperiment.nc14;


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
save([DataFolder, filesep, 'FrameInfo.mat'], 'FrameInfo', '-v6');

%Get the actual time corresponding to each frame, in minutes
ElapsedTime = getFrameElapsedTime(FrameInfo, nFrames);

ncFrames = [nc9, nc10, nc11, nc12, nc13, nc14]; %anaphases
ncFramesFull = [zeros(1,8), ncFrames, nFrames]; %more useful for some things
[anaphaseInMins, prophaseInMins, metaphaseInMins] =...
    getPhasesDurationsInMins(ncFrames, prophase, metaphase, ElapsedTime);

try correspondingNCInfo = [FrameInfo.nc]; end

[schnitzcells] = sortNuclei(sortByFrames, sortByLength, schnitzcells);

%Some flags and initial parameters
NucleusToFollow = [];
CurrentFrameWithinNucleus = 1;

%% 

cntState = CNTState(liveExperiment, schnitzcells, Ellipses,...
                FrameInfo, UseHistoneOverlay, nWorkers, projectionMode);

try
    inputChannels = liveExperiment.inputChannels; 
    cntState.CurrentChannel = spotChannels(1);
    cntState.CurrentChannelIndex= 1;
end


if ~isempty(cntState.schnitzcells)
    cntState.CurrentFrame =...
        cntState.schnitzcells(cntState.CurrentNucleus).frames(CurrentFrameWithinNucleus);

else, error('Looks like the schnitz cells structure is empty. There''s nothing to check.'); end

%load the movies
movieMat = getMovieMat(liveExperiment);  


hisMat = getHisMat(liveExperiment);
maxMat = getMaxMat(liveExperiment);

ZoomRange = 50;
snipImageHandle = [];
CurrentSnippet = [];
ImageHandle = [];
spotHandles = {};
ellipseHisHandles = {};



% Changing the intial frames and particle if justNC13
if ncRange
    
    if strcmpi('nc15', endNC), lastNCFrame = nFrames;      
    else lastNCFrame = eval(endNC) - 1; end
    
    firstNCFrame = eval(startNC);
    schnitzcellsInRange = schnitzcellsWithinFrames(Prefix, firstNCFrame, lastNCFrame);
    cntState.CurrentNucleus = schnitzcellsInRange(1);
    cntState.CurrentFrame = cntState.schnitzcells(cntState.CurrentNucleus).frames(1);
    disp(['nc range: ' num2str(NC)])
    disp(['start frame: ' num2str(firstNCFrame)])
    disp(['end frame: ' num2str(lastNCFrame)])
    disp(['Nuclei in range: ' num2str(schnitzcellsInRange)])
    disp(['Number of Nuclei: ' num2str(length(schnitzcellsInRange))])
    
end

%Define user interface
[Overlay, overlayAxes, snippetFigAxes,...
    rawDataAxes, gaussianAxes, traceFig, traceFigAxes, zProfileFigAxes,...
    zTraceAxes,HisOverlayFig,HisOverlayFigAxes, multiFig]...
    ...
    = checkNuclearTracking_drawGUI(...
    cntState.UseHistoneOverlay, ...
     ExperimentType, xSize, ySize);


blankImage = []; % dummy, to enable calling CTPState for this later



set(0, 'CurrentFigure', Overlay);

frameChangeKeyInput = NuclearFrameChangeEventHandler(cntState);
zSliceChangeKeyInput = NuclearZSliceChangeEventHandler(cntState);
particleChangeKeyInput = NuclearChangeEventHandler(cntState);
%channelSwitchKeyInput = ChannelSwitchEventHandler(cptState, numSpotChannels, cptState.UseHistoneOverlay);
%zoomParticleToggleKeyInput = ZoomParticleToggleEventHandler(cptState);
%zoomAnywhereKeyInput = ZoomAnywhereEventHandler(cptState);
%histoneContrastKeyInput = HistoneContrastChangeEventHandler(cptState);
%addSpotKeyInput = AddSpotEventHandler(cptState, Prefix);
%deleteSpotKeyInput = DeleteSpotEventHandler(cptState);
ellipsesKeyInput = NuclearEllipsesEventHandler(cntState);
tracesKeyInput = NuclearTracesEventHandler(cntState);
%nuclearTrackingKeyInput = NuclearTrackingEventHandler(cptState);
generalKeyInput = GeneralNuclerEventHandler(cntState, DataFolder, DropboxFolder, FilePrefix);

% Create the approved field if it does not exist

if ~isfield(cntState.schnitzcells, 'Approved')
    for i = 1:length(cntState.schnitzcells)
        cntState.schnitzcells(i).Approved = 0;
    end

end


currentCharacter = 1;

if ForCompileAll
    %we just want to save the data, we set 'x' as user command
    currentCharacter = 'x';
end

plotTraceSettings = PlotTraceSettings();

%% Main loop - start
while (currentCharacter ~= 'x')
    
    cntState.inputChannel = getInputChannel(Channel1, Channel2, Channel3);
   
    cntState.nameSuffix = ['_ch', iIndex(cntState.CurrentChannel, 2)];
    
    inds = find(cntState.CurrentFrame > ncFramesFull);
    currentNC = inds(end);
    
    %Get the coordinates of all the nuclei in this frame
    [x, y, z] = NucleiXYZ(cntState.schnitzcells, cntState);
    
    %If the approved field does not exist create it
    if ~isfield(cntState.schnitzcells, 'Approved')
        
        for i = 1:cntState.numNuclei()
            cntState.schnitzcells(i).Approved = 0;
        end
        
    end
    
    %Pull out the right particle if it exists in this frame
    cntState.updateCurrentNucleusCellNo();

    %This is the position of the current particle
    [xTrace, yTrace] = cntState.getXYTraces(x, y);
    
    cntState.updateZIndex(x, y, z);
    
    % GM: STOPPED HERE
    cntState.processImageMatrices(nFrames, nSlices,...
         blankImage, currentNC,...
            ncFramesFull, movieMat,...
            maxMat);
    
    if multiView && ~exist('subAx', 'var')
        tiles = tiledlayout(multiFig, 3, 3, 'TileSpacing', 'none', 'Padding', 'none');
        subAx = cell(3);
        n = 0;
        for z = 1:3
            for f = 1:3
                n = n + 1;
                subAx{z, f} = nexttile(tiles, n);
                title(subAx{z,f},['z: ', num2str(cptState.CurrentZ + z - 2),...
                    ' frame: ', num2str(cptState.CurrentFrame + f - 2)])
            end
        end
    end
    
    if ~multiView
        subAx = {};
    end
    
    
    set(0, 'CurrentFigure', Overlay);
    %     if isempty(hIm)
    ImageHandle = imshow(cptState.ImageMat,...
        cptState.DisplayRangeSpot, 'Border', 'Tight', 'Parent', overlayAxes, ...
        'InitialMagnification', 'fit');
    if multiView
        displayRangeSpot = cptState.DisplayRangeSpot;
        if isempty(displayRangeSpot)
            displayRangeSpot = [median(cptState.ImageMat(:)), max(cptState.ImageMat(:))];
            cptState.DisplayRangeSpot = displayRangeSpot;
        end
        for z = 1:size(cptState.multiImage, 1)
            for f= 1:size(cptState.multiImage, 2)
                if ~isempty(subAx{z,f}.Children)
                    subAx{z,f}.Children.CData = cptState.multiImage{z, f};
                    subAx{z, f}.CLim = displayRangeSpot;                   
                else
                    imshow(cptState.multiImage{z, f},...
                        displayRangeSpot, 'Border', 'Tight', 'Parent', subAx{z,f},...
                        'InitialMagnification', 'fit');
                end
                title(subAx{z,f},['z: ',...
                    num2str(cptState.CurrentZ + z - 2), ' frame: ',...
                    num2str(cptState.CurrentFrame + f - 2)])
            end
        end
    end
    hold(overlayAxes, 'on')
    
    if cptState.UseHistoneOverlay
        hisImage = hisMat(:, :, cptState.CurrentFrame);

        displayOverlays(overlayAxes, cptState, SpeedMode, ShowThreshold2, Overlay, nFrames, UseSchnitz, ZoomRange, fish, subAx,...
            HisOverlayFigAxes, hisImage);
    else
        displayOverlays(overlayAxes, cptState, SpeedMode, ShowThreshold2, Overlay, nFrames, UseSchnitz, ZoomRange, fish, subAx);
    end   

    if ~isempty(xTrace)
        MaxZIndex = cptState.getMaxZIndex();
        cptState.updateCurrentZIndex()
    end
    
     
    % PLOT SNIPPET
    [CurrentSnippet, snipImageHandle] = plotSnippet(snippetFigAxes, rawDataAxes, gaussianAxes, xTrace, ...
        cptState, ExperimentType, snippetSize_px, xSize, ySize, SnippetEdge, CurrentSnippet, snipImageHandle, pixelSize_nm);
    
    % PLOTS TRACE OF CURRENT PARTICLE
    if ~fish
        plotTrace(traceFigAxes, cptState, anaphaseInMins, ElapsedTime,...
            ncFrames, prophase, metaphase, prophaseInMins,...
            metaphaseInMins, Prefix, nFrames, correspondingNCInfo,...
            ExperimentType, Channels, PreProcPath, DropboxFolder,...
            plotTraceSettings);
    end
    
%%
% %AR- disabled until it's fast enough to be useful.
% %too slow at present.
% 
%     % PLOT Z SLICE RELATED FIGURES
%     plotzvars = {zProfileFigAxes, zTraceAxes, ExperimentType,...
%         xTrace, cptState, plotTraceSettings, fish};
%     if exist('MaxZProfile', 'var')
%         plotzvars = [plotzvars, MaxZProfile];
%     end
%    MaxZProfile = plotZFigures(plotzvars{:});
%%
    set(0, 'CurrentFigure', Overlay);
    
    currentCharacter = getUserKeyInput(Overlay);
    
    frameChangeKeyInput(currentCharacter);
    zSliceChangeKeyInput(currentCharacter);
    particleChangeKeyInput(currentCharacter);
    channelSwitchKeyInput(currentCharacter);
    zoomParticleToggleKeyInput(currentCharacter);
    zoomAnywhereKeyInput(currentCharacter);
    histoneContrastKeyInput(currentCharacter);
    addSpotKeyInput(currentCharacter);
    deleteSpotKeyInput(currentCharacter);
    ellipsesKeyInput(currentCharacter);
    tracesKeyInput(currentCharacter);
    nuclearTrackingKeyInput(currentCharacter);
    generalKeyInput(currentCharacter);
    
end
%% Main loop - end

% save after exiting the main loop - the user pressed 'x'

saveChanges(numSpotChannels, cptState, DataFolder, FilePrefix, DropboxFolder);

disp(['(Left off at Particle #', num2str(cptState.CurrentParticle), ')'])

CheckNucleiModified(cptState, DropboxFolder, Prefix, fish)

end