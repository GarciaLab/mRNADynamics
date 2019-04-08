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
  % plot3DGauss: plot 3D gaussian fitting intensities in tracefig
  % sortByLength: sort particles by duration instead of time observed
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

  [Prefix, Sort, sortByLength, ForCompileAll, SpeedMode, ~, ...
  ncRange, projectionMode, plot3DGauss, intScale, NC, ...
  startNC, endNC, optionalResults] = determineCheckParticleTrackingOptions(varargin);

  
  %% Information about about folders

  % Get the folders
  [~, ~, DefaultDropboxFolder, ~, PreProcPath] = DetermineLocalFolders;
  [~, ~, DropboxFolder, ~, ~] = DetermineLocalFolders(varargin{1}, optionalResults);
 
  DataFolder = [DropboxFolder, filesep, varargin{1}];


  FilePrefix = [DataFolder(length(DropboxFolder) + 2:end), '_'];

  % Now get the actual folders
  [~, ~, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(FilePrefix(1:end - 1));

  [Particles, SpotFilter, Spots, FrameInfo] = loadCheckParticleTrackingMats(DataFolder, PreProcPath);

  [xSize, ySize, pixelSize, zStep, snippet_size, LinesPerFrame, PixelsPerLine,...
    numFrames, NDigits, NChannels, Particles, Spots, SpotFilter] = getFrameInfoParams(FrameInfo, Particles, Spots, SpotFilter);

  %Add FramesApproved where necessary
  [Particles] = addFrameApproved(NChannels, Particles);

  [Ellipses, UseHistoneOverlay, UseSchnitz] = checkHistoneAndNuclearSegmentation(PreProcPath, FilePrefix, NDigits, DropboxFolder);

  % we name the variable DataFolderColumnValue to avoid shadowing previously defined DataFolder var, which is actually a subfolder inside dropbox
  [Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution, ...
    Channel1, Channel2, Objective, Power, DataFolderColumnValue, DropboxFolderName, Comments, ...
    nc9, nc10, nc11, nc12, nc13, nc14, CF, Channel3, prophase, metaphase] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);

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
  end

  % this save() is here so that a user will still get an updated
  % frameinfo.mat even if they abort checkparticletracking without saving
  % (to prevent issues with compileparticles)
  save([DataFolder, filesep, 'FrameInfo.mat'], 'FrameInfo');

  %Check if we have already determined nc
  if (~isfield(FrameInfo, 'nc')) && (~UseHistoneOverlay)
    %FrameInfo=DetermineNC(fad,Particles,FrameInfo);  AR 3/14/16: This
    %script seems to have disappeared.

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

  [Particles] = sortParticles(Sort, sortByLength, NChannels, Particles);
  
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
  CurrentFrame = Particles{1}(1).Frame(1);
  DisplayRange = [];
  DisplayRangeSpot = [];
  ZoomMode = 0;
  GlobalZoomMode = 0;
  ZoomRange = 50;
  nameSuffix = '';

  %Set up the default contrast settings for the MCP channel depending on the
  %microscope that was used used
  if strcmpi(FrameInfo(1).FileMode, 'dspin')
    %For spinning disk, we set the contrast to the maximum and minimum
    minContrast = [];
    maxContrast = [];
  else
    %For all other microscopes, we have a default. HG is not sure this will
    %actually work well beyond Leica SP8.
    minContrast = 0; % Default contrast settings for gfp channel
    maxContrast = 80;
  end

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
    [numParticles, SpotFilter, Particles, Spots, PreviousParticle] = ...
      addSpot(ZoomMode, GlobalZoomMode, Particles, CurrentChannel, ...
      CurrentParticle, CurrentFrame, CurrentZ, Overlay, snippet_size, PixelsPerLine, ...
      LinesPerFrame, Spots, ZSlices, PathPart1, PathPart2, Path3, FrameInfo, pixelSize, ...
      SpotFilter, numParticles, smart_add, xSize, ySize, NDigits, intScale);

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

  %Figure out channel-specific information
  coatChannels = []
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
    coatChannel = coatChannels(1);
  else
    error('Experiment type not recognized')
  end
  
  %Update the name suffix
  if strcmpi(ExperimentType, '2spot2color')
    nameSuffix = ['_ch', iIndex(coatChannel, 2)];
  end

  cc = 1;
  
  if ForCompileAll
    %we just want to save the data, we set 'x' as user command
    cc = 'x';
  end
  
  while (cc ~= 'x')
    % !!!!! JP: 2019-03-22
    % re-check that non-returning variables are NOT being modified inside
    % mainLoop function
    [cc, Particles, CurrentChannel, Spots, CurrentFrame, CurrentParticle, ManualZFlag, projectionMode, DisplayRangeSpot,...
      ShowThreshold2, schnitzcells, DisplayRange, SpotFilter, ZoomMode, GlobalZoomMode, xForZoom, yForZoom,...
      coatChannel, PreviousChannel, PreviousParticle, lastParticle, HideApprovedFlag, lineFitted, CurrentZ, no_clicking] =...
      checkParticleTracking_mainLoop(cc, Particles, CurrentChannel, Spots, CurrentFrame, CurrentParticle, ManualZFlag,...
  NChannels, projectionMode, PreProcPath, FilePrefix, NDigits, nameSuffix, Overlay, DisplayRangeSpot, overlayAxes,...
  UseHistoneOverlay, SpeedMode, FrameInfo, ShowThreshold2, ZSlices, numFrames, schnitzcells, UseSchnitz,...
  DisplayRange, Ellipses, SpotFilter, ZoomMode, GlobalZoomMode, ZoomRange, xForZoom, yForZoom, Prefix, coatChannel,...
  snippetFigAxes, rawDataAxes, gaussianAxes, ExperimentType, intScale, snippet_size, xSize, ySize,...
  SnippetEdge, traceFigAxes, PreviousChannel, PreviousParticle, lastParticle, HideApprovedFlag, lineFitted,...
  anaphaseInMins, ElapsedTime, plot3DGauss, anaphase, prophase, metaphase, prophaseInMins, metaphaseInMins, DefaultDropboxFolder,...
  correspondingNCInfo, Coefficients, zProfileFigAxes, zTraceAxes, frame_num, z_num, particle_num,...
  CurrentZ, no_clicking, DataFolder, DropboxFolder,HisOverlayFig,HisOverlayFigAxes, coatChannels);
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
