function [cc, Particles, CurrentChannel, Spots, CurrentFrame, CurrentParticle, ManualZFlag, projectionMode, DisplayRangeSpot,...
      ShowThreshold2, schnitzcells, DisplayRange, SpotFilter, ZoomMode, GlobalZoomMode, xForZoom, yForZoom,...
      coatChannel, PreviousChannel, PreviousParticle, lastParticle, HideApprovedFlag, lineFitted, CurrentZ, no_clicking] =...
      checkParticleTracking_mainLoop(cc, Particles, CurrentChannel, Spots, CurrentFrame, CurrentParticle, ManualZFlag,...
  NChannels, projectionMode, PreProcPath, FilePrefix, NDigits, nameSuffix, Overlay, DisplayRangeSpot, overlayAxes,...
  UseHistoneOverlay, SpeedMode, FrameInfo, ShowThreshold2, ZSlices, numFrames, schnitzcells, UseSchnitz,...
  DisplayRange, Ellipses, SpotFilter, ZoomMode, GlobalZoomMode, ZoomRange, xForZoom, yForZoom, Prefix, coatChannel,...
  snippetFigAxes, rawDataAxes, gaussianAxes, ExperimentType, intScale, snippet_size, xSize, ySize, SnippetEdge,...
  traceFigAxes, PreviousChannel, PreviousParticle, lastParticle, HideApprovedFlag, lineFitted, anaphaseInMins,...
  ElapsedTime, plot3DGauss, anaphase, prophase, metaphase, prophaseInMins, metaphaseInMins, DefaultDropboxFolder,...
  correspondingNCInfo, Coefficients, zProfileFigAxes, zTraceAxes, frame_num, z_num, particle_num, CurrentZ, no_clicking, DataFolder, DropboxFolder,HisOverlayFig,HisOverlayFigAxes)

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

  persistent CurrentZIndex
  
  if (~isempty(xTrace)) && (~ManualZFlag)
    CurrentZ = z(CurrentParticleIndex);
    CurrentZIndex = find(...
      Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).z == ...
      CurrentZ);
    ManualZFlag = 0;
  end

  if NChannels == 1% inputoutput mode can also be in this case, changed CurrentChannel to the coatChannel (YJK : 1/15/2018)

    if strcmpi(projectionMode, 'None (Default)')
      Image = imread([PreProcPath, filesep, FilePrefix(1:end - 1), filesep, ...
                    FilePrefix, iIndex(CurrentFrame, NDigits), '_z', iIndex(CurrentZ, 2), nameSuffix, '.tif']);
    elseif strcmpi(projectionMode, 'Max Z')
      [Image, ~] = zProjections(Prefix, coatChannel, CurrentFrame, ZSlices, NDigits, DropboxFolder, PreProcPath);
    elseif strcmpi(projectionMode, 'Median Z')
      [~, Image] = zProjections(Prefix, coatChannel, CurrentFrame, ZSlices, NDigits, DropboxFolder, PreProcPath);
    elseif strcmpi(projectionMode, 'Max Z and Time')

      if isempty(storedTimeProjection)

        if ncRange
          Image = timeProjection(Prefix, coatChannel, 'nc', NC);
          storedTimeProjection = Image;
        else
          Image = timeProjection(Prefix, CurrentChannel);
          storedTimeProjection = Image;
        end

      else
        Image = storedTimeProjection;
      end

    end

  elseif NChannels > 1
    Image = imread([PreProcPath, filesep, FilePrefix(1:end - 1), filesep, ...
                  FilePrefix, iIndex(CurrentFrame, NDigits), '_z', iIndex(CurrentZ, 2), ...
                  nameSuffix, '.tif']);
  else
    error('ExperimentType and/or channel not supported.')
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
    [ImageHis, xForZoom, yForZoom] = plotFrame(overlayAxes, Image, SpeedMode, FrameInfo, Particles, ...
      Spots, CurrentFrame, ShowThreshold2, ...
      Overlay, CurrentChannel, CurrentParticle, ZSlices, CurrentZ, numFrames, ...
      schnitzcells, UseSchnitz, DisplayRange, Ellipses, SpotFilter, ZoomMode, GlobalZoomMode, ...
      ZoomRange, xForZoom, yForZoom, UseHistoneOverlay, HisOverlayFigAxes, HisPath1, HisPath2);

  else
    plotFrame(overlayAxes, Image, SpeedMode, ...
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

  persistent CurrentSnippet
  if ~isempty(CurrentSnippet)
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
    [Frames, AmpIntegral, GaussIntegral, AmpIntegral3, AmpIntegral5, ...
      ErrorIntegral, ErrorIntegral3, ErrorIntegral5, backGround3, ...
      AmpIntegralGauss3D, ErrorIntegralGauss3D, PreviousParticle] = plotTrace(traceFigAxes, ...
      FrameInfo, CurrentChannel, PreviousChannel, ...
      CurrentParticle, PreviousParticle, lastParticle, HideApprovedFlag, lineFitted, anaphaseInMins, ...
      ElapsedTime, schnitzcells, Particles, plot3DGauss, anaphase, prophase, metaphase, prophaseInMins, metaphaseInMins, Prefix, ...
      DefaultDropboxFolder, numFrames, CurrentFrame, ZSlices, CurrentZ, Spots, ...
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
      DefaultDropboxFolder, numFrames, CurrentFrame, ZSlices, CurrentZ, Spots, ...
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

    if isempty(DisplayRange)
      DisplayRange = [min(min(ImageHis)), max(max(ImageHis)) / 1.5];
    else
      DisplayRange = [DisplayRange(1), DisplayRange(2) / 1.5];
    end

  elseif cc == 'b' & UseHistoneOverlay%Decrease histone channel contrast
    DisplayRange = [min(min(ImageHis)), max(max(ImageHis)) * 1.5];

  elseif cc == '#'%remove a spot from Spots and erase its frame in Particles
    [Spots, SpotFilter, ZoomMode, GlobalZoomMode, CurrentFrame, ...
      CurrentParticle, Particles, ManualZFlag, DisplayRange, lastParticle, PreviousParticle] = ...
      removeSpot(ZoomMode, GlobalZoomMode, Frames, CurrentFrame, ...
      CurrentChannel, CurrentParticle, CurrentParticleIndex, Particles, Spots, SpotFilter, ...
      numParticles, ManualZFlag, DisplayRange, lastParticle, PreviousParticle);
    elseif cc == '[' | cc == '{'%#ok<*OR2> %Add particle and all of its shadows to Spots.
              PathPart1 = [PreProcPath, filesep, FilePrefix(1:end - 1), filesep, FilePrefix];
    PathPart2 = [nameSuffix, '.tif'];
    Path3 = [PreProcPath, filesep, Prefix, filesep, Prefix];
    [numParticles, SpotFilter, Particles, Spots, PreviousParticle] = ...
      addSpot(ZoomMode, GlobalZoomMode, Particles, CurrentChannel, ...
      CurrentParticle, CurrentFrame, CurrentZ, Overlay, snippet_size, PixelsPerLine, ...
      LinesPerFrame, Spots, ZSlices, PathPart1, PathPart2, Path3, FrameInfo, pixelSize, ...
      SpotFilter, numParticles, cc, xSize, ySize, NDigits, intScale);
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
  elseif cc == 't'
    ShowThreshold2 = ~ShowThreshold2;

    %     elseif (cc=='y')&(~UseHistoneOverlay)
    %             FrameInfo=DetermineNC(fad,Particles{CurrentChannel},FrameInfo);
    %
    %AR 9/5/18- this button is deprecated. leaving this comment in
    %case we want to replace the functionality.
    %
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
    else
      disp('Try again after exiting global zoom mode by hitting ''+ ''')
    end

  elseif cc == '+'

    if ~ZoomMode

      if ~GlobalZoomMode
        %AR 12/14/2018 ginputc doesn't work at the moment. not sure
        %why. maybe related to the fact that this is an image
        %within a figure that has multiple axes.
        %                 [ConnectPositionx,ConnectPositiony]=ginputc(1,'color', 'r', 'linewidth',1);
        [ConnectPositionx, ConnectPositiony] = ginput(1);
        xForZoom = round(ConnectPositionx);
        yForZoom = round(ConnectPositiony);
      else
      end

      GlobalZoomMode = ~GlobalZoomMode;
    else
      disp('Try again after exiting zoom mode by hitting ''o''')
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
    DisplayRangeSpot = [min(min(Image)), max(max(Image)) * 1.5];
  elseif cc == '0'%Debugging mode
    keyboard;

  end

end