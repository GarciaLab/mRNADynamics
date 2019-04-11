function [Particles, SpotFilter] = performTracking(Particles, schnitzcells, NCh, Spots, app, SpotFilter, PreProcPath, Prefix, UseHistone, ParticlesFig, SpotsChannel, NDigits, NucleiFig, particlesAxes, nucAxes, Ellipses, PixelSize, SearchRadius, ExperimentType, FrameInfo, retrack, displayFigures)
  % Iterate over all channels
  for Channel = 1:NCh

    % Iterate over all frames
    for CurrentFrame = 1:length(Spots{Channel})

      if isempty(app) && displayFigures
        figure(ParticlesFig)
        set(ParticlesFig, 'units', 'normalized', 'position', [0.01, .55, .33, .33]);          
      end

      % Get the filter for this frame
      CurrentFrameFilter = logical(SpotFilter{Channel}(CurrentFrame, ~isnan(SpotFilter{Channel}(CurrentFrame, :))));

      xPos = displayParticlesFigure(app, particlesAxes, ParticlesFig, Spots, Channel, CurrentFrame, ...
        CurrentFrameFilter, PreProcPath, Prefix, SpotsChannel, FrameInfo, displayFigures);

      if UseHistone
        [Particles, SpotFilter] = trackParticlesBasedOnNuclei(PreProcPath, Prefix, CurrentFrame, NDigits, app, nucAxes, Ellipses, ...
          ExperimentType, Channel, schnitzcells, Particles, Spots, SpotFilter, PixelSize, SearchRadius, retrack, displayFigures);
      else
        [Particles] = trackParticlesBasedOnProximity(Particles, Spots, xPos, SpotFilter, Channel, CurrentFrame, PixelSize, SearchRadius, retrack);
      end

    end

  end

  if isempty(app) && displayFigures
    close(ParticlesFig)
    if UseHistone
      close(NucleiFig)
    end
  end

  
  for currentChannel = 1:NCh
    
    if ~isfield(Particles{currentChannel}, 'FrameApproved')
      
      for i = 1:length(Particles{currentChannel})
        Particles{currentChannel}(i).FrameApproved = true(size(Particles{currentChannel}(i).Frame));
      end
      
    else
      
      for i = 1:length(Particles{currentChannel})
        
        if isempty(Particles{currentChannel}(i).FrameApproved)
          Particles{currentChannel}(i).FrameApproved = true(size(Particles{currentChannel}(i).Frame));
        end
        
      end
      
    end
    
  end
  
  % If we only have one channel, then convert SpotFilter and Particles to a standard structure.
  if NCh == 1
    SpotFilter = SpotFilter{1};
    Particles = Particles{1};
  end

end

function x = displayParticlesFigure(app, particlesAxes, ParticlesFig, Spots, Channel, CurrentFrame, CurrentFrameFilter, PreProcPath, Prefix, SpotsChannel, FrameInfo, displayFigures)
  
  % Get the positions of the spots in this frame
  [x, y, ~] = SpotsXYZ(Spots{Channel}(CurrentFrame));

  if displayFigures
      % Z plane to be displayed. We use the median of all particles found
      % in this frame
      if ~isempty(Spots{Channel}(CurrentFrame).Fits)
        CurrentZ = round(median([Spots{Channel}(CurrentFrame).Fits.brightestZ]));
      else
        CurrentZ = round(FrameInfo(1).NumberSlices / 2);
      end

      % Load the corresponding mRNA image. Check whether we have multiple
      % channels saved or not.
      FileNamePrefix = [PreProcPath, filesep, Prefix, filesep, Prefix, '_', iIndex(CurrentFrame, 3), '_z', iIndex(CurrentZ, 2)];
      particleImage = imread([FileNamePrefix, '_ch', iIndex(SpotsChannel(Channel), 2), '.tif']);

      FigureName = ['Ch', num2str(Channel), '  Frame: ', num2str(CurrentFrame), '/', num2str(length(Spots{Channel}))];

      if ~isempty(app)
        ax1 = app{1};
        title(ax1, FigureName)
      else
        ax1 = particlesAxes;
        set(ParticlesFig, 'Name', FigureName);
        title(ax1, CurrentFrame)
      end

      imshow(particleImage, [], 'Parent', ax1, 'InitialMagnification', 'fit')
      hold(ax1, 'on')
      plot(ax1, x(CurrentFrameFilter), y(CurrentFrameFilter), 'or', 'MarkerSize', 10)
      plot(ax1, x(~CurrentFrameFilter), y(~CurrentFrameFilter), 'ow', 'MarkerSize', 10)
      hold(ax1, 'off')
  end
end

function [Particles, SpotFilter] = trackParticlesBasedOnNuclei(PreProcPath, Prefix, CurrentFrame, NDigits, app, nucAxes, Ellipses, ...
    ExperimentType, Channel, schnitzcells, Particles, Spots, SpotFilter, PixelSize, SearchRadius, retrack, displayFigures)

if displayFigures
  hisImage = openHistoneImage(Prefix, PreProcPath, CurrentFrame, NDigits);

  if ~isempty(app)
    ax2 = app{2};
  else
    ax2 = nucAxes;
  end

  imshow(hisImage, [], 'Border', 'Tight', 'Parent', ax2, 'InitialMagnification', 'fit')
  hold(ax2, 'on')
  PlotHandle = [];
  [NEllipses, ~] = size(Ellipses{CurrentFrame});

  for EllipsesIndex = 1:NEllipses
    PlotHandle = [PlotHandle, ellipse(...
      Ellipses{CurrentFrame}(EllipsesIndex, 3), ...
      Ellipses{CurrentFrame}(EllipsesIndex, 4), ...
      Ellipses{CurrentFrame}(EllipsesIndex, 5), ...
      Ellipses{CurrentFrame}(EllipsesIndex, 1) + 1, ...
      Ellipses{CurrentFrame}(EllipsesIndex, 2) + 1, ...
      [], [], ax2)];

    text(ax2, Ellipses{CurrentFrame}(EllipsesIndex, 1) + 1, Ellipses{CurrentFrame}(EllipsesIndex, 2) + 1, ...
      num2str(EllipsesIndex), 'BackgroundColor', [.7 .9 .7]);
  end

  set(PlotHandle, 'Color', 'r')
  hold(ax2, 'off')
  title(ax2, CurrentFrame)
  drawnow
end

  if strcmp(ExperimentType, '2spot')
    SpotsPerNucleus = 2;
  else
     SpotsPerNucleus = 1;
  end

  [Particles{Channel}, SpotFilter{Channel}] = AssignParticle2Nucleus(schnitzcells, Ellipses, ...
    Particles{Channel}, Spots{Channel}, SpotFilter{Channel}, CurrentFrame, PixelSize, SpotsPerNucleus, retrack);

end

function hisImage = openHistoneImage(Prefix, PreProcPath, CurrentFrame, NDigits)
  
  HistoneImageFileNamePrefix = [PreProcPath, filesep, Prefix, filesep, Prefix];
  HistoneImageFileNameSuffix = [iIndex(CurrentFrame, NDigits), '.tif'];

  try
    hisImage = imread([HistoneImageFileNamePrefix, '-His_', HistoneImageFileNameSuffix]);
  catch

    try
      hisImage = imread([HistoneImageFileNamePrefix, '_His_', HistoneImageFileNameSuffix]);
    catch
      hisImage = 0;
    end

  end

end

function [Particles] = trackParticlesBasedOnProximity(Particles, Spots, xPos, SpotFilter, Channel, CurrentFrame, PixelSize, SearchRadius, retrack, displayFigures)
  
  %This function is used by the performTracking subfunction of
  %trackmRNADynamics to track particles in the event there's no nuclear
  %channel. 
  
  if displayFigures
    drawnow
  end
  % Get the particles detected in the frame

  if ~isempty(xPos)

    % Find the approved spots in this frame
    %AR 12/29/2018: since spotFilter isn't really used in this context, all
    %spots are approved. 
    ApprovedSpots = find(SpotFilter{Channel}(CurrentFrame, ~isnan(SpotFilter{Channel}(CurrentFrame, :))));

    % Get the positions of ALL spots (approved and disapproved)
    [NewSpotsX, NewSpotsY] = SpotsXYZ(Spots{Channel}(CurrentFrame));

    if isempty(Particles{Channel})
      %Initialize the Particles structure if it doesn't exist yet
      for j = 1:length(ApprovedSpots)
        Particles{Channel}(j).Frame = CurrentFrame;
        Particles{Channel}(j).Index = ApprovedSpots(j);
        Particles{Channel}(j).Approved = 0;
      end

    else
      %If we already have recorded particles, we need to compare
      %them to the new ones found and try to assign them.

      %Get a list of the particles that were present in
      %the previous frame and of their positions.
      PreviousFrameParticles = [];
      xPreviousFrameParticles = [];
      yPreviousFrameParticles = [];
      [PreviousSpotsX, PreviousSpotsY] = SpotsXYZ(Spots{Channel}(CurrentFrame - 1));

      for j = 1:length(Particles{Channel})

        if Particles{Channel}(j).Frame(end) == (CurrentFrame - 1)
          PreviousFrameParticles = [PreviousFrameParticles, j];
          PreviousFrameIndex = Particles{Channel}(j).Index(end);
          xPreviousFrameParticles = [xPreviousFrameParticles, PreviousSpotsX(PreviousFrameIndex)];
          yPreviousFrameParticles = [yPreviousFrameParticles, PreviousSpotsY(PreviousFrameIndex)];
        end

      end

      % If there were particles present in the previous frame,
      % then we find their distances to the spots present in
      % the current frame. Otherwise, we create new particles
      % for each new spot.
      % We keep track of which spots goes to a new or old
      % particle using the NewParticle array.
      NewParticleFlag = true(size(ApprovedSpots));

      if ~isempty(PreviousFrameParticles)
        % Get the distances between the spots in this frame
        % and those within the particles in the previous
        % frame
        clear Distance

        % The rows of Distance correspond to the new spots.
        % The columns correspond to the particles present in
        % the previous frame. Each element is the distance.
        for j = 1:length(NewSpotsX)
          Distance(j, :) = sqrt((NewSpotsX(j) * PixelSize - xPreviousFrameParticles * PixelSize).^2 + ...
            (NewSpotsY(j) * PixelSize - yPreviousFrameParticles * PixelSize).^2);
        end

        % We want to make sure there is only one match of
        % new spot to previous particle. To make this
        % possible, we'll go through each column in the
        % matrix Distance and set all pairwise distances
        % that are not the minimum one to infinity.
        for j = 1:length(PreviousFrameParticles)
          DistanceFilter = false(size(Distance(:, j)));
          [DistMinValue, DistMindIndex] = min(Distance(:, j));
          DistanceFilter(DistMindIndex) = true;
          Distance(~DistanceFilter, j) = inf;
        end

        % The rows of distance correspond to the new particles.
        % The columns correspond to their distance to the old particles.

        % The followign tracking works well if we have more than one previous particle. If not, we need to be
        % more careful.
        if (size(Distance, 2) > 1)
          % MinIndex is a row vector. The element position correspond to the new spot and the value within it
          % correspond to the previous particle that it's closest to.
          [MinValues, MinIndex] = min(Distance');
          % Note that inf can be a distance as well. In those cases, turn MinIndex to 0.
          MinIndex(MinValues == inf) = 0;
          % Now, check that the distances are smaller than SearchRadius
          MinIndex(~(MinValues < SearchRadius)) = 0;

          % Assign the new spots to their corresponding particles.
          if sum(MinIndex)

            for j = 1:length(MinIndex)

              if MinIndex(j) > 0
                Particles{Channel}(PreviousFrameParticles(MinIndex(j))).Frame(end + 1) = CurrentFrame;
                Particles{Channel}(PreviousFrameParticles(MinIndex(j))).Index(end + 1) = ApprovedSpots(j);

                %We don't want this new spot to generate a
                %new particle further below
                NewParticleFlag(j) = false;
              end

            end

          end

        else
          % Find the new spot that is closest to the one previous particle
          [MinValues, MinIndex] = min(Distance);
          % Note that inf can be a distance as well. In those cases, turn MinIndex to 0.
          MinIndex(MinValues == inf) = 0;
          % Now, check that the distances are smaller than SearchRadius
          MinIndex(~(MinValues < SearchRadius)) = 0;

          if sum(MinIndex)
            Particles{Channel}(PreviousFrameParticles).Frame(end + 1) = CurrentFrame;
            Particles{Channel}(PreviousFrameParticles).Index(end + 1) = MinIndex;
            % We don't want this new spot to generate a new particle further below
            NewParticleFlag(MinIndex) = false;
          end

        end

      end

      % See which spots weren't assigned and add them to the structure as new particles
      NewParticles = find(NewParticleFlag);

      for j = 1:length(NewParticles)
        TotalParticles = length(Particles{Channel});
        Particles{Channel}(TotalParticles + 1).Frame = CurrentFrame;
        Particles{Channel}(TotalParticles + 1).Index = ApprovedSpots(NewParticles(j));
        Particles{Channel}(TotalParticles + 1).Approved = 0;
      end

    end

  end

end
