function [Particles] = trackParticlesBasedOnProximity(...
    Particles, Spots, xPos, SpotFilter, Channel, CurrentFrame,...
    PixelSize, SearchRadius, retrack, displayFigures)
  
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