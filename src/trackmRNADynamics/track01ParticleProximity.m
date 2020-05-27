function [Particles] = track01ParticleProximity(...
    FrameInfo, Spots, schnitzcells, NCh, PixelSize, SearchRadiusMicrons, retrack, displayFigures)
  
  %This function is the first stage of tracking performed by the 
  %performTracking subfunction. It initializes temporary data structures
  %used throughout tracking process and takes a first pass at linking spots
  %into particles based on proximity alone
  
  %NOTE: currently not supporting a retracking option or displayFigures option.
  % Need to think about this
  
  % Extract Time Vector
  TimeVec = [FrameInfo.Time];
  
  % Check to see if we have nucleus tracking info
  UseNuclei = ~isempty(schnitzcells);
  SlidingWindowSize = 2; % size of window used for time averaging
  
  for Channel = 1:NCh
    
    for CurrentFrame = 1:length(Spots{Channel})
      
      % Get the positions of ALL spots (approved and disapproved)
      [NewSpotsX, NewSpotsY, NewSpotsZ] = SpotsXYZ(Spots{Channel}(CurrentFrame));            
      
      if CurrentFrame==1 && ~isempty(NewSpotsX)
        %Initialize the Particles structure 
        for j = 1:length(NewSpotsX)
          Particles{Channel}(j).Frame = CurrentFrame;
          Particles{Channel}(j).Index = j;
          Particles{Channel}(j).Approved = 0;
          Particles{Channel}(j).FirstFrame = CurrentFrame; % not sure I'll use these
          Particles{Channel}(j).LastFrame = CurrentFrame; % not sure I'll use these
          Particles{Channel}(j).xPos = NewSpotsX(j);
          Particles{Channel}(j).yPos = NewSpotsY(j);
          Particles{Channel}(j).zPos = NewSpotsZ(j);
        end

      elseif ~isempty(NewSpotsX)
        %If particle structure is not empty, then we must link current particles to
        %prev particles wherever possible
        dT = TimeVec(CurrentFrame)-TimeVec(CurrentFrame-1);
        SearchRadius = SearchRadiusMicrons * sqrt(dT);
        
        % if we have nulceus tracking info, calculate average
        % frame-over-frame shift
        if UseNuclei
          % use sliding window to estimate average nucleus movement
          NucleiDxVec = [];
          NucleiDyVec = [];
          NewNucleiX = [];
          NewNucleiY = [];
          for i = 1:length(schnitzcells)
            StopFrame = min([length(Spots{Channel}),CurrentFrame+SlidingWindowSize]);
            StartFrame = max([1,CurrentFrame-SlidingWindowSize]);
            FrameIndices = find(ismember(schnitzcells(i).frames,[StartFrame,StopFrame]));
            if length(FrameIndices)==2
              NucleiDxVec = [NucleiDxVec diff(schnitzcells(i).cenx(FrameIndices))];
              NucleiDyVec = [NucleiDyVec diff(schnitzcells(i).ceny(FrameIndices))];
              NewNucleiX = [NewNucleiX schnitzcells(i).cenx(FrameIndices(1))];
              NewNucleiY = [NewNucleiY schnitzcells(i).ceny(FrameIndices(1))];
            end
          end

          % calculate distance to each nucleus 
          NucleusWeightMat = NaN(length(NewSpotsX),length(NewNucleiX));
          for i = 1:length(NewSpotsX)
            NucleusWeightMat(i,:) = (rand(numel(NewNucleiX),1)*0.05+vecnorm([NewSpotsX(i) NewSpotsY(i)]  - [NewNucleiX' NewNucleiY'], 2, 2)).^-2;
          end
          % assign weighted mean bulk displacement to particles
          SpotBulkDxVec = (sum(repmat(NucleiDxVec,length(NewSpotsX),1).*NucleusWeightMat,2) ./ sum(NucleusWeightMat,2) / (StopFrame-StartFrame+1))';
          SpotBulkDyVec = (sum(repmat(NucleiDyVec,length(NewSpotsX),1).*NucleusWeightMat,2) ./ sum(NucleusWeightMat,2) / (StopFrame-StartFrame+1))';        
        else
          SpotBulkDxVec = zeros(size(NewSpotsX));
          SpotBulkDyVec = zeros(size(NewSpotsY));
        end
        
        %Get a list of the particles that were present in
        %the previous frame and of their positions.
        ExtantParticles = [];
        PrevSpotsX = [];
        PrevSpotsY = [];        

        for j = 1:length(Particles{Channel})

          if Particles{Channel}(j).Frame(end) == (CurrentFrame - 1)
            ExtantParticles = [ExtantParticles, j];            
            PrevSpotsX = [PrevSpotsX, Particles{Channel}(j).xPos(end)];
            PrevSpotsY = [PrevSpotsY, Particles{Channel}(j).yPos(end)];
          end

        end
 
        % If there were particles present in the previous frame,
        % then we find their distances to the spots present in
        % the current frame. Otherwise, we create new particles
        % for each new spot.
        
        % We keep track of which spots go to a new or old
        % particle using the NewParticle array.
        NewParticleFlag = true(size(NewSpotsX));

        if ~isempty(ExtantParticles)
          % Calculate the distances between the spots in this frame
          % and those within the particles in the previous
          % frame          
          
          % Adjust for nuclear movements
          DistanceMat = sqrt((NewSpotsX'-SpotBulkDxVec'- PrevSpotsX).^2 + (NewSpotsY'-SpotBulkDyVec' - PrevSpotsY).^2)*PixelSize;

          % Find existing particles and new spots are close enough to be 
          % linked. In cases of degenerate assignemnt, take pairt that
          % minimizes jump distance
          
          [MatchIndices,~,~] = matchpairs(DistanceMat,0.5*SearchRadius);
          NewParticleFlag(MatchIndices(:,1)) = false;         
          
          % Assign matchesd spots to existing particles
          for j = 1:size(MatchIndices)
            ParticleIndex = ExtantParticles(MatchIndices(j,2));
            NewSpotIndex = MatchIndices(j,1);
            Particles{Channel}(ParticleIndex).Frame(end + 1) = CurrentFrame;
            Particles{Channel}(ParticleIndex).Index(end + 1) = NewSpotIndex;            
            Particles{Channel}(ParticleIndex).LastFrame = CurrentFrame; % not sure I'll use these
            Particles{Channel}(ParticleIndex).xPos(end + 1) = NewSpotsX(NewSpotIndex);
            Particles{Channel}(ParticleIndex).yPos(end + 1) = NewSpotsY(NewSpotIndex);
            Particles{Channel}(ParticleIndex).zPos(end + 1) = NewSpotsZ(NewSpotIndex);                           
          end                        

        end

        % See which spots weren't assigned and add them to the structure as new particles
        NewParticles = find(NewParticleFlag);

        for j = NewParticles
          TotalParticles = length(Particles{Channel});
          Particles{Channel}(TotalParticles + 1).Frame = CurrentFrame;
          Particles{Channel}(TotalParticles + 1).Index = j;
          Particles{Channel}(TotalParticles + 1).Approved = 0;
          Particles{Channel}(TotalParticles + 1).FirstFrame = CurrentFrame; % not sure I'll use these
          Particles{Channel}(TotalParticles + 1).LastFrame = CurrentFrame; % not sure I'll use these
          Particles{Channel}(TotalParticles + 1).xPos = NewSpotsX(j);
          Particles{Channel}(TotalParticles + 1).yPos = NewSpotsY(j);
          Particles{Channel}(TotalParticles + 1).zPos = NewSpotsZ(j);          
        end
      end

    end
    
    % adjust for z stack shifts
    if isfield(FrameInfo,'zPosition')      
      % need to make sure this field is now a permanent feature
      zPosVec = [FrameInfo.zPosition]*1e6 / FrameInfo(1).ZStep;
      zPosVec = zPosVec - zPosVec(1);
      frameIndex = 1:length(zPosVec);
      % generate new det-trended z variable
      for p = 1:length(Particles{Channel})
        fVec = Particles{Channel}(p).Frame;
        Particles{Channel}(p).zPosDetrended = Particles{Channel}(p).zPos - zPosVec(ismember(frameIndex,fVec));
      end
    else
      % get list of all frames and corresponding z positions
      frameVec = [Particles{Channel}.Frame];
      zPosVec = [Particles{Channel}.zPos];
      
      % get iteratable frame  list
      frameIndex = unique(frameVec);
      avgZProfile = NaN(size(frameIndex));
      
      for f = 1:length(frameIndex)
        avgZProfile(f) = nanmean(zPosVec(frameVec==frameIndex(f)));
      end
      
      % generate new det-trended z variable
      for p = 1:length(Particles{Channel})
        fVec = Particles{Channel}(p).Frame;
        Particles{Channel}(p).zPosDetrended = Particles{Channel}(p).zPos - avgZProfile(ismember(frameIndex,fVec));
      end
      
    end
  end
end