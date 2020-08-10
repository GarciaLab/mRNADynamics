function [Particles] = track01ParticleProximity(...
    FrameInfo, Spots, schnitzcells, Prefix, PixelSize, MaxSearchRadiusMicrons, ...
    UseHistone, retrack, displayFigures)
  
  %This function is the first stage of tracking performed by the 
  %performTracking subfunction. It initializes temporary data structures
  %used throughout tracking process and takes a first pass at linking spots
  %into particles based on proximity alone
  
  %NOTE: currently not supporting a retracking option or displayFigures option.
  % Need to think about this
  
  % get number of channels
  NCh = length(Spots);
  
  % get experiment type
  liveExperiment = LiveExperiment(Prefix);
  ExperimentType = liveExperiment.experimentType;
  % load particles and link info if we're retracking
  if retrack
    [Particles, SpotFilter] = getParticles(liveExperiment);
    ParticleStitchInfo = getParticleStitchInfo(liveExperiment);
  else
    %Initialize Particles for the number of spot channels we have
    Particles = cell(1,NCh);
    ParticleStitchInfo = [];
  end
  
  % Extract vector indicating nuclear cleavage cycle for each frame
  ncVec = [FrameInfo.nc]; 
  SlidingWindowSize = 2; % size of window used for time averaging of nuclear movements
    
  if ismember(ExperimentType,{'inputoutput'}) && NCh == 2    
    % Figure out which channel is the cluster channel
    spotsChannels = liveExperiment.spotChannels;
    inputChannels = liveExperiment.inputChannels;
    matchPercent(spotsChannels == intersect(inputChannels, spotsChannels)) = 50;
    matchPercent(spotsChannels ~= intersect(inputChannels, spotsChannels)) = 75;    
  else  % No clusters:
    matchPercent(1:NCh) = 75;
  end
  
  for Channel = 1:NCh    
    f = waitbar(0,['Stitching particle tracks (channel ' num2str(Channel) ')']);
    % calculate jump dist threshold
    SearchRadiusMicrons = estimateSearchRadius(Spots,ncVec,MaxSearchRadiusMicrons,PixelSize,matchPercent,Channel);
    SearchRadiusMicrons = max([0.5,SearchRadiusMicrons]); % should be at least 0.5um
    
    for CurrentFrame = 1:length(Spots{Channel})
      waitbar(CurrentFrame/length(Spots{Channel}),f);
      % Get the positions of ALL spots (approved and disapproved)
      [NewSpotsX, NewSpotsY, NewSpotsZ] = SpotsXYZ(Spots{Channel}(CurrentFrame));
      
      % Check if we're at the start of a new nuclear cycle
      NewNCFlag = false;
      if CurrentFrame > 1
        NewNCFlag = ncVec(CurrentFrame-1)~=ncVec(CurrentFrame);
      end
      
      if ~isempty(NewSpotsX)
        %If spots were detected for this frame we need to add them to the
        %Particles structure. Link them to existing particles whenever
        %possible                        
        particlesFlag = ~isempty(Particles{Channel});
        % if we have nucleus info, assign each spot to a nucleus
        [NewSpotNuclei, NewSpotDistances] = getNuclearAssigments(NewSpotsX,NewSpotsY,...
              schnitzcells,UseHistone);
        %Get a list of the particles that were present in
        %the previous frame and of their positions.
        ExtantParticles = [];
        PrevSpotsX = [];
        PrevSpotsY = [];   
        PrevNuclei = [];
        if particlesFlag
          for j = 1:length(Particles{Channel})
            if Particles{Channel}(j).Frame(end) == (CurrentFrame - 1)
              ExtantParticles = [ExtantParticles, j];            
              PrevSpotsX = [PrevSpotsX, Particles{Channel}(j).xPos(end)];
              PrevSpotsY = [PrevSpotsY, Particles{Channel}(j).yPos(end)];
              PrevNuclei = [PrevNuclei, Particles{Channel}(j).Nucleus];
            end
          end
        end
        % If there were particles present in the previous frame,
        % then we find their distances to the spots present in
        % the current frame. Otherwise, we create new particles
        % for each new spot.
        
        % We keep track of which spots go to a new or old
        % particle using the NewParticle array.
        NewParticleFlag = true(size(NewSpotsX));
        
        if ~isempty(ExtantParticles) && ~NewNCFlag
          % Generate maximum allowed jump radius between spots          
          SearchRadius = SearchRadiusMicrons;
          
          % if we have nulceus tracking info, calculate average
          % frame-over-frame shift
          [SpotBulkDxVec,SpotBulkDyVec] = getNuclearShifts(schnitzcells,...
                    CurrentFrame,NewSpotsX,NewSpotsY,UseHistone);
          
          % Calculate the distances between the spots in this frame
          % and those within the particles in the previous
          % frame          
          
          % Adjust for nuclear movements
          DistanceMat = sqrt((NewSpotsX'-SpotBulkDxVec'- PrevSpotsX).^2 + (NewSpotsY'-SpotBulkDyVec' - PrevSpotsY).^2)*PixelSize;
          % enforce consistent nucleus IDs
          mismatchMat = NewSpotNuclei'~=PrevNuclei;
          DistanceMat(mismatchMat) = Inf;
          % check for user-assigned links
          [PrevIndices, NewIndices] = findPersistentLinks(ParticleStitchInfo,CurrentFrame);
          persistentIndices = sub2ind(size(DistanceMat),NewIndices,PrevIndices);
          DistanceMat(persistentIndices) = 0;
          
          % Find existing particles and new spots are close enough to be 
          % linked. In cases of degenerate assignemnt, take pairt that
          % minimizes jump distance
          
          [MatchIndices,~,~] = matchpairs(DistanceMat,0.5*SearchRadius);
          NewParticleFlag(MatchIndices(:,1)) = false;         

          % Assign matched spots to existing particles
          for j = 1:size(MatchIndices)
            ParticleIndex = ExtantParticles(MatchIndices(j,2));
            NewSpotIndex = MatchIndices(j,1);
            Particles{Channel}(ParticleIndex).Frame(end + 1) = CurrentFrame;
            Particles{Channel}(ParticleIndex).Index(end + 1) = NewSpotIndex;            
            Particles{Channel}(ParticleIndex).LastFrame = CurrentFrame; % not sure I'll use these
            Particles{Channel}(ParticleIndex).xPos(end + 1) = NewSpotsX(NewSpotIndex);
            Particles{Channel}(ParticleIndex).yPos(end + 1) = NewSpotsY(NewSpotIndex);
            Particles{Channel}(ParticleIndex).zPos(end + 1) = NewSpotsZ(NewSpotIndex);      
            Particles{Channel}(ParticleIndex).NucleusDist(end + 1) = NewSpotDistances(NewSpotIndex);
          end                        

        end
   
        % See which spots weren't assigned and add them to the structure as new particles
        NewParticles = find(NewParticleFlag);
        TotalParticles = 0;
        if particlesFlag
          TotalParticles = length(Particles{Channel});
        end
        for j = NewParticles    
          TotalParticles = TotalParticles + 1;
          Particles{Channel}(TotalParticles).Frame = CurrentFrame;
          Particles{Channel}(TotalParticles).Index = j;
          Particles{Channel}(TotalParticles).Approved = 0;
          Particles{Channel}(TotalParticles).FirstFrame = CurrentFrame; % not sure I'll use these
          Particles{Channel}(TotalParticles).LastFrame = CurrentFrame; % not sure I'll use these
          Particles{Channel}(TotalParticles).xPos = NewSpotsX(j);
          Particles{Channel}(TotalParticles).yPos = NewSpotsY(j);
          Particles{Channel}(TotalParticles).zPos = NewSpotsZ(j);          
          Particles{Channel}(TotalParticles).Nucleus = NewSpotNuclei(j); 
          Particles{Channel}(TotalParticles).NucleusDist = NewSpotDistances(j); 
        end
      end

    end
    
    % adjust for z stack shifts
    Particles{Channel} = detrendZ(Particles{Channel},FrameInfo);
      
    close(f);
  end
end