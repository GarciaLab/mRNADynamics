function [RawParticles,SpotFilter,ParticleStitchInfo] = track01ParticleProximity(...
    FrameInfo, Spots, schnitzcells, liveExperiment, PixelSize, MaxSearchRadiusMicrons, ...
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
  ExperimentType = liveExperiment.experimentType;
  % load particles and link info if we're retracking
  RawParticles = cell(1,NCh);
  if retrack
    [Particles, SpotFilter] = getParticles(liveExperiment);
    if ~iscell(SpotFilter)
      SpotFilter = {SpotFilter};
    end
    if ~iscell(Particles)
      Particles = {Particles};
    end
    ParticleStitchInfo = getParticleStitchInfo(liveExperiment);
  else
    %Initialize Particles for the number of spot channels we have
    SpotFilter = createSpotFilter(Spots);         
    ParticleStitchInfo = cell(1,NCh);
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
    
    % If we're retracking, we need to (a) break up un-approved particles
    % and (b) make a temporary spotFilter that indicates which Spots are a
    % part of approved particles
    SpotFilterTemp = ones(size(SpotFilter{Channel}));
    if retrack
      % keep only full particles that were approved. Individual frames that
      % were linked are taken care of later on
      Particles{Channel} = Particles{Channel}([Particles{Channel}.Approved]==1);              
      % Create Temp Spot Filter to indicate which spots are fair game and
      % which are bound up in approved particles      
      for p = 1:length(Particles{Channel})
        appFrames = Particles{Channel}(p).Frame;
        appIndices = Particles{Channel}(p).Index;
        appLinIndices = sub2ind(size(SpotFilterTemp),appFrames,appIndices);
        SpotFilterTemp(appLinIndices) = 0;
      end      
    else % initialize stitch info structure
      ParticleStitchInfo{Channel}.persistentLinkIndexCell = {};
      ParticleStitchInfo{Channel}.persistentLinkFrameCell= {};
      ParticleStitchInfo{Channel}.forbiddenLinkIndexCell = {};
      ParticleStitchInfo{Channel}.forbiddenLinkFrameCell = {};
    end
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
        particlesFlag = ~isempty(RawParticles{Channel});
        % if we have nucleus info, assign each spot to a nucleus
        [NewSpotNuclei, NewSpotDistances] = getNuclearAssigments(NewSpotsX,NewSpotsY,...
              schnitzcells,CurrentFrame,UseHistone);
        %Get a list of the particles that were present in
        %the previous frame and of their positions.
        ExtantParticles = [];
        PrevSpotsX = [];
        PrevSpotsY = [];   
        PrevNuclei = [];
        if particlesFlag
          for j = 1:length(RawParticles{Channel})
            if RawParticles{Channel}(j).Frame(end) == (CurrentFrame - 1)
              ExtantParticles = [ExtantParticles, j];            
              PrevSpotsX = [PrevSpotsX, RawParticles{Channel}(j).xPos(end)];
              PrevSpotsY = [PrevSpotsY, RawParticles{Channel}(j).yPos(end)];
              PrevNuclei = [PrevNuclei, RawParticles{Channel}(j).Nucleus];
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
          % frame-over-frame shift;
          maxFrame = find(ncVec==ncVec(CurrentFrame),1,'last');
          [SpotBulkDxVec,SpotBulkDyVec] = getNuclearShifts(schnitzcells,...
                    CurrentFrame,maxFrame,NewSpotsX,NewSpotsY,SlidingWindowSize,UseHistone);
          
          % Calculate the distances between the spots in this frame
          % and those within the particles in the previous
          % frame. Adjust for nuclear movements                       
          DistanceMat = sqrt((NewSpotsX'-SpotBulkDxVec'- PrevSpotsX).^2 + (NewSpotsY'-SpotBulkDyVec' - PrevSpotsY).^2)*PixelSize;
          % enforce consistent nucleus IDs
          mismatchMat = NewSpotNuclei'~=PrevNuclei;
          DistanceMat(mismatchMat) = Inf;
          
          % apply user-assigned links that are within a single frame
          [PrevJoin, NewJoin] = findPersistentLinks(ParticleStitchInfo{Channel},CurrentFrame);
          persistentIndices = sub2ind(size(DistanceMat),NewJoin,PrevJoin);
          DistanceMat(persistentIndices) = 0;
          
          % apply user-assigned splits that are within a single frame
          [PrevSplit, NewSplit] = findForbiddenLinks(ParticleStitchInfo{Channel},CurrentFrame);
          forbiddenIndices = sub2ind(size(DistanceMat),NewSplit,PrevSplit);
          DistanceMat(forbiddenIndices) = Inf;
          
          % remove existing spots that are in approved particles 
          prevAppIndices = SpotFilterTemp(CurrentFrame-1,:)==0;
          nextAppIndices = SpotFilterTemp(CurrentFrame,:)==0;
          DistanceMat(:,prevAppIndices) = Inf;
          DistanceMat(nextAppIndices,:) = Inf;
          
          % Find existing particles and new spots are close enough to be 
          % linked. In cases of degenerate assignemnt, take pairt that
          % minimizes jump distance          
          [MatchIndices,~,~] = matchpairs(DistanceMat,0.5*SearchRadius);
          NewParticleFlag(MatchIndices(:,1)) = false;         

          % Assign matched spots to existing particles
          for j = 1:size(MatchIndices)
            ParticleIndex = ExtantParticles(MatchIndices(j,2));
            NewSpotIndex = MatchIndices(j,1);
            RawParticles{Channel}(ParticleIndex).Frame(end + 1) = CurrentFrame;
            RawParticles{Channel}(ParticleIndex).Index(end + 1) = NewSpotIndex;            
            RawParticles{Channel}(ParticleIndex).LastFrame = CurrentFrame; % not sure I'll use these
            RawParticles{Channel}(ParticleIndex).xPos(end + 1) = NewSpotsX(NewSpotIndex);
            RawParticles{Channel}(ParticleIndex).yPos(end + 1) = NewSpotsY(NewSpotIndex);
            RawParticles{Channel}(ParticleIndex).zPos(end + 1) = NewSpotsZ(NewSpotIndex);      
            RawParticles{Channel}(ParticleIndex).NucleusDist(end + 1) = NewSpotDistances(NewSpotIndex);
          end                        

        end
   
        % See which spots weren't assigned and add them to the structure as new particles
        NewParticles = find(NewParticleFlag);
        TotalParticles = 0;
        if particlesFlag
          TotalParticles = length(RawParticles{Channel});
        end
        for j = NewParticles    
          TotalParticles = TotalParticles + 1;
          RawParticles{Channel}(TotalParticles).Frame = CurrentFrame;
          RawParticles{Channel}(TotalParticles).Index = j;
          RawParticles{Channel}(TotalParticles).Approved = 0;
          RawParticles{Channel}(TotalParticles).FirstFrame = CurrentFrame; % not sure I'll use these
          RawParticles{Channel}(TotalParticles).LastFrame = CurrentFrame; % not sure I'll use these
          RawParticles{Channel}(TotalParticles).xPos = NewSpotsX(j);
          RawParticles{Channel}(TotalParticles).yPos = NewSpotsY(j);
          RawParticles{Channel}(TotalParticles).zPos = NewSpotsZ(j);          
          RawParticles{Channel}(TotalParticles).Nucleus = NewSpotNuclei(j); 
          RawParticles{Channel}(TotalParticles).NucleusDist = NewSpotDistances(j); 
        end
      end

    end
    
    % adjust for z stack shifts
    RawParticles{Channel} = detrendZ(RawParticles{Channel},FrameInfo);
      
    close(f);
  end
end