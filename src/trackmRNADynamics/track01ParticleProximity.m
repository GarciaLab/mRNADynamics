function [Particles] = track01ParticleProximity(...
    FrameInfo, Spots, schnitzcells, Prefix, PixelSize, MaxSearchRadiusMicrons, UseHistone, retrack, displayFigures)
  
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
  % Extract Time Vector
  TimeVec = [FrameInfo.Time];
  % Extract vector indicating nuclear cleavage cycle for each frame
  ncVec = [FrameInfo.nc]; 
  SlidingWindowSize = 2; % size of window used for time averaging of nuclear movements
  %Initialize Particles for the number of spot channels we have
  Particles = cell(1,NCh);
  
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
        if UseHistone     
          ExtantNucleiX = [];
          ExtantNucleiY = [];
          ncIDVec = [];
          for i = 1:length(schnitzcells)             
            CurrFT = ismember(schnitzcells(i).frames,CurrentFrame);
            if sum(CurrFT) == 1
              ExtantNucleiX = [ExtantNucleiX schnitzcells(i).cenx(CurrFT)];
              ExtantNucleiY = [ExtantNucleiY schnitzcells(i).ceny(CurrFT)];
              ncIDVec = [ncIDVec i];
            end
          end
          NucleusDistMat = NaN(length(NewSpotsX),length(ExtantNucleiX));
          for i = 1:length(NewSpotsX)            
            NucleusDistMat(i,:) = vecnorm([NewSpotsX(i) NewSpotsY(i)]  - [ExtantNucleiX' ExtantNucleiY'],2,2);
          end
          % assign spots to nearest neighbors
          [NewSpotDistances, minIndices] = min(NucleusDistMat,[],2); % note that I'm not enforcing unique assignment at this stage. Will do this later on in the process
          NewSpotNucleuss = ncIDVec(minIndices); 
        else
          NewSpotNucleuss = ones(size(NewSpotsX)); % if no nucleus info, then we set all ID values to dummy val
          NewSpotDistances = zeros(size(NewSpotsX));
        end
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
          if UseHistone            
            % use sliding window to estimate average nucleus movement
            NucleiDxVec = [];
            NucleiDyVec = []; 
            NucleiPxVec = [];
            NucleiPyVec = []; 
            StopFrame = min([length(Spots{Channel}),CurrentFrame+SlidingWindowSize]);
            StartFrame = max([1,CurrentFrame-SlidingWindowSize]);
            for i = 1:length(schnitzcells)             
              FrameFT = ismember(schnitzcells(i).frames,[StartFrame,StopFrame]);
              if sum(FrameFT)==2
                NucleiDxVec = [NucleiDxVec diff(schnitzcells(i).cenx(FrameFT))];
                NucleiDyVec = [NucleiDyVec diff(schnitzcells(i).ceny(FrameFT))];
                NucleiPxVec = [NucleiPxVec mean(schnitzcells(i).cenx(FrameFT))];
                NucleiPyVec = [NucleiPyVec mean(schnitzcells(i).ceny(FrameFT))];
              end
            end

            % calculate distancen to each nucleus 
            NucleusWeightMat = NaN(length(NewSpotsX),length(NucleiPxVec));
            for i = 1:length(NewSpotsX)
              NucleusWeightMat(i,:) = (rand(numel(NucleiPxVec),1)*0.05+vecnorm([NewSpotsX(i) NewSpotsY(i)]  - [NucleiPxVec' NucleiPyVec'], 2, 2)).^-2;              
            end
            % assign weighted mean bulk displacement to particles
            SpotBulkDxVec = (sum(repmat(NucleiDxVec,length(NewSpotsX),1).*NucleusWeightMat,2) ./ sum(NucleusWeightMat,2) / (StopFrame-StartFrame+1))';
            SpotBulkDyVec = (sum(repmat(NucleiDyVec,length(NewSpotsX),1).*NucleusWeightMat,2) ./ sum(NucleusWeightMat,2) / (StopFrame-StartFrame+1))';        
          else
            SpotBulkDxVec = zeros(size(NewSpotsX));
            SpotBulkDyVec = zeros(size(NewSpotsY));
          end
          
          % Calculate the distances between the spots in this frame
          % and those within the particles in the previous
          % frame          
          
          % Adjust for nuclear movements
          DistanceMat = sqrt((NewSpotsX'-SpotBulkDxVec'- PrevSpotsX).^2 + (NewSpotsY'-SpotBulkDyVec' - PrevSpotsY).^2)*PixelSize;
          % enforce consistent nucleus IDs
          mismatchMat = NewSpotNucleuss'~=PrevNuclei;
          DistanceMat(mismatchMat) = Inf;
          
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
          Particles{Channel}(TotalParticles).Nucleus = NewSpotNucleuss(j); 
          Particles{Channel}(TotalParticles).NucleusDist = NewSpotDistances(j); 
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
      
      for i = 1:length(frameIndex)
        avgZProfile(i) = nanmean(zPosVec(frameVec==frameIndex(i)));
      end
      
      % generate new det-trended z variable
      for p = 1:length(Particles{Channel})
        fVec = Particles{Channel}(p).Frame;
        Particles{Channel}(p).zPosDetrended = Particles{Channel}(p).zPos - avgZProfile(ismember(frameIndex,fVec));
      end
      
    end
    close(f);
  end
end