function Particles = ParticleTrackingKalman(Spots, liveExperiment, ...
          UseHistone, retrack, displayFigures)
     
  % get number of channels
  trackingInfo.NCh = length(Spots);
  trackingInfo.maxCost = -2*log(1e-4);
  kfType = 'ConstantVelocity';
  trackingInfo.maxUnobservedFrames = Inf; % sets limit on how many consecutive frames a particle can be unobserved before "capping" it
  
  % get experiment type  
  ExperimentType = liveExperiment.experimentType;
  FrameInfo = getFrameInfo(liveExperiment);
  schnitzcells = getSchnitzcells(liveExperiment);
  
  % get vector indicating stage position in Z
  zPosStage = [FrameInfo.zPosition]*1e6 / FrameInfo(1).ZStep;
  
  % reset histone option to 0 if we have no nucleus data
  UseHistone = UseHistone && ~isempty(schnitzcells);
        
  % load particles and link info if we're retracking
  Particles = cell(1,trackingInfo.NCh);  
  
  % Extract vector indicating nuclear cleavage cycle for each frame  
  ncVec = [FrameInfo.nc];  
  
  for Channel = 1:trackingInfo.NCh    
    f = waitbar(0,['Stitching particle tracks (channel ' num2str(Channel) ')']);
    
    % determine kalman filter characteristics
    kalmanOptions = determineKalmanOptions(liveExperiment,kfType);       
          
    % initialize tracking structure
    particleTracks = initializeParticleTracks();
    
    % particle ID counter
    for CurrentFrame = 1:length(Spots{Channel})
      waitbar(CurrentFrame/length(Spots{Channel}),f);
      
      % Get the positions of ALL spots (approved and disapproved)
      [NewSpotsX, NewSpotsY, NewSpotsZ] = SpotsXYZ(Spots{Channel}(CurrentFrame));
      
      % Check if we're at the start of a new nuclear cycle
      continuedNCFlag = true;
      if CurrentFrame > 1
        continuedNCFlag = ncVec(CurrentFrame-1)==ncVec(CurrentFrame);
      end
      
      % adjust Z position variable for stage movements
      NewSpotsZAdjusted = NewSpotsZ - zPosStage(CurrentFrame);
      SpotMeasurements = [NewSpotsX', NewSpotsY', NewSpotsZAdjusted'];
      
      % if we have nucleus info, assign each spot to a nucleus
      [NewSpotNuclei, ~] = getNuclearAssigments(NewSpotsX,NewSpotsY,...
              schnitzcells,CurrentFrame,UseHistone);
            
      % predict positions of extant particles
      particleTracks = predictParticleTrackLocations(particleTracks);
      
      % Perform cost-based matching
      [assignments, unassignedTracks, unassignedDetections] = ...
                  makeParticleTrackAssignment(particleTracks, SpotMeasurements, ...
                  continuedNCFlag*trackingInfo.maxCost, NewSpotNuclei,trackingInfo.maxUnobservedFrames);
                                
      % make new entries for spots that were not assigned to existing
      % particles
      particleTracks = makeNewTracks(particleTracks, SpotMeasurements,...
                                       unassignedDetections, kalmanOptions,...
                                       CurrentFrame, NewSpotsZ, NewSpotNuclei);
       
      % update existing tracks that had no match this frame
      particleTracks = updateUnassignedParticleTracks(particleTracks, unassignedTracks);          
                
      % update tracks that matched with a new particle
      particleTracks = updateAssignedParticleTracks(...
                              particleTracks, assignments, SpotMeasurements, CurrentFrame,...
                              NewSpotsZ, NewSpotNuclei);
                            
%       if ~isempty(NewSpotsX)
%         %If spots were detected for this frame we need to add them to the
%         %Particles structure. Link them to existing particles whenever
%         %possible                        
%         particlesFlag = ~isempty(Particles{Channel});
%         
%         
%             
%         %Get a list of the particles that were present in
%         %the previous frame and of their positions.
%         ExtantParticles = [];
%         ExtantParticleIndices = [];
%         PrevSpotsX = [];
%         PrevSpotsY = [];   
%         PrevNuclei = [];
%         if particlesFlag
%           for j = 1:length(Particles{Channel})
%             if Particles{Channel}(j).Frame(end) == (CurrentFrame - 1)
%               ExtantParticles = [ExtantParticles, j];            
%               ExtantParticleIndices = [ExtantParticleIndices, Particles{Channel}(j).Index(end)];            
%               PrevSpotsX = [PrevSpotsX, Particles{Channel}(j).xPos(end)];
%               PrevSpotsY = [PrevSpotsY, Particles{Channel}(j).yPos(end)];
%               PrevNuclei = [PrevNuclei, Particles{Channel}(j).Nucleus];
%             end
%           end
%         end
%         
%         % If there were particles present in the previous frame,
%         % then we find their distances to the spots present in
%         % the current frame. Otherwise, we create new particles
%         % for each new spot.
%         
%         % We keep track of which spots go to a new or old
%         % particle using the NewParticle array.
%         NewParticleFlag = true(size(NewSpotsX));
%         
% 
%         % exclude Spots corresponding to existing particle
%         rmIndices = find(SpotFilterLink(CurrentFrame,:)==0);
%         NewParticleFlag(rmIndices) = false;
%         
%         if ~isempty(ExtantParticles) && ~continuedNCFlag
%           % Generate maximum allowed jump radius between spots          
%           SearchRadius = SearchRadiusMicrons;
%           
%           % if we have nulceus tracking info, calculate average
%           % frame-over-frame shift;
%           maxFrame = find(ncVec==ncVec(CurrentFrame),1,'last');
%           [SpotBulkDxVec,SpotBulkDyVec] = getNuclearShifts(schnitzcells,...
%                     CurrentFrame,maxFrame,NewSpotsX,NewSpotsY,SlidingWindowSize,UseHistone);
%           
%           % Calculate the distances between the spots in this frame
%           % and those within the particles in the previous
%           % frame. Adjust for nuclear movements                       
%           DistanceMat = sqrt((NewSpotsX'-SpotBulkDxVec'- PrevSpotsX).^2 + (NewSpotsY'-SpotBulkDyVec' - PrevSpotsY).^2)*PixelSize;
%           % enforce consistent nucleus IDs
%           mismatchMat = NewSpotNuclei'~=PrevNuclei;
%           DistanceMat(mismatchMat) = Inf;
%           
%           % apply user-assigned links that are within a single frame
%           [PrevJoin, NewJoin] = findPersistentLinks(ParticleStitchInfo{Channel},CurrentFrame);
%           PrevJoinIndices = find(ismember(ExtantParticleIndices,PrevJoin));
%           if ~isempty(PrevJoinIndices)
%             persistentIndices = sub2ind(size(DistanceMat),NewJoin,PrevJoinIndices);          
%             DistanceMat(persistentIndices) = 0;
%           end
%           
%           % apply user-assigned splits that are within a single frame
%           [PrevSplit, NewSplit] = findForbiddenLinks(ParticleStitchInfo{Channel},CurrentFrame);
%           PrevSplitIndices = find(ismember(ExtantParticleIndices,PrevSplit));
%           if ~isempty(PrevSplitIndices)
%             forbiddenIndices = sub2ind(size(DistanceMat),NewSplit,PrevSplitIndices);
%             DistanceMat(forbiddenIndices) = Inf;
%           end
%           
%           % remove existing spots that are in approved particles 
% %           prevAppIndices = SpotFilterTemp(CurrentFrame-1,:)==0;
%           nextAppIndices = find(SpotFilterLink(CurrentFrame,:)==0);
% %           DistanceMat(:,prevAppIndices) = Inf;
%           DistanceMat(nextAppIndices,:) = Inf;
%           
%           % Find existing particles and new spots are close enough to be 
%           % linked. In cases of degenerate assignemnt, take pairt that
%           % minimizes jump distance          
%           [MatchIndices,~,~] = matchpairs(DistanceMat,0.5*SearchRadius);
%           NewParticleFlag(MatchIndices(:,1)) = false;         
% 
%           % Assign matched spots to existing particles
%           for j = 1:size(MatchIndices)
%             ParticleIndex = ExtantParticles(MatchIndices(j,2));
%             NewSpotIndex = MatchIndices(j,1);
%             Particles{Channel}(ParticleIndex).Frame(end + 1) = CurrentFrame;
%             Particles{Channel}(ParticleIndex).Index(end + 1) = NewSpotIndex;            
%             Particles{Channel}(ParticleIndex).LastFrame = CurrentFrame; % not sure I'll use these
%             Particles{Channel}(ParticleIndex).xPos(end + 1) = NewSpotsX(NewSpotIndex);
%             Particles{Channel}(ParticleIndex).yPos(end + 1) = NewSpotsY(NewSpotIndex);
%             Particles{Channel}(ParticleIndex).zPos(end + 1) = NewSpotsZ(NewSpotIndex);      
%             Particles{Channel}(ParticleIndex).NucleusDist(end + 1) = NewSpotDistances(NewSpotIndex);
%           end                        
%         end
%    
%         % See which spots weren't assigned and add them to the structure as new particles
%         NewParticles = find(NewParticleFlag);
%         TotalParticles = 0;
%         if particlesFlag
%           TotalParticles = length(Particles{Channel});
%         end
%         for j = NewParticles    
%           TotalParticles = TotalParticles + 1;
%           Particles{Channel}(TotalParticles).Frame = CurrentFrame;
%           Particles{Channel}(TotalParticles).Index = j;
%           Particles{Channel}(TotalParticles).Approved = 0;
%           Particles{Channel}(TotalParticles).FirstFrame = CurrentFrame; 
%           Particles{Channel}(TotalParticles).LastFrame = CurrentFrame; 
%           Particles{Channel}(TotalParticles).xPos = NewSpotsX(j);
%           Particles{Channel}(TotalParticles).yPos = NewSpotsY(j);
%           Particles{Channel}(TotalParticles).zPos = NewSpotsZ(j);      
%           Particles{Channel}(TotalParticles).Nucleus = NewSpotNuclei(j); 
%           Particles{Channel}(TotalParticles).NucleusDist = NewSpotDistances(j);
%           % add identifier
%           while ismember(ptIDCounter,reservedFragmentIDs)
%             ptIDCounter = ptIDCounter + 1;
%           end
%           Particles{Channel}(TotalParticles).FragmentID = ptIDCounter;
%           ptIDCounter = ptIDCounter + 1;
%         end
%       end
% 
%     end
%     
%     % adjust for z stack shifts
%     Particles{Channel} = detrendZ(Particles{Channel},FrameInfo);
%     close(f);
    end
  end