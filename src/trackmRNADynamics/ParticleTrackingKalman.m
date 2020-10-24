function Particles = ParticleTrackingKalman(Spots, liveExperiment, ...
          useHistone, retrack, displayFigures)
     
  % get number of channels
  trackingInfo.NCh = length(Spots);
  trackingInfo.maxCost = -2*log(1e-7);
  trackingInfo.maxUnobservedFrames = Inf; % sets limit on how many consecutive frames a particle can be unobserved before "capping" it
  if useHistone
    trackingInfo.maxCost = realmax;
    trackingInfo.maxUnobservedFrames = Inf;
  end
  kfType = 'ConstantVelocity';
  
  
  % get experiment type  
  ExperimentType = liveExperiment.experimentType;
  FrameInfo = getFrameInfo(liveExperiment);
  schnitzcells = getSchnitzcells(liveExperiment);
  
  % get vector indicating stage position in Z
  zPosStage = [FrameInfo.zPosition]*1e6 / FrameInfo(1).ZStep;
  
  % reset histone option to 0 if we have no nucleus data
  useHistone = useHistone && ~isempty(schnitzcells);
        
  % load particles and link info if we're retracking
  Particles = cell(1,trackingInfo.NCh);  
  
  % Extract vector indicating nuclear cleavage cycle for each frame  
  ncVec = [FrameInfo.nc];  
  nFrames = length(ncVec);
  
  for Channel = 1:trackingInfo.NCh        
    
    % determine kalman filter characteristics
    kalmanOptions = determineKalmanOptions(liveExperiment,kfType);       
          
    %% %%%%%%%%%%%%%% Stage 1: Forward tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%
    wb = waitbar(0,['Performing tracking stage 1 (channel ' num2str(Channel) ')']);
    % initialize tracking structure
    forwardTracks = initializeParticleTracks();
    
    % particle ID counter
    for CurrentFrame = 1:nFrames
      waitbar(CurrentFrame/nFrames,wb);
      
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
              schnitzcells,CurrentFrame,useHistone);
         
      % predict positions of extant particles
      forwardTracks = predictParticleTrackLocations(forwardTracks);
      
      earlyFlags = false(size(forwardTracks));
      % Perform cost-based matching
      [assignments, unassignedTracks, unassignedDetections] = ...
                  makeParticleTrackAssignment(forwardTracks, SpotMeasurements, ...
                  continuedNCFlag*trackingInfo.maxCost, NewSpotNuclei,trackingInfo.maxUnobservedFrames,...
                  [], [], earlyFlags);
                                
      % make new entries for spots that were not assigned to existing
      % particles
      forwardTracks = makeNewTracks(forwardTracks, SpotMeasurements,...
                                       unassignedDetections, kalmanOptions,...
                                       CurrentFrame, NewSpotsZ, NewSpotNuclei);
       
      % update existing tracks that had no match this frame
      forwardTracks = updateUnassignedParticleTracks(forwardTracks, unassignedTracks);          
                
      % update tracks that matched with a new particle
      forwardTracks = updateAssignedParticleTracks(...
                              forwardTracks, assignments, SpotMeasurements, CurrentFrame,...
                              NewSpotsZ);                                  
    end
    close(wb);
    %% %%%%%%%%%%%%%% Stage 2: Backwards tracking %%%%%%%%%%%%%%%%%%%%%%%%%
    wb = waitbar(0,['Performing tracking stage 2 (channel ' num2str(Channel) ')']);
    % generate array to track  active status of all particles through time
    activeArray = zeros(nFrames,length(forwardTracks));
    for p = 1:length(forwardTracks)
        activeArray(forwardTracks(p).Frame,p) = 1;
    end
    % generate array with spot indices 
    indexArray = NaN(nFrames,length(forwardTracks));
    for p = 1:length(forwardTracks)
        indexArray(forwardTracks(p).Frame,p) = forwardTracks(p).Index;
    end    
    % generate array to store assignment info 
    assignmentArray = NaN(size(indexArray));
      
    % initialize tracking structure
    backwardTracks = initializeParticleTracks();
    
    % get list of the first detection frames for each track
    firstFrameVec = [forwardTracks.firstFrame];
    
    % particle ID counter
    for CurrentFrame = nFrames:-1:1
      waitbar((nFrames-CurrentFrame+1)/nFrames,wb);
      
      % Get the positions of ALL spots (approved and disapproved)
      [NewSpotsX, NewSpotsY, NewSpotsZ] = SpotsXYZ(Spots{Channel}(CurrentFrame));
      
      % get list of foward tracks that ended in this frame
      activeArrayIndices = find(activeArray(CurrentFrame,:)==1);
      activeSpotIndices = indexArray(CurrentFrame,activeArrayIndices);
      activeParticleIndices = assignmentArray(CurrentFrame,activeArrayIndices);
      assignedParticles = assignmentArray(CurrentFrame,:);      
      % Check if we're at the start of a new nuclear cycle
      continuedNCFlag = true;
      if CurrentFrame < nFrames
        continuedNCFlag = ncVec(CurrentFrame+1)==ncVec(CurrentFrame);
      end
          
      % adjust Z position variable for stage movements
      NewSpotsZAdjusted = NewSpotsZ - zPosStage(CurrentFrame);
      SpotMeasurements = [NewSpotsX', NewSpotsY', NewSpotsZAdjusted'];
      
      % if we have nucleus info, assign each spot to a nucleus      
      [NewSpotNuclei, ~] = getNuclearAssigments(NewSpotsX,NewSpotsY,...
              schnitzcells,CurrentFrame,useHistone);
         
      % predict positions of extant particles
      backwardTracks = predictParticleTrackLocations(backwardTracks);
   
      % get list of tracks that have upcoming assigned aprticles
      ffVec = NaN(size(backwardTracks));
      for p = 1:length(ffVec)
          ffVec(p) = min(firstFrameVec(assignedParticles==p));
      end
      earlyFlags = ffVec<CurrentFrame;
      % Perform cost-based matching
      [assignments, unassignedTracks, unassignedDetections] = ...
                  makeParticleTrackAssignment(backwardTracks, SpotMeasurements, ...
                  continuedNCFlag*trackingInfo.maxCost, NewSpotNuclei,trackingInfo.maxUnobservedFrames,...
                  activeSpotIndices, activeParticleIndices,earlyFlags);
                      
      % update mapping vec     
      if ~isempty(assignments)
        for i = 1:size(assignments,1)          
            colIndex = indexArray(CurrentFrame,:)==assignments(i,2);
            if all(isnan(assignmentArray(:,colIndex)))
                assignmentArray(:,colIndex) = assignments(i,1);
            end
        end     
      end
      newIndVec = length(backwardTracks)+1:length(backwardTracks)+length(unassignedDetections);
      for i = 1:length(unassignedDetections)    
          colIndex = indexArray(CurrentFrame,:)==unassignedDetections(i);         
          assignmentArray(:,colIndex) = newIndVec(i);
      end
      
      % make new entries for spots that were not assigned to existing
      % particles      
      backwardTracks = makeNewTracks(backwardTracks, SpotMeasurements,...
                                       unassignedDetections, kalmanOptions,...
                                       CurrentFrame, NewSpotsZ, NewSpotNuclei);
       
      % update existing tracks that had no match this frame
      backwardTracks = updateUnassignedParticleTracks(backwardTracks, unassignedTracks);          
                
      % update tracks that matched with a new particle
      backwardTracks = updateAssignedParticleTracks(...
                              backwardTracks, assignments, SpotMeasurements, CurrentFrame,...
                              NewSpotsZ);   
      
                            
     
    end
    close(wb);
    
    %% %%%%%%%%%%%%% Step 3: Infer particle trajectory %%%%%%%%%%%%%%%%%%%%
    ParticlesTemp = struct;
    for b = 1:length(backwardTracks)        
      
        [frameVec, frameOrder] = sort(backwardTracks(b).Frame);
        
        %%%%% Add fields to Particles structure %%%%%%%%%%%%%%%%%%%%%%%%%%%
        ParticlesTemp(b).Frame = frameVec;
        ParticlesTemp(b).Index = backwardTracks(b).Index(frameOrder);
        ParticlesTemp(b).xPos = backwardTracks(b).MeasurementVec(frameOrder,1);
        ParticlesTemp(b).yPos = backwardTracks(b).MeasurementVec(frameOrder,2);
        ParticlesTemp(b).zPosDetrended = backwardTracks(b).MeasurementVec(frameOrder,3);
        ParticlesTemp(b).zPos = backwardTracks(b).zPos(frameOrder);
        ParticlesTemp(b).Nucleus = backwardTracks(b).Nucleus;
        ParticlesTemp(b).firstFrame = backwardTracks(b).firstFrame;
        ParticlesTemp(b).lastFrame = backwardTracks(b).lastFrame;
        
        %%%% Infer particle position %%%%%%%%%
        framesFull = frameVec(1):frameVec(end);
        posData = NaN(length(framesFull), size(backwardTracks(b).MeasurementVec,2));
        posData(ismember(framesFull,frameVec),:) = backwardTracks(b).MeasurementVec(frameOrder,:);       

        % call forward filtering function        
        KFTrack = kalmanFilterFwd(posData,kalmanOptions);
        KFTrack = kalmanFilterBkd(KFTrack);        
        
        % Add inferred position info to structure
        ParticlesTemp(b).xPosInf = KFTrack.smoothedTrack(:,1);
        ParticlesTemp(b).yPosInf = KFTrack.smoothedTrack(:,2);
        ParticlesTemp(b).zPosDetrendedInf = KFTrack.smoothedTrack(:,3);
        ParticlesTemp(b).zPosInf = KFTrack.smoothedTrack(:,3) + zPosStage(framesFull);
        ParticlesTemp(b).xPosSEInf = sqrt(KFTrack.smoothedTrackSE(:,1));
        ParticlesTemp(b).yPosSEInf = sqrt(KFTrack.smoothedTrackSE(:,2));
        ParticlesTemp(b).zPosSEInf = sqrt(KFTrack.smoothedTrackSE(:,3));        
    end
    Particles{Channel} = ParticlesTemp;
  end