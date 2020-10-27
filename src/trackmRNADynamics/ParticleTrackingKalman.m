function [Particles, trackingOptions] = ParticleTrackingKalman(Spots, liveExperiment, trackingOptions)
         
  % get experiment type  
  FrameInfo = getFrameInfo(liveExperiment);
  schnitzcells = getSchnitzcells(liveExperiment); % this will be empty and ignored if no nucleus info    
  
  % get vector indicating stage position in Z
  zPosStage = [FrameInfo.zPosition]*1e6 / FrameInfo(1).ZStep;    
  trackingOptions.zPosStage = zPosStage;
  
  % load particles and link info if we're retracking
  Particles = cell(1,trackingOptions.NCh);  
  
  % Extract vector indicating nuclear cleavage cycle for each frame  
  trackingOptions.ncVec = [FrameInfo.nc];  
  
  for Channel = 1:trackingOptions.NCh        
    
    % determine kalman filter characteristics
    kalmanOptions = determineKalmanOptions(liveExperiment,trackingOptions,Spots{Channel});       
    
    %% %%%%%%%%%%%%%% Stage 0: Calibrate matching cost %%%%%%%%%%%%%%%%%%%%           
    % Determine number of detections per frame
%     nDetections = zeros(1,length(Spots{Channel}));
%     for f = 1:length(Spots{Channel})
%       nDetections(f) = length(Spots{Channel}(f).Fits);
%     end
%     [maxN, maxF] = max(nDetections);
%     maxF = 30;
%     testTracks = initializeParticleTracks();
%     for CurrentFrame = 1:maxF+1       
% 
%         [testTracks, trackingInfo] = forwardTrackingLoop(testTracks, trackingInfo, ...
%             kalmanOptions, Spots, Channel, CurrentFrame, schnitzcells, CurrentFrame==maxF+1);  
%                                   
%     end        
    %% %%%%%%%%%%%%%% Stage 1: Forward tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%           
    
    % initialize tracking structure
    forwardTracks = initializeParticleTracks();
    
    % track fraction of unassigned detections over time
    trackingOptions.fractionUnassigned = NaN(1,trackingOptions.nFrames);
    
    wb = waitbar(0,['Stage 1: Forward particle linking (channel ' num2str(Channel) ')']);
    % iterate throug frames
    for CurrentFrame = 1:trackingOptions.nFrames
        waitbar(CurrentFrame/trackingOptions.nFrames,wb);

        forwardTracks = forwardTrackingLoop(forwardTracks, trackingOptions, ...
            kalmanOptions, Spots, Channel, CurrentFrame, schnitzcells, 0);  
                                  
    end
    close(wb);
    
    %% %%%%%%%%%%%%%% Stage 2: Backwards tracking %%%%%%%%%%%%%%%%%%%%%%%%%
    wb = waitbar(0,['Stage 2: Backwards particle linking (channel ' num2str(Channel) ')']);
    % generate array to track  active status of all particles through time
    trackingOptions.activeArray = zeros(trackingOptions.nFrames,length(forwardTracks));
    for p = 1:length(forwardTracks)
        trackingOptions.activeArray(forwardTracks(p).Frame,p) = 1;
    end
    % generate array with spot indices 
    trackingOptions.indexArray = NaN(trackingOptions.nFrames,length(forwardTracks));
    for p = 1:length(forwardTracks)
        trackingOptions.indexArray(forwardTracks(p).Frame,p) = forwardTracks(p).Index;
    end    
    % generate array to store assignment info 
    trackingOptions.assignmentArray = NaN(size(trackingOptions.indexArray));
      
    % initialize tracking structure
    backwardTracks = initializeParticleTracks();
    
    % get list of the first detection frames for each track
    trackingOptions.firstFrameVec = [forwardTracks.firstFrame];
    
    % particle ID counter
    for CurrentFrame = trackingOptions.nFrames:-1:1
        waitbar((trackingOptions.nFrames-CurrentFrame+1)/trackingOptions.nFrames,wb);

        [backwardTracks, trackingOptions] = backwardTrackingLoop(backwardTracks, trackingOptions, ...
            kalmanOptions, Spots, Channel, CurrentFrame, schnitzcells);
      
    end
    close(wb);
    
    % remove particle tracks that were flagged as duplicates
    backwardTracks = backwardTracks(~[backwardTracks.duplicateFlag]);
    
    %% %%%%%%%%%%%%% Step 3: Infer particle trajectory %%%%%%%%%%%%%%%%%%%%
    wb = waitbar(0,['Stage 3: Performing path prediction (channel ' num2str(Channel) ')']);
    ParticlesTemp = initializeParticles();
    for b = 1:length(backwardTracks)        
        waitbar(b/length(backwardTracks),wb);
        [frameVec, frameOrder] = sort(backwardTracks(b).Frame);
        
        %%%%% Add  standard fields to Particles structure %%%%%%%%%%%%%%%%%
        ParticlesTemp(b).Frame = frameVec;
        ParticlesTemp(b).Index = backwardTracks(b).Index(frameOrder);
        ParticlesTemp(b).xPos = backwardTracks(b).MeasurementVec(frameOrder,1);
        ParticlesTemp(b).yPos = backwardTracks(b).MeasurementVec(frameOrder,2);
        ParticlesTemp(b).zPosDetrended = backwardTracks(b).MeasurementVec(frameOrder,3);
        ParticlesTemp(b).zPos = backwardTracks(b).zPos(frameOrder);
        ParticlesTemp(b).Nucleus = backwardTracks(b).Nucleus;
        ParticlesTemp(b).FirstFrame = backwardTracks(b).firstFrame;
        ParticlesTemp(b).LastFrame = backwardTracks(b).lastFrame;
        ParticlesTemp(b).Approved = 0;
        ParticlesTemp(b).FrameApproved = zeros(size(ParticlesTemp(b).Frame));
        
        % make particle path predictions
        ParticlesTemp(b) = pathPrediction(ParticlesTemp(b), backwardTracks(b), trackingOptions, kalmanOptions);
             
    end    
    Particles{Channel} = ParticlesTemp;
    kalmanOptions.Channel = Channel;
    trackingOptions.kalmanOptions(Channel) = kalmanOptions;
    close(wb);
  end