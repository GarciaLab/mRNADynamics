function Particles = addQCFields(Particles,useHistone,FrameInfo,retrack,liveExperiment)

  disp('Adding QC fields...')
  % Iterate over all channels and generate additional QC flags
  % I'm employing a tiered system. Especially egregious cases will be flagge
  % with 2's and will be excluded absent user in put ("opt-in"). Less severe
  % cases will be flaged with 1's and left in unless user removes ("opt-out")
  
  
  NCh = length(Particles);
  
  %%% flag large jumps
  Time = [FrameInfo.Time];
  dT = median(diff(Time));
  distThresh1 = 0.75/20; % um/sec 
  distThresh2 = 1.1/20;
  PixelSize = FrameInfo(1).PixelSize;
  % zSize = FrameInfo(1).ZStep;

  %%% flag unlikely linkages 
  % costThresh = 3;
    
  %%% flag spots that are far from their assigned nuclei
  ncDistPrctile = 99.5;

  %%% flag isolated spots
  sizeNeighborhood = 1+2*ceil(120/dT/2);
  slidingIncrement = floor(sizeNeighborhood/2);
  mvWindow = ones(1,sizeNeighborhood);
  frameIndex = 1:length(Time);
  mvThresh = max([1, sum(mvWindow(1:slidingIncrement))-1]);

  %%% flag spots that occur too early after start of nuclar cycle 
  ncVec = [FrameInfo.nc];
  ncIndex = unique(ncVec);
  hasNCStart = [0 ones(1,length(ncIndex)-1)];
  earlyThresh1 = 270;
  earlyThresh2 = 180;
  
  %%% flag implausible frame-over-frame shifts in fluoresnce
  % NL: not currently supported. Need to implement

    
  for Channel = 1:NCh
    if useHistone
      allDistanceVec = [Particles{Channel}.NucleusDist];  
      threshDist1 = prctile(allDistanceVec,ncDistPrctile);
    end
    for p = 1:length(Particles{Channel})
      % flag cases when particle is far away from nearest nucleus
      if useHistone
        Particles{Channel}(p).ncDistFlags = int8(Particles{Channel}(p).NucleusDist>threshDist1);
        Particles{Channel}(p).ncDistFlags(Particles{Channel}(p).NucleusDist>1.25*threshDist1) = 2; % this should almost never happen
      else
        Particles{Channel}(p).ncDistFlags = false(size(Particles{Channel}(p).Frame));
      end
      Particles{Channel}(p).FrameApproved = ones(size(Particles{Channel}(p).Frame));
      %%% flag big jumps
      dt = diff(Time(Particles{Channel}(p).Frame));
      dx = diff(Particles{Channel}(p).xPos)*PixelSize;
      dy = diff(Particles{Channel}(p).yPos)*PixelSize;
  %     dz = diff(Particles{Channel}(p).zPosDetrended)*zSize; % exclude Z for now
      dr = sqrt(dx.^2+dy.^2);%+dz.^2);
      drdt1 = dr./dt>distThresh1 & dt < 120;
      drdt2 = dr./dt>distThresh2 & dt < 120;
      Particles{Channel}(p).distShiftFlags = int8([false drdt1] | [drdt1 false]);
      Particles{Channel}(p).distShiftFlags([false drdt2] | [drdt2 false]) = 2;
      Particles{Channel}(p).distShiftVec = [0 dr];

  %     %%% flag unlikely linkages
  %     if length(Particles{Channel}(p).Frame) > 30
  %       error('afsa')
  %     end
  %     Particles{Channel}(p).linkCostFlags = Particles{Channel}(p).linkCostCell>costThresh;    
  %     Particles{Channel}(p).costThresh = costThresh;
  %     

      %%% flag isolated points    
      FrameVec = Particles{Channel}(p).Frame;    
      frameFlags = zeros(size(frameIndex));
      frameFlags(FrameVec) = 1;
      cvFrames = conv(frameFlags,mvWindow);
      cvFrames = cvFrames(slidingIncrement+1:end-slidingIncrement);
      Particles{Channel}(p).fragmentFlags = int8(cvFrames(FrameVec)<=mvThresh);         
%       Particles{Channel}(p).fragmentFlags(cvFrames(FrameVec)<=1) = 2;
      Particles{Channel}(p).numNeighbors = cvFrames(FrameVec);
      % automatically remove isolated starts and ends
      Particles{Channel}(p).FrameApproved(1) = Particles{Channel}(p).FrameApproved(1) && cvFrames(FrameVec(1))>1;
      Particles{Channel}(p).FrameApproved(1) = Particles{Channel}(p).FrameApproved(end) && cvFrames(FrameVec(1))>1;
      
      %%% flag early points
      nc = ncVec(Particles{Channel}(p).Frame(1));
      ncStart = Time(find(ncVec==nc,1));
      Particles{Channel}(p).earlyFlags = int8(1*(Time(Particles{Channel}(p).Frame)-ncStart)<=earlyThresh1 & hasNCStart(ncIndex==nc))+...
                                         int8(1*(Time(Particles{Channel}(p).Frame)-ncStart)<=earlyThresh2 & hasNCStart(ncIndex==nc));

      %%% define flags-per-frame metric for use in CheckParticleTracking
      Particles{Channel}(p).flagsPerFrame = ...mean( Particles{Channel}(p).linkCostFlags) + ...
        mean(Particles{Channel}(p).distShiftFlags>0) + mean(Particles{Channel}(p).fragmentFlags>0) + ...
        mean(Particles{Channel}(p).ncDistFlags>0)+ mean(Particles{Channel}(p).earlyFlags>0);

      Particles{Channel}(p).urgentFlagsPerFrame = ...mean( Particles{Channel}(p).linkCostFlags) + ...
        mean(Particles{Channel}(p).distShiftFlags==2) + mean(Particles{Channel}(p).fragmentFlags==2) + ...
        mean(Particles{Channel}(p).ncDistFlags==2) + mean(Particles{Channel}(p).earlyFlags==2);

      %%% automatically disapprove of frames with at least one "2" flag
      Particles{Channel}(p).FrameApproved = Particles{Channel}(p).FrameApproved & ...
        Particles{Channel}(p).distShiftFlags~=2 & Particles{Channel}(p).fragmentFlags~=2 & ...
        Particles{Channel}(p).ncDistFlags~=2 & Particles{Channel}(p).earlyFlags~=2;        

      %%% automatically disapprove of particles with more than half "2" frames
      if Particles{Channel}(p).urgentFlagsPerFrame > 0.5
        Particles{Channel}(p).Approved = -1;
%         Particles{Channel}(p).FrameApproved = zeros(size(Particles{Channel}(p).FrameApproved));
      end
      
      %%% also disapprove particles that are entirely fragmented points
      if all(Particles{Channel}(p).fragmentFlags)
        Particles{Channel}(p).Approved = -1;
%         Particles{Channel}(p).FrameApproved = zeros(size(Particles{Channel}(p).FrameApproved));
      end
      
      
      % define static flag vectors to keep track of original states
      Particles{Channel}(p).distShiftFlagsOrig = Particles{Channel}(p).distShiftFlags;
      Particles{Channel}(p).ncDistFlagsOrig = Particles{Channel}(p).ncDistFlags;
      Particles{Channel}(p).fragmentFlagsOrig = Particles{Channel}(p).fragmentFlags;
      Particles{Channel}(p).earlyFlagsOrig = Particles{Channel}(p).earlyFlags;
      Particles{Channel}(p).FrameApprovedOrig = Particles{Channel}(p).FrameApproved;

    end 
    
    
  end
  
  % if we're retracking then we need to cross-reference previously
  % approved qc vectors to see whether we need to transfer any manual
  % user overrides etc. from CheckParticleTracking
  if retrack
    % get list of fields to update
    varNames = fieldnames(Particles{Channel})';
    flagFields = [{'FrameApproved'},varNames(contains(varNames,'Flags')&~contains(varNames,'Orig')&~contains(varNames,'Per'))];
    % get previous version of Particles
    PrevParticles = getParticles(liveExperiment);
    if ~iscell(PrevParticles)
      PrevParticles = {PrevParticles};
    end
    for Channel = 1:NCh
      % first find subset of particles that were changed
      alteredParticles = [];
      for p = 1:length(PrevParticles{Channel})
        if any(PrevParticles{Channel}(p).FrameApprovedOrig~=PrevParticles{Channel}(p).FrameApproved)
          alteredParticles(end+1) = p;
        end
      end
      
      % now search for matches in new version using x position vectors
      matchedParticles = NaN(size(alteredParticles));
      for a = 1:length(alteredParticles)
        xRef = PrevParticles{Channel}(alteredParticles(a)).xPos;
        for p = 1:length(Particles{Channel})
          if all(isequal(Particles{Channel}(p).xPos,xRef))
            if ~isnan(matchedParticles(a))
              error('degenerate matching variable')
            end
            matchedParticles(a) = p;
          end
        end
      end
      
      % update flag fields for matches
      updateIndices = find(~isnan(matchedParticles(a)));
      for u = 1:length(updateIndices)
        for f = 1:length(flagFields)
          Particles{Channel}(matchedParticles(updateIndices(u))).(flagFields{f}) = ...
            PrevParticles{Channel}(alteredParticles(updateIndices(u))).(flagFields{f});
        end      
      
        % update summary metrics
        p = updateIndices(u);
        %%% define flags-per-frame metric for use in CheckParticleTracking
        Particles{Channel}(p).flagsPerFrame = ...mean( Particles{Channel}(p).linkCostFlags) + ...
          mean(Particles{Channel}(p).distShiftFlags>0) + mean(Particles{Channel}(p).fragmentFlags>0) + ...
          mean(Particles{Channel}(p).ncDistFlags>0)+ mean(Particles{Channel}(p).earlyFlags>0);

        Particles{Channel}(p).urgentFlagsPerFrame = ...mean( Particles{Channel}(p).linkCostFlags) + ...
          mean(Particles{Channel}(p).distShiftFlags==2) + mean(Particles{Channel}(p).fragmentFlags==2) + ...
          mean(Particles{Channel}(p).ncDistFlags==2) + mean(Particles{Channel}(p).earlyFlags==2);

        %%% automatically disapprove of frames with at least one "2" flag
        Particles{Channel}(p).FrameApproved = Particles{Channel}(p).FrameApproved & ...
          Particles{Channel}(p).distShiftFlags~=2 & Particles{Channel}(p).fragmentFlags~=2 & ...
          Particles{Channel}(p).ncDistFlags~=2 & Particles{Channel}(p).earlyFlags~=2;
      end
    end
  end