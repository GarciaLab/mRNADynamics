function [Particles, trackingOptions] = addQCFlags(Particles, liveExperiment, trackingOptions)

  disp('Adding QC fields...')
  % Iterate over all channels and generate additional QC flags
  % I'm employing a tiered system. Especially egregious cases will be flagge
  % with 2's and will be excluded absent user in put ("opt-in"). Less severe
  % cases will be flaged with 1's and left in unless user removes ("opt-out")
  
  % basic frame info
  FrameInfo = getFrameInfo(liveExperiment);    
  Time = [FrameInfo.Time]; 

  %%% flag spots that occur too early after start of nuclar cycle 
  ncIndex = unique(trackingOptions.ncVec);
  hasNCStart = [0 ones(1,length(ncIndex)-1)];  
  trackingOptions.earlyThresh = 3*60;  
  
  %%% set default thresholds used to flag unlikely points and traces
  trackingOptions.SpotlogLThreshold = repelem(-trackingOptions.matchCostDefault,trackingOptions.NCh);
  trackingOptions.TracelogLThreshold = repelem(-0.75*trackingOptions.matchCostDefault,trackingOptions.NCh);
  
  % Iterate through channels
  for Channel = 1:trackingOptions.NCh
    if ~trackingOptions.useHistone
      trackingOptions.SpotlogLThreshold(Channel) = prctile(vertcat(Particles{Channel}.logL),99);
      trackingOptions.TracelogLThreshold(Channel) = prctile([Particles{Channel}.logLMean],99);
    end
    for p = 1:length(Particles{Channel})      
      
      %%% flag expecially unlikely points according to the motion model
      Particles{Channel}(p).SpotlogLFlags = Particles{Channel}(p).logL(Particles{Channel}(p).obsFrameFilter) < ...
                                                                           trackingOptions.SpotlogLThreshold(Channel);
      Particles{Channel}(p).TracelogLFlag = Particles{Channel}(p).logLMean < trackingOptions.TracelogLThreshold(Channel);
      
      %%% flag early points
      nc = trackingOptions.ncVec(Particles{Channel}(p).Frame(1));
      ncStart = Time(find(trackingOptions.ncVec==nc,1));
      Particles{Channel}(p).earlyFlags = int8(1*(Time(Particles{Channel}(p).Frame)-ncStart)<=...
                    trackingOptions.earlyThresh & hasNCStart(ncIndex==nc));

      %%% automatically disapprove of frames with at least one "2" flag
      Particles{Channel}(p).FrameApproved = Particles{Channel}(p).FrameApproved & ...
        ~Particles{Channel}(p).SpotlogLFlags & ~Particles{Channel}(p).earlyFlags;                   

    end 
    
    
  end
  
  % if we're retracking then we need to cross-reference previously
  % approved qc vectors to see whether we need to transfer any manual
  % user overrides etc. from CheckParticleTracking
%   if retrack
%     % get list of fields to update
%     varNames = fieldnames(Particles{Channel})';
%     flagFields = [{'FrameApproved'},varNames(contains(varNames,'Flags')&~contains(varNames,'Orig')&~contains(varNames,'Per'))];
%     % get previous version of Particles
%     PrevParticles = getParticles(liveExperiment);
%     if ~iscell(PrevParticles)
%       PrevParticles = {PrevParticles};
%     end
%     for Channel = 1:NCh
%       % first find subset of particles that were changed
%       alteredParticles = [];
%       for p = 1:length(PrevParticles{Channel})
%         if any(PrevParticles{Channel}(p).FrameApprovedOrig~=PrevParticles{Channel}(p).FrameApproved)
%           alteredParticles(end+1) = p;
%         end
%       end
%       
%       % now search for matches in new version using x position vectors
%       matchedParticles = NaN(size(alteredParticles));
%       for a = 1:length(alteredParticles)
%         xRef = PrevParticles{Channel}(alteredParticles(a)).xPos;
%         for p = 1:length(Particles{Channel})
%           if all(isequal(Particles{Channel}(p).xPos,xRef))
%             if ~isnan(matchedParticles(a))
%               error('degenerate matching variable')
%             end
%             matchedParticles(a) = p;
%           end
%         end
%       end
%       
%       % update flag fields for matches
%       updateIndices = find(~isnan(matchedParticles(a)));
%       for u = 1:length(updateIndices)
%         for f = 1:length(flagFields)
%           Particles{Channel}(matchedParticles(updateIndices(u))).(flagFields{f}) = ...
%             PrevParticles{Channel}(alteredParticles(updateIndices(u))).(flagFields{f});
%         end      
%       
%         % update summary metrics
%         p = updateIndices(u);
%         %%% define flags-per-frame metric for use in CheckParticleTracking
%         Particles{Channel}(p).flagsPerFrame = ...mean( Particles{Channel}(p).linkCostFlags) + ...
%           mean(Particles{Channel}(p).distShiftFlags>0) + mean(Particles{Channel}(p).fragmentFlags>0) + ...
%           mean(Particles{Channel}(p).ncDistFlags>0)+ mean(Particles{Channel}(p).earlyFlags>0);
% 
%         Particles{Channel}(p).urgentFlagsPerFrame = ...mean( Particles{Channel}(p).linkCostFlags) + ...
%           mean(Particles{Channel}(p).distShiftFlags==2) + mean(Particles{Channel}(p).fragmentFlags==2) + ...
%           mean(Particles{Channel}(p).ncDistFlags==2) + mean(Particles{Channel}(p).earlyFlags==2);
% 
%         %%% automatically disapprove of frames with at least one "2" flag
%         Particles{Channel}(p).FrameApproved = Particles{Channel}(p).FrameApproved & ...
%           Particles{Channel}(p).distShiftFlags~=2 & Particles{Channel}(p).fragmentFlags~=2 & ...
%           Particles{Channel}(p).ncDistFlags~=2 & Particles{Channel}(p).earlyFlags~=2;
%       end
%     end
  end