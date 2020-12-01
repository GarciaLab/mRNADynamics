function [Particles, trackingOptions] = addQCFlags(Particles, liveExperiment, trackingOptions)

  % specify which variables to use for likelihood score calculations
  % just use position variables for now
  logLIndices = 1:3;
  
  startParallelPool(8, 0, 1);
  % Iterate over all channels and generate additional QC flags  
  FrameInfo = getFrameInfo(liveExperiment);
  
  %%% load nucleus probability maps if they exist
  nucleusProbDirRaw = [liveExperiment.procFolder 'nucleusProbabilityMaps' filesep];
  nucleusProbDirFinal = [liveExperiment.procFolder 'nucleusProbabilityMapsFull' filesep];
  probFiles = dir([nucleusProbDirFinal '*.tif']);
  probFilesRaw = dir([nucleusProbDirRaw '*.tif']);
  hasFullProbFiles = length(probFiles) == length(FrameInfo);
  
  if hasFullProbFiles
      disp('Incorporating nucleus probabilities...')
      Particles = addNucleusProbabilities(liveExperiment, trackingOptions, FrameInfo, Particles);      
  else
      if exist([nucleusProbDirRaw 'nucleusInfo.mat'], 'file')
          load([nucleusProbDirRaw 'nucleusInfo.mat'], 'nucleusInfo')
          if length(nucleusInfo.originalFileNames) == length(probFilesRaw)
              compileNuclearTiffStacks(liveExperiment.Prefix);
              disp('Incorporating nucleus probabilities...')
              Particles = addNucleusProbabilities(liveExperiment, trackingOptions, FrameInfo, Particles);      
          end
      else
          warning('No nucleus probabilities found. Nucleus boundaries will not be used to assess particle quality.')
          for Channel = 1:trackingOptions.NCh
              for  p = 1:length(Particles{Channel})
                  Particles{Channel}(p).nucleusProbability = ones(size(Particles{Channel}(p).Frame));
              end
          end
      end
  end
    
  
  % basic frame info
  Time = [FrameInfo.Time]; 

  %%% flag spots that occur too early after start of nuclar cycle 
  ncIndex = unique(trackingOptions.ncVec);
  hasNCStart = [0 ones(1,length(ncIndex)-1)];  
  trackingOptions.earlyThresh = 3*60;  
  
  %%% set default thresholds used to flag unlikely points and traces
  trackingOptions.SpotlogLThreshold = repelem(-trackingOptions.matchCostDefault,trackingOptions.NCh);
  trackingOptions.SpotlogLThreshholdSingle = trackingOptions.SpotlogLThreshold/length(logLIndices);
  trackingOptions.TracelogLThreshold = repelem(-0.75*trackingOptions.matchCostDefault,trackingOptions.NCh);   
  trackingOptions.qcFieldNames = {'xPos','yPos','zPos','nucleusMask','earlyflag'};
  trackingOptions.nucleusProbThresh = 0.40; %NL: we can afford to be a little generous here
  
  disp('Adding QC fields...')
  % Iterate through channels
  for Channel = 1:trackingOptions.NCh
    
    % calcualte average likelihood quantities
%     for p = 1:length(Particles{Channel})
%         Particles{Channel}(p).logL = nansum(Particles{Channel}(p).logLDistance(:,logLIndices),2);
%         Particles{Channel}(p).logLMean = nanmean(Particles{Channel}(p).logL);
%     end
    
%     if true%~trackingOptions.useHistone
    trackingOptions.SpotlogLThreshold(Channel) = 75;%prctile(-vertcat(Particles{Channel}.logL),99);
    trackingOptions.SpotlogLThresholdSingle(Channel) = trackingOptions.SpotlogLThreshold(Channel)/length(logLIndices);
    trackingOptions.TracelogLThreshold(Channel) = prctile(-[Particles{Channel}.logLMean],99);
%     end
    
      for p = 1:length(Particles{Channel})      

          %%% flag particles that appear to be outside of nuclear envelope
          Particles{Channel}(p).NucleusBoundaryFlags = Particles{Channel}(p).nucleusProbability < trackingOptions.nucleusProbThresh;

          %%% flag expecially unlikely points according to the motion model
          nFrames = length(Particles{Channel}(p).Frame);
          Particles{Channel}(p).SpotlogLFlags = (-Particles{Channel}(p).logL(Particles{Channel}(p).obsFrameFilter)) > ...
                                                                               trackingOptions.SpotlogLThreshold(Channel);

          Particles{Channel}(p).SpotlogLFlags = Particles{Channel}(p).SpotlogLFlags&(nFrames>2);

          Particles{Channel}(p).TracelogLFlag = (-Particles{Channel}(p).logLMean) > trackingOptions.TracelogLThreshold(Channel);
          Particles{Channel}(p).TracelogLFlag = Particles{Channel}(p).TracelogLFlag&&(nFrames>2);

          %%% flag early points
          nc = trackingOptions.ncVec(Particles{Channel}(p).Frame(1));
          ncStart = Time(find(trackingOptions.ncVec==nc,1));     

          Particles{Channel}(p).earlyFlags = int8(1*(Time(Particles{Channel}(p).Frame)-ncStart)<=...
                        trackingOptions.earlyThresh & hasNCStart(ncIndex==nc));

          %%% automatically disapprove of frames   
          Particles{Channel}(p).FrameApproved = Particles{Channel}(p).FrameApproved & ...
            ~Particles{Channel}(p).SpotlogLFlags' & ~Particles{Channel}(p).earlyFlags & ...
            ~Particles{Channel}(p).NucleusBoundaryFlags;                   

          %%% calculatemean number of disapproved frames
          Particles{Channel}(p).FlaggedFraction = mean(~Particles{Channel}(p).FrameApproved);

          %%% geenrate array to use for QC plot
          timeDeltaVec = 2*trackingOptions.earlyThresh - (Time(Particles{Channel}(p).Frame)'-ncStart);
          timeDeltaVec(timeDeltaVec<0) = 0;
          Particles{Channel}(p).qcScoreArray = [-Particles{Channel}(p).logLArray(Particles{Channel}(p).obsFrameFilter,logLIndices)...
                                                1-Particles{Channel}(p).nucleusProbability'...
                                                timeDeltaVec*hasNCStart(ncIndex==nc)];   
          Particles{Channel}(p).qcScoreArray(isnan(Particles{Channel}(p).qcScoreArray)) = 0;                                  
          Particles{Channel}(p).qcThreshVec = [repelem(trackingOptions.SpotlogLThresholdSingle(Channel),length(logLIndices)),1-trackingOptions.nucleusProbThresh,trackingOptions.earlyThresh];

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