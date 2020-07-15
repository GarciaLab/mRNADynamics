function StitchedParticles = track04StitchTracks(...
                          RawParticles, Prefix, UseHistone, retrack, displayFigures)
                        
%   UseHistone = false;
  % grab useful info for experiment
  liveExperiment = LiveExperiment(Prefix);
  ExperimentType = liveExperiment.experimentType;
  FrameInfo = getFrameInfo(liveExperiment);
  % set useful parameters
  ncVec = [FrameInfo.nc];
  frameIndex = 1:length(ncVec);
  NCh = length(RawParticles);
  matchCostMaxDefault = 3; % maximum number of sigmas away (this is reset to Inf if we have nuclei)
  matchCostMax = matchCostMaxDefault * ones(1,NCh); % can be different between channels
  
  % Set max spots per nucleus per frame, can be different between channels
  spotsPerNucleus = inf(1,NCh); 
  if ismember(ExperimentType,{'1spot'}) && UseHistone
    spotsPerNucleus = 1;
    matchCostMax = realmax;
  elseif ismember(ExperimentType,{'2spot'}) && UseHistone
    spotsPerNucleus = 2;
    matchCostMax = realmax;
  elseif ismember(ExperimentType,{'2spot2color'}) && UseHistone
    spotsPerNucleus = [1,1];
    matchCostMax = [realmax,realmax];
  elseif ismember(ExperimentType,{'inputoutput'}) && UseHistone
    % Do we have TF clusters in the input channel?
    % No clusters:
    if NCh == 1
        spotsPerNucleus = 1;
        matchCostMax = realmax;
    % Clusters:
    elseif NCh == 2
        % Figure out which channel is the cluster channel
        spotsChannels = liveExperiment.spotChannels;
        inputChannels = liveExperiment.inputChannels;
        clusterChannel = find(spotsChannels == intersect(inputChannels, spotsChannels));
        outputChannel = find(spotsChannels ~= intersect(inputChannels, spotsChannels));
        % Cluster channel not limited in spotsPerNucleus
        spotsPerNucleus(clusterChannel) = Inf; %#ok<FNDSB>
        spotsPerNucleus(outputChannel) = 1; %#ok<FNDSB>
        % Cluster channel should still have a sigma limit on fragment
        % matching
        matchCostMax(clusterChannel) = matchCostMaxDefault;
        matchCostMax(outputChannel) = realmax;
    else
        error('No spot channels, or too many, detected.')
    end
  else
      error(['''',ExperimentType,''' ExperimentType not supported by track04StitchTracks'])
  end
  
  % initialize data structure
  StitchedParticles = cell(1,NCh);
  for Channel = 1:NCh
    
    % number of distinct parameters we're using for linking
    nParams = length(RawParticles{Channel}(1).hmmModel);
    
    % determine nucleus each fragment corresponds to
    nucleusIDVec = NaN(1,length(RawParticles{Channel})); 
    if UseHistone
      for p = 1:length(RawParticles{Channel})   
        nucleusIDVec(p) = RawParticles{Channel}(p).NucleusID(1);
      end
    else
      nucleusIDVec(:) = 1;
    end
    % see how many unique nucleus groups we have
    nucleusIDIndex = unique(nucleusIDVec);    
    % initialize cell structure to temporatily store results for each
    % assignment group
    tempParticles = struct;
    nIter = 1;
    % we only need to perform cost-based tracking within each nucleus group
    f = waitbar(0,'Stitching particle fragments');
    for n = 1:length(nucleusIDIndex)
      waitbar(n/length(nucleusIDIndex),f);
      NucleusID = nucleusIDIndex(n);
      
      [pathArray, sigmaArray, extantFrameArray, particleIDArray, linkIDArray, ...
          linkCostArray, mDistanceMatrix] = performParticleStitching(...
          NucleusID, nucleusIDVec, frameIndex, RawParticles, Channel, ncVec, matchCostMax(Channel));

      % generate local structure to store results       
      assigmentFlags = UseHistone & (sum(extantFrameArray,2)>spotsPerNucleus(Channel))';
      rmVec = [];
      % check for degenerate particle-nucleus assignments
      if any(assigmentFlags)        
        localKernel = 10; % number of leading and trailing frames to examine
        % find problematic frames
        errorIndices = find(assigmentFlags);
        clusterIndices = find([1 diff(errorIndices)>1 1]);
        % initialize vector to track particle IDs to remove        
        % iterate through these and guess which spots are anamolous based
        % on local connectivity
        for e = 1:length(clusterIndices)-1
          % get problematic frame list
          cFrames = errorIndices(clusterIndices(e):clusterIndices(e+1)-1);          
          % get list of correspnding particles
          ptList = nanmax(particleIDArray(cFrames,:),[],1);
          
          ff = max([1,cFrames(1)-localKernel]);
          lf = min([length(frameIndex),cFrames(end)+localKernel]);          
          lcVec = ff:lf;
          lcVec = lcVec(~ismember(lcVec,cFrames));
          localCounts = sum(extantFrameArray(lcVec,:));
          [~,rankVec] = sort(localCounts);
          % add lowest ranking indices
          rmVec = [rmVec ptList(rankVec(1:end-spotsPerNucleus(Channel)))];
          rmVec = rmVec(~isnan(rmVec));
        end
        % reset nucleus ID values for these particles to NaN
        nucleusIDVecNew = nucleusIDVec;
        nucleusIDVecNew(rmVec) = NaN;           
                
        % reset values to originals
        for p = 1:length(rmVec)
          % approval 
          tempParticles(nIter).Approved = false;
          % extant frames
          tempParticles(nIter).Frame = RawParticles{Channel}(rmVec(p)).Frame;
          % position info
          tempParticles(nIter).xPos = RawParticles{Channel}(rmVec(p)).xPos;
          tempParticles(nIter).yPos = RawParticles{Channel}(rmVec(p)).yPos;
          tempParticles(nIter).zPosDetrended = RawParticles{Channel}(rmVec(p)).zPosDetrended;
          % full projected path and error          
          tempParticles(nIter).pathArray = NaN(length(frameIndex),nParams);
          tempParticles(nIter).sigmaArray = NaN(length(frameIndex),nParams);
          nc_ft = ismember(ncVec,ncVec(RawParticles{Channel}(rmVec(p)).FirstFrame));
          for np = 1:nParams
            tempParticles(nIter).pathArray(nc_ft,np) = RawParticles{Channel}(rmVec(p)).hmmModel(np).pathVec;
            tempParticles(nIter).sigmaArray(nc_ft,np) = RawParticles{Channel}(rmVec(p)).hmmModel(np).sigmaVec;        
          end 
          % record info vectors
          tempParticles(nIter).origIDs = rmVec(p);
          tempParticles(nIter).NucleusID = NaN;
          tempParticles(nIter).NucleusIDOrig = NucleusID;
          tempParticles(nIter).linkIDs = zeros(1,length(RawParticles{Channel}(rmVec(p)).Frame));
          tempParticles(nIter).linkCosts = zeros(1,length(RawParticles{Channel}(rmVec(p)).Frame));
          tempParticles(nIter).assigmentFlags = assigmentFlags;
          tempParticles(nIter).NucleusDist = RawParticles{Channel}(rmVec(p)).NucleusDist;
          % increment
          nIter = nIter + 1;
        end 
        
        % aaaaaaand rerun the assignment steps
        [pathArray, sigmaArray, extantFrameArray, particleIDArray, linkIDArray, linkCostArray, mDistanceMatrix]  = performParticleStitching(...
              NucleusID,nucleusIDVecNew,frameIndex,RawParticles,Channel,ncVec,matchCostMax(Channel));
         if size(extantFrameArray,2) ~= spotsPerNucleus(Channel)
           error('problem with spot-nucleus reassigment')
         end
      end
      
      % add particles to structure
      for p = 1:size(extantFrameArray,2)
        % approval 
        tempParticles(nIter).Approved = false;
        % extant frames
        tempParticles(nIter).Frame = find(extantFrameArray(:,p)');
        % position info
        tempParticles(nIter).xPos = pathArray(tempParticles(nIter).Frame,p,1)';
        tempParticles(nIter).yPos = pathArray(tempParticles(nIter).Frame,p,2)';
        tempParticles(nIter).zPosDetrended = pathArray(tempParticles(p).Frame,p,3)';
        % full projected path and error
        tempParticles(nIter).pathArray = reshape(pathArray(:,p,:),[],3);
        tempParticles(nIter).sigmaArray = reshape(sigmaArray(:,p,:),[],3);
        % record info vectors
        tempParticles(nIter).origIDs = particleIDArray(tempParticles(nIter).Frame,p)';
        tempParticles(nIter).linkIDs = linkIDArray(tempParticles(nIter).Frame,p)';
        tempParticles(nIter).linkCosts = linkCostArray(tempParticles(nIter).Frame,p)';
        tempParticles(nIter).NucleusID = NucleusID;
        tempParticles(nIter).NucleusIDOrig = NucleusID;
        tempParticles(nIter).assigmentFlags = assigmentFlags;
        % add other info from original particles
        particleVec = unique(tempParticles(nIter).origIDs,'stable');
        particleVec = particleVec(~isnan(particleVec));  
        ncDist = [];
        for o = particleVec
          ncDist = [ncDist RawParticles{Channel}(o).NucleusDist];
        end
        tempParticles(nIter).NucleusDist = ncDist;
        % increment
        nIter = nIter + 1;
      end   
    end
    close(f)       
    StitchedParticles{Channel} = tempParticles;
  end