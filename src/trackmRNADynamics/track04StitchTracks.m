function [StitchedParticles,ParticleStitchInfo] = track04StitchTracks(...
                          SimParticles, SpotFilter, ApprovedParticlesFull, ParticleStitchInfo,...
                          Prefix, useHistone, retrack, displayFigures)
                        
%   useHistone = false;
  % grab useful info for experiment
  liveExperiment = LiveExperiment(Prefix);
  ExperimentType = liveExperiment.experimentType;
  FrameInfo = getFrameInfo(liveExperiment);
  % set useful parameters
  ncVec = [FrameInfo.nc];
  frameIndex = 1:length(ncVec);
  NCh = length(SimParticles);
  
  % initialize data structure
  StitchedParticles = cell(1,NCh);
  
  % Set max spots per nucleus per frame, can be different between channels
  matchCostMax = repelem(realmax,NCh);
  if ismember(ExperimentType,{'1spot'}) 
    spotsPerNucleus = 1;    
  elseif ismember(ExperimentType,{'2spot'}) 
    spotsPerNucleus = 2;    
  elseif ismember(ExperimentType,{'2spot2color'}) 
    spotsPerNucleus = [1,1];    
  elseif ismember(ExperimentType,{'inputoutput'}) 
    % Do we have TF clusters in the input channel?
    % No clusters:
    if NCh == 1
        spotsPerNucleus = 1;        
    % Clusters:
    elseif NCh == 2
        % Figure out which channel is the cluster channel
        spotsChannels = liveExperiment.spotChannels;
        inputChannels = liveExperiment.inputChannels;     
        % Cluster channel not limited in spotsPerNucleus
        spotsPerNucleus(spotsChannels == intersect(inputChannels, spotsChannels)) = Inf; 
        spotsPerNucleus(spotsChannels ~= intersect(inputChannels, spotsChannels)) = 1;
    else
        error('No spot channels, or too many, detected.')
    end
  else
      error(['''',ExperimentType,''' ExperimentType not supported by track04StitchTracks'])
  end
  
  for Channel = 1:NCh
    %first, reintegrate approved particles (if they exist)
    idsToExclude = [];
    if retrack 
      % make list of IDs to exclude      
      for r = 1:length(ApprovedParticlesFull.rawParticleIDs{Channel})
        appID = ApprovedParticlesFull.rawParticleIDs{Channel}(r);
        if appID <= length(SimParticles{Channel})
          SimParticles{Channel} = [SimParticles{Channel}(1:appID-1) ...
            ApprovedParticlesFull.SimParticles{Channel}(r) SimParticles{Channel}(appID+1:end)];
          idsToExclude(r) = appID;
        else
          SimParticles{Channel} = [SimParticles{Channel} ApprovedParticlesFull.SimParticles{Channel}(r)];
          idsToExclude(r) = length(SimParticles{Channel});
        end
      end      
    end
      
      
    % get full list of pre-assigned links and breaks (will be empty unless
    % retracking)    
    linkStruct = generateLinkStructure(ParticleStitchInfo{Channel},SpotFilter{Channel});

    % number of distinct parameters we're using for linking
    nParams = length(SimParticles{Channel}(1).hmmModel);
    
    % determine which  nucleus each fragment corresponds to
    nucleusIDVec = NaN(1,length(SimParticles{Channel})); 
    if useHistone
      for p = 1:length(SimParticles{Channel})   
        nucleusIDVec(p) = SimParticles{Channel}(p).Nucleus(1);
      end
    else
      nucleusIDVec(:) = 1;
    end
    nucleusIDVec(idsToExclude) = NaN;
    
    % see how many unique nucleus groups we have
    nucleusIDIndex = unique(nucleusIDVec);    
    nucleusIDIndex = nucleusIDIndex(~isnan(nucleusIDIndex));
    % initialize cell structure to temporarily store results for each
    % assignment group
    tempParticles = struct;    
    nIter = 1;
    % we only need to perform cost-based tracking within each nucleus group
    f = waitbar(0,'Stitching particle fragments');
    for n = 1:length(nucleusIDIndex)
      waitbar(n/length(nucleusIDIndex),f);
      Nucleus = nucleusIDIndex(n);
      
      ncIndices = find(nucleusIDVec==nucleusIDIndex(n));
      [ForceSplitCell, ForceMatchCell] = ...
        checkForAssignedLinkInfo(ncIndices,SimParticles{Channel},linkStruct,SpotFilter{Channel});      
  
      [pathArray, sigmaArray, extantFrameArray, particleIDArray, linkIDCell, ...
              linkCostVec, linkAdditionCell,linkCostCell, linkFrameCell, linkParticleCell] = ...
              performParticleStitching(...
              Nucleus, nucleusIDVec, frameIndex, SimParticles{Channel},  ncVec, matchCostMax(Channel),...
              ForceMatchCell,ForceSplitCell); 

      % check for conflicts (cases where there are more detections per frame than ins permitted)    
      assignmentFlags = useHistone & (sum(extantFrameArray,2)>(spotsPerNucleus(Channel)+length(ForceSplitCell)))';
      rmVec = [];
      % check for degenerate particle-nucleus assignments
      if any(assignmentFlags)        
        localKernel = 10; % number of leading and trailing frames to examine
        % find problematic frames
        errorIndices = find(assignmentFlags);
        clusterIndices = find([1 diff(errorIndices)>1 1]);
        % initialize vector to track particle IDs to remove        
        % iterate through these and guess which spots are anamolous based
        % on local connectivity
        for e = 1:length(clusterIndices)-1
          % get problematic frame list
          cFrames = errorIndices(clusterIndices(e):clusterIndices(e+1)-1); 

          % calculate first and last frames over which to conduct
          % comparison on connectivity (i.e. number of detections)
          ff = max([1,cFrames(1)-localKernel]);
          lf = min([length(frameIndex),cFrames(end)+localKernel]);          
          lcVec = ff:lf;
          lcVec = lcVec(~ismember(lcVec,cFrames));
          % get counts of linked particles for each conflicting detection
          localCounts = sum(extantFrameArray(lcVec,:));
          [~,rankVec] = sort(localCounts);
          % flag lowest ranking particles for removal
          ptList = reshape(unique(particleIDArray(cFrames,rankVec(1:end-spotsPerNucleus(Channel)))),1,[]);          
          rmVec = [rmVec ptList];
          rmVec = rmVec(~isnan(rmVec));
        end
        % reset nucleus ID values for these particles to NaN
%         nucleusIDVecNew = nucleusIDVec;
        nucleusIDVec(rmVec) = NaN;           
                
        % reset values to originals
        for p = 1:length(rmVec)          
          % extant frames
          tempParticles(nIter).Frame = SimParticles{Channel}(rmVec(p)).Frame;
          tempParticles(nIter).FirstFrame = tempParticles(nIter).Frame(1);
          tempParticles(nIter).LastFrame = tempParticles(nIter).Frame(end);
          % approval 
          tempParticles(nIter).Approved = false;
          tempParticles(nIter).FrameApproved = true(size(tempParticles(nIter).Frame));
          % position info
          tempParticles(nIter).xPos = SimParticles{Channel}(rmVec(p)).xPos;
          tempParticles(nIter).yPos = SimParticles{Channel}(rmVec(p)).yPos;
          tempParticles(nIter).zPos = SimParticles{Channel}(rmVec(p)).zPos;
          tempParticles(nIter).zPosDetrended = SimParticles{Channel}(rmVec(p)).zPosDetrended;
          % full projected path and error          
          tempParticles(nIter).pathArray = NaN(length(frameIndex),nParams);
          tempParticles(nIter).sigmaArray = NaN(length(frameIndex),nParams);
          nc_ft = ismember(ncVec,ncVec(SimParticles{Channel}(rmVec(p)).FirstFrame));
          for np = 1:nParams
            tempParticles(nIter).pathArray(nc_ft,np) = SimParticles{Channel}(rmVec(p)).hmmModel(np).pathVec;
            tempParticles(nIter).sigmaArray(nc_ft,np) = SimParticles{Channel}(rmVec(p)).hmmModel(np).sigmaVec;        
          end 
          % record info vectors
          tempParticles(nIter).idVec = NaN(1,size(particleIDArray,1));
          tempParticles(nIter).idVec(tempParticles(nIter).Frame) = rmVec(p);   
          tempParticles(nIter).linkCostCell = [0];
          tempParticles(nIter).linkFrameCell = {unique([SimParticles{Channel}(rmVec(p)).Frame(1) SimParticles{Channel}(rmVec(p)).Frame(end)])};
          tempParticles(nIter).linkParticleCell = {rmVec(p)};
          tempParticles(nIter).stitchInfoPointer = n;
          tempParticles(nIter).Nucleus = NaN;
          tempParticles(nIter).NucleusOrig = Nucleus;
          tempParticles(nIter).linkStateString = num2str(rmVec(p));          
%           tempParticles(nIter).assignmentFlags = assignmentFlags;
          tempParticles(nIter).NucleusDist = SimParticles{Channel}(rmVec(p)).NucleusDist;
          tempParticles(nIter).Index = SimParticles{Channel}(rmVec(p)).Index;
          % increment
          nIter = nIter + 1;
        end 
        
        % aaaaaaand rerun the assignment steps
        [pathArray, sigmaArray, extantFrameArray, particleIDArray, linkIDCell, ...
              linkCostVec, linkAdditionCell, linkCostCell, linkFrameCell, linkParticleCell] = performParticleStitching(...
              Nucleus, nucleusIDVec, frameIndex, SimParticles{Channel},  ncVec, matchCostMax(Channel),...
              ForceMatchCell,ForceSplitCell); 
         if size(extantFrameArray,2) ~= (spotsPerNucleus(Channel) + length(ForceSplitCell))
           error('problem with spot-nucleus reassignment')
         end
      end

      % identify elements in particle array with multiple fragments
      linkAdditionIDCell = {};
      for p = 1:length(linkAdditionCell)
        fragments = strsplit(linkAdditionCell{p},'|');
        linkAdditionIDCell{p} = cellfun(@str2num,fragments);
      end
      
      % update stitch info
      ParticleStitchInfo{Channel}.linkCostVec = [ParticleStitchInfo{Channel}.linkCostVec linkCostVec];
      ParticleStitchInfo{Channel}.linkAdditionCell = [ParticleStitchInfo{Channel}.linkAdditionCell linkAdditionCell];
      ParticleStitchInfo{Channel}.linkAdditionIDCell = [ParticleStitchInfo{Channel}.linkAdditionIDCell linkAdditionIDCell];
      ParticleStitchInfo{Channel}.linkApprovedVec = [ParticleStitchInfo{Channel}.linkApprovedVec repelem(0,length(linkAdditionCell))];

      % add particles to structure
      for p = 1:size(extantFrameArray,2)
        % extant frames
        tempParticles(nIter).Frame = find(extantFrameArray(:,p)');
        tempParticles(nIter).FirstFrame = tempParticles(nIter).Frame(1);
        tempParticles(nIter).LastFrame = tempParticles(nIter).Frame(end);
        % approval 
        tempParticles(nIter).Approved = false;
        tempParticles(nIter).FrameApproved = true(size(tempParticles(nIter).Frame));
        % position info
        tempParticles(nIter).xPos = pathArray(tempParticles(nIter).Frame,p,1)';
        tempParticles(nIter).yPos = pathArray(tempParticles(nIter).Frame,p,2)';        
        tempParticles(nIter).zPosDetrended = pathArray(tempParticles(nIter).Frame,p,3)';
        % full projected path and error
        tempParticles(nIter).pathArray = reshape(pathArray(:,p,:),[],3);
        tempParticles(nIter).sigmaArray = reshape(sigmaArray(:,p,:),[],3);
        % record info vectors
        tempParticles(nIter).idVec = particleIDArray(:,p)';
        tempParticles(nIter).linkCostCell = linkCostCell{p};
        tempParticles(nIter).linkFrameCell = linkFrameCell{p};
        tempParticles(nIter).linkParticleCell = linkParticleCell{p};
        tempParticles(nIter).linkStateString = linkIDCell{p};
        tempParticles(nIter).stitchInfoPointer = n;
        tempParticles(nIter).Nucleus = Nucleus;
        tempParticles(nIter).NucleusOrig = Nucleus;
%         tempParticles(nIter).assignmentFlags = assignmentFlags;        
        % add other info from original particles
        particleVec = unique(tempParticles(nIter).idVec(~isnan(tempParticles(nIter).idVec)));
        ncDist = [];
        zOrig = [];
        indexVec = [];
        for o = particleVec
          ncDist = [ncDist SimParticles{Channel}(o).NucleusDist];
          zOrig = [zOrig SimParticles{Channel}(o).zPos];
          indexVec = [indexVec SimParticles{Channel}(o).Index];
        end
        tempParticles(nIter).NucleusDist = ncDist;
        tempParticles(nIter).zPos = zOrig;
        tempParticles(nIter).Index = indexVec;
        % increment
        nIter = nIter + 1;
      end   
    end
    close(f)       
    % Now, add in approved particles
    if retrack
      tempFields = fieldnames(tempParticles);
      appFields = fieldnames(ApprovedParticlesFull.Particles{Channel});
      fieldsToRemove = ~ismember(appFields,tempFields);
      mergeStruct = rmfield(ApprovedParticlesFull.Particles{Channel},appFields(fieldsToRemove));
      tempParticles = [tempParticles mergeStruct];
    end
    StitchedParticles{Channel} = tempParticles;
  end