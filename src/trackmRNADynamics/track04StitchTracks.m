function [StitchedParticles,ParticleStitchInfo] = track04StitchTracks(...
                          SimParticles, Prefix, useHistone, retrack, displayFigures)
                        
%   useHistone = false;
  % grab useful info for experiment
  liveExperiment = LiveExperiment(Prefix);
  ExperimentType = liveExperiment.experimentType;
  FrameInfo = getFrameInfo(liveExperiment);
  % set useful parameters
  ncVec = [FrameInfo.nc];
  frameIndex = 1:length(ncVec);
  NCh = length(SimParticles);
%   matchCostMaxDefault = 3; % maximum number of sigmas away (this is reset to Inf if we have nuclei)
%   matchCostMax = matchCostMaxDefault * ones(1,NCh); % can be different between channels
  
  % simplify the logic a bit
%   traceCost = matchCostMaxDefault;
%   if useHistone
%     traceCost = realmax;
%   end
  
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
  
  % initialize data structure
  StitchedParticles = cell(1,NCh);
  ParticleStitchInfo = cell(1,NCh);
  for Channel = 1:NCh
    
    % number of distinct parameters we're using for linking
    nParams = length(SimParticles{Channel}(1).hmmModel);
    
    % determine which  nucleus each fragment corresponds to
    nucleusIDVec = NaN(1,length(SimParticles{Channel})); 
    if useHistone
      for p = 1:length(SimParticles{Channel})   
        nucleusIDVec(p) = SimParticles{Channel}(p).NucleusID(1);
      end
    else
      nucleusIDVec(:) = 1;
    end
    % see how many unique nucleus groups we have
    nucleusIDIndex = unique(nucleusIDVec);    
    % initialize cell structure to temporarily store results for each
    % assignment group
    tempParticles = struct;    
    nIter = 1;
    % we only need to perform cost-based tracking within each nucleus group
    f = waitbar(0,'Stitching particle fragments');
    for n = 1:length(nucleusIDIndex)
      waitbar(n/length(nucleusIDIndex),f);
      NucleusID = nucleusIDIndex(n);
      
      [pathArray, sigmaArray, extantFrameArray, particleIDArray, linkIDCell, ...
              linkCostVec, linkAdditionCell,linkCostCell, linkFrameCell, linkParticleCell] = performParticleStitching(...
          NucleusID, nucleusIDVec, frameIndex, SimParticles, Channel, ncVec, matchCostMax(Channel));

      % check for conflicts (cases where there are more detections per frame than ins permitted)    
      assignmentFlags = useHistone & (sum(extantFrameArray,2)>spotsPerNucleus(Channel))';
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
        nucleusIDVecNew = nucleusIDVec;
        nucleusIDVecNew(rmVec) = NaN;           
                
        % reset values to originals
        for p = 1:length(rmVec)
          % approval 
          tempParticles(nIter).Approved = false;
          % extant frames
          tempParticles(nIter).Frame = SimParticles{Channel}(rmVec(p)).Frame;
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
          tempParticles(nIter).linkCostCell = {0};
          tempParticles(nIter).linkFrameCell = {unique([SimParticles{Channel}(rmVec(p)).Frame(1) SimParticles{Channel}(rmVec(p)).Frame(end)])};
          tempParticles(nIter).linkParticleCell = {rmVec(p)};
          tempParticles(nIter).NucleusID = NaN;
          tempParticles(nIter).NucleusIDOrig = NucleusID;
          tempParticles(nIter).linkStateString = num2str(rmVec(p));          
          tempParticles(nIter).assignmentFlags = assignmentFlags;
          tempParticles(nIter).NucleusDist = SimParticles{Channel}(rmVec(p)).NucleusDist;
          % increment
          nIter = nIter + 1;
        end 
        
        % aaaaaaand rerun the assignment steps
        [pathArray, sigmaArray, extantFrameArray, particleIDArray, linkIDCell, ...
              linkCostVec, linkAdditionCell, linkCostCell, linkFrameCell, linkParticleCell] = performParticleStitching(...
              NucleusID,nucleusIDVecNew,frameIndex,SimParticles,Channel,ncVec,matchCostMax(Channel));
         if size(extantFrameArray,2) ~= spotsPerNucleus(Channel)
           error('problem with spot-nucleus reassignment')
         end
      end
      % update stitch info
      ParticleStitchInfo{Channel}(n).linkCostVec = linkCostVec;
      ParticleStitchInfo{Channel}(n).linkAdditionCell = linkAdditionCell;
      % add particles to structure
      for p = 1:size(extantFrameArray,2)
        % approval 
        tempParticles(nIter).Approved = false;
        % extant frames
        tempParticles(nIter).Frame = find(extantFrameArray(:,p)');
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
        tempParticles(nIter).NucleusID = NucleusID;
        tempParticles(nIter).NucleusIDOrig = NucleusID;
        tempParticles(nIter).assignmentFlags = assignmentFlags;        
        % add other info from original particles
        particleVec = unique(tempParticles(nIter).idVec(~isnan(tempParticles(nIter).idVec)));
        ncDist = [];
        zOrig = [];
        for o = particleVec
          ncDist = [ncDist SimParticles{Channel}(o).NucleusDist];
          zOrig = [zOrig SimParticles{Channel}(o).zPos];
        end
        tempParticles(nIter).NucleusDist = ncDist;
        tempParticles(nIter).zPos = zOrig;
        % increment
        nIter = nIter + 1;
      end   
    end
    close(f)       
    StitchedParticles{Channel} = tempParticles;
  end