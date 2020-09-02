function [ParticlesFull,ParticleStitchInfo,SpotFilter] = track04StitchTracks(...
                          ParticlesFull, SpotFilter, ReviewedParticlesFull, ParticleStitchInfo,...
                          Prefix, useHistone, retrack, displayFigures)
                        
  % grab useful info for experiment
  SimParticles = ParticlesFull.SimParticles;
  liveExperiment = LiveExperiment(Prefix);
  ExperimentType = liveExperiment.experimentType;
  FrameInfo = getFrameInfo(liveExperiment);
  % set useful parameters
  ncVec = [FrameInfo.nc];
  frameIndex = 1:length(ncVec);
  NCh = length(SimParticles);
  
  % initialize data structure
  ParticlesFull.StitchedParticles = cell(1,NCh);
  
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
%         matchCostMax(spotsChannels == intersect(inputChannels, spotsChannels)) = 3; % NL: this is a temporary fix for cluster tracking
        spotsPerNucleus(spotsChannels ~= intersect(inputChannels, spotsChannels)) = 1;
    else
        error('No spot channels, or too many, detected.')
    end
  else
      error(['''',ExperimentType,''' ExperimentType not supported by track04StitchTracks'])
  end
  
  for Channel = 1:NCh
    reservedFragmentIDs = ParticleStitchInfo{Channel}.reservedFragmentIDs;
              
    % get full list of pre-assigned links and breaks (will be empty unless
    % retracking)    
    linkStruct = generateLinkStructure(ParticleStitchInfo{Channel},SpotFilter{Channel});
    
    % determine which  nucleus each fragment corresponds to
    NucleusIDVec = NaN(1,length(SimParticles{Channel})); 
    if useHistone
      for p = 1:length(SimParticles{Channel})   
        NucleusIDVec(p) = SimParticles{Channel}(p).Nucleus(1);
      end
    else
      NucleusIDVec(:) = 1;
    end
    FragmentIDVec = [SimParticles{Channel}.FragmentID];    

    % see how many unique nucleus groups we have
    nucleusIDIndex = unique(NucleusIDVec);    
    nucleusIDIndex = nucleusIDIndex(~isnan(nucleusIDIndex));
    
    % initialize cell structure to temporarily store results for each
    % assignment group
    tempParticles = struct;    
    nIter = 1;
    % we only need to perform cost-based tracking within each nucleus group
    wb = waitbar(0,'Stitching particle fragments');
    for n = 1:length(nucleusIDIndex)
      waitbar(n/length(nucleusIDIndex),wb);
      Nucleus = nucleusIDIndex(n);
      
      ncIndices = find(NucleusIDVec==nucleusIDIndex(n));
      [ForceSplitCell, ForceMatchCell] = ...
        checkForAssignedLinkInfo(ncIndices,SimParticles{Channel},linkStruct,SpotFilter{Channel});      
  
      [pathArray, sigmaArray, extantFrameArray, particleIDArray, linkIDCell, ...
              linkCostVec, linkAdditionCell,linkCostCell, linkFrameCell, linkParticleCell] = ...
              performParticleStitching(...
              Nucleus, NucleusIDVec, frameIndex, SimParticles{Channel},  ncVec, matchCostMax(Channel),...
              ForceMatchCell,ForceSplitCell,FragmentIDVec); 
      
      % check for conflicts (cases where there are more detections per frame than ins permitted)     
      nExtantVec = sum(extantFrameArray,2)';
      assignmentFlags = useHistone & (nExtantVec>spotsPerNucleus(Channel));%+length(ForceSplitCell)))';
     
      rmVec = [];
      rmFrameCell = {};
      % check for degenerate particle-nucleus assignments
      if any(assignmentFlags)        
        localKernel1 = 5; % number of leading and trailing frames to examine
        localKernel2 = 15; % number of leading and trailing frames to examine
        % find problematic frames
        errorIndices = find(assignmentFlags);
%         clusterIndices = find([1 diff(errorIndices)>1 1]);
        % initialize vector to track particle IDs to remove        
        % iterate through these and guess which spots are anamolous based
        % on local connectivity
        for e = 1:length(errorIndices)%1:length(clusterIndices)-1
          
          % get problematic frame list
          cFrame = errorIndices(e);%errorIndices(clusterIndices(e):clusterIndices(e+1)-1); 
          activeIndices = find(~isnan(particleIDArray(cFrame,:)));
          % calculate first and last frames over which to conduct
          % comparison on connectivity (i.e. number of detections)
          ff1 = max([1,cFrame(1)-localKernel1]);
          lf1 = min([length(frameIndex),cFrame(end)+localKernel1]);          
          lcVec1 = ff1:lf1;
          lcVec1 = lcVec1(~ismember(lcVec1,cFrame));
          
          ff2 = max([1,cFrame(1)-localKernel2]);
          lf2 = min([length(frameIndex),cFrame(end)+localKernel2]);          
          lcVec2 = ff2:lf2;
          lcVec2 = lcVec2(~ismember(lcVec2,cFrame));
          
          % get counts of linked particles for each conflicting detection
          localCounts1 = sum(extantFrameArray(lcVec1,activeIndices));
          localCounts2 = sum(extantFrameArray(lcVec2,activeIndices));
          [~,rankVec] = sortrows([localCounts1' localCounts2']);
          
          % flag particles with fewest connections within time window and
          % remove
          ptList = particleIDArray(cFrame,activeIndices(rankVec(1:end-spotsPerNucleus(Channel))));          
%           ptList = reshape(unique(particleIDArray(cFrame,rankVec(1:end-spotsPerNucleus(Channel)))),1,[]);          
          rmVec = [rmVec ptList];
          rmFrameCell = [rmFrameCell repelem({cFrame},length(ptList))];
%           rmFrameCell = rmFrameCell(~isnan(rmVec));
%           rmVec = rmVec(~isnan(rmVec));
        end             
        % consolidate rmVec 
        rmIndex = unique(rmVec);
        for r = 1:length(rmIndex)
          rmFrameCellComp(r) = {[rmFrameCell{rmVec==rmIndex(r)}]};
        end
    
        % adjust SimParticles to account for
        for p = 1:length(rmIndex)          
          ptIndex = find(FragmentIDVec==rmIndex(p));
          rmFrames = rmFrameCellComp{p};
          % first check to see how many fragmetn frames fall inside problem
          % region
          overlapFilter = ismember(SimParticles{Channel}(ptIndex).Frame,rmFrames);
          % get frames and indices of offending spots
          overlapFrames = SimParticles{Channel}(ptIndex).Frame(overlapFilter);
          overlapIndices = SimParticles{Channel}(ptIndex).Index(overlapFilter);                    
          
          if all(overlapFilter)
            % Case 1: if the full fragment falls inside region, then remove entirely
            updateIndex = [1:ptIndex-1 ptIndex+1:length(SimParticles{Channel})];
            SimParticles{Channel} = SimParticles{Channel}(updateIndex); 
            NucleusIDVec = NucleusIDVec(updateIndex); 
            FragmentIDVec = FragmentIDVec(updateIndex); 
            % other particle structures
            ParticlesFull.RawParticles{Channel} = ParticlesFull.RawParticles{Channel}(updateIndex);
            ParticlesFull.HMMParticles{Channel} = ParticlesFull.HMMParticles{Channel}(updateIndex);
                          
          elseif sum(~overlapFilter)==1 || max(diff(find(~overlapFilter))) == 1
            % case 2: if a contiguous fragment remains, update entry            
            ParticlesFull.RawParticles{Channel}(ptIndex) = updateVectorFields(ParticlesFull.RawParticles{Channel}(ptIndex),~overlapFilter);

            % now we need to re-generate motion model and path predictions
            globalMotionModel = getGlobalMotionModel(liveExperiment);
            [HMMTemp, ~] = track02TrainGHMM(...
              {ParticlesFull.RawParticles{Channel}(ptIndex)}, globalMotionModel, false);
            
            SimTemp = track03PredictParticlePaths(HMMTemp, FrameInfo, false);
            
            % update structures
            ParticlesFull.HMMParticles{Channel}(ptIndex) = HMMTemp{1};
            SimParticles{Channel}(ptIndex) = SimTemp{1};
          else            
            % case 3: in the unlikely event that a hole is created in the 
            % middle of a contiguous fragment, then we must create a new
            % entry for the trailing fragment
            tempIDVec = bwlabel(~overlapFilter);
            tempIDIndex = unique(tempIDVec(tempIDVec~=0));
            newIndices = [ptIndex length(NucleusIDVec)+(1:max(tempIDVec-1))];
            RawOrig = ParticlesFull.RawParticles{Channel}(ptIndex);
            for id = 1:length(tempIDIndex)
              subFilter = tempIDVec==tempIDIndex(id);
              % first update raw particles               
              ParticlesFull.RawParticles{Channel}(newIndices(id)) = ...
                updateVectorFields(RawOrig,subFilter);

              % change fragment ID if necessary
              ptID = FragmentIDVec(ptIndex);
              if id > 1
                ptID = nanmax(FragmentIDVec)+1;
                while ismember(ptID,reservedFragmentIDs)
                  ptID = ptID + 1;
                end
                FragmentIDVec(newIndices(id)) = ptID;
                NucleusIDVec(newIndices(id)) = Nucleus;
              end
              ParticlesFull.RawParticles{Channel}(newIndices(id)).FragmentID = ptID;
              % then retrain motion model
              globalMotionModel = getGlobalMotionModel(liveExperiment);
              [HMMTemp, ~] = track02TrainGHMM(...
                {ParticlesFull.RawParticles{Channel}(newIndices(id))}, globalMotionModel, false);
              % then make path predictions
              SimTemp = track03PredictParticlePaths(HMMTemp, FrameInfo, false);
              % update structures
              ParticlesFull.HMMParticles{Channel}(newIndices(id)) = HMMTemp{1};
              SimParticles{Channel}(newIndices(id)) = SimTemp{1};              
            end
            
          end  
          % Update SpotFilter
          for f = 1:length(overlapFrames)
            SpotFilter{Channel}(overlapFrames(f),overlapIndices(f)) = 0;
          end  
        end 
  
        % aaaaaaand rerun the assignment steps
        [pathArray, sigmaArray, extantFrameArray, particleIDArray, linkIDCell, ...
              linkCostVec, linkAdditionCell, linkCostCell, linkFrameCell, linkParticleCell] = performParticleStitching(...
              Nucleus, NucleusIDVec, frameIndex, SimParticles{Channel},  ncVec, matchCostMax(Channel),...
              ForceMatchCell,ForceSplitCell,FragmentIDVec); 
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
      
      if length(ParticleStitchInfo{Channel}.linkAdditionIDCell)~=length(ParticleStitchInfo{Channel}.linkAdditionCell)
        error('goddammit')
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
        tempParticles(nIter).Nucleus = Nucleus;
        tempParticles(nIter).NucleusOrig = Nucleus;
        % add other info from original particles
        particleIndexVec = find(ismember(FragmentIDVec,unique(tempParticles(nIter).idVec(~isnan(tempParticles(nIter).idVec)))));
        ncDist = [];
        zOrig = [];
        xOrig = [];
        fOrig = [];
        indexVec = [];
        for o = particleIndexVec
          ncDist = [ncDist SimParticles{Channel}(o).NucleusDist];
          zOrig = [zOrig SimParticles{Channel}(o).zPos];
          xOrig = [xOrig SimParticles{Channel}(o).xPos];
          fOrig = [fOrig SimParticles{Channel}(o).Frame];
          indexVec = [indexVec SimParticles{Channel}(o).Index];
        end        
        [~,si] = sort(fOrig);
        tempParticles(nIter).NucleusDist = ncDist(si);
        tempParticles(nIter).zPos = zOrig(si);
        tempParticles(nIter).Index = indexVec(si);
        % increment
        nIter = nIter + 1;
      end   
    end

    close(wb)       
    % Now, add in approved particles
    if retrack && ~isempty(ReviewedParticlesFull.Particles{Channel})
      tempFields = fieldnames(tempParticles);
      appFields = fieldnames(ReviewedParticlesFull.Particles{Channel});      
      fieldsToAdd = ~ismember(appFields,tempFields);
      mergeStruct = rmfield(ReviewedParticlesFull.Particles{Channel},appFields(fieldsToAdd));   
      tempParticles = [tempParticles mergeStruct];
    end
    ParticlesFull.FullParticles{Channel} = tempParticles;
    ParticlesFull.SimParticles{Channel} = SimParticles{Channel};
  end