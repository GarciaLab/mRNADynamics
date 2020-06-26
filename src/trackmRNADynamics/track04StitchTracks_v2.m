function StitchedParticles = track04StitchTracks(...
                          RawParticles, FrameInfo, ExperimentType, retrack, displayFigures)

  % set useful parameters
  NCh = length(RawParticles);
  ncVec = [FrameInfo.nc];
  frameIndex = 1:length(ncVec);
  maxDist = sum(ncVec==mode(ncVec)); % max time gap between linkable points (in frames)
  matchCostMax = 3; % maximum number of sigmas away (the 0.5 factor is an adjustment for 
  matchCostMin = 0.5; % cost for initial pair matching (pairs must be within 1/2 sigma to be matched)
  nMatchRounds = sum(ncVec==mode(ncVec)); % number of distinct match costs to iterate through
  costVec = linspace(matchCostMin,matchCostMax,nMatchRounds);  
  
  % initialize data structure
  StitchedParticles = cell(1,NCh);
  for Channel = 1:NCh 
    % check to see if we have nucleus info
    UseNuclei = all([RawParticles{Channel}.NucleusID]==1);
    % record extant frames
    rightPointVec = NaN(length(RawParticles{Channel}),1);    
    leftPointVec = NaN(length(RawParticles{Channel}),1);
    segmentIDVec = NaN(length(RawParticles{Channel}),1);
    for p = 1:length(RawParticles{Channel})
      rightPointVec(p) = RawParticles{Channel}(p).LastFrame;
      leftPointVec(p) = RawParticles{Channel}(p).FirstFrame;
      segmentIDVec(p) = p;
    end
    
    % now store particle positions and positional errors in arrays
    nDims = length(RawParticles{Channel}(1).hmmModel);
    pathArray = Inf(length(frameIndex),length(RawParticles{Channel}),nDims);   
    sigmaArray = Inf(length(frameIndex),length(RawParticles{Channel}),nDims); 
    NucleusIDVec = NaN(1,length(RawParticles{Channel}));
    for p = 1:length(RawParticles{Channel})      
      nc_ft = ismember(ncVec,ncVec(RawParticles{Channel}(p).FirstFrame));
      for n = 1:nDims
        pathArray(nc_ft,p,n) = RawParticles{Channel}(p).hmmModel(n).pathVec;
        sigmaArray(nc_ft,p,n) = RawParticles{Channel}(p).hmmModel(n).sigmaVec;        
      end
      NucleusIDVec(p) = RawParticles{Channel}(p).NucleusID(1);
    end
    extantArray = sigmaArray(:,:,1)==0;

    % for each particle generate list of overlapping (incompatible)
    % particles
    avoidanceCell = cell(1,length(RawParticles{Channel}));
    nucleusIDCell = cell(1,length(RawParticles{Channel}));
    activeFrameCell = cell(1,length(RawParticles{Channel}));
    origIDCell = cell(1,length(RawParticles{Channel}));
    linkIDCell = cell(1,length(RawParticles{Channel}));
    costCell = cell(1,length(RawParticles{Channel}));
    particleIndexCell = cell(1,length(RawParticles{Channel}));
    for p = 1:length(RawParticles{Channel})
      overlapVec = max(extantArray(extantArray(:,p),:),[],1)==1;
      ncMismatchVec = NucleusIDVec~=NucleusIDVec(p);
      nucleusIDCell{p} = RawParticles{Channel}(p).NucleusID;
      avoidanceCell{p} = find(overlapVec);
      activeFrameCell{p} = RawParticles{Channel}(p).Frame;
      origIDCell{p} = repelem(p,length(activeFrameCell{p}));
      linkIDCell{p} = repelem(1,length(activeFrameCell{p}));
      costCell{p} = repelem(0,length(activeFrameCell{p}));
      particleIndexCell{p} = RawParticles{Channel}(p).Index;
    end
    %%% calculate Mahalanobis Distance between particles in likelihood space
    forwardDistanceMat = Inf(length(RawParticles{Channel}),length(RawParticles{Channel}));   
    backwardDistanceMat = Inf(length(RawParticles{Channel}),length(RawParticles{Channel}));   
    for p = 1:length(RawParticles{Channel})
      % calculate forward distances
      rpFrame = rightPointVec(p);
      rpDelta = pathArray(rpFrame,p,:)-pathArray(rpFrame,:,:);      
      lpSigVec = sigmaArray(rpFrame,:,:);
      forwardDistanceMat(:,p) = mean((rpDelta ./ lpSigVec).^2,3);
      forwardDistanceMat(leftPointVec-rpFrame<=0|leftPointVec-rpFrame>maxDist,p) = Inf;
      forwardDistanceMat(avoidanceCell{p},p) = Inf;
      % calculate backwards distances 
      lpFrame = leftPointVec(p);
      lpDelta = pathArray(lpFrame,p,:) - pathArray(lpFrame,:,:);      
      rpSigVec = sigmaArray(lpFrame,:,:);
      backwardDistanceMat(p,:) = mean((lpDelta ./ rpSigVec).^2,3);
      backwardDistanceMat(p,rightPointVec-lpFrame>=0|rightPointVec-lpFrame<-maxDist) = Inf;
      backwardDistanceMat(p,avoidanceCell{p}) = Inf;
    end    
    % if there are a fixed number of spots per nucleus, increase the max
    mDistanceMat = sqrt((backwardDistanceMat + forwardDistanceMat)/2);
    mDistanceMat(isnan(mDistanceMat)) = Inf;
    
    % cost such that all fragments are matched
    if ismember(ExperimentType,{'inputoutput','1spot','2spot'}) %&& UseNuclei
      costVec = linspace(matchCostMin,1.1*max(mDistanceMat(~isinf(mDistanceMat))),nMatchRounds);
    end
    
    % iterate through cost layers
    f = waitbar(0,['Stitching particle tracks (channel ' num2str(Channel) ')']);
    for c = 1:length(costVec)
      % udate waitbar
      waitbar(c/nMatchRounds,f);
      
      % perform pair matching using optimal linear assignment      
      [matchArray,~,~] = matchpairs(mDistanceMat,0.5*costVec(c)); % note that the 0.5 factor is intentional here

      % update path arrays and info cells  
      alteredParticles = [];
      for m = 1:size(matchArray,1)
        % extract match indices
        pKeep = matchArray(m,2);
        pDrop = matchArray(m,1);
        cost = mDistanceMat(pDrop,pKeep);
        if ismember(pKeep,avoidanceCell{pDrop})
          error('wtf')
        end
        % update ID vec
        toID = segmentIDVec(pKeep);
        fromID = segmentIDVec(pDrop);
        segmentIDVec(segmentIDVec==fromID) = toID;
        toFilter = segmentIDVec==toID;
        alteredParticles = [alteredParticles find(toFilter)'];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % update info tracking cells
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        avoidanceCell(toFilter) = {unique([avoidanceCell{[pKeep pDrop]}])}; % tracks incompatible segments
       
        [activeFrames,si] = unique([activeFrameCell{[pKeep pDrop]}]); % tracks active frames             
        
        ptIndices = [particleIndexCell{[pKeep pDrop]}]; % tracks index of particles within Spots
        particleIndexCell(toFilter) = {ptIndices(si)};
        
        nucleusIDs = [nucleusIDCell{[pKeep pDrop]}];        
        nucleusIDCell(toFilter) = {nucleusIDs(si)};

        currentLinks = [linkIDCell{pKeep} repelem(max(linkIDCell{pKeep})+1,length(activeFrameCell{pDrop}))]; % track ordering of links
        linkIDCell(toFilter) = {currentLinks(si)};

        currentCosts = [costCell{pKeep} repelem(cost,length(activeFrameCell{pDrop}))]; % linking cosst
        costCell(toFilter) = {currentCosts(si)};
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate updated combined path using variance-weighted average 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        firstFrame = activeFrames(1);
        lastFrame = activeFrames(end);
        midVec = firstFrame:lastFrame;
        gapVec = midVec(~ismember(midVec,activeFrames));
        
        % calculate variance weights
        sg1 = sigmaArray(:,pKeep,:).^-2;
        sg2 = sigmaArray(:,pDrop,:).^-2;
        
        % first assign known points
        newPath = NaN(size(sigmaArray(:,1,:)));
        newPath(activeFrameCell{pDrop},:,:) = pathArray(activeFrameCell{pDrop},pDrop,:);
        newPath(activeFrameCell{pKeep},:,:) = pathArray(activeFrameCell{pKeep},pKeep,:);
        
        pathArray(activeFrameCell{pDrop},toFilter,:) = repmat(pathArray(activeFrameCell{pDrop},pDrop,:),1,sum(toFilter),1);
        pathArray(activeFrameCell{pKeep},toFilter,:) = repmat(pathArray(activeFrameCell{pKeep},pKeep,:),1,sum(toFilter),1);
                        
        %%% next calculate projections for missing points that are between
        
        % existing points
        newPath(gapVec,:,:) = (pathArray(gapVec,pKeep,:).*sg1(gapVec') ...
                                      + pathArray(gapVec,pDrop,:).*sg2(gapVec')) ...
                                      ./ (sg1(gapVec') + sg2(gapVec'));
        
        %%% calculate projections for backwards and forwards trajectories
        
        % forward
        mFrame = max(frameIndex);
        df1 = pathArray(lastFrame+1:mFrame,pKeep,:)-pathArray(lastFrame:mFrame-1,pKeep,:);
        df2 = pathArray(lastFrame+1:mFrame,pDrop,:)-pathArray(lastFrame:mFrame-1,pDrop,:);
        dfMean = (df1.*sg1(lastFrame+1:mFrame)'+df2.*sg2(lastFrame+1:mFrame)')./(sg1(lastFrame+1:mFrame)' + sg2(lastFrame+1:mFrame)');
        newPath(lastFrame+1:mFrame,:,:) = newPath(lastFrame,:,:)+cumsum(dfMean,1);
        
        % backward
        db1 = pathArray(1:firstFrame-1,pKeep,:)-pathArray(2:firstFrame,pKeep,:);
        db2 = pathArray(1:firstFrame-1,pDrop,:)-pathArray(2:firstFrame,pDrop,:);
        dbMean = (db1.*sg1(1:firstFrame-1)'+db2.*sg2(1:firstFrame-1)')./(sg1(1:firstFrame-1)' + sg2(1:firstFrame-1)');
        newPath(1:firstFrame-1,:,:) = newPath(firstFrame,:,:)+flipud(cumsum(flipud(dbMean),1));
        
        % update variance and path arrays
        sigmaArray(:,toFilter,:) = repmat(sqrt(1 ./ (sg1 + sg2)),1,sum(toFilter),1);
        pathArray(:,toFilter,:) = repmat(newPath,1,sum(toFilter),1);
        
        % update list of active frames
        activeFrameCell(toFilter) = {activeFrames};
      end

      % update distance array       
      for p = 1:length(RawParticles{Channel})%alteredParticles
        % calculate forward distances
        rpFrame = rightPointVec(p);
        rpDelta = pathArray(rpFrame,p,:)-pathArray(rpFrame,:,:);      
        lpSigVec = sigmaArray(rpFrame,:,:);
        forwardDistanceMat(:,p) = mean((rpDelta ./ lpSigVec).^2,3);
        forwardDistanceMat(leftPointVec-rpFrame<=0|leftPointVec-rpFrame>maxDist,p) = Inf;
        forwardDistanceMat(avoidanceCell{p},p) = Inf;

        % calculate backwards distances 
        lpFrame = leftPointVec(p);
        lpDelta = pathArray(lpFrame,p,:) - pathArray(lpFrame,:,:);      
        rpSigVec = sigmaArray(lpFrame,:,:);
        backwardDistanceMat(p,:) = mean((lpDelta ./ rpSigVec).^2,3);
        backwardDistanceMat(p,rightPointVec-lpFrame>=0|rightPointVec-lpFrame<-maxDist) = Inf;
        backwardDistanceMat(p,avoidanceCell{p}) = Inf;
      end
      
      mDistanceMat = sqrt((backwardDistanceMat + forwardDistanceMat)/2);
      mDistanceMat(isnan(mDistanceMat)) = Inf;
    end    
    close(f);
    % generate new stitched particles struct
    newParticleIndex = unique(segmentIDVec);    
    
    StitchedParticlesSub = struct;
    for p = 1:length(newParticleIndex)
      matchIndices = find(segmentIDVec==newParticleIndex(p));
      % approval 
      StitchedParticlesSub(p).Approved = false;
      % extant frames
      StitchedParticlesSub(p).Frame = activeFrameCell{matchIndices(1)};
      % position info
      StitchedParticlesSub(p).xPos = pathArray(StitchedParticlesSub(p).Frame,matchIndices(1),1);
      StitchedParticlesSub(p).yPos = pathArray(StitchedParticlesSub(p).Frame,matchIndices(1),2);
      StitchedParticlesSub(p).zPosDetrended = pathArray(StitchedParticlesSub(p).Frame,matchIndices(1),3);
      % full projected path and error
      StitchedParticlesSub(p).pathArray = reshape(pathArray(:,matchIndices(1),:),[],3);
      StitchedParticlesSub(p).sigmaArray = reshape(sigmaArray(:,matchIndices(1),:),[],3);
      % record info vectors
      StitchedParticlesSub(p).origIDs = origIDCell{matchIndices(1)};
      StitchedParticlesSub(p).linkIDs = linkIDCell{matchIndices(1)};
      StitchedParticlesSub(p).linkCosts = costCell{matchIndices(1)};
      StitchedParticlesSub(p).particleIndices = particleIndexCell{matchIndices(1)};
      StitchedParticlesSub(p).NucleusID = nucleusIDCell{matchIndices(1)};
    end    
    StitchedParticles{Channel} = StitchedParticlesSub;
  end