function StitchedParticles = track04StitchTracks(...
                          RawParticles, schnitzcells, FrameInfo, retrack, displayFigures)

  % set useful parameters
  NCh = length(RawParticles);
  ncVec = [FrameInfo.nc];
  frameIndex = 1:length(ncVec);
  maxDist = sum(ncVec==mode(ncVec))/2; % max time gap between linkable points (in frames)
  matchCostMax = 3; % maximum number of sigmas away (the 0.5 factor is an adjustment for 
  matchCostMin = 0.5; % cost for initial pair matching (pairs must be within 1/2 sigma to be matched)
  granularity = 50; % number of distinct match costs to iterate through
  costVec = linspace(matchCostMin,matchCostMax,granularity);
  
  % Check to see if we have nucleus tracking info
  useNuclei = ~isempty(schnitzcells);
  
  % initialize data structure
  StitchedParticles = cell(1,NCh);
  for Channel = 1:NCh              
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
    for p = 1:length(RawParticles{Channel})      
      nc_ft = ismember(ncVec,ncVec(RawParticles{Channel}(p).FirstFrame));
      for n = 1:nDims
        pathArray(nc_ft,p,n) = RawParticles{Channel}(p).hmmModel(n).pathVec;
        sigmaArray(nc_ft,p,n) = RawParticles{Channel}(p).hmmModel(n).sigmaVec;        
      end
    end
    extantArray = sigmaArray(:,:,1)==0;
    % set sigmas to be nonzero for weighted averaging step
%     sigmaArray(sigmaArray==0) = sqrt(realmin) * 1e8; 
    % for each particle generate list of overlapping (incompatible)
    % particles
    overlapCell = cell(1,length(RawParticles{Channel}));
    activeFrameCell = cell(1,length(RawParticles{Channel}));
    origIDCell = cell(1,length(RawParticles{Channel}));
    linkIDCell = cell(1,length(RawParticles{Channel}));
    costCell = cell(1,length(RawParticles{Channel}));
    particleIndexCell = cell(1,length(RawParticles{Channel}));
    for p = 1:length(RawParticles{Channel})
      overlapCell{p} = find(max(extantArray(extantArray(:,p),:),[],1)==1);
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
      forwardDistanceMat(overlapCell{p},p) = Inf;
      % calculate backwards distances 
      lpFrame = leftPointVec(p);
      lpDelta = pathArray(lpFrame,p,:) - pathArray(lpFrame,:,:);      
      rpSigVec = sigmaArray(lpFrame,:,:);
      backwardDistanceMat(p,:) = mean((lpDelta ./ rpSigVec).^2,3);
      backwardDistanceMat(rightPointVec-lpFrame>=0|rightPointVec-lpFrame<-maxDist,p) = Inf;
      backwardDistanceMat(overlapCell{p},p) = Inf;
    end      
    for c = 1:length(costVec)
      % perform pair matching using optimal linear assignment
      mDistanceMat = sqrt((backwardDistanceMat + forwardDistanceMat)/2);
      mDistanceMat(isnan(mDistanceMat)) = Inf;
      [matchArray,~,~] = matchpairs(mDistanceMat,0.5*costVec(c)); % note that the 0.5 factor is intentional here

      % update path arrays and info cells      
      for m = 1:size(matchArray,1)
        % extract match indices
        pKeep = matchArray(m,2);
        pDrop = matchArray(m,1);
        cost = mDistanceMat(pDrop,pKeep);
        
        % update ID vec
        toID = segmentIDVec(pKeep);
        fromID = segmentIDVec(pDrop);
        segmentIDVec(segmentIDVec==fromID) = toID;
        toFilter = segmentIDVec==toID;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % update info tracking cells
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        overlapCell(toFilter) = {unique([overlapCell{[pKeep pDrop]}])}; % tracks incompatible segments
       
        [activeFrames,si] = unique([activeFrameCell{[pKeep pDrop]}]); % tracks active frames             
        
        ptIndices = [particleIndexCell{[pKeep pDrop]}]; % tracks index of particles within Spots
        particleIndexCell(toFilter) = {ptIndices(si)};
        
        origIDs = [origIDCell{[pKeep pDrop]}]; % tracks original particle identities (before stitching)
        origIDCell(toFilter) = {origIDs(si)};

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
      alteredParticles = reshape(matchArray,1,[]);
      for p = alteredParticles
        % calculate forward distances
        rpFrame = rightPointVec(p);
        rpDelta = pathArray(rpFrame,p,:)-pathArray(rpFrame,:,:);      
        lpSigVec = sigmaArray(rpFrame,:,:);
        forwardDistanceMat(:,p) = mean((rpDelta ./ lpSigVec).^2,3);
        forwardDistanceMat(leftPointVec-rpFrame<=0|leftPointVec-rpFrame>maxDist,p) = Inf;
        forwardDistanceMat(overlapCell{p},p) = Inf;

        % calculate backwards distances 
        lpFrame = leftPointVec(p);
        lpDelta = pathArray(lpFrame,p,:) - pathArray(lpFrame,:,:);      
        rpSigVec = sigmaArray(lpFrame,:,:);
        backwardDistanceMat(p,:) = mean((lpDelta ./ rpSigVec).^2,3);
        backwardDistanceMat(rightPointVec-lpFrame>=0|rightPointVec-lpFrame<-maxDist,p) = Inf;
        backwardDistanceMat(p,overlapCell{p}) = Inf;
      end   
    end    
    
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
    end    
    StitchedParticles{Channel} = StitchedParticlesSub;
  end