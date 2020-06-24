function StitchedParticles = track04StitchTrackNoNuclei(...
                          RawParticles, FrameInfo, retrack, displayFigures)
                        
                        

  NCh = length(RawParticles);
  ncVec = [FrameInfo.nc];
  frameIndex = 1:length(ncVec);
  maxDist = 10; % max time gap between linkable points (in frames)
  gateValue = 3; % maximum number of sigmas away
  matchCost = 0.5*.5; % cost for initial pair matching (pairs must be within 1/2 sigma to be matched)
  
  for Channel = 1:NCh          
    tic
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
    sigmaArray(sigmaArray==0) = sqrt(realmin) * 1e8; 
    % for each particle generate list of overlapping (incompatible)
    % particles
    overlapCell = cell(1,length(RawParticles{Channel}));
    activeFrameCell = cell(1,length(RawParticles{Channel}));
    origIDCell = cell(1,length(RawParticles{Channel}));
    linkIDCell = cell(1,length(RawParticles{Channel}));
    costCell = cell(1,length(RawParticles{Channel}));
    for p = 1:length(RawParticles{Channel})
      overlapCell{p} = find(max(extantArray(extantArray(:,p),:),[],1)==1);
      activeFrameCell{p} = find(extantArray(:,p))';
      origIDCell{p} = repelem(p,length(activeFrameCell{p}));
      linkIDCell{p} = repelem(1,length(activeFrameCell{p}));
      costCell{p} = repelem(0,length(activeFrameCell{p}));
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
    toc
    tic
    % perform initial pair matching using optimal linear assignment
    mDistanceMat = sqrt((backwardDistanceMat + forwardDistanceMat)/2);
    mDistanceMat(isnan(mDistanceMat)) = Inf;
    [matchArray,~,~] = matchpairs(mDistanceMat,matchCost);
    
    % update path arrays and info cells
    for m = 1:size(matchArray,1)
      % extract match indices
      pKeep = matchArray(m,2);
      pDrop = matchArray(m,1);
      minimum = mDistanceMat(pDrop,pKeep);
      % update ID vec
      toID = segmentIDVec(pKeep);
      fromID = segmentIDVec(pDrop);
      segmentIDVec(segmentIDVec==fromID) = toID;
      toFilter = segmentIDVec==toID;
      
      % update info tracking cells
      overlapCell(toFilter) = {unique([overlapCell{[pKeep pDrop]}])};
      
      [activeFrames,si] = sort([activeFrameCell{[pKeep pDrop]}]);
      activeFrameCell(toFilter) = {activeFrames};
      
      origIDs = [origIDCell{[pKeep pDrop]}];
      origIDCell(toFilter) = {origIDs(si)};
      
      currentLinks = [linkIDCell{pKeep} repelem(max(linkIDCell{pKeep})+1,length(activeFrameCell{pDrop}))];
      linkIDCell(toFilter) = {currentLinks(si)};
      
      currentCosts = [costCell{pKeep} repelem(minimum,length(activeFrameCell{pDrop}))];
      costCell(toFilter) = {currentCosts(si)};
      % generate updated combined path using variance-weighted average
      sg1 = sigmaArray(:,pKeep,:).^-2;
      sg2 = sigmaArray(:,pDrop,:).^-2;
      pathArray(:,toFilter,:) = repmat((pathArray(:,pKeep,:).*sg1 ...
                                    + pathArray(:,pDrop,:).*sg2) ...
                                    ./ (sg1 + sg2),1,sum(toFilter),1);
        
      % do the same for the variance
      sigmaArray(:,toFilter,:) = repmat(sqrt(1 ./ (sg1 + sg2)),1,sum(toFilter),1);
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
      backwardDistanceMat(overlapCell{p},p) = Inf;
    end   
    
    % use greedy algorithm to link additional trace fragments
    mDistanceMat = sqrt((backwardDistanceMat + forwardDistanceMat)/2);
    minimum = min(mDistanceMat(:));
    iter = 1;
    while minimum < gateValue
      % get coordinates of minimum
      [pDrop,pKeep]=find(mDistanceMat==minimum,1);
      
      % update ID vec
      toID = segmentIDVec(pKeep);
      fromID = segmentIDVec(pDrop);
      segmentIDVec(segmentIDVec==fromID) = toID;
      toFilter = segmentIDVec==toID;
      
      % update info tracking cells
      overlapCell(toFilter) = {unique([overlapCell{[pKeep pDrop]}])};
      
      [activeFrames,si] = sort([activeFrameCell{[pKeep pDrop]}]);
      activeFrameCell(toFilter) = {activeFrames};
      
      origIDs = [origIDCell{[pKeep pDrop]}];
      origIDCell(toFilter) = {origIDs(si)};
      
      currentLinks = [linkIDCell{pKeep} repelem(max(linkIDCell{pKeep})+1,length(activeFrameCell{pDrop}))];
      linkIDCell(toFilter) = {currentLinks(si)};
      
      currentCosts = [costCell{pKeep} repelem(minimum,length(activeFrameCell{pDrop}))];
      costCell(toFilter) = {currentCosts(si)};
      % generate updated combined path using variance-weighted average
      sg1 = sigmaArray(:,pKeep,:).^-2;
      sg2 = sigmaArray(:,pDrop,:).^-2;
      pathArray(:,toFilter,:) = repmat((pathArray(:,pKeep,:).*sg1 ...
                                    + pathArray(:,pDrop,:).*sg2) ...
                                    ./ (sg1 + sg2),1,sum(toFilter),1);
        
      % do the same for the variance
      sigmaArray(:,toFilter,:) = repmat(sqrt(1 ./ (sg1 + sg2)),1,sum(toFilter),1);             
      
      % update distance calculations
      for p = find(toFilter)'
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
      % calculate distances
      mDistanceMat = sqrt((backwardDistanceMat + forwardDistanceMat)/2);
      minimum = min(mDistanceMat(:));
      iter = iter + 1;
    end
    
    toc
  end