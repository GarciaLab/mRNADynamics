function [pathArray, sigmaArray, extantFrameArray, particleIDArray, linkIDCell, ...
              linkCostCell, linkLevelCell, mDistanceMat] = performParticleStitching(...
              NucleusID,nucleusIDVec,frameIndex,RawParticles,Channel,ncVec,matchCostMax)
  
  nParams = length(RawParticles{Channel}(1).hmmModel);
  % see how many fragments match this nucleus
  nucleusFilter = nucleusIDVec==NucleusID;
  nucleusIndices = find(nucleusFilter);
  nFragments = length(nucleusIndices);

  % arrays to track active frames (dynamically updated throughout)
  endpointFrameArray = false(length(frameIndex),nFragments);
  extantFrameArray = false(length(frameIndex),nFragments);

  % arrays to track projected particle positions
  pathArray = Inf(length(frameIndex),nFragments,nParams);   
  sigmaArray = Inf(length(frameIndex),nFragments,nParams); 

  % particle ID tracker
  particleIDArray = NaN(length(frameIndex),nFragments);
  linkIDCell = cell(1,nFragments);
  linkCostCell = cell(1,nFragments);
  linkLevelCell = cell(1,nFragments);
  
  % add fragment info
  i_pass = 1;
  for p = nucleusIndices      
    endpointFrameArray([RawParticles{Channel}(p).FirstFrame,RawParticles{Channel}(p).LastFrame],i_pass) = true;
    extantFrameArray(RawParticles{Channel}(p).FirstFrame:RawParticles{Channel}(p).LastFrame,i_pass) = true;      
    particleIDArray(RawParticles{Channel}(p).FirstFrame:RawParticles{Channel}(p).LastFrame,i_pass) = p;
    linkIDCell{i_pass} = {p};
    linkCostCell{i_pass} = [0];
    linkLevelCell{i_pass} = [0];
    nc_ft = ismember(ncVec,ncVec(RawParticles{Channel}(p).FirstFrame));
    for np = 1:nParams
      pathArray(nc_ft,i_pass,np) = RawParticles{Channel}(p).hmmModel(np).pathVec;
      sigmaArray(nc_ft,i_pass,np) = RawParticles{Channel}(p).hmmModel(np).sigmaVec;        
    end     
    i_pass = i_pass + 1;
  end    

  %%% calculate Mahalanobis Distance between particles in likelihood space
  cumActivityArray = cumsum(extantFrameArray);
  mDistanceMatRaw = Inf(nFragments,nFragments);     
  maxFrame = max(frameIndex); 
  for m = 1:nFragments
    % pull activity indicators
    calcFrames = find(endpointFrameArray(:,m)'); 
    epFrameVec = ([1 calcFrames maxFrame]);   
      
    % flag overlaps
    optionVec = max(extantFrameArray(extantFrameArray(:,m),:),[],1)==0; 
    if any(optionVec)
      % we want to exclude points that are not at an interface with another
      % fragment            
      includeMat = diff(cumActivityArray(epFrameVec,optionVec))>0;
      includeMat = repmat(includeMat(1:end-1,:) | includeMat(2:end,:),1,1,nParams);          
      
      % calculate distances      
      deltaMat = (pathArray(calcFrames,m,:) - pathArray(calcFrames,optionVec,:)).^2;        
      sigmaMat = sigmaArray(calcFrames,optionVec,:).^2;
      distanceMat = deltaMat./sigmaMat;
      distanceMat(~includeMat) = NaN;        
      mDistanceMatRaw(m,optionVec) = nanmean(nanmean(distanceMat,1),3);
    end
  end

  % calculate cost matrix
  mDistanceMat = tril(sqrt((mDistanceMatRaw+mDistanceMatRaw')/2));
  mDistanceMat(mDistanceMat==0) = Inf;    
  %       maxStitchNum = size(mDistanceMat,1)-1; % maximum number of stitches possible    
  % calculate  current lowest cost
  minCost = min(mDistanceMat(:));  
  % iterate through cost layers
  %       f = waitbar(0,['Stitching particle tracks (channel ' num2str(Channel) ')']);
  iter = 0;
  while minCost < matchCostMax      

    % get coordinates of minimum
    [pDrop,pKeep]=find(mDistanceMat==minCost,1);                    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update tracking info arrays
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    % pull active frames
    afKeep = find(extantFrameArray(:,pKeep))';
    afDrop = find(extantFrameArray(:,pDrop))';

    % condense active frame trackers
    endpointFrameArray(:,pKeep) = endpointFrameArray(:,pKeep) | endpointFrameArray(:,pDrop);
    endpointFrameArray = endpointFrameArray(:,[1:pDrop-1 pDrop+1:end]);

    extantFrameArray(:,pKeep) = extantFrameArray(:,pKeep) | extantFrameArray(:,pDrop);
    extantFrameArray = extantFrameArray(:,[1:pDrop-1 pDrop+1:end]);

    % condense particle and nucleus trackers
    particleIDArray(afDrop,pKeep) = particleIDArray(afDrop,pDrop);
    particleIDArray = particleIDArray(:,[1:pDrop-1 pDrop+1:end]);

    nucleusIDVec = nucleusIDVec([1:pDrop-1 pDrop+1:end]);
    
    % condense link ID tracker
    newLinkEntry = {[linkIDCell{pKeep}{end} linkIDCell{pDrop}{end}]};   
    linkIDCell{pKeep} = [linkIDCell{pKeep} linkIDCell{pDrop} newLinkEntry];
    linkIDCell = linkIDCell(1,[1:pDrop-1 pDrop+1:end]);
    
    % condense link cost trackers
    linkCostCell{pKeep} = [linkCostCell{pKeep} linkCostCell{pDrop} minCost];
    linkCostCell = linkCostCell(1,[1:pDrop-1 pDrop+1:end]);
    
    
    % condense node hierarchy tracker
    newLevel = max([linkLevelCell{pKeep} linkLevelCell{pDrop}])+1;
    linkLevelCell{pKeep} = [linkLevelCell{pKeep} linkLevelCell{pDrop} newLevel];
    linkLevelCell = linkLevelCell(1,[1:pDrop-1 pDrop+1:end]);
    
    % condense distance array
    mDistanceMat = mDistanceMat([1:pDrop-1 pDrop+1:end],[1:pDrop-1 pDrop+1:end]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate updated combined path using variance-weighted average 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % identify gaps where different fragments meet
    [activeFrames, si] = sort([afKeep afDrop]);
    cbIDVec = [ones(size(afKeep)) 2*ones(size(afDrop))];
    cbIDVec = cbIDVec(si);         
    
    activeFrames = [0 activeFrames maxFrame+1];
    cbIDVec = [cbIDVec(1) cbIDVec cbIDVec(end)];
    
    dfGapFlags = find(diff(cbIDVec)~=0);    
    smGapFlags = find(diff(cbIDVec)==0&diff(activeFrames)>1);

    % calculate variance weights
    nK = 1;%length(afKeep);
    nD = 1;%length(afDrop);
    sg1 = nK*sigmaArray(:,pKeep,:).^-2;
    sg2 = nD*sigmaArray(:,pDrop,:).^-2;

    % initialize arrays
    newPath = NaN(size(sigmaArray(:,1,:)));
    newSigma = NaN(size(sigmaArray(:,1,:)));
    
    % assign detected points        
    newPath(afDrop,:,:) = pathArray(afDrop,pDrop,:);
    newPath(afKeep,:,:) = pathArray(afKeep,pKeep,:);
    newSigma(afDrop,:,:) = sigmaArray(afDrop,pDrop,:);
    newSigma(afKeep,:,:) = sigmaArray(afKeep,pKeep,:);
    
    %%% next calculate updated projections for gaps between interfaces of
    %%% the two joined fragments (use variance-weighted mean) 
    cFrames = [];
    for f = 1:length(dfGapFlags)
      cFrames = [cFrames activeFrames(dfGapFlags(f))+1:activeFrames(dfGapFlags(f)+1)-1]; 
    end
    % update path info
    newPath(cFrames,:,:) = (pathArray(cFrames,pKeep,:).*sg1(cFrames,:,:) ...
                                  + pathArray(cFrames,pDrop,:).*sg2(cFrames,:,:)) ...
                                  ./ (sg1(cFrames,:,:) + sg2(cFrames,:,:));
    % update error info
    newSigma(cFrames,:,:) = sqrt(1 ./ (sg1(cFrames,:,:) + sg2(cFrames,:,:)));
    % determine which indices to use from each parent
    kFrames = [];
    dFrames = [];
    for f = 1:length(smGapFlags)
      % extract relevant frames
      frameIndices = activeFrames(smGapFlags(f))+1:activeFrames(smGapFlags(f)+1)-1;      
      if cbIDVec(smGapFlags(f)) == 1
        kFrames = [kFrames frameIndices];
      else
        dFrames = [dFrames frameIndices];
      end      
    end
    
    %%% update 
    newPath(kFrames,:,:) = pathArray(kFrames,pKeep,:);
    newPath(dFrames,:,:) = pathArray(dFrames,pDrop,:);
    newSigma(kFrames,:,:) = sigmaArray(kFrames,pKeep,:);
    newSigma(dFrames,:,:) = sigmaArray(dFrames,pDrop,:);   
    
    tVec = [afKeep afDrop kFrames dFrames cFrames];
    if length(unique(tVec))~=length(frameIndex) || length(unique(tVec))~=length(tVec)
      error('inconsistent update indices')
    end
    % update variance and path arrays
    sigmaArray(:,pKeep,:) = newSigma;
    pathArray(:,pKeep,:) = newPath;

    % condense path arrays
    sigmaArray = sigmaArray(:,[1:pDrop-1 pDrop+1:end],:);
    pathArray = pathArray(:,[1:pDrop-1 pDrop+1:end],:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update cost matrix 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    validIndices = max(extantFrameArray(extantFrameArray(:,pKeep),:),[],1)==0;% & nucleusIDVec(pKeep)==nucleusIDVec;
    newDistCol = Inf(size(extantFrameArray,2),1);
    %%% paths FROM new particle first
    for p = find(validIndices)
      % pull activity indicators
      calcFrames = find(endpointFrameArray(:,p)'); 
      epFrameVec = [1 find(endpointFrameArray(:,p)') maxFrame];      
      
      % we want to exclude points that are not at an interface with another
      % fragment            
      includeMat = diff(cumActivityArray(epFrameVec,pKeep))>0;
      includeMat = repmat(includeMat(1:end-1,:) | includeMat(2:end,:),1,1,nParams);

      % calculate distances
      deltaMat = (pathArray(calcFrames,pKeep,:) - pathArray(calcFrames,p,:)).^2;        
      sigmaMat = sigmaArray(calcFrames,pKeep,:).^2;
      distanceMat = deltaMat./sigmaMat;
      distanceMat(~includeMat) = NaN;        
      newDistCol(p) = nanmean(nanmean(distanceMat,1),3); 
    end       
    
    %%% next paths TO new particle from existing particles      
    % pull activity indicators
    calcFrames = find(endpointFrameArray(:,pKeep)');    
    epFrameVec = [1 calcFrames maxFrame];    
    newDistRow = Inf(1,size(extantFrameArray,2));
    % flag overlaps
    optionVec = max(extantFrameArray(extantFrameArray(:,pKeep),:),[],1)==0;    
    if any(optionVec)
      % we want to exclude exterior-most points from intersection of
      % segments. perform calculations to determine which points to
      % exclude in each case
      includeMat = diff(cumActivityArray(epFrameVec,optionVec))>0;
      includeMat = repmat(includeMat(1:end-1,:) | includeMat(2:end,:),1,1,nParams);
    
      % calculate distances
      deltaMat = (pathArray(calcFrames,pKeep,:) - pathArray(calcFrames,optionVec,:)).^2;
      sigmaMat = sigmaArray(calcFrames,optionVec,:).^2;
      distanceMat = deltaMat./sigmaMat;
      distanceMat(~includeMat) = NaN;        
      newDistRow(optionVec) = nanmean(nanmean(distanceMat,1),3);
    end
    % incorporate into existing distance array
    newDistVec = sqrt((newDistRow+newDistCol')/2);      
    mDistanceMat(pKeep,:) = newDistVec;
    mDistanceMat(:,pKeep) = newDistVec;
    mDistanceMat = tril(mDistanceMat);
    mDistanceMat(mDistanceMat==0) = Inf;
    % calculate  current lowest cost
    minCost = min(mDistanceMat(:));
    % increment
    iter = iter + 1;        
  end  