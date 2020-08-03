function [pathArray, sigmaArray, extantFrameArray, particleIDArray, linkIDCell, ...
              linkCostVec, LinkAdditionCell, mDistanceMat] = performParticleStitching(...
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
  pathArray = NaN(length(frameIndex),nFragments,nParams);   
  sigmaArray = NaN(length(frameIndex),nFragments,nParams); 

  % particle ID tracker
  particleIDArray = NaN(length(frameIndex),nFragments);
  linkIDCell = cell(1,nFragments);
  LinkAdditionCell = {};
  linkCostVec = [];
  
  % add fragment info
  i_pass = 1;
  for p = nucleusIndices      
    endpointFrameArray([RawParticles{Channel}(p).FirstFrame,RawParticles{Channel}(p).LastFrame],i_pass) = true;
    extantFrameArray(RawParticles{Channel}(p).FirstFrame:RawParticles{Channel}(p).LastFrame,i_pass) = true;      
    particleIDArray(RawParticles{Channel}(p).FirstFrame:RawParticles{Channel}(p).LastFrame,i_pass) = p;
    linkIDCell{i_pass} = num2str(p);
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
    
    % condense particle and nucleus trackers
    particleIDArray(afDrop,pKeep) = particleIDArray(afDrop,pDrop);
    particleIDArray = particleIDArray(:,[1:pDrop-1 pDrop+1:end]);

    nucleusIDVec = nucleusIDVec([1:pDrop-1 pDrop+1:end]);                          
    
    % condense link ID tracker
    index_vec = 1:length(linkIDCell);
    nKeep = max(diff(find(linkIDCell{pKeep}~='|')))-1;
    nDrop = max(diff(find(linkIDCell{pDrop}~='|')))-1;
    divNum = max([nKeep nDrop]);
    if isempty(divNum)
      divNum = 0;
    end
    linkIDCell{pKeep} = [ linkIDCell{pKeep} repelem('|',divNum+1) linkIDCell{pDrop} ];
    linkIDCell = linkIDCell(~ismember(index_vec,pDrop));
    
    % same for link cost cell (use nesting logic here as well
    linkCostVec = [linkCostVec minCost];    
    
    % add new link 
    LinkAdditionCell{length(linkCostVec)} = linkIDCell{pKeep};
    
    % condense distance array
    mDistanceMat = mDistanceMat(~ismember(index_vec,pDrop),~ismember(index_vec,pDrop));   
        
    % call path update function
    [pathArray(:,pKeep,:), sigmaArray(:,pKeep,:)] = updatePaths(...
      pathArray(:,[pKeep pDrop],:),sigmaArray(:,[pKeep pDrop],:),extantFrameArray(:,[pKeep pDrop]));        
            
    % condense path arrays
    sigmaArray = sigmaArray(:,[1:pDrop-1 pDrop+1:end],:);
    pathArray = pathArray(:,[1:pDrop-1 pDrop+1:end],:);
    
    % update active frame info
    extantFrameArray(:,pKeep) = extantFrameArray(:,pKeep) | extantFrameArray(:,pDrop);
    extantFrameArray = extantFrameArray(:,[1:pDrop-1 pDrop+1:end]);
    
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