function [pathArray, sigmaArray, extantFrameArray, particleIDArray, linkIDArray, linkCostArray, mDistanceMat] = performParticleStitching(...
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
  linkIDArray = NaN(length(frameIndex),nFragments);
  linkCostArray = NaN(length(frameIndex),nFragments);

  % add fragment info
  i_pass = 1;
  for p = nucleusIndices      
    endpointFrameArray([RawParticles{Channel}(p).FirstFrame,RawParticles{Channel}(p).LastFrame],i_pass) = true;
    extantFrameArray(RawParticles{Channel}(p).FirstFrame:RawParticles{Channel}(p).LastFrame,i_pass) = true;      
    particleIDArray(RawParticles{Channel}(p).FirstFrame:RawParticles{Channel}(p).LastFrame,i_pass) = p;
    linkIDArray(RawParticles{Channel}(p).FirstFrame:RawParticles{Channel}(p).LastFrame,i_pass) = 0;
    linkCostArray(RawParticles{Channel}(p).FirstFrame:RawParticles{Channel}(p).LastFrame,i_pass) = 0;
    nc_ft = ismember(ncVec,ncVec(RawParticles{Channel}(p).FirstFrame));
    for np = 1:nParams
      pathArray(nc_ft,i_pass,np) = RawParticles{Channel}(p).hmmModel(np).pathVec;
      sigmaArray(nc_ft,i_pass,np) = RawParticles{Channel}(p).hmmModel(np).sigmaVec;        
    end     
    i_pass = i_pass + 1;
  end    

  %%% calculate Mahalanobis Distance between particles in likelihood space
  mDistanceMatRaw = Inf(nFragments,nFragments);      
  for m = 1:nFragments
    % pull activity indicators
    epFrameVec1 = find(endpointFrameArray(:,m));    
    f1 = epFrameVec1(1);
    f2 = epFrameVec1(end);
    % flag overlaps
    optionVec = max(extantFrameArray(extantFrameArray(:,m),:),[],1)==0; 
    if any(optionVec)
      % we want to exclude exterior-most points from intersection of
      % segments. perform calculations to determine which points to
      % exclude in each case
      excludeMat = false(length(epFrameVec1),size(extantFrameArray,2));
      if f1<f2        
        excludeMat(1,:) = max(extantFrameArray(1:f1,:),[],1)~=1;
        excludeMat(end,:) = max(extantFrameArray(f2:end,:),[],1)~=1;
      end
      % calculate distances
      deltaMat = (pathArray(epFrameVec1,m,:) - pathArray(epFrameVec1,optionVec,:)).^2;        
      sigmaMat = sigmaArray(epFrameVec1,optionVec,:).^2;
      distanceMat = deltaMat./sigmaMat;
      distanceMat(excludeMat(:,optionVec)) = NaN;        
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
    linkIDArray(afDrop,pKeep) = nanmax(linkIDArray(:,pKeep))+1;
    linkIDArray = linkIDArray(:,[1:pDrop-1 pDrop+1:end]);
    
    % condense link cost trackers
    linkCostArray(afDrop,pKeep) = minCost;
    linkCostArray = linkCostArray(:,[1:pDrop-1 pDrop+1:end]);
    
    % condense distance array
    mDistanceMat = mDistanceMat([1:pDrop-1 pDrop+1:end],[1:pDrop-1 pDrop+1:end]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate updated combined path using variance-weighted average 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    activeFrames = sort([afKeep afDrop]);
    firstFrame = activeFrames(1);
    lastFrame = activeFrames(end);
    midVec = firstFrame:lastFrame;
    gapVec = midVec(~ismember(midVec,activeFrames));

    % calculate variance weights
    sg1 = sigmaArray(:,pKeep,:).^-2;
    sg2 = sigmaArray(:,pDrop,:).^-2;

    % first assign known points
    newPath = NaN(size(sigmaArray(:,1,:)));
    newPath(afDrop,:,:) = pathArray(afDrop,pDrop,:);
    newPath(afKeep,:,:) = pathArray(afKeep,pKeep,:);

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
    sigmaArray(:,pKeep,:) = sqrt(1 ./ (sg1 + sg2));
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
      epFrameVec1 = find(endpointFrameArray(:,p));    
      f1 = epFrameVec1(1);
      f2 = epFrameVec1(end);

      % we want to exclude exterior-most points from intersection of
      % segments. perform calculations to determine which points to
      % exclude in each case
      excludeVec = false(length(epFrameVec1),1);
      if f1<f2        
        excludeVec(1,:) = max(extantFrameArray(1:f1,pKeep),[],1)~=1;
        excludeVec(end,:) = max(extantFrameArray(f2:end,pKeep),[],1)~=1;
      end
      % calculate distances
      deltaMat = (pathArray(epFrameVec1,pKeep,:) - pathArray(epFrameVec1,p,:)).^2;
      sigmaMat = sigmaArray(epFrameVec1,pKeep,:).^2;
      distanceMat = deltaMat./sigmaMat;
      distanceMat(excludeVec) = NaN;        
      newDistCol(p) = nanmean(nanmean(distanceMat,1),3);                               
    end

    %%% next paths TO new particle from existing particles      
    % pull activity indicators
    epFrameVec1 = find(endpointFrameArray(:,pKeep));    
    f1 = epFrameVec1(1);
    f2 = epFrameVec1(end);
    newDistRow = Inf(1,size(extantFrameArray,2));
    % flag overlaps
    optionVec = max(extantFrameArray(extantFrameArray(:,pKeep),:),[],1)==0;    
    if any(optionVec)
      % we want to exclude exterior-most points from intersection of
      % segments. perform calculations to determine which points to
      % exclude in each case
      excludeMat = false(length(epFrameVec1),size(extantFrameArray,2));
      if f1<f2        
        excludeMat(1,:) = max(extantFrameArray(1:f1,:),[],1)~=1;
        excludeMat(end,:) = max(extantFrameArray(f2:end,:),[],1)~=1;
      end
      % calculate distances
      deltaMat = (pathArray(epFrameVec1,pKeep,:) - pathArray(epFrameVec1,optionVec,:)).^2;
      sigmaMat = sigmaArray(epFrameVec1,optionVec,:).^2;
      distanceMat = deltaMat./sigmaMat;
      distanceMat(excludeMat(:,optionVec)) = NaN;        
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