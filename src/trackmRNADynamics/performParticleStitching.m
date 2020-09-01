function [pathArray, sigmaArray, extantFrameArray, particleIDArray, linkIDCell, ...
              linkCostVec, LinkAdditionCell, linkCostCell, linkFrameCell, linkParticleCell] = ...
              performParticleStitching(...
              NucleusID,nucleusIDVec,frameIndex,SimParticles,ncVec,matchCostMax,...
              ForceMatchCell,ForceSplitCell,FragmentIDVec)

  nParams = length(SimParticles(1).hmmModel);
  % see how many fragments match this nucleus
  nucleusFilter = nucleusIDVec==NucleusID;
  nucleusIndices = find(nucleusFilter);
  fragmentIDs = FragmentIDVec(nucleusIndices);
  nFragments = length(nucleusIndices);

  % arrays to track active frames (dynamically updated throughout)
  endpointFrameArray = false(length(frameIndex),nFragments);
  extantFrameArray = false(length(frameIndex),nFragments);
  linkCostCell = cell(1,size(extantFrameArray,2));
  linkFrameCell = cell(1,size(extantFrameArray,2));
  linkParticleCell = cell(1,size(extantFrameArray,2));
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
    endpointFrameArray([SimParticles(p).FirstFrame,SimParticles(p).LastFrame],i_pass) = true;
    extantFrameArray(SimParticles(p).Frame,i_pass) = true;      
    particleIDArray(SimParticles(p).Frame,i_pass) = fragmentIDs(i_pass);
    linkIDCell{i_pass} = num2str(fragmentIDs(i_pass));
    linkCostCell{i_pass} = [0];
    linkFrameCell{i_pass} = {unique([SimParticles(p).FirstFrame,SimParticles(p).LastFrame])};
    linkParticleCell{i_pass} = {fragmentIDs(i_pass)};
%     nc_ft = ismember(ncVec,ncVec(SimParticles(p).FirstFrame));  
    for np = 1:nParams
      pathArray(:,i_pass,np) = SimParticles(p).hmmModel(np).pathVec;
      sigmaArray(:,i_pass,np) = SimParticles(p).hmmModel(np).sigmaVec;        
    end     
    i_pass = i_pass + 1;
  end    
  
  %%% calculate Mahalanobis Distance between particles in likelihood space

  cumActivityArrayDown = cumsum(extantFrameArray);
  cumActivityArrayUp = flipud(cumsum(flipud(extantFrameArray)));
  mDistanceMatRaw = Inf(nFragments,nFragments);     
  maxFrame = max(frameIndex); 
  for m = 1:nFragments
    % pull activity indicators
    calcFrames = find(endpointFrameArray(:,m)'); 
    epFrameVec = ([1 calcFrames maxFrame]);   
      
    % find non-overlapping segments
    optionVec = max(extantFrameArray(extantFrameArray(:,m),:),[],1)==0; 
    if any(optionVec)
      % we want to exclude points that are not at an interface with another
      % fragment            
      includeMat = diff(cumActivityArrayUp(epFrameVec,optionVec))<0 | diff(cumActivityArrayDown(epFrameVec,optionVec))>0;
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
  
  % enforce user-assigned matches 
  [mDistanceMat, ForceMatchCell] = imposeLinkAssigments(mDistanceMat,ForceMatchCell,ForceSplitCell,particleIDArray);                                    

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
    afDrop = find(extantFrameArray(:,pDrop))';
    
    % condense particle and nucleus trackers
    particleIDArray(afDrop,pKeep) = particleIDArray(afDrop,pDrop);
    particleIDArray = particleIDArray(:,[1:pDrop-1 pDrop+1:end]);

    nucleusIDVec = nucleusIDVec([1:pDrop-1 pDrop+1:end]);
    
    
    %%%%%%%%%%%%%%%%%%
    % update link info array
    activeFrames = [find(endpointFrameArray(:,pKeep)') find(endpointFrameArray(:,pDrop)')];
    keepFrames = ([1 find(endpointFrameArray(:,pKeep)') maxFrame]);
    includeVecKeep = diff(cumActivityArrayUp(keepFrames,pDrop)')<0|diff(cumActivityArrayDown(keepFrames,pDrop)')>0;
    includeVecKeep = includeVecKeep(1:end-1) | includeVecKeep(2:end);
    
    dropFrames = ([1 find(endpointFrameArray(:,pDrop)') maxFrame]);
    includeVecDrop = diff(cumActivityArrayUp(dropFrames,pDrop)')<0|diff(cumActivityArrayDown(dropFrames,pDrop)')>0;
    includeVecDrop = includeVecDrop(1:end-1) | includeVecDrop(2:end);
    
    % frames
    linkFrameCell{pKeep} = [linkFrameCell{pKeep} linkFrameCell{pDrop} {activeFrames([includeVecKeep includeVecDrop])}]; 
    linkFrameCell = linkFrameCell (:,[1:pDrop-1 pDrop+1:end]);
    % cost
    linkCostCell{pKeep} = [linkCostCell{pKeep} linkCostCell{pDrop} minCost];
    linkCostCell = linkCostCell(:,[1:pDrop-1 pDrop+1:end]);
    % particles 
    ptIDs = unique(particleIDArray(~isnan(particleIDArray(:,pKeep)),pKeep))';
    linkParticleCell{pKeep} = [linkParticleCell{pKeep} linkParticleCell{pDrop} {ptIDs}];
    linkParticleCell = linkParticleCell(:,[1:pDrop-1 pDrop+1:end]);

    % condense active frame trackers
    endpointFrameArray(:,pKeep) = endpointFrameArray(:,pKeep) | endpointFrameArray(:,pDrop);
    endpointFrameArray = endpointFrameArray(:,[1:pDrop-1 pDrop+1:end]);                                  
    
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
    
    % add link cost
    linkCostVec = [linkCostVec minCost];    
    
    % add new link info
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
    
    % update cumulative activity arrays
    cumActivityArrayDown = cumsum(extantFrameArray);
    cumActivityArrayUp = flipud(cumsum(flipud(extantFrameArray)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update cost matrix 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    
    optionVec = max(extantFrameArray(extantFrameArray(:,pKeep),:),[],1)==0;%
    newDistCol = Inf(size(extantFrameArray,2),1);
    %%% paths FROM new particle first
    for p = find(optionVec)
      % pull activity indicators
      calcFrames = find(endpointFrameArray(:,p)'); 
      epFrameVec = [1 find(endpointFrameArray(:,p)') maxFrame];      
      
      % we want to exclude points that are not at an interface with another
      % fragment            
      includeMat = diff(cumActivityArrayUp(epFrameVec,pKeep))<0 | diff(cumActivityArrayDown(epFrameVec,pKeep))>0;%diff(cumActivityArray(epFrameVec,pKeep))>0;
      includeMat = repmat(includeMat(1:end-1,:) | includeMat(2:end,:),1,1,nParams);

      % calculate distances
      deltaMat = (pathArray(calcFrames,pKeep,:) - pathArray(calcFrames,p,:)).^2;        
      sigmaMat = sigmaArray(calcFrames,pKeep,:).^2;
      distanceMat = deltaMat./sigmaMat;
      distanceMat(~includeMat) = NaN;        
      newDistCol(p) = nanmean(nanmean(distanceMat,1),3); 
    end       
    
    %%% next paths TO new particle from existing particles      
    % here the endpoint frames are fixed
    calcFrames = find(endpointFrameArray(:,pKeep)');    
    epFrameVec = [1 calcFrames maxFrame];    
    newDistRow = Inf(1,size(extantFrameArray,2));
   
    if any(optionVec)
      % we want to exclude endpoints that are not at interface with other
      % fragments
      includeMat = diff(cumActivityArrayUp(epFrameVec,optionVec))<0 | diff(cumActivityArrayDown(epFrameVec,optionVec))>0;%diff(cumActivityArray(epFrameVec,optionVec))>0;
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
    
    % impose user-assigned links/splits
    try
      [mDistanceMat, ForceMatchCell] = imposeLinkAssigments(mDistanceMat,ForceMatchCell,ForceSplitCell,particleIDArray);
    catch
      error('wtf')
    end
    % calculate  current lowest cost
    minCost = min(mDistanceMat(:));
    % increment
    iter = iter + 1;        
  end  