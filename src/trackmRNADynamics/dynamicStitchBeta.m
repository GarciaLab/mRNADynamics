function NewParticles = dynamicStitchBeta(FullParticles,SimParticles,ParticleStitchInfo,Prefix,matchCostMax,Channel)
  
  NewParticles = FullParticles;
  liveExperiment = LiveExperiment(Prefix);
  FrameInfo = getFrameInfo(liveExperiment);
  ncVec = [FrameInfo.nc];
  nFrames = length(ncVec);
 
  % cobine linkage costs across nuclei
  costVecOrig = [ParticleStitchInfo{Channel}.linkCostVec];
  % this contains the info about order in which links were added
  additionCellOrig = [ParticleStitchInfo{Channel}.linkAdditionCell];
  levelVec = NaN(1,length(additionCellOrig));
  for a = 1:length(additionCellOrig)
    levelVec(a) = max(diff(find(additionCellOrig{a}~='|')))-1;
  end
  % order additions by level of linkage (max # of consecutive | characters)
  [~, sortOrder] = sort(levelVec,'descend');
  additionCellOrig = additionCellOrig(sortOrder);
  costVecOrig = costVecOrig(sortOrder);
  % identify links that are above desired threshold
  linkFlags = costVecOrig>matchCostMax;
  linksToBreak = find(linkFlags);
  % check for cases where smaller child is costlier than parent. If this
  % happens, we must also break the parent
  refVec = 1:length(additionCellOrig);
  for i = fliplr(linksToBreak)
    linkedTrace = additionCellOrig{i};
    addFlags = contains(additionCellOrig,linkedTrace) & ~linkFlags & refVec<i;
    linkFlags(addFlags) = true;
  end
  % update ID list
  linksToBreak = find(linkFlags);
  
  for i = linksToBreak
    linkedTrace = additionCellOrig{i};
    % find max number of dividers
    maxDivFrag = max(diff(find(linkedTrace~='|')))-1;    
    % find old entry in particles and update
    currStateCell = {NewParticles{Channel}.linkStateString};
    oldInd = find(contains(currStateCell,linkedTrace)); 
    if length(oldInd) ~= 1
      error('uh oh')
    end
    % split
    fullTrace = NewParticles{Channel}(oldInd).linkStateString;
    if ~strcmp(fullTrace,linkedTrace)
      error('wtf')
    end
    maxDivFull = max(diff(find(fullTrace~='|')))-1;  
    newTraces = {fullTrace};
    for m = maxDivFull:-1:maxDivFrag
      temp = {};
      for n = 1:length(newTraces)
        temp = [temp{:} strsplit(newTraces{n},repelem('|',m))];
      end
      newTraces = temp;
    end
%     newTraces = strsplit(fullTrace,repelem('|',min(maxDivFull,maxDivFrag)));
    % extract linkage info cells
    linkFrameCell = NewParticles{Channel}(oldInd).linkFrameCell;
    linkCostCell = NewParticles{Channel}(oldInd).linkCostCell;
    linkParticleCell = NewParticles{Channel}(oldInd).linkParticleCell;
    % extract particle vec
    particleVec = NewParticles{Channel}(oldInd).idVec;
    particleVecNN = particleVec(~isnan(particleVec));
    ptIDsFull = unique(particleVecNN);
    % get atom ids for each child particle
    for j = 1:length(newTraces)

      % initialize structure
      temp = struct;

      % add id fields
      temp.idVec =  NaN(size(particleVec));    
      ptIDs = cellfun(@str2num,strsplit(newTraces{j},'|')); % get particleIDs
      
      % determine subset of link cells to bring along
      ptIDsAlt = ptIDsFull(~ismember(ptIDsFull,ptIDs));
      inVec = cellfun(@(x) any(intersect(ptIDs,x)), linkParticleCell , 'UniformOutput', true);
      exVec = cellfun(@(x) any(intersect(ptIDsAlt,x)), linkParticleCell , 'UniformOutput', true);
      keepFilter = inVec&~exVec;
      
      temp.linkFrameCell = linkFrameCell(keepFilter);
      temp.linkCostCell = linkCostCell(keepFilter);
      temp.linkParticleCell = linkParticleCell(keepFilter);
      
      temp.idVec(ismember(particleVec,ptIDs)) = particleVec(ismember(particleVec,ptIDs));
      temp.Frame = find(~isnan(temp.idVec));

      % add position info 
      temp.xPos = NewParticles{Channel}(oldInd).xPos(ismember(particleVecNN,ptIDs));
      temp.yPos = NewParticles{Channel}(oldInd).yPos(ismember(particleVecNN,ptIDs));
      temp.zPos = NewParticles{Channel}(oldInd).zPos(ismember(particleVecNN,ptIDs));
      temp.zPosDetrended = NewParticles{Channel}(oldInd).zPosDetrended(ismember(particleVecNN,ptIDs));
      % add link info
      temp.linkStateString = newTraces{j};

      %%%%%%%%%%%%%%%%%%%%%%%%%
      % update predicted paths 
      %%%%%%%%%%%%%%%%%%%%%%%%%    
    
      % generate useful arrays 
      nc_ft = ismember(ncVec,ncVec(temp.Frame(1))); 
      % initialize 
      nParams = length(SimParticles{Channel}(1).hmmModel);
      pathArray = NaN(nFrames,length(ptIDs),nParams);
      sigmaArray = NaN(nFrames,length(ptIDs),nParams);
      extantFrameArray = NaN(nFrames,length(ptIDs));
      extantFrameArray(nc_ft,:) = 0;
      for p = 1:length(ptIDs)
        frames = SimParticles{Channel}(ptIDs(p)).Frame;     
        extantFrameArray(frames,p) = 1;
        pathArray(nc_ft,p,:) = cat(1,SimParticles{Channel}(ptIDs(p)).hmmModel.pathVec)';
        sigmaArray(nc_ft,p,:) = cat(1,SimParticles{Channel}(ptIDs(p)).hmmModel.sigmaVec)';
      end
      % call path update function
      [temp.pathArray, temp.sigmaArray] = updatePaths(pathArray,sigmaArray,extantFrameArray);        

      % add other info 
      temp.NucleusID = NewParticles{Channel}(oldInd).NucleusID;
      temp.NucleusIDOrig = NewParticles{Channel}(oldInd).NucleusIDOrig;
      temp.NucleusDist = NewParticles{Channel}(oldInd).NucleusDist(ismember(particleVecNN,ptIDs));
%       temp.assignmentFlags = zeros(size(NewParticles{Channel}(oldInd).assignmentFlags));
%       temp.assignmentFlags(ismember(particleVec,ptIDs)) = NewParticles{Channel}(oldInd).assignmentFlags(ismember(particleVec,ptIDs));
      temp.Approved = 0;
      
      % incorporate into structure
      NewParticles{Channel}(length(NewParticles{Channel})+1) = temp;
    end
    NewParticles{Channel} = NewParticles{Channel}([1:oldInd-1 oldInd+1:end]);   
  end
  for i = 1:length(NewParticles{Channel})
    NewParticles{Channel}(i).matchCost = matchCostMax;
  end
  