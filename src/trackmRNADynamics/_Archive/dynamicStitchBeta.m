function NewParticles = dynamicStitchBeta(FullParticles,SimParticles,...
                  ParticleStitchInfo,Prefix,matchCostMax)
  
  NewParticles = FullParticles;
  liveExperiment = LiveExperiment(Prefix);
  FrameInfo = getFrameInfo(liveExperiment);
  ncVec = [FrameInfo.nc];
  nFrames = length(ncVec);
 
  % cobine linkage costs across nuclei
  costVecOrig = ParticleStitchInfo.linkCostVec;
  % this contains the info about order in which links were added
  additionCellOrig = ParticleStitchInfo.linkAdditionCell;
  % link approval info
  approvalVec = ParticleStitchInfo.linkApprovedVec;
  % set all approved links to have negative costs 
  costVecOrig(approvalVec==1) = -Inf;
  
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
    currStateCell = {NewParticles.linkStateString};
    oldInd = find(contains(currStateCell,linkedTrace)); 
    if length(oldInd) ~= 1
      error('uh oh')
    end
    % split
    fullTrace = NewParticles(oldInd).linkStateString;
 
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
    linkFrameCell = NewParticles(oldInd).linkFrameCell;
    linkCostCell = NewParticles(oldInd).linkCostCell;
    linkParticleCell = NewParticles(oldInd).linkParticleCell;
    % extract particle vec
    particleVec = NewParticles(oldInd).idVec;
    particleVecNN = particleVec(~isnan(particleVec));
    ptIDsFull = unique(particleVecNN);
    % get indexing vector for SimParticles
    FragmentIDVec = [SimParticles.FragmentID];
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
      temp.xPos = NewParticles(oldInd).xPos(ismember(particleVecNN,ptIDs));
      temp.yPos = NewParticles(oldInd).yPos(ismember(particleVecNN,ptIDs));
      temp.zPos = NewParticles(oldInd).zPos(ismember(particleVecNN,ptIDs));
      temp.zPosDetrended = NewParticles(oldInd).zPosDetrended(ismember(particleVecNN,ptIDs));
      % add link info
      temp.linkStateString = newTraces{j};

      %%%%%%%%%%%%%%%%%%%%%%%%%
      % update predicted paths 
      %%%%%%%%%%%%%%%%%%%%%%%%%    
    
      % generate useful arrays 
      nc_ft = ismember(ncVec,ncVec(temp.Frame(1))); 
      % initialize 
      nParams = length(SimParticles(1).hmmModel);
      pathArray = NaN(nFrames,length(ptIDs),nParams);
      sigmaArray = NaN(nFrames,length(ptIDs),nParams);
      extantFrameArray = NaN(nFrames,length(ptIDs));
      extantFrameArray(nc_ft,:) = 0;
      for p = 1:length(ptIDs)
        ptFilter = ismember(FragmentIDVec,ptIDs(p));
        frames = SimParticles(ptFilter).Frame;     
        extantFrameArray(frames,p) = 1;
        pathArray(:,p,:) = cat(3,SimParticles(ptFilter).hmmModel.pathVec);
        sigmaArray(:,p,:) = cat(3,SimParticles(ptFilter).hmmModel.sigmaVec);
      end
      % call path update function
      [temp.pathArray, temp.sigmaArray] = updatePaths(pathArray,sigmaArray,extantFrameArray);        

      % add other info 
      temp.Nucleus = NewParticles(oldInd).Nucleus;
      temp.NucleusOrig = NewParticles(oldInd).NucleusOrig;
      temp.NucleusDist = NewParticles(oldInd).NucleusDist(ismember(particleVecNN,ptIDs));
      temp.FrameApproved = NewParticles(oldInd).FrameApproved(ismember(particleVecNN,ptIDs));
      temp.Index = NewParticles(oldInd).Index(ismember(particleVecNN,ptIDs));
      temp.FirstFrame = temp.Frame(1);
      temp.LastFrame = temp.Frame(end);
      temp.Approved = 0;

      % incorporate into structure
      NewParticles(length(NewParticles)+1) = temp;
    end
    NewParticles = NewParticles([1:oldInd-1 oldInd+1:end]);   
  end
  for i = 1:length(NewParticles)
    NewParticles(i).matchCost = matchCostMax;
  end
  