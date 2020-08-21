function cptState = JoinParticleTraces(cptState,ClickedParticle)

%This function joins two particle traces and renumbers all particles in the
%cptState.Particles structure accordingly
%Check to make sure nucleus IDs are consistent
%Transfer the information to the original particle
Ch = cptState.CurrentChannelIndex;
if isnan(cptState.Particles{Ch}(ClickedParticle).Nucleus) ...
    || cptState.Particles{Ch}(ClickedParticle).Nucleus==...
    cptState.Particles{Ch}(cptState.CurrentParticle).Nucleus
  
  % initialize temporary structure 
  newParticleTemp = cptState.Particles{Ch}(cptState.CurrentParticle);  

  % find closest frames
  distMat = abs(cptState.Particles{Ch}(cptState.CurrentParticle).Frame' - ...
            cptState.Particles{Ch}(ClickedParticle).Frame);
  [~,minOrig] = min(distMat,[],1);
  [~,minNew] = min(distMat,[],2);  
  
  % get frame values at join
  joinFrames = [cptState.Particles{Ch}(cptState.CurrentParticle).Frame(minOrig)...
    cptState.Particles{Ch}(ClickedParticle).Frame(minNew(1))];
  % get fragment ids at join
  joinParticles = [cptState.Particles{Ch}(cptState.CurrentParticle).idVec(joinFrames(1))...
    cptState.Particles{Ch}(ClickedParticle).idVec(joinFrames(2))];
  
  % get list of all fragment IDs involved
  allIDs = [cptState.Particles{Ch}(cptState.CurrentParticle).idVec cptState.Particles{Ch}(ClickedParticle).idVec];
  allIDs = unique(allIDs(~isnan(allIDs)));
  
  % generate other info needed to re-run stitching
  FragmentIDVec = [cptState.SimParticles{Ch}.FragmentID];
  FragmentIDVec = [FragmentIDVec nanmax(FragmentIDVec)+1];
  nucleusIDVec = false(1,length(cptState.SimParticles{Ch}));
  nucleusIDVec(allIDs) = 1;
  frameIndex = 1:length(cptState.Particles{Ch}(cptState.CurrentParticle).idVec);
  ncVec = [cptState.FrameInfo.nc];
  
  % Re-run stitching
  [pathArray, sigmaArray, ~, particleIDArray, newParticleTemp.linkStateString, ...
              linkCostVec, linkAdditionCell, newParticleTemp.linkCostCell, ...
              newParticleTemp.linkFrameCell, newParticleTemp.linkParticleCell] = performParticleStitching(...
              1, nucleusIDVec, frameIndex, cptState.SimParticles{Ch}, ncVec, ...
              cptState.Particles{Ch}(cptState.CurrentParticle).matchCost,{joinParticles},{},FragmentIDVec);
            
  % incorporate path and link-related info
  newParticleTemp.pathArray = squeeze(pathArray);
  newParticleTemp.sigmaArray = squeeze(sigmaArray);
  newParticleTemp.linkCostCell = newParticleTemp.linkCostCell{1};
  newParticleTemp.linkCostFlags = newParticleTemp.linkCostCell > cptState.Particles{Ch}(cptState.CurrentParticle).costThresh;
  newParticleTemp.idVec = particleIDArray(:,p)';
  
  %Concatentate vector quantities with one entry per frame
  varNames = fieldnames(cptState.Particles{Ch}(cptState.CurrentParticle))';
  catIndices = find(ismember(varNames,cptState.frameLevelFields));
  
  % get frame ordering
  [newParticleTemp.Frame, sortIndices] = sort([cptState.Particles{Ch}(cptState.CurrentParticle).Frame,cptState.Particles{Ch}(ClickedParticle).Frame]);

  for c = 1:length(catIndices)
    newVec = [cptState.Particles{Ch}(cptState.CurrentParticle).(varNames{catIndices(c)}),cptState.Particles{Ch}(ClickedParticle).(varNames{catIndices(c)})];
    newParticleTemp.(varNames{catIndices(c)}) = newVec(sortIndices);
  end
  newParticleTemp.Approved=0;
  newParticleTemp.FrameApproved=logical(newParticleTemp.FrameApproved);
  newParticleTemp.FirstFrame = newParticleTemp.Frame(1);
  newParticleTemp.LastFrame = newParticleTemp.Frame(end);


  %Now, assign temp and get rid of the clicked particle
  cptState.Particles{Ch}(cptState.CurrentParticle) = newParticleTemp;
  cptState.Particles{Ch} = cptState.Particles{Ch}([1:ClickedParticle-1,ClickedParticle+1:end]);
  
  
  %Lastly, remove outdated stitch info and replace with new 
  
  % identify elements in particle array with multiple fragments
  linkAdditionIDCell = {};
  for p = 1:length(linkAdditionCell)
    fragments = strsplit(linkAdditionCell{p},'|');
    linkAdditionIDCell{p} = cellfun(@str2num,fragments);
  end
  
  linksToRemove = cellfun(@(x) all(ismember(x,allIDs)),cptState.ParticleStitchInfo{Ch}.linkAdditionIDCell);
  cptState.ParticleStitchInfo{Ch}.linkAdditionCell = [cptState.ParticleStitchInfo{Ch}.linkAdditionCell(~linksToRemove) linkAdditionCell];
  cptState.ParticleStitchInfo{Ch}.linkAdditionIDCell = [cptState.ParticleStitchInfo{Ch}.linkAdditionIDCell(~linksToRemove) linkAdditionIDCell];
  cptState.ParticleStitchInfo{Ch}.linkCostVec = [cptState.ParticleStitchInfo{Ch}.linkCostVec(~linksToRemove) linkCostVec];
  cptState.ParticleStitchInfo{Ch}.linkApprovedVec = [cptState.ParticleStitchInfo{Ch}.linkApprovedVec(~linksToRemove) repelem(0,length(linkCostVec))];
  
  % add persistent frames
  joinFilter = ismember(cptState.Particles{Ch}(cptState.CurrentParticle).Frame,joinFrames);
  cptState.ParticleStitchInfo{Ch}.persistentLinkFrameCell(end+1) = {joinFrames};
  cptState.ParticleStitchInfo{Ch}.persistentLinkIndexCell(end+1) = {cptState.Particles{Ch}(cptState.CurrentParticle).Index(joinFilter)};
else
  warning("Mismatching nucleus IDs. Aborting linkage")
end