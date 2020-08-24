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
  origParticleTemp = cptState.Particles{Ch}(cptState.CurrentParticle);  
  addParticleTemp = cptState.Particles{Ch}(ClickedParticle);  

  % generate list of frames to join
  joinLinkFrameCell = {};
  joinLinkIndexCell = {};
  joinLinkLinCell = {};
  joinParticles = {};
  for i = 1:length(addParticleTemp.Frame)
    for j = 1:length(origParticleTemp.Frame)
      joinFrames = [addParticleTemp.Frame(i) origParticleTemp.Frame(j)];
      joinLinkFrameCell(end+1) = {joinFrames};
      joinIndices = [addParticleTemp.Index(addParticleTemp.Frame==joinFrames(1)) origParticleTemp.Index(origParticleTemp.Frame==joinFrames(2))];
      joinLinkIndexCell(end+1) = {joinIndices};
      joinLinkLinCell(end+1) = {sub2ind(size(cptState.SpotFilter{Ch}),joinFrames,joinIndices)};
      joinParticles(end+1) = {[addParticleTemp.idVec(joinFrames(1)) origParticleTemp.idVec(joinFrames(2))]};
    end
  end
  
  % get list of all fragment IDs involved
  allIDs = [cptState.Particles{Ch}(cptState.CurrentParticle).idVec cptState.Particles{Ch}(ClickedParticle).idVec];
  allIDs = unique(allIDs(~isnan(allIDs)));
  
  % generate other info needed to re-run stitching
  FragmentIDVec = [cptState.SimParticles{Ch}.FragmentID];
  FragmentIDVec = [FragmentIDVec nanmax(FragmentIDVec)+1];
  nucleusIDVec = false(1,length(cptState.SimParticles{Ch}));
  nucleusIDVec(allIDs) = 1; % using "nucleus" as a dummy grouper variable here
  frameIndex = 1:length(cptState.Particles{Ch}(cptState.CurrentParticle).idVec);
  ncVec = [cptState.FrameInfo.nc];
  
  % Re-run stitching
  [pathArray, sigmaArray, ~, particleIDArray, origParticleTemp.linkStateString, ...
              linkCostVec, linkAdditionCell, origParticleTemp.linkCostCell, ...
              origParticleTemp.linkFrameCell, origParticleTemp.linkParticleCell] = performParticleStitching(...
              1, nucleusIDVec, frameIndex, cptState.SimParticles{Ch}, ncVec, ...
              cptState.Particles{Ch}(cptState.CurrentParticle).matchCost,joinParticles,{},FragmentIDVec);
            
  % incorporate path and link-related info
  origParticleTemp.pathArray = squeeze(pathArray);
  origParticleTemp.sigmaArray = squeeze(sigmaArray);
  origParticleTemp.linkCostCell = origParticleTemp.linkCostCell{1};
  origParticleTemp.linkCostFlags = origParticleTemp.linkCostCell > cptState.Particles{Ch}(cptState.CurrentParticle).costThresh;
  origParticleTemp.idVec = particleIDArray(:,1)';
  
  %Concatentate vector quantities with one entry per frame
  varNames = fieldnames(cptState.Particles{Ch}(cptState.CurrentParticle))';
  catIndices = find(ismember(varNames,cptState.frameLevelFields));
  
  % get frame ordering
  [origParticleTemp.Frame, sortIndices] = sort([cptState.Particles{Ch}(cptState.CurrentParticle).Frame,cptState.Particles{Ch}(ClickedParticle).Frame]);

  for c = 1:length(catIndices)
    newVec = [cptState.Particles{Ch}(cptState.CurrentParticle).(varNames{catIndices(c)}),cptState.Particles{Ch}(ClickedParticle).(varNames{catIndices(c)})];
    origParticleTemp.(varNames{catIndices(c)}) = newVec(sortIndices);
  end
  origParticleTemp.Approved=0;
  origParticleTemp.FrameApproved=logical(origParticleTemp.FrameApproved);
  origParticleTemp.FirstFrame = origParticleTemp.Frame(1);
  origParticleTemp.LastFrame = origParticleTemp.Frame(end);


  %Now, assign temp and get rid of the clicked particle
  cptState.Particles{Ch}(cptState.CurrentParticle) = origParticleTemp;
  cptState.Particles{Ch} = cptState.Particles{Ch}([1:ClickedParticle-1,ClickedParticle+1:end]);  
  
  % generate numeric join tracker structure
  linkAdditionIDCell = {};
  for p = 1:length(linkAdditionCell)
    fragments = strsplit(linkAdditionCell{p},'|');
    linkAdditionIDCell{p} = cellfun(@str2num,fragments);
  end
  
  % remove old links and replace with new ones
  linksToRemove = cellfun(@(x) all(ismember(x,allIDs)),cptState.ParticleStitchInfo{Ch}.linkAdditionIDCell);
  cptState.ParticleStitchInfo{Ch}.linkAdditionCell = [cptState.ParticleStitchInfo{Ch}.linkAdditionCell(~linksToRemove) linkAdditionCell];
  cptState.ParticleStitchInfo{Ch}.linkAdditionIDCell = [cptState.ParticleStitchInfo{Ch}.linkAdditionIDCell(~linksToRemove) linkAdditionIDCell];
  cptState.ParticleStitchInfo{Ch}.linkCostVec = [cptState.ParticleStitchInfo{Ch}.linkCostVec(~linksToRemove) linkCostVec];
  cptState.ParticleStitchInfo{Ch}.linkApprovedVec = [cptState.ParticleStitchInfo{Ch}.linkApprovedVec(~linksToRemove) repelem(0,length(linkCostVec))];
  
  % check to see if prior user-assigned link info exists for this pair of
  % spots
  
  % get exisiting index pairs
  linkStruct = generateLinkStructure(cptState.ParticleStitchInfo{Ch},cptState.SpotFilter{Ch});
  % look for overlaps
  oldPLinks = false(size(linkStruct.persistentLinIndices));
  oldFLinks = false(size(linkStruct.forbiddenLinIndices));
  for f = 1:length(joinLinkLinCell)
    oldPLinks = cellfun(@(x) all(ismember(x,joinLinkLinCell{f})),linkStruct.persistentLinIndices) | oldPLinks;
    oldFLinks = cellfun(@(x) all(ismember(x,joinLinkLinCell{f})),linkStruct.forbiddenLinIndices) | oldFLinks;
  end

  % override previous spearation if it exists
  cptState.ParticleStitchInfo{Ch}.forbiddenLinkFrameCell = cptState.ParticleStitchInfo{Ch}.forbiddenLinkFrameCell(~oldFLinks);
  cptState.ParticleStitchInfo{Ch}.forbiddenLinkIndexCell = cptState.ParticleStitchInfo{Ch}.forbiddenLinkIndexCell(~oldFLinks);
  
  % add persistent info if not already present      
  cptState.ParticleStitchInfo{Ch}.persistentLinkFrameCell = [cptState.ParticleStitchInfo{Ch}.persistentLinkFrameCell(~oldPLinks) joinLinkFrameCell];
  cptState.ParticleStitchInfo{Ch}.persistentLinkIndexCell = [cptState.ParticleStitchInfo{Ch}.persistentLinkIndexCell(~oldPLinks) joinLinkIndexCell];  
else
  warning("Mismatching nucleus IDs. Aborting linkage")
end