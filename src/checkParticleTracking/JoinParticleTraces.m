function cptState = JoinParticleTraces(cptState,ClickedParticle)

%This function joins two particle traces and renumbers all particles in the
%cptState.Particles structure accordingly
%Check to make sure nucleus IDs are consistent
%Transfer the information to the original particle
Channel = cptState.CurrentChannelIndex;
if isnan(cptState.Particles{Channel}(ClickedParticle).Nucleus) ...
    || cptState.Particles{Channel}(ClickedParticle).Nucleus==...
    cptState.Particles{Channel}(cptState.CurrentParticle).Nucleus
  
  % initialize temporary structure 
  temp = cptState.Particles{Channel}(cptState.CurrentParticle);  

  %Concatentate vector quantities with one entry per frame
  varNames = fieldnames(cptState.Particles{Channel}(cptState.CurrentParticle))';
  catIndices = find(ismember(varNames,cptState.frameLevelFields));
  
  % get frame ordering
  [temp.Frame, sortIndices] = sort([cptState.Particles{Channel}(cptState.CurrentParticle).Frame,cptState.Particles{Channel}(ClickedParticle).Frame]);

  for c = 1:length(catIndices)
    newVec = [cptState.Particles{Channel}(cptState.CurrentParticle).(varNames{catIndices(c)}),cptState.Particles{Channel}(ClickedParticle).(varNames{catIndices(c)})];
    temp.(varNames{catIndices(c)}) = newVec(sortIndices);
  end
  temp.Approved=0;
  temp.FrameApproved=logical(temp.FrameApproved);

  %Next update projected paths
  ptIDs = [cptState.CurrentParticle ClickedParticle];
  nParams = size(cptState.Particles{Channel}(cptState.CurrentParticle).pathArray,2);
  nFrames = size(cptState.Particles{Channel}(cptState.CurrentParticle).pathArray,1);
  pathArray = NaN(nFrames,2,nParams);
  sigmaArray = NaN(nFrames,2,nParams);
  extantFrameArray = zeros(nFrames,2);  
  for p = 1:length(ptIDs)
    frames = cptState.Particles{Channel}(ptIDs(p)).Frame;     
    extantFrameArray(frames,p) = 1;
    pathArray(:,p,:) = cptState.Particles{Channel}(ptIDs(p)).pathArray;
    sigmaArray(:,p,:) = cptState.Particles{Channel}(ptIDs(p)).sigmaArray;
  end
  
  % call path update function
  [temp.pathArray, temp.sigmaArray] = updatePaths(pathArray,sigmaArray,extantFrameArray);

  %Now update link info 
  % frames
  binActVec = max(extantFrameArray,[],2);
  binActVec(extantFrameArray(:,2)==1)=2;  % set new frames to 2
  
  newDist = bwdist(binActVec==2)';
  newDist(~extantFrameArray(:,1)) = Inf;
  [minDist,baseFrame] = min(newDist);
%   newFrames = find(newDist==minDist);
  
  baseDist = bwdist(binActVec==1)';
  baseDist(~extantFrameArray(:,2)) = Inf;    
  [minDist,newFrame] = min(baseDist);
%   baseFrames = find(baseDist==minDist);
    
  joinFrames = [baseFrame newFrame];
  temp.linkFrameCell = [cptState.Particles{Channel}(cptState.CurrentParticle).linkFrameCell cptState.Particles{Channel}(ClickedParticle).linkFrameCell {joinFrames}]; 

  % cost
  temp.linkCostCell = [cptState.Particles{Channel}(cptState.CurrentParticle).linkCostCell cptState.Particles{Channel}(ClickedParticle).linkCostCell -1]; % flag cost with negative val

  % particle vec
  temp.idVec = NaN(size(cptState.Particles{Channel}(cptState.CurrentParticle).idVec));
  ptIDsOrig = cptState.Particles{Channel}(cptState.CurrentParticle).idVec(~isnan(cptState.Particles{Channel}(cptState.CurrentParticle).idVec));
  ptIDsClick = cptState.Particles{Channel}(ClickedParticle).idVec(~isnan(cptState.Particles{Channel}(ClickedParticle).idVec));
  temp.idVec(~isnan(cptState.Particles{Channel}(cptState.CurrentParticle).idVec)) = ptIDsOrig;
  temp.idVec(~isnan(cptState.Particles{Channel}(ClickedParticle).idVec)) = ptIDsClick;

  % particle cell
  temp.linkParticleCell = [cptState.Particles{Channel}(cptState.CurrentParticle).linkParticleCell cptState.Particles{Channel}(ClickedParticle).linkParticleCell {[unique(ptIDsOrig) unique(ptIDsClick)]}];

  % update link string
  nKeep = max(diff(find(cptState.Particles{Channel}(cptState.CurrentParticle).linkStateString~='|')))-1;
  nDrop = max(diff(find(cptState.Particles{Channel}(ClickedParticle).linkStateString~='|')))-1;
  divNum = max([nKeep nDrop]);
  if isempty(divNum)
    divNum = 0;
  end
  temp.linkStateString = [cptState.Particles{Channel}(cptState.CurrentParticle).linkStateString repelem('|',divNum+1) cptState.Particles{Channel}(ClickedParticle).linkStateString];

  %Now, assign temp and get rid of the clicked particle
  cptState.Particles{Channel}(cptState.CurrentParticle) = temp;
  cptState.Particles{Channel} = cptState.Particles{Channel}([1:ClickedParticle-1,ClickedParticle+1:end]);
  
  %Lastly, add assigned link to stitch info struct 
  joinFilter = ismember(cptState.Particles{Channel}(cptState.CurrentParticle).Frame,joinFrames);
  cptState.ParticleStitchInfo{Channel}(cptState.Particles{Channel}(cptState.CurrentParticle).stitchInfoPointer).persistentLinkFrameCell(end+1) = {joinFrames};
  cptState.ParticleStitchInfo{Channel}(cptState.Particles{Channel}(cptState.CurrentParticle).stitchInfoPointer).persistentLinkIndexCell(end+1) = {cptState.Particles{Channel}(cptState.CurrentParticle).Index(joinFilter)};
else
  warning("Mismatching nucleus IDs. Aborting linkage")
end