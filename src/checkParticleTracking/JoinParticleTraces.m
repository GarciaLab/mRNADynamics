function [Particles, ParticleStitchInfo] = ...
  JoinParticleTraces(OriginalParticle,ClickedParticle,Particles,ParticleStitchInfo)

%This function joins two particle traces and renumbers all particles in the
%Particles structure accordingly
%Check to make sure nucleus IDs are consistent
%Transfer the information to the original particle

if isnan(Particles(ClickedParticle).Nucleus) || Particles(ClickedParticle).Nucleus==Particles(OriginalParticle).Nucleus
  % initialize temporary structure 
  temp = Particles(OriginalParticle);  

  %Concatentate vector quantities with one entry per frame
  varNames = fieldnames(Particles(OriginalParticle))';
  catIndices = [find(strcmp(varNames,'FrameApproved')) find(strcmp(varNames,'nc')) find(contains(varNames,{'Pos','Dist','Shift'})) find(strcmp(varNames,'Index'))];
  [temp.Frame, sortIndices] = sort([Particles(OriginalParticle).Frame,Particles(ClickedParticle).Frame]);

  for c = 1:length(catIndices)
    newVec = [Particles(OriginalParticle).(varNames{catIndices(c)}),Particles(ClickedParticle).(varNames{catIndices(c)})];
    temp.(varNames{catIndices(c)}) = newVec(sortIndices);
  end
  temp.Approved=0;
  temp.FrameApproved=logical(temp.FrameApproved);

  %Next update projected paths
  ptIDs = [OriginalParticle ClickedParticle];
  nParams = size(Particles(OriginalParticle).pathArray,2);
  nFrames = size(Particles(OriginalParticle).pathArray,1);
  pathArray = NaN(nFrames,2,nParams);
  sigmaArray = NaN(nFrames,2,nParams);
  extantFrameArray = zeros(nFrames,2);  
  for p = 1:length(ptIDs)
    frames = Particles(ptIDs(p)).Frame;     
    extantFrameArray(frames,p) = 1;
    pathArray(:,p,:) = Particles(ptIDs(p)).pathArray;
    sigmaArray(:,p,:) = Particles(ptIDs(p)).sigmaArray;
  end
  
  % call path update function
  [temp.pathArray, temp.sigmaArray] = updatePaths(pathArray,sigmaArray,extantFrameArray);

  %Now update link info 
  % frames
  binActVec = max(extantFrameArray,[],2);
  binActVec(extantFrameArray(:,2)==1)=2;  % set new frames to 2
  
  newDist = bwdist(binActVec==2)';
  newDist(~extantFrameArray(:,1)) = Inf;
  [minDist,newFrame] = min(newDist);
%   newFrames = find(newDist==minDist);
  
  baseDist = bwdist(binActVec==1)';
  baseDist(~extantFrameArray(:,2)) = Inf;    
  [minDist,baseFrame] = min(baseDist);
%   baseFrames = find(baseDist==minDist);
    
  joinFrames = [baseFrame newFrame];
  temp.linkFrameCell = [Particles(OriginalParticle).linkFrameCell Particles(ClickedParticle).linkFrameCell {joinFrames}]; 

  % cost
  temp.linkCostCell = [Particles(OriginalParticle).linkCostCell Particles(ClickedParticle).linkCostCell -1]; % flag cost with negative val

  % particle vec
  temp.idVec = NaN(size(Particles(OriginalParticle).idVec));
  ptIDsOrig = Particles(OriginalParticle).idVec(~isnan(Particles(OriginalParticle).idVec));
  ptIDsClick = Particles(ClickedParticle).idVec(~isnan(Particles(ClickedParticle).idVec));
  temp.idVec(~isnan(Particles(OriginalParticle).idVec)) = ptIDsOrig;
  temp.idVec(~isnan(Particles(ClickedParticle).idVec)) = ptIDsClick;

  % particle cell
  temp.linkParticleCell = [Particles(OriginalParticle).linkParticleCell Particles(ClickedParticle).linkParticleCell {[unique(ptIDsOrig) unique(ptIDsClick)]}];

  % update link string
  nKeep = max(diff(find(Particles(OriginalParticle).linkStateString~='|')))-1;
  nDrop = max(diff(find(Particles(ClickedParticle).linkStateString~='|')))-1;
  divNum = max([nKeep nDrop]);
  if isempty(divNum)
    divNum = 0;
  end
  temp.linkStateString = [Particles(OriginalParticle).linkStateString repelem('|',divNum+1) Particles(ClickedParticle).linkStateString];

  %Now, assign temp and get rid of the clicked particle
  Particles(OriginalParticle) = temp;
  Particles = Particles([1:ClickedParticle-1,ClickedParticle+1:end]);
  
  %Lastly, add assigned link to stitch info struct 
  joinFilter = ismember(Particles(OriginalParticle).Frame,joinFrames);
  ParticleStitchInfo(Particles(OriginalParticle).Nucleus).persistentLinkFrames(end+1) = {joinFrames};
  ParticleStitchInfo(Particles(OriginalParticle).Nucleus).persistentLinkIndices(end+1) = {Particles(OriginalParticle).Index(joinFilter)};
else
  warning("Mismatching nucleus IDs. Aborting linkage")
end