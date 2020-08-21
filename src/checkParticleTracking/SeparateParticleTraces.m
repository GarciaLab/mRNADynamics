function cptState = SeparateParticleTraces(cptState)

%Separate the particle trace at the specified position
%%
Ch = cptState.CurrentChannelIndex;
%get list of particle fragments to exclude
ptList = cptState.Particles{Ch}(cptState.CurrentParticle).idVec;
ptList = unique(ptList(~isnan(ptList)));

%List of cptState.Particles{Ch} variables
fieldNameParticles = fieldnames(cptState.Particles{Ch})';

%Generate temporary particle structures
for f = 1:length(fieldNameParticles)
  if contains(fieldNameParticles{f},'Cell')
    currParticleTemp.(fieldNameParticles{f}) = {};
    newParticleTemp.(fieldNameParticles{f}) = {};
  else
    currParticleTemp.(fieldNameParticles{f}) = [];
    newParticleTemp.(fieldNameParticles{f}) = [];
  end
end

%identify vector fields to split
splitIndices = find(ismember(fieldNameParticles,cptState.frameLevelFields));
  
frameFilter = cptState.Particles{Ch}(cptState.CurrentParticle).Frame < cptState.CurrentFrame;

for c = 1:length(splitIndices)  
  origVec = cptState.Particles{Ch}(cptState.CurrentParticle).(fieldNameParticles{splitIndices(c)});
  newParticleTemp.(fieldNameParticles{splitIndices(c)}) = origVec(~frameFilter);
  currParticleTemp.(fieldNameParticles{splitIndices(c)}) = origVec(frameFilter);
end

%Reset overall trace approval flags
currParticleTemp.Approved=0;
newParticleTemp.Approved=0;

%Set Nucleus ID Info
currParticleTemp.matchCost = cptState.Particles{Ch}(cptState.CurrentParticle).matchCost;
currParticleTemp.Nucleus = cptState.Particles{Ch}(cptState.CurrentParticle).Nucleus;
currParticleTemp.NucleusOrig = cptState.Particles{Ch}(cptState.CurrentParticle).NucleusOrig;
currParticleTemp.stitchInfoPointer = cptState.Particles{Ch}(cptState.CurrentParticle).stitchInfoPointer;

newParticleTemp.matchCost = cptState.Particles{Ch}(cptState.CurrentParticle).matchCost;
newParticleTemp.Nucleus = NaN; % Not sure how to handle this yet...
newParticleTemp.NucleusOrig = cptState.Particles{Ch}(cptState.CurrentParticle).NucleusOrig;
newParticleTemp.stitchInfoPointer = currParticleTemp.stitchInfoPointer+1;

%Next update particle ID info
ptIDsFull = cptState.Particles{Ch}(cptState.CurrentParticle).idVec;
currParticleTemp.idVec = NaN(size(ptIDsFull));
newParticleTemp.idVec = NaN(size(ptIDsFull));
currParticleTemp.idVec(1:cptState.CurrentFrame-1) = ptIDsFull(1:cptState.CurrentFrame-1);
newParticleTemp.idVec(cptState.CurrentFrame:end) = ptIDsFull(cptState.CurrentFrame:end);

% check to see if part of same particle fragment appears in both new traces
overlapFilter = ismember(newParticleTemp.idVec,currParticleTemp.idVec);
overlapID = unique(newParticleTemp.idVec(overlapFilter)); 
if length(overlapID) > 1
  error('Problem with particle assignment');
end
  
prevIDs = unique(currParticleTemp.idVec(~isnan(currParticleTemp.idVec)));
postIDs = unique(newParticleTemp.idVec(~isnan(newParticleTemp.idVec)));

%Generate ID vectors
frameIndex = 1:length(ptIDsFull);
ncVec = [cptState.FrameInfo.nc];

%Re-run tracking for particles affiliated with early fragment
FragmentIDVec = [cptState.SimParticles{Ch}.FragmentID];
nucleusIDVecPrev = false(1,length(cptState.SimParticles{Ch}));
nucleusIDVecPrev(prevIDs) = 1;
 
[pathArrayPrev, sigmaArrayPrev, ~,~, currParticleTemp.linkStateString, ...
              linkCostVecOrig,...
              linkAdditionCellOrig,...
              currParticleTemp.linkCostCell, ...
              currParticleTemp.linkFrameCell, currParticleTemp.linkParticleCell] = performParticleStitching(...
              1, nucleusIDVecPrev, frameIndex, cptState.SimParticles{Ch}, ncVec, ...
              cptState.Particles{Ch}(cptState.CurrentParticle).matchCost,{},{},FragmentIDVec);
 
%Reshape and assign path info        
currParticleTemp.pathArray = squeeze(pathArrayPrev);
currParticleTemp.sigmaArray = squeeze(sigmaArrayPrev);
currParticleTemp.linkCostCell = currParticleTemp.linkCostCell{1};
currParticleTemp.linkCostFlags = currParticleTemp.linkCostCell > cptState.Particles{Ch}(cptState.CurrentParticle).costThresh;

%Same for second fragment
nucleusIDVecPost = false(1,length(cptState.SimParticles{Ch}));
nucleusIDVecPost(postIDs) = 1;


[pathArrayPost, sigmaArrayPost, ~,~, newParticleTemp.linkStateString, ...
              linkCostVecNew, ...
              linkAdditionCellNew,...
              newParticleTemp.linkCostCell, ...
              newParticleTemp.linkFrameCell, newParticleTemp.linkParticleCell] = performParticleStitching(...
              1, nucleusIDVecPost, frameIndex, cptState.SimParticles{Ch}, ncVec, ...
              cptState.Particles{Ch}(cptState.CurrentParticle).matchCost,{},{},FragmentIDVec);
            
%Reshape and assign path info        
newParticleTemp.pathArray = squeeze(pathArrayPost);
newParticleTemp.sigmaArray = squeeze(sigmaArrayPost);
newParticleTemp.linkCostCell = newParticleTemp.linkCostCell{1};
newParticleTemp.linkCostFlags = newParticleTemp.linkCostCell > cptState.Particles{Ch}(cptState.CurrentParticle).costThresh;

%Now update particles structure
cptState.Particles{Ch} = [cptState.Particles{Ch}(1:cptState.CurrentParticle-1) currParticleTemp newParticleTemp cptState.Particles{Ch}(cptState.CurrentParticle+1:end)];

%And update pointer variable for later entries
for c = cptState.CurrentParticle+2:length(cptState.Particles{Ch})
  cptState.Particles{Ch}(c).stitchInfoPointer = cptState.Particles{Ch}(c).stitchInfoPointer+1;
end

%Lastly, add split info. This is inelegant, but we need to make sure every
%pair of points before and after split are forbidden
forbiddenLinkFrameCell = {};
forbiddenLinkIndexCell = {};
for i = 1:length(currParticleTemp.Frame)
  for j = 1:length(newParticleTemp.Frame)
    splitFrames = [currParticleTemp.Frame(i) newParticleTemp.Frame(j)];
    forbiddenLinkFrameCell(end+1) = {splitFrames};
    forbiddenLinkIndexCell(end+1) = ...
      {[currParticleTemp.Index(currParticleTemp.Frame==splitFrames(1)) newParticleTemp.Index(newParticleTemp.Frame==splitFrames(2))]};
  end
end

%%%%%%%%%%%%%%%
%Update Stitch structure
cptState.ParticleStitchInfo{Ch}.forbiddenLinkFrameCell = [cptState.ParticleStitchInfo{Ch}.forbiddenLinkFrameCell forbiddenLinkFrameCell];
cptState.ParticleStitchInfo{Ch}.forbiddenLinkIndexCell = [cptState.ParticleStitchInfo{Ch}.forbiddenLinkIndexCell forbiddenLinkIndexCell];

% genernate addition ID cells 
linkAdditionIDCellOrig = {};
for p = 1:length(linkAdditionCellOrig)
  fragments = strsplit(linkAdditionCellOrig{p},'|');
  linkAdditionIDCellOrig{p} = cellfun(@str2num,fragments);
end

linkAdditionIDCellNew = {};
for p = 1:length(linkAdditionCellNew)
  fragments = strsplit(linkAdditionCellNew{p},'|');
  linkAdditionIDCellNew{p} = cellfun(@str2num,fragments);
end

%We need to identify and remove info corresponding to old particle and ad
%new inf
linksToRemove = cellfun(@(x) all(ismember(x,ptList)),cptState.ParticleStitchInfo{Ch}.linkAdditionIDCell);
cptState.ParticleStitchInfo{Ch}.linkAdditionCell = [cptState.ParticleStitchInfo{Ch}.linkAdditionCell(~linksToRemove) linkAdditionCellOrig linkAdditionCellNew];
cptState.ParticleStitchInfo{Ch}.linkAdditionIDCell = [cptState.ParticleStitchInfo{Ch}.linkAdditionIDCell(~linksToRemove) linkAdditionIDCellOrig linkAdditionIDCellNew];
cptState.ParticleStitchInfo{Ch}.linkCostVec = [cptState.ParticleStitchInfo{Ch}.linkCostVec(~linksToRemove) linkCostVecOrig linkCostVecNew];
cptState.ParticleStitchInfo{Ch}.linkApprovedVec = [cptState.ParticleStitchInfo{Ch}.linkApprovedVec(~linksToRemove) repelem(0,length(linkCostVec))];
