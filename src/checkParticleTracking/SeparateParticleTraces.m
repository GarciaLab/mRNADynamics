function [Particles,ParticleStitchInfo]=SeparateParticleTraces(CurrentParticle,CurrentFrame,Particles,SimParticles,ParticleStitchInfo,FrameInfo)

%Separate the particle trace at the specified position
%%
%List of stitchInfovariables
fieldNamesStitch = fieldnames(ParticleStitchInfo)';

%Generate temporary particle structures
for f = 1:length(fieldNamesStitch)
  if contains(fieldNamesStitch{f},'Cell')
%     currStitchTemp.(fieldNamesStitch{f}) = {};
    newStitchTemp.(fieldNamesStitch{f}) = {};
  else
%     currStitchTemp.(fieldNamesStitch{f}) = [];
    newStitchTemp.(fieldNamesStitch{f}) = [];
  end
end


%List of Particles variables
fieldNamesParticles = fieldnames(Particles)';

%Generate temporary particle structures
for f = 1:length(fieldNamesParticles)
  if contains(fieldNamesParticles{f},'Cell')
    currParticleTemp.(fieldNamesParticles{f}) = {};
    newParticleTemp.(fieldNamesParticles{f}) = {};
  else
    currParticleTemp.(fieldNamesParticles{f}) = [];
    newParticleTemp.(fieldNamesParticles{f}) = [];
  end
end

%identify vector fields to split
splitIndicesParticles = [find(strcmp(fieldNamesParticles,'nc')) ...
    find(~contains(fieldNamesParticles,{'Cell'})&contains(fieldNamesParticles,{'Frame','Pos','Dist','Shift'}))...
    find(strcmp(fieldNamesParticles,'Index'))];
  
frameFilter = Particles(CurrentParticle).Frame < CurrentFrame;

for c = 1:length(splitIndicesParticles)  
  origVec = Particles(CurrentParticle).(fieldNamesParticles{splitIndicesParticles(c)});
  newParticleTemp.(fieldNamesParticles{splitIndicesParticles(c)}) = origVec(~frameFilter);
  currParticleTemp.(fieldNamesParticles{splitIndicesParticles(c)}) = origVec(frameFilter);
end

%Reset overall trace approval flags
currParticleTemp.Approved=0;
newParticleTemp.Approved=0;

%Set Nucleus ID Info
currParticleTemp.matchCost = Particles(CurrentParticle).matchCost;
currParticleTemp.Nucleus = Particles(CurrentParticle).Nucleus;
currParticleTemp.NucleusOrig = Particles(CurrentParticle).NucleusOrig;
currParticleTemp.stitchInfoPointer = Particles(CurrentParticle).stitchInfoPointer;

newParticleTemp.matchCost = Particles(CurrentParticle).matchCost;
newParticleTemp.Nucleus = NaN; % Not sure how to handle this yet...
newParticleTemp.NucleusOrig = Particles(CurrentParticle).NucleusOrig;
newParticleTemp.stitchInfoPointer = currParticleTemp.stitchInfoPointer+1;

%Next update particle ID info
ptIDsFull = Particles(CurrentParticle).idVec;
currParticleTemp.idVec = NaN(size(ptIDsFull));
newParticleTemp.idVec = NaN(size(ptIDsFull));
currParticleTemp.idVec(1:CurrentFrame-1) = ptIDsFull(1:CurrentFrame-1);
newParticleTemp.idVec(CurrentFrame:end) = ptIDsFull(CurrentFrame:end);

%Assign Stitch info 
currStitchTemp = ParticleStitchInfo(currParticleTemp.stitchInfoPointer);



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
ncVec = [FrameInfo.nc];

%Re-run tracking for particles affiliated with early fragment
nucleusIDVecPrev = false(1,length(SimParticles));
nucleusIDVecPrev(prevIDs) = 1;
 
[pathArrayPrev, sigmaArrayPrev, ~,~, currParticleTemp.linkStateString, ...
              currStitchTemp.linkCostVec,...
              currStitchTemp.linkAdditionCell,...
              currParticleTemp.linkCostCell, ...
              currParticleTemp.linkFrameCell, currParticleTemp.linkParticleCell] = performParticleStitching(...
              1, nucleusIDVecPrev, frameIndex, {SimParticles}, 1, ncVec, Particles(CurrentParticle).matchCost);
        
%Reshape and assign path info        
currParticleTemp.pathArray = squeeze(pathArrayPrev);
currParticleTemp.sigmaArray = squeeze(sigmaArrayPrev);
currParticleTemp.linkCostCell = currParticleTemp.linkCostCell{1};
currParticleTemp.linkCostFlags = currParticleTemp.linkCostCell > Particles(CurrentParticle).costThresh;

%Same for second fragment
nucleusIDVecPost = false(1,length(SimParticles));
nucleusIDVecPost(postIDs) = 1;


[pathArrayPost, sigmaArrayPost, ~,~, newParticleTemp.linkStateString, ...
              newStitchTemp.linkCostVec, ...
              newStitchTemp.linkAdditionCell,...
              newParticleTemp.linkCostCell, ...
              newParticleTemp.linkFrameCell, newParticleTemp.linkParticleCell] = performParticleStitching(...
              1, nucleusIDVecPost, frameIndex, {SimParticles}, 1, ncVec, Particles(CurrentParticle).matchCost);
%Reshape and assign path info        
newParticleTemp.pathArray = squeeze(pathArrayPost);
newParticleTemp.sigmaArray = squeeze(sigmaArrayPost);
newParticleTemp.linkCostCell = newParticleTemp.linkCostCell{1};
newParticleTemp.linkCostFlags = newParticleTemp.linkCostCell > Particles(CurrentParticle).costThresh;

%Now update particles structure
Particles = [Particles(1:CurrentParticle-1) currParticleTemp newParticleTemp Particles(CurrentParticle+1:end)];

%And update pointer variable for later entries
for c = CurrentParticle+2:length(Particles)
  Particles(c).stitchInfoPointer = Particles(c).stitchInfoPointer+1;
end

%Lastly, add split info...
splitFrames = [max(currParticleTemp.Frame) min(newParticleTemp.Frame)];
currStitchTemp.forbiddenLinkFrameCell(end+1) = {splitFrames};
currStitchTemp.forbiddenLinkIndexCell(end+1) = ...
  {[currParticleTemp.Index(currParticleTemp.Frame==splitFrames(1)) newParticleTemp.Index(newParticleTemp.Frame==splitFrames(2))]};

%...and update stitch structure
ParticleStitchInfo = [ParticleStitchInfo(1:newParticleTemp.stitchInfoPointer-1) currStitchTemp ...
  newStitchTemp ParticleStitchInfo(newParticleTemp.stitchInfoPointer+1:end)];