function [SpotFilter,Particles]=...
    TransferParticle(Spots,SpotFilter,Particles,CurrentFrame,IndexOutput,Prefix)

liveExperiment = LiveExperiment(Prefix);
FrameInfo = getFrameInfo(liveExperiment);
globalMotionModel = getGlobalMotionModel(liveExperiment);

%First, approve the particle
SpotFilter(CurrentFrame,IndexOutput)=1;

% set all fields to NaN by default 
varNames = fieldnames(Particles(end));
addIndex = length(Particles) + 1;
for v = 1:length(varNames)
  Particles(addIndex).(varNames{v}) = NaN;
end

%Add this spot as a new particle to the end of the Particles structure
Particles(addIndex).Frame=CurrentFrame;
Particles(addIndex).Index=IndexOutput;
Particles(addIndex).Approved=false;
Particles(addIndex).FrameApproved=true;

%Get the position of the this particle
[x,y,z]=SpotsXYZ(Spots(CurrentFrame));
Particles(addIndex).xPos=x(IndexOutput);
Particles(addIndex).yPos=y(IndexOutput);
Particles(addIndex).zPos=z(IndexOutput);
Particles(addIndex) = detrendZ(Particles(addIndex),FrameInfo);
% error('add z detrend function')

%%%%%
%Next update projected paths
anaphaseFrames = [liveExperiment.anaphaseFrames' size(SpotFilter,1)+1];
currentCycleStart = anaphaseFrames(find(CurrentFrame>anaphaseFrames,1,'last'));
nextCycleStart = anaphaseFrames(find(CurrentFrame<=anaphaseFrames,1));
ncFrameFilter = 1:size(SpotFilter,1) >= currentCycleStart & 1:size(SpotFilter,1) < nextCycleStart;

[~, Particles(addIndex).pathArray, Particles(addIndex).sigmaArray] = simulatePathsWrapper(Particles(addIndex),globalMotionModel,ncFrameFilter);

%%%%%
%Generate path info 

%Now update link info   
Particles(addIndex).linkFrameCell = {CurrentFrame}; 
Particles(addIndex).linkCostCell = [0];
% find largest link number
newID = nanmax([Particles.idVec])+1;
Particles(addIndex).idVec = NaN(size(Particles(addIndex-1).idVec));
Particles(addIndex).idVec(CurrentFrame) = newID;
Particles(addIndex).linkParticleCell = {newID}; 
Particles(addIndex).linkStateString = num2str(newID);

%%%%%
%Other fields
Particles(addIndex).FirstFrame = CurrentFrame;
Particles(addIndex).LastFrame = CurrentFrame;