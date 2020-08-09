function [SpotFilter,Particles]=...
    TransferParticle(Spots,SpotFilter,Particles,CurrentFrame,IndexOutput,Prefix)

%First, approve the particle
SpotFilter(CurrentFrame,IndexOutput)=1;

%HG + AR: We deleted this because we don't assign xPos and yPos for the
%Particles that are detected in the first place. %AR 8/16/2018- reinstated
%this because it was causing bugs. 
% %Get the position of the this particle
[x,y,z]=SpotsXYZ(Spots(CurrentFrame));
x=x(IndexOutput);
y=y(IndexOutput);
z=z(IndexOutput);

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

%HG + AR: We deleted this because we don't assign xPos and yPos for the
%Particles that are detected in the first place.%AR 8/16/2018- reinstated
%this because it was causing bugs. 
Particles(addIndex).xPos=x;
Particles(addIndex).yPos=y;
Particles(addIndex).zPos=z;
% error('add z detrend function')

%Next update projected paths
liveExperiment = LiveExperiment(Prefix);
globalMotionModel = getGlobalMotionModel(liveExperiment);
anaphaseFrames = [liveExperiment.anaphaseFrames size(SpotFilte.r,1)];
currentCycleStart = find(CurrentFrame>anaphaseFrames,1,'last');
nextCycleStart = find(CurrentFrame<=anaphaseFrames,1);
ncFrameFilter = 1:size(SpotFilter,1) <= nextCycleStart & 1:size(SpotFilter,1) > currentCycleStart;

[~, Particles(addIndex).pathArray, Particles(addIndex).sigmaArray] = simulatePathsWrapper(Particles(addIndex),globalMotionModel,ncFrameFilter);

%Generate path info 

%Now update link info   
Particles(addIndex).linkFrameCell = {}; 
Particles(addIndex).linkCostCell = [0];
Particles(addIndex).ptIDVec = NaN(size(Particles(end-1).idVec));
Particles(addIndex).ptIDVec(CurrentFrame) = -1;  
Particles(addIndex).linkParticleCell = {-1}; 
Particles(addIndex).linkStateString = '-1';