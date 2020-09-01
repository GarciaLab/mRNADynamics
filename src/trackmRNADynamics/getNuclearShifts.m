function [SpotBulkDxVec,SpotBulkDyVec] = getNuclearShifts(schnitzcells,...
                    CurrentFrame,maxFrame,NewSpotsX,NewSpotsY,SlidingWindowSize,UseHistone)

if UseHistone            
  % use sliding window to estimate average nucleus movement
  NucleiDxVec = [];
  NucleiDyVec = []; 
  NucleiPxVec = [];
  NucleiPyVec = []; 
  StopFrame = min([maxFrame,CurrentFrame+SlidingWindowSize]);
  StartFrame = max([1,CurrentFrame-SlidingWindowSize]);
  for i = 1:length(schnitzcells)             
    FrameFT = ismember(schnitzcells(i).frames,[StartFrame,StopFrame]);
    if sum(FrameFT)==2
      NucleiDxVec = [NucleiDxVec diff(schnitzcells(i).cenx(FrameFT))];
      NucleiDyVec = [NucleiDyVec diff(schnitzcells(i).ceny(FrameFT))];
      NucleiPxVec = [NucleiPxVec mean(schnitzcells(i).cenx(FrameFT))];
      NucleiPyVec = [NucleiPyVec mean(schnitzcells(i).ceny(FrameFT))];
    end
  end
  
  if ~isempty(NucleiPxVec)
    % calculate distance to each nucleus 
    NucleusWeightMat = NaN(length(NewSpotsX),length(NucleiPxVec));
    for i = 1:length(NewSpotsX)
      NucleusWeightMat(i,:) = (rand(numel(NucleiPxVec),1)*0.05+vecnorm([NewSpotsX(i) NewSpotsY(i)]  - [NucleiPxVec' NucleiPyVec'], 2, 2)).^-2;              
    end
    % assign weighted mean bulk displacement to particles
    SpotBulkDxVec = (sum(repmat(NucleiDxVec,length(NewSpotsX),1).*NucleusWeightMat,2) ./ sum(NucleusWeightMat,2) / (StopFrame-StartFrame+1))';
    SpotBulkDyVec = (sum(repmat(NucleiDyVec,length(NewSpotsX),1).*NucleusWeightMat,2) ./ sum(NucleusWeightMat,2) / (StopFrame-StartFrame+1))';     
  else
    SpotBulkDxVec = zeros(size(NewSpotsX));
    SpotBulkDyVec = zeros(size(NewSpotsY));
  end
else
  SpotBulkDxVec = zeros(size(NewSpotsX));
  SpotBulkDyVec = zeros(size(NewSpotsY));
end