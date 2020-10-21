function [NewSpotNuclei, NewSpotDistances] = getNuclearAssigments(NewSpotsX,NewSpotsY,...
              schnitzcells,CurrentFrame,UseHistone)
if UseHistone     
  
  ExtantNucleiX = [];
  ExtantNucleiY = [];
  ncIDVec = [];
  for i = 1:length(schnitzcells)             
    CurrFT = ismember(schnitzcells(i).frames,CurrentFrame);
    if sum(CurrFT) == 1
      ExtantNucleiX = [ExtantNucleiX schnitzcells(i).cenx(CurrFT)];
      ExtantNucleiY = [ExtantNucleiY schnitzcells(i).ceny(CurrFT)];
      ncIDVec = [ncIDVec i];
    end
  end
  NucleusDistMat = NaN(length(NewSpotsX),length(ExtantNucleiX));
  for i = 1:length(NewSpotsX)            
    NucleusDistMat(i,:) = vecnorm([NewSpotsX(i) NewSpotsY(i)]  - [ExtantNucleiX' ExtantNucleiY'],2,2);
  end
  
  % assign spots to nearest neighbors
  [NewSpotDistances, minIndices] = min(NucleusDistMat,[],2); % note that I'm not enforcing unique assignment at this stage. Will do this later on in the process
  NewSpotNuclei = ncIDVec(minIndices); 
  
else
  
  NewSpotNuclei = ones(size(NewSpotsX)); % if no nucleus info, then we set all ID values to dummy val
  NewSpotDistances = zeros(size(NewSpotsX));
end