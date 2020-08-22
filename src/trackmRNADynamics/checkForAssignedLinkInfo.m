function [ForceSplitCell, ForceMatchCell] = checkForAssignedLinkInfo(...
          fragmentIndices,SimParticles,linkStruct,SpotFilter)

  % Check frames and indices agains link lists  
  FrameVec = [];
  IndexVec = [];
  IDVec = [];
  for i = 1:length(fragmentIndices)
    FrameVec = [FrameVec SimParticles(fragmentIndices(i)).Frame];
    IndexVec = [IndexVec SimParticles(fragmentIndices(i)).Index];
    IDVec = [IDVec repelem(fragmentIndices(i),length(SimParticles(fragmentIndices(i)).Frame))];
  end
  LinVec = sub2ind(size(SpotFilter),FrameVec,IndexVec);
  
  % Check for forced joins first

  % iterate through dupVals and generate lists of IDs to match
  ForceSplitCell = {};  
  for d = 1:length(linkStruct.forbiddenLinIndices)    
    overlaps = ismember(LinVec,linkStruct.forbiddenLinIndices{d});
    splitIDs = IDVec(overlaps);    
    splitIDs = unique(splitIDs(~isnan(splitIDs)));
    if length(splitIDs)==2
      ForceSplitCell{end+1} = splitIDs;
    end
  end

  % Check for forced splits 

  % find elements that are in the lists      
  ForceMatchCell = {};  
  for d = 1:length(linkStruct.persistentLinIndices)    
    overlaps = ismember(LinVec,linkStruct.persistentLinIndices{d});
    joinIDs = IDVec(overlaps);    
    joinIDs = unique(joinIDs(~isnan(joinIDs)));
    if length(joinIDs)==2
      ForceMatchCell{end+1} = joinIDs;
    end
  end