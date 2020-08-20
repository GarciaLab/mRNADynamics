function [ForceSplitCell, ForceMatchCell] = checkForAssignedLinkInfo(ncIndices,SimParticles,linkStruct,SpotFilter)

  % Check frames and indices agains link lists
  
  FrameVec = [];
  IndexVec = [];
  IDVec = [];
  for i = 1:length(ncIndices)
    FrameVec = [FrameVec SimParticles(ncIndices(i)).Frame];
    IndexVec = [IndexVec SimParticles(ncIndices(i)).Index];
    IDVec = [IDVec repelem(ncIndices(i),length(SimParticles(ncIndices(i)).Frame))];
  end
  LinVec = sub2ind(size(SpotFilter),FrameVec,IndexVec);
  % Check for forced joins first

  % find elements that are in the lists      
  ft1 = ismember(linkStruct.persistentLinIndices,LinVec);
  % find eleements with two or more members
  subIDs = linkStruct.persistentLinkSubIDVec(ft1);
  [~,ia] = unique(subIDs);
  dupVals = subIDs(~ismember(1:length(subIDs),ia));
  
  % iterate through dupVals and generate lists of IDs to match
  ForceMatchCell = cell(1,length(dupVals));
  
  for d = 1:length(dupVals)
    ft2 = ft1&linkStruct.persistentLinkSubIDVec==dupVals;
    ft3 = ismember(LinVec,linkStruct.persistentLinIndices(ft2));
    joinIds = IDVec(ft3);    
    joinIds = unique(joinIds(~isnan(joinIDs)));
    if length(joinIds)<2
      error('sigh')
    end
    ForceMatchCell{d} = joinIds;
  end


  % Check for forced splits 

  % find elements that are in the lists      
  ft1 = ismember(linkStruct.forbiddenLinIndices,LinVec);
  % find elements with two or more members
  subIDs = linkStruct.forbiddenLinkSubIDVec(ft1);
  [~,ia] = unique(subIDs);
  dupVals = subIDs(~ismember(1:length(subIDs),ia));
  
  % iterate through dupVals and generate lists of IDs to match
  ForceSplitCell = cell(1,length(dupVals));
  
  for d = 1:length(dupVals)
    ft2 = ft1&linkStruct.forbiddenLinkSubIDVec==dupVals;
    ft3 = ismember(LinVec,linkStruct.forbiddenLinIndices(ft2));
    splitIds = IDVec(ft3);    
    splitIds = unique(splitIds(~isnan(splitIds)));
    if length(splitIds)<2
      error('sigh')
    end
    ForceSplitCell{d} = splitIds;
  end