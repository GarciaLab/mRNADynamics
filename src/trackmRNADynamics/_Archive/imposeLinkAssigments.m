function [mDistanceMat, ForceMatchCell] = imposeLinkAssigments(mDistanceMat,ForceMatchCell,ForceSplitCell,particleIDArray)

  for m = 1:length(ForceMatchCell)
    FragmentIndices = ForceMatchCell{m};
    StitchIndices = [find(nanmax(particleIDArray==FragmentIndices(1))==1) find(nanmax(particleIDArray==FragmentIndices(2))==1)];
    % check to see if match has been made 
%     if StitchIndices(1)==StitchIndices(2)
%       ForceMatchCell = ForceMatchCell([1:m-1 m+1:end]);
    if StitchIndices(1)~=StitchIndices(2)
      mi = max(StitchIndices);
      mDistanceMat(StitchIndices(StitchIndices==mi),StitchIndices(StitchIndices~=mi)) = -Inf;
    end
  end
  % enforce user-assigned splits
  for m = 1:length(ForceSplitCell)
    FragmentIndices = ForceSplitCell{m};
    StitchIndices = [find(nanmax(particleIDArray==FragmentIndices(1))==1) find(nanmax(particleIDArray==FragmentIndices(2))==1)];
    % check to see if match has been made 
    if StitchIndices(1)==StitchIndices(2)      
      error('Issue with imposing user-assigned splits')
    else
      mi = max(StitchIndices);
      mDistanceMat(StitchIndices(StitchIndices==mi),StitchIndices(StitchIndices~=mi)) = Inf;
    end
  end