function [pathNew, sigmaNew] = updatePaths(pathArray,sigmaArray,extantFrameArray)

refVec = 1:size(pathArray,1);
extantFramesBin = max(extantFrameArray,[],2)==1;

% calculate which frames to count for each 
firstFrame = find(~isnan(extantFrameArray(:,1)),1);
lastFrame = find(~isnan(extantFrameArray(:,1)),1,'last');
ffVec = NaN(1,size(extantFrameArray,2));
lfVec = NaN(1,size(extantFrameArray,2));
for p = 1:size(extantFrameArray,2)
  ffVec(p) = find(extantFrameArray(:,p)==1,1);
  lfVec(p) = find(extantFrameArray(:,p)==1,1,'last');
end
[ffVec, sortIndices] = sort(ffVec);
lfVec = lfVec(sortIndices);

startList = [firstFrame ffVec-1 lastFrame];
endList = [firstFrame lfVec+1 lastFrame];

% initialize stacks for path update calculations
weightArray = NaN(size(sigmaArray));

% pull data
for p = 2:length(ffVec)+1
  % get indices for frames to update
  updateFilter = ismember(refVec,endList(p-1):startList(p+1))';
  % remove frames where a different particle is active
  gapFilter = updateFilter&~extantFramesBin;
  activeFilter = extantFrameArray(:,sortIndices(p-1))==1;
  % add
  weightArray(gapFilter,p-1,:) = sigmaArray(gapFilter,sortIndices(p-1),:).^-2; 
  weightArray(activeFilter,p-1,:) = 1; 
end

% calculate final path and error info
pathNew = nansum(pathArray.*weightArray,2)./nansum(weightArray,2);
sigmaNew = nansum(weightArray,3).^-0.5;
sigmaNew(extantFramesBin,:) = 0;

