function [pathNew, sigmaNew] = updatePaths(pathArray,sigmaArray,extantFrameArray)

  refVec = 1:size(pathArray,1);
  extantFramesBin = max(extantFrameArray,[],2)==1;

  % calculate which frames to count for each 
  firstFrame = find(~isnan(pathArray(:,1,1)),1);
  lastFrame = find(~isnan(pathArray(:,1,1)),1,'last');
  ffVec = [];%NaN(1,size(extantFrameArray,2));
  lfVec = [];%NaN(1,size(extantFrameArray,2));
  idVec = [];
  for p = 1:size(extantFrameArray,2)
    starts = find([0 diff(extantFrameArray(:,p))'==1]);
    if extantFrameArray(1,p) == 1
      starts = [1 starts];
    end
    ends = find(diff(extantFrameArray(:,p))==-1)';
    if extantFrameArray(end,p) == 1
      ends = [ends size(extantFrameArray,1)];
    end
    ffVec = [ffVec starts];
    lfVec = [lfVec ends];
    idVec = [idVec repelem(p,length(starts))];
  end
  [ffVec, sortIndices] = sort(ffVec);
  lfVec = lfVec(sortIndices);
  idVec = idVec(sortIndices);
  
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
    activeFilter = updateFilter&extantFrameArray(:,idVec(p-1))==1;
    % add
    weightArray(gapFilter,idVec(p-1),:) = sigmaArray(gapFilter,idVec(p-1),:).^-2; 
    weightArray(activeFilter,idVec(p-1),:) = 1; 
  end

  % calculate final path and error info
  pathNew = nansum(pathArray.*weightArray,2)./nansum(weightArray,2);
  sigmaNew = nansum(weightArray,2).^-0.5;
  sigmaNew(extantFramesBin,:,:) = 0;
  sigmaNew([1:firstFrame-1 lastFrame+1:end],:,:) = NaN;

