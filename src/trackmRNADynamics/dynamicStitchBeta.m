function newParticles = dynamicStitchBeta(Particles,SimParticles,FrameInfo,newCost)

% newCost = 2;
ncVec = [FrameInfo.nc];
newParticles = Particles;

Channel = 1;

% get final state of system (always start here and work backwards)
fullStateCell = {newParticles{Channel}.linkIDs};
% cobine linkage costs across nuclei
[costVecFull, si] = sort([newParticles{Channel}.linkCosts],'descend');
% this contains the info about order in which links were added
idCellFull = [newParticles{Channel}.linkAdditionCell];
idCellFull = idCellFull(si);
% identify links that are above desired threshold
linksToBreak = find(costVecFull>newCost);
% get list of fieldnames
fnames = fieldnames(Particles{Channel});
% iterate through these links and break the corresponding particles
newStateCell = fullStateCell;

tic
for i = linksToBreak
  origTrace = idCellFull{i};
  % find max number of dividers
  maxDiv = max(diff(find(diff(origTrace~='|')~=0)));
  % split
  newTraces = strsplit(origTrace,repelem('|',maxDiv));
  
  % find old entry in particles and update
  currStateCell = {newParticles{Channel}.linkIDs};
  oldInd = strcmp(fullStateCell,origTrace);
  if ~any(oldInd)
    oldInd = contains(fullStateCell,origTrace);
  end
  
  % extract particle vec
  particleVec = newParticles{Channel}(oldInd).idVec;
  particleVecNN = particleVec(~isnan(particleVec));
  
  % get atom ids for each child particle
  for j = 1:length(newTraces)
    
    % initialize structure
    temp = struct;
    for f = 1:length(fnames)
      temp.(fnames{f}) = [];
    end
    
    % add frame fields
    temp.idVec =  NaN(size(particleVec));    
    ptIDs = cellfun(@str2num,strsplit(newTraces{j},'|'));  
    temp.idVec(ismember(particleVec,ptIDs)) = particleVec(ismember(particleVec,ptIDs));
    temp.Frame = find(~isnan(temp.idVec));
    
    % add position info 
    temp.xPos = newParticles{Channel}(oldInd).xPos(ismember(particleVecNN,ptIDs));
    temp.yPos = newParticles{Channel}(oldInd).yPos(ismember(particleVecNN,ptIDs));
    
    % add link info
    temp.linkIDs = newTraces{j};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % update predicted paths 
    %%%%%%%%%%%%%%%%%%%%%%%%%    
    origPaths = newParticles{Channel}(oldInd).pathArray;   
    % generate useful arrays 
    nc_ft = ismember(ncVec,ncVec(temp.Frame(1))); 
    % initialize 
    pathArray = NaN(size(origPaths,1),length(ptIDs),size(origPaths,2));
    sigmaArray = NaN(size(origPaths,1),length(ptIDs),size(origPaths,2));
    extantFrameArray = NaN(size(origPaths,1),length(ptIDs));
    extantFrameArray(nc_ft,:) = 0;
    for p = 1:length(ptIDs)
      frames = SimParticles{Channel}(ptIDs(p)).Frame;     
      extantFrameArray(frames,p) = 1;
      pathArray(nc_ft,p,:) = cat(1,SimParticles{Channel}(ptIDs(p)).hmmModel.pathVec)';
      sigmaArray(nc_ft,p,:) = cat(1,SimParticles{Channel}(ptIDs(p)).hmmModel.sigmaVec)';
    end
    % call path update function
    [temp.pathNew, temp.sigmaNew] = updatePaths(pathArray,sigmaArray,extantFrameArray);        
    
  end
  newParticles{Channel} = newParticles{Channel}([1:oldInd-1 oldInd+1:end]);
end
toc