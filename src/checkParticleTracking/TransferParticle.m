function cptState = TransferParticle(cptState,NewSpotIndex,Prefix)

  CF = cptState.CurrentFrame;
  CC = cptState.CurrentChannelIndex;
  
  liveExperiment = LiveExperiment(Prefix);
  FrameInfo = cptState.FrameInfo;
  globalMotionModel = getGlobalMotionModel(LiveExperiment(Prefix)); 
  
  %First, approve the particle
  cptState.SpotFilter{CC}(CF,NewSpotIndex)=1; % set to 2 to distinguish from original particles

  % set all fields to NaN by default 
  varNames = fieldnames(cptState.Particles{CC}(end));
  addIndex = length(cptState.Particles{CC}) + 1;
  for v = 1:length(varNames)
    cptState.Particles{CC}(addIndex).(varNames{v}) = NaN;
  end

  %Add this spot as a new particle to the end of the cptState.Particles{CC} structure
  cptState.Particles{CC}(addIndex).Frame=CF;
  cptState.Particles{CC}(addIndex).Index=NewSpotIndex;
  cptState.Particles{CC}(addIndex).Approved=false;
  cptState.Particles{CC}(addIndex).FrameApproved=true;

  %Get the position of the this particle
  [x,y,z]=SpotsXYZ(cptState.Spots{CC}(CF));
  cptState.Particles{CC}(addIndex).xPos=x(NewSpotIndex);
  cptState.Particles{CC}(addIndex).yPos=y(NewSpotIndex);
  cptState.Particles{CC}(addIndex).zPos=z(NewSpotIndex);
  cptState.Particles{CC}(addIndex) = detrendZ(cptState.Particles{CC}(addIndex),FrameInfo);

  %%%%%
  %Next update projected paths
  anaphaseFrames = [liveExperiment.anaphaseFrames' size(cptState.SpotFilter{CC},1)+1];
  currentCycleStart = anaphaseFrames(find(CF>anaphaseFrames,1,'last'));
  nextCycleStart = anaphaseFrames(find(CF<=anaphaseFrames,1));
  ncFrameFilter = 1:size(cptState.SpotFilter{CC},1) >= currentCycleStart & 1:size(cptState.SpotFilter{CC},1) < nextCycleStart;

  [~, cptState.Particles{CC}(addIndex).pathArray, cptState.Particles{CC}(addIndex).sigmaArray] = ...
    simulatePathsWrapper(cptState.Particles{CC}(addIndex),globalMotionModel,ncFrameFilter);

  %%%%%
  %Generate path info 

  %Now update link info   
  cptState.Particles{CC}(addIndex).linkFrameCell = {CF}; 
  cptState.Particles{CC}(addIndex).linkCostCell = [0];
  
  % find largest link number
  newID = nanmax([cptState.Particles{CC}.idVec])+1;
  cptState.Particles{CC}(addIndex).idVec = NaN(size(cptState.Particles{CC}(addIndex-1).idVec));
  cptState.Particles{CC}(addIndex).idVec(CF) = newID;
  cptState.Particles{CC}(addIndex).linkParticleCell = {newID}; 
  cptState.Particles{CC}(addIndex).linkStateString = num2str(newID);

  %%%%%
  %Other fields
  cptState.Particles{CC}(addIndex).FirstFrame = CF;
  cptState.Particles{CC}(addIndex).LastFrame = CF;
  
  
  %Update auxiliary particle structures...
  auxNames = {'RawParticles','HMMParticles','SimParticles'};
  commonVarNames = fieldnames(cptState.RawParticles{CC});
  for n = 1:length(auxNames)
    auxIndex = length(cptState.(auxNames{n}){CC})+1;
    for v = 1:length(commonVarNames)
      if ~strcmp(commonVarNames{v},'FragmentID')
        cptState.(auxNames{n}){CC}(auxIndex).(commonVarNames{v}) = cptState.Particles{CC}(addIndex).(commonVarNames{v});
      else
        cptState.(auxNames{n}){CC}(auxIndex).FragmentID = newID;
      end
    end
  end
  % add additional fields
  cptState.HMMParticles{CC}(end).hmmModel = globalMotionModel;
  cptState.SimParticles{CC}(end).hmmModel = globalMotionModel;
  for i = 1:size(cptState.Particles{CC}(addIndex).pathArray,2)
    cptState.SimParticles{CC}(end).hmmModel(i).pathVec = cptState.Particles{CC}(addIndex).pathArray(:,i);
    cptState.SimParticles{CC}(end).hmmModel(i).sigmaVec = cptState.Particles{CC}(addIndex).sigmaArray(:,i);
  end