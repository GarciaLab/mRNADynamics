function linkStruct = generateLinkStructure(ParticleStitchInfo,SpotFilter)

  fLinkIndexCell = ParticleStitchInfo.forbiddenLinkIndexCell;
  fLinkFrameCell = ParticleStitchInfo.forbiddenLinkFrameCell;
  
  linkStruct.forbiddenLinIndices = cell(1,length(ParticleStitchInfo.forbiddenLinkIndexCell));
  szVec = size(SpotFilter);
  for p = 1:length(fLinkIndexCell)
    % convert to linear indices    
    linkStruct.forbiddenLinIndices{p} = sub2ind(szVec,fLinkFrameCell{p},fLinkIndexCell{p});
  end
  
  pLinkIndexCell = ParticleStitchInfo.persistentLinkIndexCell;
  pLinkFrameCell = ParticleStitchInfo.persistentLinkFrameCell;
  
  linkStruct.persistentLinIndices = cell(1,length(ParticleStitchInfo.persistentLinkIndexCell)); 
  for p = 1:length(pLinkIndexCell)
    % convert to linear indices    
    linkStruct.persistentLinIndices{p} = sub2ind(szVec,pLinkFrameCell{p},pLinkIndexCell{p});
  end
