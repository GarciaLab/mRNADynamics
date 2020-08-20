function linkStruct = generateLinkStructure(ParticleStitchInfo,SpotFilter)

  persistentLinkIndexVec = [];
  persistentLinkFrameVec = [];
  linkStruct.persistentLinkIDVec = [];
  linkStruct.persistentLinkSubIDVec = [];
  forbiddenLinkFrameVec = [];
  forbiddenLinkIndexVec = [];    
  linkStruct.forbiddenLinkIDVec = [];
  linkStruct.forbiddenLinkSubIDVec = [];

  f_iter = 0;
  p_iter = 0;
  for p = 1:length(ParticleStitchInfo)
    pLinks = ParticleStitchInfo(p).persistentLinkIndexCell;
    pSubIDVec = [];
    for i = 1:length(pLinks)
      pSubIDVec = [pSubIDVec repelem(i,length(pLinks{i}))];
    end
    pIDVec = repelem(p,length(pSubIDVec));
    if ~isempty(pLinks)
      persistentLinkIndexVec = [persistentLinkIndexVec pLinks{:}];
      persistentLinkFrameVec = [persistentLinkFrameVec ParticleStitchInfo(p).persistentLinkFrameCell{:}];
      linkStruct.persistentLinkIDVec = [linkStruct.persistentLinkIDVec pIDVec];
      linkStruct.persistentLinkSubIDVec = [linkStruct.persistentLinkSubIDVec p_iter+pSubIDVec];
      p_iter = max(linkStruct.persistentLinkSubIDVec);
    end

    fLinks = ParticleStitchInfo(p).forbiddenLinkIndexCell;
    fSubIDVec = [];
    for i = 1:length(fLinks)
      fSubIDVec = [fSubIDVec repelem(i,length(fLinks{i}))];
    end
    fIDVec = repelem(p,length(fSubIDVec));
    if ~isempty(fLinks)
      forbiddenLinkIndexVec = [forbiddenLinkIndexVec fLinks{:}];
      forbiddenLinkFrameVec = [forbiddenLinkFrameVec ParticleStitchInfo(p).forbiddenLinkFrameCell{:}];      
      linkStruct.forbiddenLinkIDVec = [linkStruct.forbiddenLinkIDVec fIDVec];
      linkStruct.forbiddenLinkSubIDVec = [linkStruct.forbiddenLinkSubIDVec f_iter+fSubIDVec];
      f_iter = max(linkStruct.forbiddenLinkSubIDVec);
    end
  end
  
  % convert to linear indices
  linkStruct.persistentLinIndices = sub2ind(size(SpotFilter),persistentLinkFrameVec,persistentLinkIndexVec);
  linkStruct.forbiddenLinIndices = sub2ind(size(SpotFilter),forbiddenLinkFrameVec,forbiddenLinkIndexVec);