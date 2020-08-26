function SimParticles = track03PredictParticlePaths(...
                          HMMParticles, FrameInfo, displayFigures)
                        
  % This script uses inferred GHMM models to (a) infer "true" motion states
  % for observed particle frames and to predict particles positions before
  % and after first and last observed frames. These hypothetical tracks can
  % then be used to stitch together particle fragments
  
  SimParticles = HMMParticles;
  NCh = length(SimParticles);
  ncVec = [FrameInfo.nc];
  frameIndex = 1:length(ncVec);
  
  for Channel = 1:NCh
    if length(length(SimParticles{Channel}))>1
      wb = waitbar(0,['Simulating particle paths (channel ' num2str(Channel) ')']);
    end
    for p = 1:length(SimParticles{Channel})
      if length(length(SimParticles{Channel}))>1
        waitbar(p/length(SimParticles{Channel}),wb);      
      end
      
      nc = ncVec(SimParticles{Channel}(p).Frame(1)==frameIndex);
      ncFrameFilter = ncVec==nc;
      hmmModel = simulatePathsWrapper(SimParticles{Channel}(p),SimParticles{Channel}(p).hmmModel,ncFrameFilter);
      
      SimParticles{Channel}(p).hmmModel = hmmModel;
    end
    if length(length(SimParticles{Channel}))>1
      close(wb);
    end
  end
  