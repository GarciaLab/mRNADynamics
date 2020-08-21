function SimParticles = track03PredictParticlePaths(...
                          RawParticles, FrameInfo, displayFigures)
                        
  % This script uses inferred GHMM models to (a) infer "true" motion states
  % for observed particle frames and to predict particles positions before
  % and after first and last observed frames. These hypothetical tracks can
  % then be used to stitch together particle fragments
  
  SimParticles = RawParticles;
  NCh = length(SimParticles);
  ncVec = [FrameInfo.nc];
  frameIndex = 1:length(ncVec);
  
  for Channel = 1:NCh
    wb = waitbar(0,['Simulating particle paths (channel ' num2str(Channel) ')']);
    for p = 1:length(SimParticles{Channel})
      waitbar(p/length(SimParticles{Channel}),wb);      
      
      nc = ncVec(SimParticles{Channel}(p).Frame(1)==frameIndex);
      ncFrameFilter = ncVec==nc;
      hmmModel = simulatePathsWrapper(SimParticles{Channel}(p),SimParticles{Channel}(p).hmmModel,ncFrameFilter);
      
      SimParticles{Channel}(p).hmmModel = hmmModel;
    end
    close(wb);
  end
  