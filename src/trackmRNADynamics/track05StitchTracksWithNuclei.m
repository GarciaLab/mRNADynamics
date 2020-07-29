function Particles = track05StitchTracksWithNuclei(...
                          StitchedParticles, RawParticles, schnitzcells, FrameInfo, ExperimentType, retrack, displayFigures)
                        
  % This function uses proximity to nucleus centers as a way to stitch
  % together particle fragments that could not otherwise be linked, for
  % instance, if they're too far seperated in the temporal domain. If the
  % number of spots per nucleus is well defined (often we expect either 1
  % or 2 spots per nucleus) this information is used to further esnure
  % correct particle linking
  
  
  % set useful parameters
  NCh = length(RawParticles);
  ncVec = [FrameInfo.nc];
  frameIndex = 1:length(ncVec);
  
  % inialize particles structure
  Particles = cell(1,NCh);
  
  for Channel = 1:NCh
    
    
    
    %
  end
