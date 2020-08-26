function particleStructure = updateVectorFields(particleStructure,frameFilter)
  vecFields = {'Frame','Index','xPos','yPos','zPos','zPosDetrended','NucleusDist'};
  
  for v = 1:length(vecFields)
    vec = particleStructure.(vecFields{v});
    particleStructure.(vecFields{v}) = vec(frameFilter);
  end
  particleStructure.FirstFrame = particleStructure.Frame(1);
  particleStructure.LastFrame = particleStructure.Frame(end);