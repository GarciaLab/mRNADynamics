function cptState = separateTraces(cptState)
  
  
%SEPARATETRACES Summary of this function goes here
%   Detailed explanation goes here

%The separated particle (the trace following current frame) won't have a nucleus assigned!
cptState.PreviousParticle = 0;

%Check that the particle does actually exist in this frame
if ~(cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle).Frame(1)==cptState.CurrentFrame) && ...
  any(cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle).Frame==cptState.CurrentFrame)

  % call splitting function
  [cptState.Particles{cptState.CurrentChannelIndex},cptState.ParticleStitchInfo{cptState.CurrentChannelIndex}]=...
    SeparateParticleTraces(cptState.CurrentParticle,cptState.CurrentFrame,...
    cptState.Particles{cptState.CurrentChannelIndex},cptState.SimParticles{cptState.CurrentChannelIndex},...
    cptState.ParticleStitchInfo{cptState.CurrentChannelIndex},cptState.FrameInfo); 
  
  % save
      
elseif length(cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle).Frame)==1
  
  cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle).Nucleus=[];
  
else
  
  disp('Cannot divide a trace at the first time point')
end
end

