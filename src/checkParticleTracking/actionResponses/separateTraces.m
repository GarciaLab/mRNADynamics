function [Particles, PreviousParticle, ParticleStitchInfo] = separateTraces(Prefix,Particles, ...
    CurrentChannel, CurrentFrame, CurrentParticle)
%SEPARATETRACES Summary of this function goes here
%   Detailed explanation goes here

%The separated particle (the trace following current frame) won't have a nucleus assigned!
PreviousParticle=0;

%Check that the particle does actually exist in this frame
if ~(Particles{CurrentChannel}(CurrentParticle).Frame(1)==CurrentFrame) && ...
  any(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame)

  % load Particle and Frame information
  liveExperiment = LiveExperiment(Prefix);
  SimParticles = getSimParticles(liveExperiment);
  ParticleStitchInfo = getParticleStitchInfo(liveExperiment);
  FrameInfo = getFrameInfo(liveExperiment);
  
  % call splitting function
  [Particles{CurrentChannel},ParticleStitchInfo{CurrentChannel}]=SeparateParticleTraces(CurrentParticle,CurrentFrame,...
    Particles{CurrentChannel},SimParticles{CurrentChannel},ParticleStitchInfo{CurrentChannel},FrameInfo); 
  
  % save
      
elseif length(Particles{CurrentChannel}(CurrentParticle).Frame)==1
  
  Particles{CurrentChannel}(CurrentParticle).Nucleus=[];
  
else
  
  disp('Cannot divide a trace at the first time point')
end
end

