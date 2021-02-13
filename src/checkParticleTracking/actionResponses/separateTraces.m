function [Particles, PreviousParticle] = separateTraces(Particles, ...
    CurrentChannelIndex, CurrentFrame, CurrentParticle)
%SEPARATETRACES Summary of this function goes here
%   Detailed explanation goes here

%The separated particle (the trace following current frame) won't have a nucleus assigned!
PreviousParticle=0;
%Check that the particle does actually exist in this frame
if ~(Particles{CurrentChannelIndex}(CurrentParticle).Frame(1)==CurrentFrame)
    if sum(Particles{CurrentChannelIndex}(CurrentParticle).Frame==CurrentFrame)
        Particles{CurrentChannelIndex}=SeparateParticleTraces(CurrentParticle,CurrentFrame,Particles{CurrentChannelIndex});
    end
elseif length(Particles{CurrentChannelIndex}(CurrentParticle).Frame)==1
    Particles{CurrentChannelIndex}(CurrentParticle).Nucleus=[];
else
    disp('Cannot divide a trace at the first time point')
end
end

