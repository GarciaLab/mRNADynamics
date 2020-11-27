function [Particles, PreviousParticle] = separateTraces(Particles, ...
    CurrentChannel, CurrentFrame, CurrentParticle)
%SEPARATETRACES Summary of this function goes here
%   Detailed explanation goes here

%The separated particle (the trace following current frame) won't have a nucleus assigned!
PreviousParticle=0;
%Check that the particle does actually exist in this frame
if ~(Particles{CurrentChannel}(CurrentParticle).Frame(1)==CurrentFrame)
    if sum(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame)
        Particles{CurrentChannel}=SeparateParticleTraces(CurrentParticle,CurrentFrame,Particles{CurrentChannel});
    end
elseif length(Particles{CurrentChannel}(CurrentParticle).Frame)==1
    Particles{CurrentChannel}(CurrentParticle).Nucleus=[];
else
    disp('Cannot divide a trace at the first time point')
end
end

