function Particles=SeparateParticleTraces(CurrentParticle,CurrentFrame,Particles)

%Separate the particle trace at the specified position


fieldNames = fields(Particles);
fieldNames(cellfun(@(x) strcmpi(x, 'Approved'), fieldNames) |  cellfun(@(x) strcmpi(x, 'Nucleus'), fieldNames) ) = [];

FrameFilter=Particles(CurrentParticle).Frame < CurrentFrame;

%Create a gap for a new particle
NewParticles(1:CurrentParticle)=Particles(1:CurrentParticle);
NewParticles(CurrentParticle+2:length(Particles)+1)=Particles(CurrentParticle+1:end);

for f = 1:length(fieldNames)
    
    %Delete the information in the current particle
    
    NewParticles(CurrentParticle).(fieldNames{f}) = NewParticles(CurrentParticle).(fieldNames{f})(FrameFilter);
    NewParticles(CurrentParticle).Approved=0;

    %Move the information to the new particle
    
    NewParticles(CurrentParticle+1).(fieldNames{f})=Particles(CurrentParticle).(fieldNames{f})(~FrameFilter);
    NewParticles(CurrentParticle+1).Approved=0;
    NewParticles(CurrentParticle+1).Nucleus=[];


end

% NewParticles(CurrentParticle).Frame=NewParticles(CurrentParticle).Frame(FrameFilter);
% NewParticles(CurrentParticle).Index=NewParticles(CurrentParticle).Index(FrameFilter);
% NewParticles(CurrentParticle).Approved=0;
% NewParticles(CurrentParticle).FrameApproved=NewParticles(CurrentParticle).FrameApproved(FrameFilter);
% 
% 
% %Move the information to the new particle
% NewParticles(CurrentParticle+1).Frame=Particles(CurrentParticle).Frame(~FrameFilter);
% NewParticles(CurrentParticle+1).Index=Particles(CurrentParticle).Index(~FrameFilter);
% NewParticles(CurrentParticle+1).Approved=0;
% NewParticles(CurrentParticle+1).FrameApproved=Particles(CurrentParticle).FrameApproved(~FrameFilter);
% NewParticles(CurrentParticle+1).Nucleus=[];


Particles=NewParticles;