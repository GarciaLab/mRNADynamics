function schnitzcells=SeparateNuclearTraces(CurrentNucleus,CurrentFrame,schnitzcells)

%Separate the particle trace at the specified position


fieldNames = fields(schnitzcells);
fieldNames(cellfun(@(x) strcmpi(x, 'Approved'), fieldNames)) = [];

FrameFilter=schnitzcells.Frame < CurrentFrame;

%Create a gap for a new particle
NewSchnitzCells(1:CurrentNucleus)=schnitzcells(1:CurrentNucleus);
NewSchnitzCells(CurrentNucleus+2:length(schnitzcells)+1)=schnitzcells(CurrentNucleus+1:end);

for f = 1:length(fieldNames)
    
    %Delete the information in the current particle
    
    NewSchnitzCells(CurrentNucleus).(fieldNames{f}) = NewSchnitzCells(CurrentNucleus).(fieldNames{f})(FrameFilter);
    NewSchnitzCells(CurrentNucleus).Approved=0;

    %Move the information to the new particle
    
    NewSchnitzCells(CurrentNucleus+1).(fieldNames{f})=schnitzcells(CurrentNucleus).(fieldNames{f})(~FrameFilter);
    NewSchnitzCells(CurrentNucleus+1).Approved=0;
    NewSchnitzCells(CurrentNucleus+1).Nucleus=[];


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


schnitzcells=NewSchnitzCells;