function schnitzcells=SeparateNuclearTraces(CurrentNucleus,CurrentFrame,schnitzcells, FrameInfo, ncFrames)

%Separate the particle trace at the specified position


fieldNames = fields(schnitzcells);
fieldNames(cellfun(@(x) strcmpi(x, 'Approved'), fieldNames)) = [];

FrameFilter=schnitzcells(CurrentNucleus).frames < CurrentFrame;

%Create a gap for a new particle
NewSchnitzCells(1:CurrentNucleus)=schnitzcells(1:CurrentNucleus);
NewSchnitzCells(CurrentNucleus+2:length(schnitzcells)+1)=schnitzcells(CurrentNucleus+1:end);

for f = 1:length(fieldNames)
    
    %Delete the information in the current particle
    NewSchnitzCells(CurrentNucleus).P = schnitzcells(CurrentNucleus).P;
    NewSchnitzCells(CurrentNucleus).E = schnitzcells(CurrentNucleus).E;
    NewSchnitzCells(CurrentNucleus).D = schnitzcells(CurrentNucleus).D;
    NewSchnitzCells(CurrentNucleus).frames = schnitzcells(CurrentNucleus).frames(FrameFilter);
    NewSchnitzCells(CurrentNucleus).cenx = schnitzcells(CurrentNucleus).cenx(FrameFilter);
    NewSchnitzCells(CurrentNucleus).ceny = schnitzcells(CurrentNucleus).ceny(FrameFilter);
    NewSchnitzCells(CurrentNucleus).len = schnitzcells(CurrentNucleus).len(FrameFilter);
    NewSchnitzCells(CurrentNucleus).cellno = schnitzcells(CurrentNucleus).cellno(FrameFilter);
    NewSchnitzCells(CurrentNucleus).AlreadyUsed= schnitzcells(CurrentNucleus).AlreadyUsed;
    NewSchnitzCells(CurrentNucleus).ExtendedIntoFutureAlready= schnitzcells(CurrentNucleus).ExtendedIntoFutureAlready;
    NewSchnitzCells(CurrentNucleus).StitchedTo= schnitzcells(CurrentNucleus).StitchedTo;
    NewSchnitzCells(CurrentNucleus).StitchedFrom= schnitzcells(CurrentNucleus).StitchedFrom;
    NewSchnitzCells(CurrentNucleus).Fluo = schnitzcells(CurrentNucleus).Fluo(FrameFilter, :);
    minFrame = min(NewSchnitzCells(CurrentNucleus).frames );
    NewSchnitzCells(CurrentNucleus).cycle = 14-length(ncFrames(ncFrames>minFrame));
    if isfield(schnitzcells,'timeSinceAnaphase')
        ncFrames(ncFrames==0) = 1;
        ind = find(isnan(ncFrames));
        ncFrames(ind) = ncFrames(ind-1);
        time = [FrameInfo.Time]/60; %frame times in minutes 
        ncTimes = time(ncFrames);

        NewSchnitzCells(CurrentNucleus).timeSinceAnaphase = time(NewSchnitzCells(CurrentNucleus).frames) - ncTimes(NewSchnitzCells(CurrentNucleus).cycle-8);

    end
    NewSchnitzCells(CurrentNucleus).FrameApproved = schnitzcells(CurrentNucleus).FrameApproved(FrameFilter);
    if isfield(schnitzcells,'APpos')
        NewSchnitzCells(CurrentNucleus).APpos = schnitzcells(CurrentNucleus).APpos(FrameFilter);
    end
    if isfield(schnitzcells,'DVpos')
        NewSchnitzCells(CurrentNucleus).DVpos = schnitzcells(CurrentNucleus).DVpos(FrameFilter);
    end
    NewSchnitzCells(CurrentNucleus).Approved=0;
    NewSchnitzCells(CurrentNucleus+1).Checked=0;
    NewSchnitzCells(CurrentNucleus+1).Flag=0;

    %Move the information to the new particle
    NewSchnitzCells(CurrentNucleus+1).P = schnitzcells(CurrentNucleus).P;
    NewSchnitzCells(CurrentNucleus+1).E = schnitzcells(CurrentNucleus).E;
    NewSchnitzCells(CurrentNucleus+1).D = schnitzcells(CurrentNucleus).D;
    NewSchnitzCells(CurrentNucleus+1).frames = schnitzcells(CurrentNucleus).frames(~FrameFilter);
    NewSchnitzCells(CurrentNucleus+1).cenx = schnitzcells(CurrentNucleus).cenx(~FrameFilter);
    NewSchnitzCells(CurrentNucleus+1).ceny = schnitzcells(CurrentNucleus).ceny(~FrameFilter);
    NewSchnitzCells(CurrentNucleus+1).len = schnitzcells(CurrentNucleus).len(~FrameFilter);
    NewSchnitzCells(CurrentNucleus+1).cellno = schnitzcells(CurrentNucleus).cellno(~FrameFilter);
    NewSchnitzCells(CurrentNucleus+1).AlreadyUsed= schnitzcells(CurrentNucleus).AlreadyUsed;
    NewSchnitzCells(CurrentNucleus+1).ExtendedIntoFutureAlready= schnitzcells(CurrentNucleus).ExtendedIntoFutureAlready;
    NewSchnitzCells(CurrentNucleus+1).StitchedTo= schnitzcells(CurrentNucleus).StitchedTo;
    NewSchnitzCells(CurrentNucleus+1).StitchedFrom= schnitzcells(CurrentNucleus).StitchedFrom;
    NewSchnitzCells(CurrentNucleus+1).Fluo = schnitzcells(CurrentNucleus).Fluo(~FrameFilter, :);
    minFrame = min(NewSchnitzCells(CurrentNucleus+1).frames );
    NewSchnitzCells(CurrentNucleus+1).cycle = 14-length(ncFrames(ncFrames>minFrame));
    if isfield(schnitzcells,'timeSinceAnaphase')
        ncFrames(ncFrames==0) = 1;
        ind = find(isnan(ncFrames));
        ncFrames(ind) = ncFrames(ind-1);
        time = [FrameInfo.Time]/60; %frame times in minutes 
        ncTimes = time(ncFrames);

        NewSchnitzCells(CurrentNucleus+1).timeSinceAnaphase = time(NewSchnitzCells(CurrentNucleus+1).frames) - ncTimes(NewSchnitzCells(CurrentNucleus+1).cycle-8);

    end
    NewSchnitzCells(CurrentNucleus+1).FrameApproved = schnitzcells(CurrentNucleus).FrameApproved(~FrameFilter);
    if isfield(schnitzcells,'APpos')
        NewSchnitzCells(CurrentNucleus+1).APpos = schnitzcells(CurrentNucleus).APpos(~FrameFilter);
    end
    if isfield(schnitzcells,'DVpos')
        NewSchnitzCells(CurrentNucleus+1).DVpos = schnitzcells(CurrentNucleus).DVpos(~FrameFilter);
    end
    NewSchnitzCells(CurrentNucleus+1).Approved=0;
    NewSchnitzCells(CurrentNucleus+1).Checked=0;
    NewSchnitzCells(CurrentNucleus+1).Flag=0;
    %NewSchnitzCells(CurrentNucleus+1).Nucleus=[];


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