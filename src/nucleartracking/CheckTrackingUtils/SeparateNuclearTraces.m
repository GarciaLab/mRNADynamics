function schnitzcells=SeparateNuclearTraces(CurrentNucleus,CurrentFrame,schnitzcells, FrameInfo, ncFrames)

%Separate the particle trace at the specified position


fieldNames = fields(schnitzcells);
fieldNames(cellfun(@(x) strcmpi(x, 'Approved'), fieldNames)) = [];

FrameFilter=schnitzcells(CurrentNucleus).frames < CurrentFrame;
%Create a gap for a new particle
NewSchnitzCells(1:CurrentNucleus)=schnitzcells(1:CurrentNucleus);
NewSchnitzCells(CurrentNucleus+2:length(schnitzcells)+1)=schnitzcells(CurrentNucleus+1:end);


%Delete the information in the current particle
if isfield(schnitzcells, 'P')
    NewSchnitzCells(CurrentNucleus).P = schnitzcells(CurrentNucleus).P;
end
if isfield(schnitzcells, 'E')
    NewSchnitzCells(CurrentNucleus).E = schnitzcells(CurrentNucleus).E;
end
if isfield(schnitzcells, 'D')
    NewSchnitzCells(CurrentNucleus).D = schnitzcells(CurrentNucleus).D;
end
if isfield(schnitzcells, 'frames')
    NewSchnitzCells(CurrentNucleus).frames = schnitzcells(CurrentNucleus).frames(FrameFilter);
end
if isfield(schnitzcells, 'cenx')
    NewSchnitzCells(CurrentNucleus).cenx = schnitzcells(CurrentNucleus).cenx(FrameFilter);
end
if isfield(schnitzcells, 'ceny')
    NewSchnitzCells(CurrentNucleus).ceny = schnitzcells(CurrentNucleus).ceny(FrameFilter);
end
if isfield(schnitzcells, 'len')
    NewSchnitzCells(CurrentNucleus).len = schnitzcells(CurrentNucleus).len(FrameFilter);
end
if isfield(schnitzcells, 'cellno')
    NewSchnitzCells(CurrentNucleus).cellno = schnitzcells(CurrentNucleus).cellno(FrameFilter);
end
if isfield(schnitzcells, 'AlreadyUsed')
    NewSchnitzCells(CurrentNucleus).AlreadyUsed= schnitzcells(CurrentNucleus).AlreadyUsed;
end
if isfield(schnitzcells, 'ExtendedIntoFutureAlready')
    NewSchnitzCells(CurrentNucleus).ExtendedIntoFutureAlready= schnitzcells(CurrentNucleus).ExtendedIntoFutureAlready;
end
if isfield(schnitzcells, 'StitchedTo')
    NewSchnitzCells(CurrentNucleus).StitchedTo= schnitzcells(CurrentNucleus).StitchedTo;
end
if isfield(schnitzcells, 'StitchedFrom')
    NewSchnitzCells(CurrentNucleus).StitchedFrom= schnitzcells(CurrentNucleus).StitchedFrom;
end
FluoVars = fieldNames(contains(fieldNames, 'Fluo'));
for i = 1:length(fieldNames(contains(fieldNames, 'Fluo')))
    FluoVar = FluoVars{i};
    NewSchnitzCells(CurrentNucleus).(FluoVar) = schnitzcells(CurrentNucleus).(FluoVar)(FrameFilter, :);
end
if isfield(schnitzcells, 'anaphaseFrame')
    NewSchnitzCells(CurrentNucleus).anaphaseFrame = schnitzcells(CurrentNucleus).anaphaseFrame;
end
if isfield(schnitzcells, 'inferredAnaphaseFrame')
    NewSchnitzCells(CurrentNucleus).inferredAnaphaseFrame = schnitzcells(CurrentNucleus).inferredAnaphaseFrame;
end
minFrame = min(NewSchnitzCells(CurrentNucleus).frames );
if isfield(schnitzcells, 'cycle')
    NewSchnitzCells(CurrentNucleus).cycle = 14-length(ncFrames(ncFrames>minFrame));
end
if isfield(schnitzcells,'timeSinceAnaphase')
    ncFrames(ncFrames==0) = 1;
    ind = find(isnan(ncFrames));
    ncFrames(ind) = ncFrames(ind-1);
    time = [FrameInfo.Time]/60; %frame times in minutes 
    ncTimes = time(ncFrames);
    if isfield(schnitzcells, 'anaphaseFrame')
        if ~isempty(NewSchnitzCells(CurrentNucleus).anaphaseFrame(NewSchnitzCells(CurrentNucleus).anaphaseFrame > 0))
            NewSchnitzCells(CurrentNucleus).timeSinceAnaphase = time(NewSchnitzCells(CurrentNucleus).frames) -time(NewSchnitzCells(CurrentNucleus).anaphaseFrame);
        else
            NewSchnitzCells(CurrentNucleus).timeSinceAnaphase = time(NewSchnitzCells(CurrentNucleus).frames) - ncTimes(NewSchnitzCells(CurrentNucleus).cycle-8);
        end
    else
        NewSchnitzCells(CurrentNucleus).timeSinceAnaphase = time(NewSchnitzCells(CurrentNucleus).frames) - ncTimes(NewSchnitzCells(CurrentNucleus).cycle-8);
    end
end
if isfield(schnitzcells, 'FrameApproved')
NewSchnitzCells(CurrentNucleus).FrameApproved = schnitzcells(CurrentNucleus).FrameApproved(FrameFilter);
end
if isfield(schnitzcells,'APpos')
    NewSchnitzCells(CurrentNucleus).APpos = schnitzcells(CurrentNucleus).APpos(FrameFilter);
end
if isfield(schnitzcells,'DVpos')
    NewSchnitzCells(CurrentNucleus).DVpos = schnitzcells(CurrentNucleus).DVpos(FrameFilter);
end
if isfield(schnitzcells,'Approved')
    NewSchnitzCells(CurrentNucleus).Approved=schnitzcells(CurrentNucleus).Approved;
end
if isfield(schnitzcells,'Checked')
    NewSchnitzCells(CurrentNucleus+1).Checked=0;
end
if isfield(schnitzcells,'Flag')
    NewSchnitzCells(CurrentNucleus+1).Flag=0;
end
%Move the information to the new particle
if isfield(schnitzcells, 'P')
    NewSchnitzCells(CurrentNucleus+1).P = schnitzcells(CurrentNucleus).P;
end
if isfield(schnitzcells, 'E')
	NewSchnitzCells(CurrentNucleus+1).E = schnitzcells(CurrentNucleus).E;
end
if isfield(schnitzcells, 'D')
    NewSchnitzCells(CurrentNucleus+1).D = schnitzcells(CurrentNucleus).D;
end
if isfield(schnitzcells, 'frames')
    NewSchnitzCells(CurrentNucleus+1).frames = schnitzcells(CurrentNucleus).frames(~FrameFilter);
end
if isfield(schnitzcells, 'cenx')
    NewSchnitzCells(CurrentNucleus+1).cenx = schnitzcells(CurrentNucleus).cenx(~FrameFilter);
end
if isfield(schnitzcells, 'ceny')
    NewSchnitzCells(CurrentNucleus+1).ceny = schnitzcells(CurrentNucleus).ceny(~FrameFilter);
end
if isfield(schnitzcells, 'len')
    NewSchnitzCells(CurrentNucleus+1).len = schnitzcells(CurrentNucleus).len(~FrameFilter);
end
if isfield(schnitzcells, 'cellno')
    NewSchnitzCells(CurrentNucleus+1).cellno = schnitzcells(CurrentNucleus).cellno(~FrameFilter);
end
if isfield(schnitzcells, 'AlreadyUsed')
    NewSchnitzCells(CurrentNucleus+1).AlreadyUsed= schnitzcells(CurrentNucleus).AlreadyUsed;
end
if isfield(schnitzcells, 'ExtendedIntoFutureAlready')
    NewSchnitzCells(CurrentNucleus+1).ExtendedIntoFutureAlready= schnitzcells(CurrentNucleus).ExtendedIntoFutureAlready;
end
if isfield(schnitzcells, 'StitchedTo')
    NewSchnitzCells(CurrentNucleus+1).StitchedTo= schnitzcells(CurrentNucleus).StitchedTo;
end
if isfield(schnitzcells, 'StitchedFrom')
    NewSchnitzCells(CurrentNucleus+1).StitchedFrom= schnitzcells(CurrentNucleus).StitchedFrom;
end
if isfield(schnitzcells, 'Fluo')
    NewSchnitzCells(CurrentNucleus+1).Fluo = schnitzcells(CurrentNucleus).Fluo(~FrameFilter, :);
end
for i = 1:length(FluoVars)
    FluoVar = FluoVars{i};
    NewSchnitzCells(CurrentNucleus+1).(FluoVar) = schnitzcells(CurrentNucleus).(FluoVar)(~FrameFilter, :);
end

if isfield(schnitzcells, 'anaphaseFrame')
    NewSchnitzCells(CurrentNucleus+1).anaphaseFrame = [];
end
if isfield(schnitzcells, 'inferredAnaphaseFrame')
    NewSchnitzCells(CurrentNucleus+1).inferredAnaphaseFrame = false;
end

minFrame = min(NewSchnitzCells(CurrentNucleus+1).frames );
if isfield(schnitzcells, 'cycle')
NewSchnitzCells(CurrentNucleus+1).cycle = 14-length(ncFrames(ncFrames>minFrame));
end

if isfield(schnitzcells,'timeSinceAnaphase')
    ncFrames(ncFrames==0) = 1;
    ind = find(isnan(ncFrames));
    ncFrames(ind) = ncFrames(ind-1);
    time = [FrameInfo.Time]/60; %frame times in minutes 
    ncTimes = time(ncFrames);
    if isfield(schnitzcells, 'anaphaseFrame')
        if ~isempty(NewSchnitzCells(CurrentNucleus+1).anaphaseFrame)
            NewSchnitzCells(CurrentNucleus+1).timeSinceAnaphase = time(NewSchnitzCells(CurrentNucleus+1).frames) -time(NewSchnitzCells(CurrentNucleus+1).anaphaseFrame);
        else
            NewSchnitzCells(CurrentNucleus+1).timeSinceAnaphase = time(NewSchnitzCells(CurrentNucleus+1).frames) - ncTimes(NewSchnitzCells(CurrentNucleus+1).cycle-8);
        end
    else
        NewSchnitzCells(CurrentNucleus+1).timeSinceAnaphase = time(NewSchnitzCells(CurrentNucleus+1).frames) - ncTimes(NewSchnitzCells(CurrentNucleus+1).cycle-8);
    end
end
if isfield(schnitzcells, 'FrameApproved')
    NewSchnitzCells(CurrentNucleus+1).FrameApproved = schnitzcells(CurrentNucleus).FrameApproved(~FrameFilter);
end
if isfield(schnitzcells,'APpos')
    NewSchnitzCells(CurrentNucleus+1).APpos = schnitzcells(CurrentNucleus).APpos(~FrameFilter);
end
if isfield(schnitzcells,'DVpos')
    NewSchnitzCells(CurrentNucleus+1).DVpos = schnitzcells(CurrentNucleus).DVpos(~FrameFilter);
end
if isfield(schnitzcells, 'Approved')
    NewSchnitzCells(CurrentNucleus+1).Approved=schnitzcells(CurrentNucleus).Approved;
end
if isfield(schnitzcells, 'Checked')
    NewSchnitzCells(CurrentNucleus+1).Checked=0;
end
if isfield(schnitzcells, 'FirstPass')
    NewSchnitzCells(CurrentNucleus+1).FirstPass=1;
end
if isfield(schnitzcells, 'Flag')
    NewSchnitzCells(CurrentNucleus+1).Flag=0;
end


schnitzcells=NewSchnitzCells;