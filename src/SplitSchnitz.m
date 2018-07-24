function [Particles,schnitzcells]=SplitSchnitz(Particles,schnitzcells,...
    CurrentFrame,CurrentParticle)

%Splits the nucleus associated with CurrentParticle at the position
%CurrentFrame. If the particle continues after CurrentFrame it splits it as
%well.

%V2: Modified to support Laurent's schnizcells

CurrentSchnitz=Particles(CurrentParticle).Nucleus;

FrameFilter=((schnitzcells(CurrentSchnitz).frames)>=CurrentFrame);


%Create the new schnitz
Nschnitz=length(schnitzcells);

%Assign the parent and daughter information
schnitzcells(Nschnitz+1).P=0;
schnitzcells(Nschnitz+1).D=schnitzcells(CurrentSchnitz).D;
schnitzcells(Nschnitz+1).E=schnitzcells(CurrentSchnitz).E;

schnitzcells(CurrentSchnitz).D=0;
schnitzcells(CurrentSchnitz).E=0;


%Now transfer the information
Fields=fieldnames(schnitzcells);
for i=1:length(Fields)
    if ~(strcmp(Fields{i},'P')|strcmp(Fields{i},'E')|strcmp(Fields{i},'D'))
        InfoAll=getfield(schnitzcells(CurrentSchnitz),Fields{i});
        if length(InfoAll)==length(FrameFilter)
            InfoKeep=InfoAll(~FrameFilter);
            InfoMove=InfoAll(FrameFilter);

            schnitzcells(Nschnitz+1)=setfield(schnitzcells(Nschnitz+1),Fields{i},InfoMove);
            schnitzcells(CurrentSchnitz)=setfield(schnitzcells(CurrentSchnitz),Fields{i},InfoKeep);
        elseif (strcmp(Fields{i},'len')|strcmp(Fields{i},'cellno')) %I had to add this because
                                                                    %some schnitzcells had their ang and length
                                                                    %incomplete. In we leave the field empty
                                                                    %that case
            schnitzcells(Nschnitz+1)=setfield(schnitzcells(Nschnitz+1),Fields{i},[]);
        end
    end
end


%Now see if we need to split the particle
if sum(ismember(schnitzcells(Nschnitz+1).frames,Particles(CurrentParticle).Frame))
    Particles(CurrentParticle).Nucleus=Nschnitz+1;
else
    %error('Do these cases')
end

