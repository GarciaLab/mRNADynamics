% JoinNuclearTraces.m
% author: Gabriella Martini
% date created: 9/13/20
% date last modified: 9/13/20
function schnitzcells=JoinNuclearTraces(OriginalNucleus,ClickedNucleus,schnitzcells,FrameInfo,ncFrames)

%This function joins two nuclear traces and renumbers all nuclei in the
%schnitzcells structure accordingly

%Transfer the information to the original particle
if length(unique([schnitzcells(OriginalNucleus).frames;schnitzcells(ClickedNucleus).frames])) ==...
        length([schnitzcells(OriginalNucleus).frames;schnitzcells(ClickedNucleus).frames])
    schnitzcells(OriginalNucleus).frames=[schnitzcells(OriginalNucleus).frames;schnitzcells(ClickedNucleus).frames];

    if isfield(schnitzcells,'cenx')
        schnitzcells(OriginalNucleus).cenx=[schnitzcells(OriginalNucleus).cenx,schnitzcells(ClickedNucleus).cenx];
    end
    if isfield(schnitzcells,'ceny')
        schnitzcells(OriginalNucleus).ceny=[schnitzcells(OriginalNucleus).ceny,schnitzcells(ClickedNucleus).ceny];
    end
    if isfield(schnitzcells,'len')
        schnitzcells(OriginalNucleus).len=[schnitzcells(OriginalNucleus).len,schnitzcells(ClickedNucleus).len];
    end
    if isfield(schnitzcells,'cellno')
        schnitzcells(OriginalNucleus).cellno=[schnitzcells(OriginalNucleus).cellno,schnitzcells(ClickedNucleus).cellno];
    end
    if isfield(schnitzcells, 'AlreadyUsed')
        schnitzcells(OriginalNucleus).AlreadyUsed = max([schnitzcells(OriginalNucleus).AlreadyUsed, schnitzcells(ClickedNucleus).AlreadyUsed]);
    end
    if isfield(schnitzcells, 'ExtendedIntoFutureAlready')
        schnitzcells(OriginalNucleus).ExtendedIntoFutureAlready = max([schnitzcells(OriginalNucleus).ExtendedIntoFutureAlready, schnitzcells(ClickedNucleus).ExtendedIntoFutureAlready]);
    end
    if isfield(schnitzcells,'StitchedTo')
        schnitzcells(OriginalNucleus).StitchedTo=[schnitzcells(OriginalNucleus).StitchedTo,...
            schnitzcells(ClickedNucleus).StitchedTo];
    end
    if isfield(schnitzcells,'StitchedFrom')
        schnitzcells(OriginalNucleus).StitchedFrom=[schnitzcells(OriginalNucleus).StitchedFrom,...
            schnitzcells(ClickedNucleus).StitchedFrom];
    end
    if isfield(schnitzcells,'Fluo')
        schnitzcells(OriginalNucleus).Fluo=[schnitzcells(OriginalNucleus).Fluo;...
            schnitzcells(ClickedNucleus).Fluo];
    end
    schnitzcells(OriginalNucleus).cycle=max([schnitzcells(OriginalNucleus).cycle, schnitzcells(ClickedNucleus).cycle]);

    if isfield(schnitzcells,'timeSinceAnaphase')
        ncFrames(ncFrames==0) = 1;
        ind = find(isnan(ncFrames));
        ncFrames(ind) = ncFrames(ind-1);
        time = [FrameInfo.Time]/60; %frame times in minutes 
        ncTimes = time(ncFrames);

        schnitzcells(OriginalNucleus).timeSinceAnaphase = time(schnitzcells(OriginalNucleus).frames) - ncTimes(schnitzcells(OriginalNucleus).cycle-8);

    end

    % Remaining fields: 
    % APpos, DVpos, 

    %Particles(OriginalParticle).nc=[Particles(OriginalParticle).nc,Particles(ClickedParticle).nc];
    if isfield(schnitzcells,'FrameApproved')
        schnitzcells(OriginalNucleus).FrameApproved=logical([schnitzcells(OriginalNucleus).FrameApproved,schnitzcells(ClickedNucleus).FrameApproved]);
    end
    if isfield(schnitzcells,'APpos')
        schnitzcells(OriginalNucleus).APpos=[schnitzcells(OriginalNucleus).APpos,...
            schnitzcells(ClickedNucleus).APpos];
    end
    if isfield(schnitzcells,'DVpos')
        schnitzcells(OriginalNucleus).DVpos=[schnitzcells(OriginalNucleus).DVpos,...
            schnitzcells(ClickedNucleus).DVpos];
    end

    schnitzcells(OriginalNucleus).Approved=0;
    schnitzcells(OriginalNucleus).Checked=0;
    schnitzcells(OriginalNucleus).Flag=0;



    %Now, get rid of the clicked particle
    schnitzcells=schnitzcells([1:ClickedNucleus-1,ClickedNucleus+1:end]);

    %Deals with the indexing changing because of the removal of
    %the old particle.
     if ClickedNucleus<OriginalNucleus
         OriginalNucleus=OriginalNucleus-1;
     end

    %Sort the frames within the particle. This is useful if we
    %connected to a particle that came before.
    [~,Permutations]=sort(schnitzcells(OriginalNucleus).frames);
    schnitzcells(OriginalNucleus).frames=schnitzcells(OriginalNucleus).frames(Permutations);

    if isfield(schnitzcells,'cenx')
        schnitzcells(OriginalNucleus).cenx=schnitzcells(OriginalNucleus).cenx(Permutations);
    end
    if isfield(schnitzcells,'ceny')
        schnitzcells(OriginalNucleus).ceny=schnitzcells(OriginalNucleus).ceny(Permutations);
    end
    if isfield(schnitzcells,'len')
        schnitzcells(OriginalNucleus).len=schnitzcells(OriginalNucleus).len(Permutations);
    end
    if isfield(schnitzcells,'cellno')
        schnitzcells(OriginalNucleus).cellno=schnitzcells(OriginalNucleus).cellno(Permutations);
    end
    if isfield(schnitzcells,'Fluo')
        schnitzcells(OriginalNucleus).Fluo=schnitzcells(OriginalNucleus).Fluo(Permutations,:);
    end
    if isfield(schnitzcells,'timeSinceAnaphase')
        schnitzcells(OriginalNucleus).timeSinceAnaphase=schnitzcells(OriginalNucleus).timeSinceAnaphase(Permutations);
    end
    if isfield(schnitzcells,'FrameApproved')
        schnitzcells(OriginalNucleus).FrameApproved=schnitzcells(OriginalNucleus).FrameApproved(Permutations);
    end
    if isfield(schnitzcells,'APpos')
        schnitzcells(OriginalNucleus).APpos=schnitzcells(OriginalNucleus).APpos(Permutations);
    end
    if isfield(schnitzcells,'DVpos')
        schnitzcells(OriginalNucleus).DVpos=schnitzcells(OriginalNucleus).DVpos(Permutations);
    end
    


    % 9/4 EL: Added the if statement to remove error of the index (Permutations)
    % exceeding matrix dimension. Did this always create an error when the
    % frames of the original particle was not approved before?
    if ~isfield(schnitzcells,'FrameApproved')
        schnitzcells(OriginalNucleus).FrameApproved=schnitzcells(OriginalNucleus).FrameApproved(Permutations);
    end
else
    disp('These two schnitz cells cannot be combined as they contain overlapping frames.')
end

end