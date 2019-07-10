function [Spots, SpotFilter, CurrentFrame, ...
    CurrentParticle, Particles, ManualZFlag, lastParticle, PreviousParticle] =...
    ...
    removeSpot(Frames, CurrentFrame, ...
    ...   
    CurrentChannel, CurrentParticle, CurrentParticleIndex, Particles, Spots, SpotFilter, ...
    numParticles)
%
%REMOVESPOT removes a spot from the spots and particles structure 
%  removes a spot from the spots and particles structure 

del = 0;
CurrentFrameWithinParticle = find(Frames==CurrentFrame);

lastParticle = CurrentParticle;
PreviousParticle = CurrentParticle;
ManualZFlag = 0;

if ~isempty(CurrentFrameWithinParticle)
    choice = questdlg('Are you sure you want to delete this spot? This can''t be undone.', ...
        '', 'Delete spot','Cancel','Cancel');
    switch choice
        case 'Delete spot'
            disp 'Deleting spot.'
            del = 1;
        case 'Cancel'
            disp 'Spot deletion cancelled.'
    end
end

if del
    ind = Particles{CurrentChannel}(CurrentParticle).Index(CurrentFrameWithinParticle);
    onlyFrame = length(Particles{CurrentChannel}(CurrentParticle).Frame) == 1;
    if onlyFrame
        Particles{CurrentChannel}(CurrentParticle) = [];
        numParticles = numParticles - 1;
    else
        particleFields = fieldnames(Particles{CurrentChannel});
        for i = 1:numel(particleFields)
            if ~strcmpi(particleFields{i},'Nucleus') && ~strcmpi(particleFields{i},'Approved')
                try
                    Particles{CurrentChannel}(CurrentParticle).(particleFields{i})(CurrentFrameWithinParticle) = [];
                end
            end
        end
    end
    %and this part changes the the index of other particles
    %in the frame.
    for i=1:length(Particles{CurrentChannel})
        for j = 1:length(Particles{CurrentChannel}(i).Frame)
            if Particles{CurrentChannel}(i).Frame(j) == CurrentFrame
                if Particles{CurrentChannel}(i).Index(j) > ind
                    Particles{CurrentChannel}(i).Index(j) = Particles{CurrentChannel}(i).Index(j) - 1;
                end
            end
        end
    end
    %and this part deletes from the spots structure.
    CurrentSpot = CurrentParticleIndex; %renaming this to make it clear what it actually is
    Spots{CurrentChannel}(CurrentFrame).Fits(CurrentSpot)= [];
    if isempty(Spots{CurrentChannel}(CurrentFrame).Fits)
        Spots{CurrentChannel}(CurrentFrame).Fits = [];
    end
    %now delete from spotfilter
    spotRow = SpotFilter{CurrentChannel}(CurrentFrame,:);
    spotRow(CurrentSpot) = [];
    spotRow(end+1) = NaN;
    SpotFilter{CurrentChannel}(CurrentFrame,:) = spotRow;
    
    if onlyFrame
        %switch to another particle just to avoid any potential weirdness with
        %checkparticletracking refreshing. simpler version of the
        %'m' button
        NextParticle = CurrentParticle+1;
        if NextParticle>numParticles
            NextParticle=NextParticle-2; %go backwards one particle if the deleted particle was the last.
        end
        if numParticles == 1
            lastParticle = 1;
        end
        CurrentParticle=NextParticle;
        CurrentFrame=Particles{CurrentChannel}(CurrentParticle).Frame(1);
        ParticleToFollow=[];
    elseif CurrentFrame > 1
        CurrentFrame=CurrentFrame-1;
        ManualZFlag=0;
        ParticleToFollow=[];
        PreviousParticle = 0; % this is done so that the trace is updated
    elseif CurrentFrame < length({Spots{1}.Fits})
        CurrentFrame=CurrentFrame+1;
        ManualZFlag=0;
        ParticleToFollow=[];
        PreviousParticle = 0; % this is done so that the trace is updated
    else
        error('something''s wrong.')
    end
elseif CurrentFrame > 1
    CurrentFrame=CurrentFrame-1;
    ManualZFlag=0;
elseif CurrentFrame < length({Spots{1}.Fits})
    CurrentFrame=CurrentFrame+1;
    ManualZFlag=0;
else
    error('something''s wrong.')
end

disp('Spot deleted successfully. Trace figures will refresh after switching particles.')

end

