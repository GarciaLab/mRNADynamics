function [Spots, SpotFilter, ZoomMode, GlobalZoomMode, CurrentFrame, ...
    CurrentParticle, Particles, ManualZFlag, DisplayRange, lastParticle, PreviousParticle] =...
    removeSpot(ZoomMode, GlobalZoomMode, Frames, CurrentFrame, ...
    CurrentChannel, CurrentParticle, CurrentParticleIndex, Particles, Spots, SpotFilter, ...
    numParticles, ManualZFlag, DisplayRange, lastParticle, PreviousParticle)
%REMOVESPOT Summary of this function goes here
%   Detailed explanation goes here

%Check that we're in zoom mode. If not, set it up.
if ~(ZoomMode || GlobalZoomMode)
    disp('You need to be in Zoom Mode to do this. You can switch using ''o'' or ''+''. Run the ''#'' command again.')
else
    %delete from particles
    del = 0;
    choice = questdlg('Are you sure you want to delete this spot? This can''t be undone.', ...
        '', 'Delete spot','Cancel','Cancel');
    switch choice
        case 'Delete spot'
            disp 'Deleting spot.'
            del = 1;
        case 'Cancel'
            disp 'Spot deletion cancelled.'
            del = 0;
    end

    if del
        CurrentFrameWithinParticle = find(Frames==CurrentFrame);
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
        try
            SpotFilter{CurrentChannel}(CurrentFrame,:) = spotRow;
        catch
            error('There probably wasn''t a spot in the frame you were trying to delete.')
        end
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
            DisplayRange=[];
        elseif CurrentFrame > 1
            CurrentFrame=CurrentFrame-1;
            ManualZFlag=0;
            ParticleToFollow=[];
            DisplayRange=[];
            PreviousParticle = 0; % this is done so that the trace is updated
        elseif CurrentFrame < length({Spots{1}.Fits})
            CurrentFrame=CurrentFrame+1;
            ManualZFlag=0;
            ParticleToFollow=[];
            DisplayRange=[];
            PreviousParticle = 0; % this is done so that the trace is updated
        else
            error('something''s wrong.')
        end
        disp 'Spot deleted successfully. Trace figures will refresh after switching particles.'
    elseif CurrentFrame > 1
        CurrentFrame=CurrentFrame-1;
        ManualZFlag=0;
    elseif CurrentFrame < length({Spots{1}.Fits})
        CurrentFrame=CurrentFrame+1;
        ManualZFlag=0;
    else
        error('something''s wrong.')
    end
    disp 'Spot deleted successfully. Trace figures will refresh after switching particles.'
end
ZoomMode=0;
GlobalZoomMode=0;
end

