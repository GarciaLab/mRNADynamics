function [Spots, SpotFilter, CurrentFrame, ...
    CurrentParticle, Particles, ManualZFlag, lastParticle, PreviousParticle] =...
    ...
    removeSpot(Frames, CurrentFrame, ...
    ...   
    CurrentChannelIndex, CurrentParticle,...
    CurrentParticleIndex, Particles, Spots, SpotFilter, shouldQueryUser)
%
%REMOVESPOT removes a spot from the spots and particles structure 
%  removes a spot from the spots and particles structure 


numParticles = length(Particles{CurrentChannelIndex});

del = false;
CurrentFrameWithinParticle = find(Frames==CurrentFrame);

lastParticle = CurrentParticle;
PreviousParticle = CurrentParticle;
ManualZFlag = false;

if shouldQueryUser
if ~isempty(CurrentFrameWithinParticle)
    choice = questdlg('Are you sure you want to delete this spot? This can''t be undone.', ...
        '', 'Delete spot','Cancel','Cancel');
    switch choice
        case 'Delete spot'
            disp 'Deleting spot.'
            del = true;
        case 'Cancel'
            disp 'Spot deletion cancelled.'
    end
end
end

if del
    
    ind = Particles{CurrentChannelIndex}(CurrentParticle).Index(CurrentFrameWithinParticle);
    onlyFrame = length(Particles{CurrentChannelIndex}(CurrentParticle).Frame) == 1;
    if onlyFrame
        Particles{CurrentChannelIndex}(CurrentParticle) = [];
        numParticles = numParticles - 1;
    else
        particleFields = fieldnames(Particles{CurrentChannelIndex});
        for i = 1:numel(particleFields)
            if ~strcmpi(particleFields{i},'Nucleus') && ~strcmpi(particleFields{i},'Approved')
                try
                    Particles{CurrentChannelIndex}...
                        (CurrentParticle).(particleFields{i})...
                        (CurrentFrameWithinParticle) = [];
                end
            end
        end
    end
    %and this part changes the the index of other particles
    %in the frame.
    for i=1:length(Particles{CurrentChannelIndex})
        for j = 1:length(Particles{CurrentChannelIndex}(i).Frame)
            if Particles{CurrentChannelIndex}(i).Frame(j) == CurrentFrame
                if Particles{CurrentChannelIndex}(i).Index(j) > ind
                    Particles{CurrentChannelIndex}(i).Index(j) = Particles{CurrentChannelIndex}(i).Index(j) - 1;
                end
            end
        end
    end
    
    %and this part deletes from the spots structure.
    CurrentSpot = CurrentParticleIndex; %renaming this to make it clear what it actually is
    Spots{CurrentChannelIndex}(CurrentFrame).Fits(CurrentSpot)= [];
    if isempty(Spots{CurrentChannelIndex}(CurrentFrame).Fits)
        Spots{CurrentChannelIndex}(CurrentFrame).Fits = [];
    end
    
    %now delete from spotfilter
    spotRow = SpotFilter{CurrentChannelIndex}(CurrentFrame,:);
    spotRow(CurrentSpot) = [];
    spotRow(end+1) = NaN;
    SpotFilter{CurrentChannelIndex}(CurrentFrame,:) = spotRow;
    disp('Spot deleted successfully.');
    

end

end

