function cptState =...
    removeSpot(cptState, shouldQueryUser)
%
%
%REMOVESPOT removes a spot from the spots and particles structure
%  removes a spot from the spots and particles structure


del = false;
CurrentFrameWithinParticle = find(cptState.Particles...
    {cptState.CurrentChannelIndex}...
    (cptState.CurrentParticle).Frame==cptState.CurrentFrame);

cptState.lastParticle = cptState.CurrentParticle;
cptState.PreviousParticle = cptState.CurrentParticle;
cptState.ManualZFlag = false;

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
else
    del = true;
end

if del
    
    ind = cptState.Particles{cptState.CurrentChannelIndex}...
        (cptState.CurrentParticle).Index(CurrentFrameWithinParticle);
    
    onlyFrame = length(cptState.Particles...
        {cptState.CurrentChannelIndex}(cptState.CurrentParticle).Frame) == 1;
    
    if onlyFrame
        cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle) = [];
    else
        particleFields = fieldnames(cptState.Particles{cptState.CurrentChannelIndex});
        for i = 1:numel(particleFields)
            if ~strcmpi(particleFields{i},'Nucleus') && ~strcmpi(particleFields{i},'Approved')
                try
                    cptState.Particles{cptState.CurrentChannelIndex}...
                        (cptState.CurrentParticle).(particleFields{i})...
                        (CurrentFrameWithinParticle) = [];
                end
            end
        end
    end
    %and this part changes the the index of other particles
    %in the frame.
    for i=1:length(cptState.Particles{cptState.CurrentChannelIndex})
        for j = 1:length(cptState.Particles{cptState.CurrentChannelIndex}(i).Frame)
            
            if cptState.Particles{cptState.CurrentChannelIndex}(i).Frame(j) == cptState.CurrentFrame
                if cptState.Particles{cptState.CurrentChannelIndex}(i).Index(j) > ind
                    cptState.Particles{cptState.CurrentChannelIndex}(i).Index(j) =...
                        cptState.Particles{cptState.CurrentChannelIndex}(i).Index(j) - 1;
                end
            end
            
        end
    end
    
    %and this part deletes from the spots structure.
    CurrentSpot = cptState.CurrentParticleIndex;
    cptState.Spots{cptState.CurrentChannelIndex}(cptState.CurrentFrame).Fits(CurrentSpot)= [];
    if isempty(cptState.Spots{cptState.CurrentChannelIndex}(cptState.CurrentFrame).Fits)
        cptState.Spots{cptState.CurrentChannelIndex}(cptState.CurrentFrame).Fits = [];
    end
    
    %now delete from spotfilter
    spotRow = cptState.SpotFilter{cptState.CurrentChannelIndex}(cptState.CurrentFrame,:);
    spotRow(CurrentSpot) = [];
    spotRow(end+1) = NaN;
    cptState.SpotFilter{cptState.CurrentChannelIndex}(cptState.CurrentFrame,:) = spotRow;
    disp('Spot deleted successfully.');
    
    
end

end

