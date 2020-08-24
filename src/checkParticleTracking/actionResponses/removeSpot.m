function cptState = removeSpot(cptState, shouldQueryUser)
%
%REMOVESPOT removes a spot from the spots and particles structure
Ch = cptState.CurrentChannelIndex;
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
    % get index and frame of spot to remove
    ind = cptState.Particles{cptState.CurrentChannelIndex}...
        (cptState.CurrentParticle).Index(CurrentFrameWithinParticle);
    frame = cptState.CurrentFrame;
    
    % get exisiting index pairs
    linkStruct = generateLinkStructure(cptState.ParticleStitchInfo{Ch},cptState.SpotFilter{Ch});    
    
    % check for overlap    
    spotLinIndex = sub2ind(size(cptState.SpotFilter{Ch}),frame,ind);
    rmPLinks = cellfun(@(x) any(ismember(x,spotLinIndex)),linkStruct.persistentLinIndices);
    rmFLinks = cellfun(@(x) any(ismember(x,spotLinIndex)),linkStruct.forbiddenLinIndices);
    
    % remove links involving this spot
    cptState.ParticleStitchInfo{Ch}.persistentLinkFrameCell = cptState.ParticleStitchInfo{Ch}.persistentLinkFrameCell(~rmPLinks);
    cptState.ParticleStitchInfo{Ch}.persistentLinkIndexCell = cptState.ParticleStitchInfo{Ch}.persistentLinkIndexCell(~rmPLinks);
    cptState.ParticleStitchInfo{Ch}.forbiddenLinkFrameCell = cptState.ParticleStitchInfo{Ch}.forbiddenLinkFrameCell(~rmFLinks);
    cptState.ParticleStitchInfo{Ch}.forbiddenLinkIndexCell = cptState.ParticleStitchInfo{Ch}.forbiddenLinkIndexCell(~rmFLinks);
    
    % remove Particles entry    
    
    onlyFrame = length(cptState.Particles...
        {cptState.CurrentChannelIndex}(cptState.CurrentParticle).Frame) == 1;
    
    if onlyFrame
        cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle) = [];
    else
        particleFields = fieldnames(cptState.Particles{cptState.CurrentChannelIndex});
        for i = 1:numel(particleFields)      
          if any(strcmpi(particleFields{i},cptState.frameLevelFields))                                    
              cptState.Particles{cptState.CurrentChannelIndex}...
                  (cptState.CurrentParticle).(particleFields{i})...
                  (CurrentFrameWithinParticle) = [];   
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
    
    % likewise adjust indices in stitch structure
    % links...
    for f = 1:length(cptState.ParticleStitchInfo{Ch}.persistentLinkFrameCell)
      frameFilter = ismember(cptState.ParticleStitchInfo{Ch}.persistentLinkFrameCell{f},frame);
      indFilter = cptState.ParticleStitchInfo{Ch}.persistentLinkIndexCell{f}>ind;
      if any(frameFilter&indFilter)
        cptState.ParticleStitchInfo{Ch}.persistentLinkIndexCell{f}(frameFilter&indFilter) = ...
          cptState.ParticleStitchInfo{Ch}.persistentLinkIndexCell{f}(frameFilter&indFilter)-1;
      end
    end
    % ...and separations
    for f = 1:length(cptState.ParticleStitchInfo{Ch}.forbiddenLinkFrameCell)
      frameFilter = ismember(cptState.ParticleStitchInfo{Ch}.forbiddenLinkFrameCell{f},frame);
      indFilter = cptState.ParticleStitchInfo{Ch}.forbiddenLinkIndexCell{f}>ind;
      if any(frameFilter&indFilter)
        cptState.ParticleStitchInfo{Ch}.forbiddenLinkIndexCell{f}(frameFilter&indFilter) = ...
          cptState.ParticleStitchInfo{Ch}.forbiddenLinkIndexCell{f}(frameFilter&indFilter)-1;
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

