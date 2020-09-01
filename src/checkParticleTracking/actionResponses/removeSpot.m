function cptState = removeSpot(cptState, shouldQueryUser, Prefix)
%
%REMOVESPOT removes a spot from the spots and particles structure
CC = cptState.CurrentChannelIndex;
CP = cptState.CurrentParticle;
CF = cptState.CurrentFrame;

del = false;
CurrentFrameWithinParticle = find(cptState.Particles...
    {CC}...
    (CP).Frame==cptState.CurrentFrame);

cptState.lastParticle = CP;
cptState.PreviousParticle = CP;
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
    CPIndex = cptState.Particles{CC}...
        (CP).Index(CurrentFrameWithinParticle);    
    
    % get exisiting index pairs
    linkStruct = generateLinkStructure(cptState.ParticleStitchInfo{CC},cptState.SpotFilter{CC});    
    
    % check for overlap    
    spotLinIndex = sub2ind(size(cptState.SpotFilter{CC}),CF,CPIndex);
    rmPLinks = cellfun(@(x) any(ismember(x,spotLinIndex)),linkStruct.persistentLinIndices);
    rmFLinks = cellfun(@(x) any(ismember(x,spotLinIndex)),linkStruct.forbiddenLinIndices);
    
    % remove links involving this spot
    cptState.ParticleStitchInfo{CC}.persistentLinkFrameCell = cptState.ParticleStitchInfo{CC}.persistentLinkFrameCell(~rmPLinks);
    cptState.ParticleStitchInfo{CC}.persistentLinkIndexCell = cptState.ParticleStitchInfo{CC}.persistentLinkIndexCell(~rmPLinks);
    cptState.ParticleStitchInfo{CC}.forbiddenLinkFrameCell = cptState.ParticleStitchInfo{CC}.forbiddenLinkFrameCell(~rmFLinks);
    cptState.ParticleStitchInfo{CC}.forbiddenLinkIndexCell = cptState.ParticleStitchInfo{CC}.forbiddenLinkIndexCell(~rmFLinks);
    
    % update auxiliary particles structures
    cptState = removeSpotFromAuxParticles(cptState,Prefix);
    
    
    % remove Particles entry    
    
    onlyFrame = length(cptState.Particles...
        {CC}(CP).Frame) == 1;
    
    if onlyFrame
      cptState.Particles{CC}(CP) = [];
    else
      particleFields = fieldnames(cptState.Particles{CC})';
      for i = 1:numel(particleFields)      
        if any(strcmpi(particleFields{i},cptState.frameLevelFields))                                    
            cptState.Particles{CC}...
                (CP).(particleFields{i})...
                (CurrentFrameWithinParticle) = [];   
        end
      end    
    
%       otherFields = [{'FrameApproved'},particleFields(contains(particleFields,...
%         'Flags')&~contains(particleFields,'Per')),{'distShiftVec'},{'numNeighbors'}];
%       for i = 1:numel(otherFields)              
%             cptState.Particles{CC}...
%                 (CP).(otherFields{i})...
%                 (CurrentFrameWithinParticle) = [];           
%       end 
    end
    
    %and this part changes the the index of other particles
    %in the frame.
    for i=1:length(cptState.Particles{CC})
        for j = 1:length(cptState.Particles{CC}(i).Frame)
            
            if cptState.Particles{CC}(i).Frame(j) == cptState.CurrentFrame
                if cptState.Particles{CC}(i).Index(j) > CPIndex
                    cptState.Particles{CC}(i).Index(j) =...
                        cptState.Particles{CC}(i).Index(j) - 1;
                end
            end
            
        end
    end
    
    % likewise adjust indices in stitch structure
    % links...
    for f = 1:length(cptState.ParticleStitchInfo{CC}.persistentLinkFrameCell)
      frameFilter = ismember(cptState.ParticleStitchInfo{CC}.persistentLinkFrameCell{f},CF);
      indFilter = cptState.ParticleStitchInfo{CC}.persistentLinkIndexCell{f}>CPIndex;
      if any(frameFilter&indFilter)
        cptState.ParticleStitchInfo{CC}.persistentLinkIndexCell{f}(frameFilter&indFilter) = ...
          cptState.ParticleStitchInfo{CC}.persistentLinkIndexCell{f}(frameFilter&indFilter)-1;
      end
    end
    % ...and separations
    for f = 1:length(cptState.ParticleStitchInfo{CC}.forbiddenLinkFrameCell)
      frameFilter = ismember(cptState.ParticleStitchInfo{CC}.forbiddenLinkFrameCell{f},CF);
      indFilter = cptState.ParticleStitchInfo{CC}.forbiddenLinkIndexCell{f}>CPIndex;
      if any(frameFilter&indFilter)
        cptState.ParticleStitchInfo{CC}.forbiddenLinkIndexCell{f}(frameFilter&indFilter) = ...
          cptState.ParticleStitchInfo{CC}.forbiddenLinkIndexCell{f}(frameFilter&indFilter)-1;
      end
    end
    
    %and this part deletes from the spots structure.    
    cptState.Spots{CC}(cptState.CurrentFrame).Fits(CPIndex)= [];
    if isempty(cptState.Spots{CC}(cptState.CurrentFrame).Fits)
        cptState.Spots{CC}(cptState.CurrentFrame).Fits = [];
    end
    
    %now delete from spotfilter
    spotRow = cptState.SpotFilter{CC}(cptState.CurrentFrame,:);
    spotRow(CPIndex) = [];
    spotRow(end+1) = NaN;
    cptState.SpotFilter{CC}(cptState.CurrentFrame,:) = spotRow;
    disp('Spot deleted successfully.');
    
    
end

end

