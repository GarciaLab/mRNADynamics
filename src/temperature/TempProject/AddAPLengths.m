function APLengths = AddAPLengths(ltmp)
    % 1. Define DropboxFolder
    [~,~,DropboxFolder,~,~]=...
        DetermineLocalFolders;
    % 2. Load APDetection and find pixel size for each embryo
    APLengths = []; % in microns
    for i = 1:length(ltmp.ExperimentPrefixes)
        if ltmp.ExperimentStatuses{i}.hasAddedParticlePosition
            load([DropboxFolder,filesep,ltmp.ExperimentPrefixes{i},filesep,'APDetection.mat']);
            % I think AP Length  is in pixels
            APLengths(i) = ltmp.Experiments{i}.pixelSize_um*APLength;
        else
            APLengths(i) = NaN;
        end
    end
end