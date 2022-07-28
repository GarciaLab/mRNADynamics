function schnitzcells = filter_nuclear_traces(Prefix)

min_frames = 15; % 5 min for my datasets, empirical

liveExperiment = LiveExperiment(Prefix);
        
dropbox_folder = liveExperiment.userResultsFolder;

schnitzcells_filename = [dropbox_folder,filesep,Prefix,filesep,Prefix,'_lin.mat']; 
schnitzcells = getSchnitzcells(liveExperiment);

 

% Automated disapproval for common trace quality issues
for sc = 1:numel(schnitzcells)
    n_frames = numel(schnitzcells(sc).frames);
    
    % QC based on overall trace length
    if n_frames < min_frames
        schnitzcells(sc).Approved = 0;
    end
    
    % QC based on the nuclear fluorescence of each trace
    if isfield(schnitzcells, 'Fluo')
        nuclear_fluo = schnitzcells(sc).Fluo;
        % Disapprove if all frames are missing a Fluo measurement
        frames_with_fluo = sum(~isnan(nuclear_fluo(:,1)));
        if frames_with_fluo == 0
            schnitzcells(sc).Approved = 0;
        % Disapprove if over half the frames are missing a Fluo measurement
        elseif frames_with_fluo < n_frames/2
            schnitzcells(sc).Approved = 0;
        end
    end

end

% Save changes
save2(schnitzcells_filename, schnitzcells);
