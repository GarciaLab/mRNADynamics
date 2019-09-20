function [Spots, SpotFilter] = loadSpotsAndCreateSpotFilter(DropboxFolder, Prefix, NCh)

disp('Loading Spots.mat...')
if ~exist('Spots', 'var')
    load([DropboxFolder, filesep, Prefix, filesep, 'Spots.mat'], 'Spots')
    
    % If there's only one channel, Particles, Spots and other structures are
    % not saved as cells. We turn them into a cell to ensure
    % compatibility.
    if NCh == 1
        Spots = {Spots};
    end
    
    MaxSpots = cell(NCh);
    
    for Channel = 1:NCh
        %Determine the maximum number of spots in a given frame for the
        %whole movie
        MaxSpots{Channel} = 0;
        
        for i = 1:length(Spots{Channel})
            MaxSpots{Channel} = max([MaxSpots{Channel}, length(Spots{Channel}(i).Fits)]);
        end
        
        %This filter tells us whether a spot is above the threshold.
        if exist('SpotFilter', 'var')
            
            if ~iscell(SpotFilter)
                SpotFilter = {SpotFilter};
            end
            
        end
        
        SpotFilter{Channel} = nan(length(Spots{Channel}), MaxSpots{Channel});
        % Populate the filter
        for i = 1:length(Spots{Channel})
            
            for j = 1:length(Spots{Channel}(i).Fits)
                
                % Initializes the filter as all 1s, since the Threshold has been removed
                SpotFilter{Channel}(i, j) = 1;
                
            end
            
        end
        
    end
    
end

end