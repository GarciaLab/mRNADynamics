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
        
    end
    
end

end