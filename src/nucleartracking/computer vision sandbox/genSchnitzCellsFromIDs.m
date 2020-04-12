    function schnitzcells = genSchnitzCellsFromIDs(...
        schnitzcells, tracks, frameIndex)

    % This is a subfunction for tracknuclei_computervision 

         
        for k = 1:length(tracks)
            if length(schnitzcells) < k
                
                schnitzcells(k).cenx(1) = tracks(k).kalmanFilter.State(1);
                schnitzcells(k).ceny(1) = tracks(k).kalmanFilter.State(4);
                schnitzcells(k).frames(1) = frameIndex;
                schnitzcells(k).len(1) = tracks(k).kalmanFilter.State(7);
                
                
            else
                
                schnitzcells(k).cenx(end+1) = tracks(k).kalmanFilter.State(1);
                schnitzcells(k).ceny(end+1) = tracks(k).kalmanFilter.State(4);
                schnitzcells(k).frames(end+1) = frameIndex;
                schnitzcells(k).len(end+1) = tracks(k).kalmanFilter.State(7);
                
                
            end
        end
        
        
    end