    function schnitzcells = genSchnitzCellsFromIDs(...
        schnitzcells, tracks, frameIndex)

    % This is a subfunction for tracknuclei_computervision 
        
        
         
        for k = 1:length(tracks)
            if length(schnitzcells) < k
               
                schnitzcells(k).frames(1) = frameIndex;
                schnitzcells(k).cenx(1) =  max(tracks(k).kalmanFilter.State(1), 0);
                schnitzcells(k).ceny(1) =  max(tracks(k).kalmanFilter.State(4), 0);
                schnitzcells(k).smaj(1) =  max(tracks(k).kalmanFilter.State(7), 0);
                schnitzcells(k).smin(1) =  max(tracks(k).kalmanFilter.State(10),0);
                schnitzcells(k).orientationAngle(1) = max(tracks(k).kalmanFilter.State(13), 0); 

                
            else
                
               schnitzcells(k).frames(end+1) = frameIndex;
                schnitzcells(k).cenx(end+1) = max(tracks(k).kalmanFilter.State(1),0);
                schnitzcells(k).ceny(end+1) = max(tracks(k).kalmanFilter.State(1+3), 0);
                schnitzcells(k).smaj(end+1) = max(tracks(k).kalmanFilter.State(1+2*3), 0); 
                schnitzcells(k).smin(end+1) = max(tracks(k).kalmanFilter.State(1+3*3), 0);
                schnitzcells(k).orientationAngle(end+1) = max(tracks(k).kalmanFilter.State(1+4*3), 0);

                
            end
        end
        
        
    end