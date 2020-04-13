    function schnitzcells = genSchnitzCellsFromIDs(...
        schnitzcells, tracks, frameIndex)

    % This is a subfunction for tracknuclei_computervision 

         
        for k = 1:length(tracks)
            if length(schnitzcells) < k
               
                schnitzcells(k).frames(1) = frameIndex;
                schnitzcells(k).cenx(1) = tracks(k).kalmanFilter.State(1);
                schnitzcells(k).ceny(1) = tracks(k).kalmanFilter.State(4);
                schnitzcells(k).smaj(1) = tracks(k).kalmanFilter.State(7);
                schnitzcells(k).smin(1) = tracks(k).kalmanFilter.State(10);
                schnitzcells(k).orientationAngle(1) = tracks(k).kalmanFilter.State(13);


                
            else
                
               schnitzcells(k).frames(end+1) = frameIndex;
                schnitzcells(k).cenx(end+1) = tracks(k).kalmanFilter.State(1);
                schnitzcells(k).ceny(end+1) = tracks(k).kalmanFilter.State(4);
                schnitzcells(k).smaj(end+1) = tracks(k).kalmanFilter.State(7);
                schnitzcells(k).smin(end+1) = tracks(k).kalmanFilter.State(10);
                schnitzcells(k).orientationAngle(end+1) = tracks(k).kalmanFilter.State(13);

                
            end
        end
        
        
    end