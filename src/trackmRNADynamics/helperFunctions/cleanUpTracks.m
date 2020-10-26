function particleTracks = cleanUpTracks(particleTracks, maxSpotsPerNucleus, lastFrameFlag)
    
    % get list of nuclei assigned to each track
    trackNuclei = [particleTracks.Nucleus];
    [ncIDs, ~, ncIDx] = unique(trackNuclei(:));
    ncCounts = accumarray(ncIDx, 1);        
    
    % flag cases with too many assignments
    degenerateIndices = ncIDs(ncCounts>maxSpotsPerNucleus);    
    for i = 1:length(degenerateIndices)
        ind = degenerateIndices(i);
        trackIndices = find(trackNuclei==ind);
        % figure out how many we need to cut
        nExtraTracks = length(trackIndices) - maxSpotsPerNucleus;
        
        % rank tracks by length and keep longest. If this is not 
        % conclusive (if both tracks have length 1) then wait
        lengthVec = NaN(size(trackIndices));        
        
        for j = 1:length(trackIndices)
            lengthVec(j) = length(particleTracks(trackIndices(j)).Frame);
        end                
        
        % flag tracks for removal
        if sum(lengthVec==1) == nExtraTracks || lastFrameFlag
            [~, lengthOrder] = sort(lengthVec, 'ascend');
            flagIDs = trackIndices(lengthOrder(1:nExtraTracks));
        else
            flagIDs = [];
        end
        for f = 1:length(flagIDs)
            particleTracks(flagIDs(f)).duplicateFlag = true;
        end
    end  
    