function Spots = segmentNanocagesZTracking(pixelSize, Spots)
% author: G. Martini
% date created: 4/19/21
%%
numFrames = length(Spots);
neighborhood = round(200 / pixelSize); %empirically, the spot should be within ~1.3um of the same spot in other z-planes
waitbarFigure = waitbar(0, 'Finding z-columns');

if ~isempty([Spots.Fits])
    fields = fieldnames([Spots.Fits]);
    for currentFrame = 1:numFrames
        waitbar(currentFrame / numFrames, waitbarFigure)
        for spot_index = 1:length(Spots(currentFrame).Fits)
            Spots(currentFrame).Fits(spot_index).LastZ = Spots(currentFrame).Fits(spot_index).z(end);
            Spots(currentFrame).Fits(spot_index).FirstZ = Spots(currentFrame).Fits(spot_index).z(1);
            Spots(currentFrame).Fits(spot_index).zCount = length(Spots(currentFrame).Fits(spot_index).z);
            Spots(currentFrame).Fits(spot_index).nearest_neighbor = NaN;
            Spots(currentFrame).Fits(spot_index).PixelList = {};
            Spots(currentFrame).Fits(spot_index).PixelList{1} = Spots(currentFrame).Fits(spot_index).bwPixelList;
            Spots(currentFrame).Fits(spot_index).zDistances = [];
        end
        if ~isempty(Spots(currentFrame).Fits)
            all_zs = unique([Spots(currentFrame).Fits(:).LastZ]);
            for j = 1:length(all_zs)
                CurrentZ = all_zs(j);
                NextZ = CurrentZ + 1;
                changes = 1;
                
                
                CurrentZSpots = find([Spots(currentFrame).Fits(:).LastZ] == CurrentZ);
                x0 = NaN(1, length(CurrentZSpots));
                y0 = NaN(1, length(CurrentZSpots));
                for spotIndex = 1:length(CurrentZSpots)
                    x0(spotIndex) = Spots(currentFrame).Fits(CurrentZSpots(spotIndex)).xDoG(end);
                    y0(spotIndex) = Spots(currentFrame).Fits(CurrentZSpots(spotIndex)).yDoG(end);
                end
            
                PointLocations = [x0.' y0.'];
                Distances = squareform(pdist(PointLocations));
                Distances(1:length(CurrentZSpots)+1:end) = max(max(Distances));
                minDistances = min(Distances);
                for spotIndex = 1:length(CurrentZSpots)
                    Spots(currentFrame).Fits(CurrentZSpots(spotIndex)).nearest_neighbor(end) = minDistances(spotIndex); % pixels
                end
                
                NextZSpots = find([Spots(currentFrame).Fits(:).FirstZ] == NextZ);
                if isempty(NextZSpots)
                    continue
                end
                AllSpots = [CurrentZSpots NextZSpots];
                x0 = NaN(1, length(AllSpots));
                y0 = NaN(1, length(AllSpots));
                
                for spotIndex2 = 1:length(AllSpots)
                    x0(spotIndex2) = Spots(currentFrame).Fits(AllSpots(spotIndex2)).xDoG(end);
                    y0(spotIndex2) = Spots(currentFrame).Fits(AllSpots(spotIndex2)).yDoG(end);
                end

                PointLocations = [x0.' y0.'];
                Distances = squareform(pdist(PointLocations));
                MaxDistance = max(max(Distances));
                Distances = Distances(length(CurrentZSpots)+1:length(CurrentZSpots)+length(NextZSpots), 1:length(CurrentZSpots));
                [MinDistances, MinIndices] = min(Distances);
                
                MinMinDistance = min(MinDistances);
                SpotsToMerge = zeros(0,2,'uint16');
                while MinMinDistance <= neighborhood
                    [sortedDistances, sortedSpotPerm] = sort(MinDistances);
                    [MinDistance, MinIndex] = min(Distances(:, sortedSpotPerm(1)));
                    SpotsToMerge = [SpotsToMerge; CurrentZSpots(sortedSpotPerm(1)) NextZSpots(MinIndex)];
                    Distances(:,sortedSpotPerm(1)) = MaxDistance;
                    Distances(MinIndex,:) = MaxDistance;
                    
                    [MinDistances, MinIndices] = min(Distances);
                    
                    MinMinDistance = min(MinDistances);
                end
%                 for rowIndex = 1:size(SpotsToMerge, 1)
%                     disp([num2str(Spots(currentFrame).Fits(SpotsToMerge(rowIndex, 1)).xDoG),', ',num2str(Spots(currentFrame).Fits(SpotsToMerge(rowIndex, 2)).xDoG),', ',...
%                         num2str(Spots(currentFrame).Fits(SpotsToMerge(rowIndex, 1)).yDoG),', ',num2str(Spots(currentFrame).Fits(SpotsToMerge(rowIndex, 2)).yDoG)]);
%                 end
                Spots = mergeNanocageSlices(Spots, currentFrame, SpotsToMerge);
            end
      
        end
        
    end
    
else
    fields = {};
end

close(waitbarFigure)

end
