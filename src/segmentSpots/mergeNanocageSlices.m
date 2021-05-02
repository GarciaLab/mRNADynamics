function Spots = mergeNanocageSlices(Spots, currentFrame, SpotsToMerge)

fields = fieldnames([Spots(currentFrame).Fits]);
for i = 1:size(SpotsToMerge, 1)
    spotIndex = SpotsToMerge(i, 1);
    spotIndex2 = SpotsToMerge(i, 2);
    if ismember('FixedAreaIntensity', fields)
        Spots(currentFrame).Fits(spotIndex).FixedAreaIntensity = ...
            [Spots(currentFrame).Fits(spotIndex).FixedAreaIntensity Spots(currentFrame).Fits(spotIndex2).FixedAreaIntensity];
    end
    if ismember('xFit', fields)
        Spots(currentFrame).Fits(spotIndex).xFit = ...
            [Spots(currentFrame).Fits(spotIndex).xFit Spots(currentFrame).Fits(spotIndex2).xFit];
        Spots(currentFrame).Fits(spotIndex).yFit = ...
            [Spots(currentFrame).Fits(spotIndex).yFit Spots(currentFrame).Fits(spotIndex2).yFit];
        try
            Spots(currentFrame).Fits(spotIndex).zDistances(end+1) =...
                sqrt((Spots(currentFrame).Fits(spotIndex).xFit(end)-Spots(currentFrame).Fits(spotIndex).xFit(end-1))^2+...
                (Spots(currentFrame).Fits(spotIndex).yFit(end)-Spots(currentFrame).Fits(spotIndex).yFit(end-1))^2);
        end
    end
    if ismember('Offset', fields)
        Spots(currentFrame).Fits(spotIndex).Offset = ...
            [Spots(currentFrame).Fits(spotIndex).Offset Spots(currentFrame).Fits(spotIndex2).Offset];
    end
    if ismember('xDoG', fields)
        Spots(currentFrame).Fits(spotIndex).xDoG = ...
            [Spots(currentFrame).Fits(spotIndex).xDoG Spots(currentFrame).Fits(spotIndex2).xDoG];
        Spots(currentFrame).Fits(spotIndex).yDoG = ...
            [Spots(currentFrame).Fits(spotIndex).yDoG Spots(currentFrame).Fits(spotIndex2).yDoG];
    end
    if ismember('xFitWidth', fields)
        Spots(currentFrame).Fits(spotIndex).xFitWidth = ...
            [Spots(currentFrame).Fits(spotIndex).xFitWidth Spots(currentFrame).Fits(spotIndex2).xFitWidth];
        Spots(currentFrame).Fits(spotIndex).yFitWidth = ...
            [Spots(currentFrame).Fits(spotIndex).yFitWidth Spots(currentFrame).Fits(spotIndex2).yFitWidth];
    end
    if ismember('GaussianIntensity', fields)
        Spots(currentFrame).Fits(spotIndex).GaussianIntensity = ...
            [Spots(currentFrame).Fits(spotIndex).GaussianIntensity Spots(currentFrame).Fits(spotIndex2).GaussianIntensity];
    end
    if ismember('GaussianInfo', fields)
        Spots(currentFrame).Fits(spotIndex).GaussianInfo = ...
            [Spots(currentFrame).Fits(spotIndex).GaussianInfo Spots(currentFrame).Fits(spotIndex2).GaussianInfo];
        Spots(currentFrame).Fits(spotIndex).GaussianError = ...
            [Spots(currentFrame).Fits(spotIndex).GaussianError Spots(currentFrame).Fits(spotIndex2).GaussianError];
    end
    if ismember('CentralIntensity', fields)
        Spots(currentFrame).Fits(spotIndex).CentralIntensity = ...
            [Spots(currentFrame).Fits(spotIndex).CentralIntensity Spots(currentFrame).Fits(spotIndex2).CentralIntensity];
    end
    if ismember('DOGIntensity', fields)
        Spots(currentFrame).Fits(spotIndex).DOGIntensity = ...
            [Spots(currentFrame).Fits(spotIndex).DOGIntensity Spots(currentFrame).Fits(spotIndex2).DOGIntensity];
    end
    if ismember('GaussianResiduals', fields)
        Spots(currentFrame).Fits(spotIndex).GaussianResiduals = ...
            [Spots(currentFrame).Fits(spotIndex).GaussianResiduals Spots(currentFrame).Fits(spotIndex2).GaussianResiduals];
    end
    if ismember('GaussianFitValues', fields)
        Spots(currentFrame).Fits(spotIndex).GaussianFitValues = ...
            [Spots(currentFrame).Fits(spotIndex).GaussianFitValues Spots(currentFrame).Fits(spotIndex2).GaussianFitValues];
    end
    if ismember('ConfidenceIntervals', fields)
        Spots(currentFrame).Fits(spotIndex).ConfidenceIntervals = ...
            [Spots(currentFrame).Fits(spotIndex).ConfidenceIntervals Spots(currentFrame).Fits(spotIndex2).ConfidenceIntervals];
    end
    if ismember('gaussParams', fields)
        Spots(currentFrame).Fits(spotIndex).gaussParams = ...
            [Spots(currentFrame).Fits(spotIndex).gaussParams Spots(currentFrame).Fits(spotIndex2).gaussParams];
    end
    if ismember('dogFixedAreaIntensity', fields)
        Spots(currentFrame).Fits(spotIndex).dogFixedAreaIntensity = ...
            [Spots(currentFrame).Fits(spotIndex).dogFixedAreaIntensity Spots(currentFrame).Fits(spotIndex2).dogFixedAreaIntensity];
    end
    if ismember('z', fields)
        Spots(currentFrame).Fits(spotIndex).z = ...
            [Spots(currentFrame).Fits(spotIndex).z Spots(currentFrame).Fits(spotIndex2).z];
    end
    if ismember('bwIntensity', fields)
        Spots(currentFrame).Fits(spotIndex).bwIntensity = ...
            [Spots(currentFrame).Fits(spotIndex).bwIntensity Spots(currentFrame).Fits(spotIndex2).bwIntensity];
    end
    if ismember('bwArea', fields)
        Spots(currentFrame).Fits(spotIndex).bwArea = ...
            [Spots(currentFrame).Fits(spotIndex).bwArea Spots(currentFrame).Fits(spotIndex2).bwArea];
    end
    if ismember('bwDiameter', fields)
        Spots(currentFrame).Fits(spotIndex).bwDiameter = ...
            [Spots(currentFrame).Fits(spotIndex).bwDiameter Spots(currentFrame).Fits(spotIndex2).bwDiameter];
    end
    if ismember('bwCircularity', fields)
        Spots(currentFrame).Fits(spotIndex).bwCircularity = ...
            [Spots(currentFrame).Fits(spotIndex).bwCircularity Spots(currentFrame).Fits(spotIndex2).bwCircularity];
    end
    if ismember('bwEccentricity', fields)
        Spots(currentFrame).Fits(spotIndex).bwEccentricity = ...
            [Spots(currentFrame).Fits(spotIndex).bwEccentricity Spots(currentFrame).Fits(spotIndex2).bwEccentricity];
    end
    if ismember('bwMajorAxisLength', fields)
        Spots(currentFrame).Fits(spotIndex).bwMajorAxisLength = ...
            [Spots(currentFrame).Fits(spotIndex).bwMajorAxisLength Spots(currentFrame).Fits(spotIndex2).bwMajorAxisLength];
    end
    if ismember('bwMinorAxisLength', fields)
        Spots(currentFrame).Fits(spotIndex).bwMinorAxisLength = ...
            [Spots(currentFrame).Fits(spotIndex).bwMinorAxisLength Spots(currentFrame).Fits(spotIndex2).bwMinorAxisLength];
    end
    if ismember('bwPixelList', fields)
        if ~ismember('PixelList', fields)
            Spots(currentFrame).Fits(spotIndex).PixelList = {};
            Spots(currentFrame).Fits(spotIndex).PixelList{1} = Spots(currentFrame).Fits(spotIndex).bwPixelList;
        end
        
        Spots(currentFrame).Fits(spotIndex).PixelList = ...
            [Spots(currentFrame).Fits(spotIndex).PixelList Spots(currentFrame).Fits(spotIndex2).PixelList];
    end
    if ismember('discardThis', fields)
        Spots(currentFrame).Fits(spotIndex).discardThis = ...
            [Spots(currentFrame).Fits(spotIndex).discardThis Spots(currentFrame).Fits(spotIndex2).discardThis];
    end
    if ismember('discardThis', fields)
        Spots(currentFrame).Fits(spotIndex).discardThis = ...
            [Spots(currentFrame).Fits(spotIndex).discardThis Spots(currentFrame).Fits(spotIndex2).discardThis];
    end
    if ismember('r', fields)
        Spots(currentFrame).Fits(spotIndex).r = 1;
    end
    if ismember('LastZ', fields)
        Spots(currentFrame).Fits(spotIndex).LastZ = ...
            Spots(currentFrame).Fits(spotIndex).z(end);
    end
    if ismember('FirstZ', fields)
        Spots(currentFrame).Fits(spotIndex).FirstZ = ...
            Spots(currentFrame).Fits(spotIndex).z(1);
    end
    if ismember('zCount', fields)
        Spots(currentFrame).Fits(spotIndex).zCount = ...
            length(Spots(currentFrame).Fits(spotIndex).z);
    end
    if ismember('nearest_neighbor', fields)
        Spots(currentFrame).Fits(spotIndex).nearest_neighbor = ...
            [Spots(currentFrame).Fits(spotIndex).nearest_neighbor Spots(currentFrame).Fits(spotIndex2).nearest_neighbor];
    end
end

MatchedNextZSpots = SpotsToMerge(:,2).';
AllSpotIndices = 1:length(Spots(currentFrame).Fits);
KeptSpots = ~ismember(AllSpotIndices, MatchedNextZSpots);
Spots(currentFrame).Fits = Spots(currentFrame).Fits(KeptSpots);
