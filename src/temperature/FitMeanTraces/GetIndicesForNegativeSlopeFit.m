function FittingIndices = GetIndicesForNegativeSlopeFit(SmoothedTrace)
LikelyMaxValues = find(SmoothedTrace >= .8*max(SmoothedTrace));

% Next find indices to use for fitting the unloading rate
ReversedSmoothedTrace = flip(SmoothedTrace);
ReversedLikelyMaxValues = find(ReversedSmoothedTrace >= .8*max(ReversedSmoothedTrace));

% First find the indices to use for fitting the loading rate
ReversedFluoDeltas = diff(ReversedSmoothedTrace).';

ReversedPositiveSlopeIndices = find(ReversedFluoDeltas > 0);
if isempty(ReversedPositiveSlopeIndices)
    FittingIndices = [];
else
    ReversedGapsBetweenPositiveSlopeFrames = diff(ReversedPositiveSlopeIndices);
    
    
    ReversedCandidatePositiveStretches = {};
    ReversedFluoGains = [];
    ReversedLargerGapsBetweenPositiveSlopes = find(ReversedGapsBetweenPositiveSlopeFrames > 1);
    ReversedLargerGapsBetweenPositiveSlopes = [0, ReversedLargerGapsBetweenPositiveSlopes];
    if ReversedLargerGapsBetweenPositiveSlopes(end) < length(ReversedPositiveSlopeIndices)
        ReversedLargerGapsBetweenPositiveSlopes = [ReversedLargerGapsBetweenPositiveSlopes, length(ReversedPositiveSlopeIndices)];
    end
    for idx = 1:(length(ReversedLargerGapsBetweenPositiveSlopes)-1)
        StartIndex = ReversedLargerGapsBetweenPositiveSlopes(idx)+1;
        EndIndex = ReversedLargerGapsBetweenPositiveSlopes(idx+1);
        ReversedPositiveSlopeIndicesSubset = ReversedPositiveSlopeIndices(uint16(StartIndex):uint16(EndIndex));
        ReversedPositiveSlopeIndicesSubset = [ReversedPositiveSlopeIndicesSubset, uint16(max(ReversedPositiveSlopeIndicesSubset)+1)];
        ReversedCandidatePositiveStretches{idx} = ReversedPositiveSlopeIndicesSubset;
        ReversedFluoGains(idx) = ReversedSmoothedTrace(ReversedPositiveSlopeIndicesSubset(end))-ReversedSmoothedTrace(ReversedPositiveSlopeIndicesSubset(1));
    end
    
    FittingIndices = sort(length(SmoothedTrace)+1-ReversedCandidatePositiveStretches{find(ReversedFluoGains == max(ReversedFluoGains), 1)});
    
%     if isempty(intersect(FittingIndices, LikelyMaxValues))
%         warning('May not be the correct set of time points to fit on loading rate.')
%     end
end