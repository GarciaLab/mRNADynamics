function FittingIndices = GetIndicesForPositiveSlopeFit(SmoothedTrace)
% Fitting Approach one-Determine range for use, t_off, and t_peak
LikelyMaxValues = find(SmoothedTrace >= .8*max(SmoothedTrace));

% First find the indices to use for fitting the loading rate
FluoDeltas = diff(SmoothedTrace).';

PositiveSlopeIndices = find(FluoDeltas > 0);
if isempty(PositiveSlopeIndices)
    FittingIndices = [];
else
    GapsBetweenPositiveSlopeFrames = diff(PositiveSlopeIndices);
    
    
    CandidatePositiveStretches = {};
    FluoGains = [];
    LargerGapsBetweenPositiveSlopes = find(GapsBetweenPositiveSlopeFrames > 1);
    LargerGapsBetweenPositiveSlopes = [0, LargerGapsBetweenPositiveSlopes];
    if LargerGapsBetweenPositiveSlopes(end) < length(PositiveSlopeIndices)
        LargerGapsBetweenPositiveSlopes = [LargerGapsBetweenPositiveSlopes, length(PositiveSlopeIndices)];
    end
    for idx = 1:(length(LargerGapsBetweenPositiveSlopes)-1)
        StartIndex = LargerGapsBetweenPositiveSlopes(idx)+1;
        EndIndex = LargerGapsBetweenPositiveSlopes(idx+1);
        PositiveSlopeIndicesSubset = PositiveSlopeIndices(uint16(StartIndex):uint16(EndIndex));
        PositiveSlopeIndicesSubset = [PositiveSlopeIndicesSubset, uint16(max(PositiveSlopeIndicesSubset)+1)];
        CandidatePositiveStretches{idx} = PositiveSlopeIndicesSubset;
        FluoGains(idx) = SmoothedTrace(PositiveSlopeIndicesSubset(end))-SmoothedTrace(PositiveSlopeIndicesSubset(1));
    end
    
    FittingIndices = CandidatePositiveStretches{find(FluoGains == max(FluoGains), 1)};
%     if isempty(intersect(FittingIndices, LikelyMaxValues))
%         warning('May not be the correct set of time points to fit on loading rate.')
%     end
end