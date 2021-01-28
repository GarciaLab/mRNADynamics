function [Ymax, Ymin] = getEmbryoStatsGlobalYlims(MeanParams, SEParams)
allowed_idx = ~isnan(MeanParams);
SEParams(isnan(SEParams)) = 0;
Ymax = max(max(max(MeanParams(allowed_idx) + SEParams(allowed_idx))));
Ymin = min(min(min(MeanParams(allowed_idx) - SEParams(allowed_idx))));
end