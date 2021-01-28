function [Ymax, Ymin, Ymax2, Ymin2] = getGlobalYlims(MeanParams, SEParams, R2s, R2bound)
if ~all(all(all(isnan(SEParams))))
    allowed_idx = (R2s >= R2bound) & (MeanParams./SEParams >= 1) & ~isnan(MeanParams);
else
    allowed_idx = (R2s >= R2bound)  & ~isnan(MeanParams);
end
Ymax = max(max(max(MeanParams(allowed_idx) + SEParams(allowed_idx))));
Ymin = min(min(min(MeanParams(allowed_idx) - SEParams(allowed_idx))));

MaxMean = max(max(max(MeanParams(allowed_idx))));
MinMean = min(min(min(MeanParams(allowed_idx))));



Ymax2 = min([MaxMean + (MaxMean-MinMean)*.5, Ymax]);
Ymin2 = max([MinMean-(MaxMean-MinMean)*.5, Ymin]);

end