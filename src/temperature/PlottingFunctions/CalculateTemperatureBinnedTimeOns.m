function [BinnedEmbryoTimeOns,SEBinnedEmbryoTimesOns,BinnedEmbryoCounts] =...
    CalculateTemperatureBinnedTimeOns(ltm)
% Hardcoded in anaphase aligned trace type
Temperatures = fliplr(unique(ltm.Temp_sps));
NumTemperatures = length(Temperatures);
EmbryoTemperatures = ltm.Temp_sps;
NumEmbryos = length(ltm.Temp_sps);
IncludeExperimentTF = ismember(1:NumEmbryos, ltm.IncludedExperiments);

NumAPbins = size(ltm.MeanProfiles{1}.AnaphaseAlignedCycleMeanTraces, 2);
EmbryoIndices = find(IncludeExperimentTF);
BinnedEmbryoTimeOns = NaN(NumTemperatures, NumAPbins, 6);
SEBinnedEmbryoTimesOns = NaN(NumTemperatures, NumAPbins, 6);
BinnedEmbryoCounts = zeros(NumTemperatures, NumAPbins, 6);

TimeOns = ltm.TimeOns.AnaphaseAligned;
SETimeOns = ltm.TimeOns.AnaphaseAlignedStdError;
R2s = ltm.MeanFitR2s.AnaphaseAligned;
R2bound = 0.8;
AllAPbins = 1:NumAPbins;
TimeOnAPbins = AllAPbins >= 5 & AllAPbins <= 25;
for idx = 1:NumTemperatures
    TempEmbryoIndices = ismember(1:NumEmbryos, EmbryoIndices) & round(ltm.Temp_sps, 1) == round(Temperatures(idx),1);
    for APindex = 1:NumAPbins
        for NC = 9:14
            EmbTimeOns = TimeOns(TempEmbryoIndices,APindex,NC-8);
            EmbSETimeOns = SETimeOns(TempEmbryoIndices,APindex,NC-8);
            EmbR2s = R2s(TempEmbryoIndices,APindex,NC-8);
            GoodTimeOns = EmbTimeOns(EmbR2s >= R2bound);
            if length(GoodTimeOns) >= 1
                BinnedEmbryoTimeOns(idx,APindex,NC-8) = mean(GoodTimeOns);
                SEBinnedEmbryoTimesOns(idx,APindex,NC-8) = std(GoodTimeOns);
                BinnedEmbryoCounts(idx, APindex, NC-8) = length(GoodTimeOns);
            end
        end
    end
end


%% First make average profiles binning anaphase aligned traces for each temperature 
BinnedEmbryoTimeOns(BinnedEmbryoTimeOns == 0) = NaN;
SEBinnedEmbryoTimesOns(isnan(SEBinnedEmbryoTimesOns)) = NaN;

