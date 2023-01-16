function [BinnedEmbryoTotalmRNA,SEBinnedEmbryoTotalmRNA,BinnedEmbryoCounts,...
    BinnedEmbryoTotalmRNAV2, SEBinnedEmbryoTotalmRNAV2,BinnedEmbryoCountsV2,...
    TempElongationTimes, SETempElongationTimes] = ...
    CalculateTemperatureBinnedIntegratedFluomRNA(ltm, IncludeFractionOn, UseRescaledFluo)
% Hardcoded in anaphase aligned trace type
if ~exist('IncludeFractionOn', 'var')
    IncludeFractionOn = false;
end
if ~exist('UseRescaledFluo', 'var')
    UseRescaledFluo = false;
end

if UseRescaledFluo 
    FluoCoeffs = ltm.FluoCoeffs;
else
    FluoCoeffs = ones(1,5);
end

Temperatures = fliplr(unique(ltm.Temp_sps));
NumTemperatures = length(Temperatures);
EmbryoTemperatures = ltm.Temp_sps;
NumEmbryos = length(ltm.Temp_sps);
IncludeExperimentTF = ismember(1:NumEmbryos, ltm.IncludedExperiments);


NumAPbins = size(ltm.MeanProfiles{1}.AnaphaseAlignedCycleMeanTraces, 2);
EmbryoIndices = find(IncludeExperimentTF);
IntegratedFluos = NaN(NumEmbryos, NumAPbins, 6); % AU*min
EmbryoTotalmRNA = NaN(NumEmbryos, NumAPbins, 6); 
BinnedEmbryoTotalmRNA = NaN(NumTemperatures, NumAPbins, 6);
SEBinnedEmbryoTotalmRNA = NaN(NumTemperatures, NumAPbins, 6);
BinnedEmbryoCounts = NaN(NumTemperatures, NumAPbins, 6);
BinnedEmbryoTotalmRNAV2 = NaN(NumTemperatures, NumAPbins, 6); 
SEBinnedEmbryoTotalmRNAV2 = NaN(NumTemperatures, NumAPbins, 6); 
BinnedEmbryoCountsV2 = NaN(NumTemperatures, NumAPbins, 6);
EmbryoElongTimes = NaN(NumEmbryos,6);
EmbryoSEElongTimes = NaN(NumEmbryos, 6);

ElongationTimes = ltm.ElongationTimes.AnaphaseAligned;
SEElongationTimes = ltm.ElongationTimes.AnaphaseAlignedStdError;
R2s = ltm.MeanFitR2s.AnaphaseAligned;
R2bound = 0.95;
AllAPbins = 1:NumAPbins;
ElongAPbins = AllAPbins >= 9 & AllAPbins <= 17;
for EmbryoIndex = 1:NumEmbryos
    if ismember(EmbryoIndex, EmbryoIndices)
        for NC = 9:14
            EmbElongTimes = ElongationTimes(EmbryoIndex,ElongAPbins,NC-8);
            EmbSEElongTimes = SEElongationTimes(EmbryoIndex,ElongAPbins,NC-8);
            EmbR2s = R2s(EmbryoIndex,ElongAPbins,NC-8);
            GoodElongTimes = EmbElongTimes(EmbR2s >= R2bound);
            if length(GoodElongTimes) >= 1
                EmbryoElongTimes(EmbryoIndex,NC-8) = mean(GoodElongTimes);
                EmbryoSEElongTimes(EmbryoIndex,NC-8) = std(GoodElongTimes);
            end
            
            nc_index = NC-8;
            FrameTimes = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleFrameTimes{nc_index}/60;
            AllAPProfiles = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleMeanTraces(:,:,nc_index);
            if IncludeFractionOn
                NumOnNuclei = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleNumOnNuclei(:,:,nc_index);
                NumNuclei = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleNumNuclei(:,:,nc_index);
                AllFractionOns = NumOnNuclei./NumNuclei;
                
            end
            for APindex = 1:NumAPbins
                FullAPProfile = AllAPProfiles(:,APindex).';
                if IncludeFractionOn
                    FullFractionOnProfile = AllFractionOns(:,APindex).';
         
                    FullAPProfile = FullAPProfile.*FullFractionOnProfile;
                end
                UsableTimes = ~isnan(FullAPProfile);
               
                APFrameTimes =  FrameTimes(UsableTimes);
                APProfile = FullAPProfile(UsableTimes);
                if sum(UsableTimes) == 0
                    continue
                end
                
                TimeBinAveragedProfile = (APProfile(1:end-1)+APProfile(2:end))/2;
                TimeBinWidths = diff(APFrameTimes);
                IntegratedFluos(EmbryoIndex, APindex, NC-8) = sum(TimeBinAveragedProfile.*TimeBinWidths);
                EmbryoTotalmRNA(EmbryoIndex, APindex, NC-8) = IntegratedFluos(EmbryoIndex, APindex, NC-8)/EmbryoElongTimes(EmbryoIndex,NC-8);
            end
        end
    end
end

TempElongationTimes = NaN(NumTemperatures, 6);
SETempElongationTimes = NaN(NumTemperatures, 6);
for idx = 1:NumTemperatures
    TempEmbryoIndices = ismember(1:NumEmbryos, EmbryoIndices) & round(ltm.Temp_sps, 1) == round(Temperatures(idx),1);
    
    for NC = 9:14
        TElongTime = mean(EmbryoElongTimes(TempEmbryoIndices,NC-8), 'omitnan');
        SETElongTime = std(EmbryoElongTimes(TempEmbryoIndices,NC-8), 'omitnan');
        TempElongationTimes(idx, NC-8) = TElongTime;
        SETempElongationTimes(idx, NC-8) = SETElongTime;
        for APindex = 1:NumAPbins
            BinnedEmbryoTotalmRNA(idx,APindex,NC-8) = FluoCoeffs(idx)*mean(EmbryoTotalmRNA(TempEmbryoIndices, APindex, NC-8), 'omitnan');
            SEBinnedEmbryoTotalmRNA(idx,APindex,NC-8) = FluoCoeffs(idx)*std(EmbryoTotalmRNA(TempEmbryoIndices, APindex, NC-8), 'omitnan');
            BinnedEmbryoCounts(idx, APindex, NC-8) = sum(~isnan(EmbryoTotalmRNA(TempEmbryoIndices, APindex, NC-8)));
            BinnedEmbryoTotalmRNAV2(idx,APindex,NC-8) = FluoCoeffs(idx)*mean(IntegratedFluos(TempEmbryoIndices, APindex, NC-8), 'omitnan')/TElongTime;
            SEBinnedEmbryoTotalmRNAV2(idx,APindex,NC-8) = FluoCoeffs(idx)*sqrt(std(IntegratedFluos(TempEmbryoIndices, APindex, NC-8), 'omitnan')^2/TElongTime^2+...
                SETElongTime^2*mean(IntegratedFluos(TempEmbryoIndices, APindex, NC-8), 'omitnan')^2/TElongTime^4);
            BinnedEmbryoCountsV2(idx, APindex, NC-8) = sum(~isnan(IntegratedFluos(TempEmbryoIndices, APindex, NC-8)/TElongTime));
        end
    end
end
%% First make average profiles binning anaphase aligned traces for each temperature 
BinnedEmbryoTotalmRNA(BinnedEmbryoTotalmRNA == 0) = NaN;
SEBinnedEmbryoTotalmRNA(isnan(BinnedEmbryoTotalmRNA)) = NaN;
BinnedEmbryoTotalmRNAV2(BinnedEmbryoTotalmRNAV2 == 0) = NaN;
SEBinnedEmbryoTotalmRNAV2(isnan(BinnedEmbryoTotalmRNAV2)) = NaN;
TempElongationTimes(TempElongationTimes == 0) = NaN;
SETempElongationTimes(isnan(TempElongationTimes)) = NaN;
