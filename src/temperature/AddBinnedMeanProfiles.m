function ltm = AddBinnedMeanProfiles(ltm, deltaTbinWidth, MinEmbryos)
%%
if ~exist('deltaTbinWidth', 'var')
    deltaTbinWidth = ltm.time_delta;
end
if ~exist('MinEmbryos', 'var')
    MinEmbryos = ltm.MinimumEmbryos;
end


Temperatures = fliplr(unique(ltm.Temp_sps));
NumTemperatures = length(Temperatures);
EmbryoTemperatures = ltm.Temp_sps;
NumEmbryos = length(ltm.Temp_sps);
IncludeExperimentTF = ismember(1:NumEmbryos, ltm.IncludedExperiments);

LengthAnaphaseAlignedVector = [];
LengthTbinnedVector = [];
for exp_index = ltm.IncludedExperiments
    LengthAnaphaseAlignedVector(end+1) = size(ltm.MeanProfiles{exp_index}.AnaphaseAlignedCycleMeanTraces, 1);
    LengthTbinnedVector(end+1) =size(ltm.MeanProfiles{exp_index}.TbinnedCycleMeanTraces, 1);
end

LengthBinnedAnaphaseAlignedProfile = max(LengthAnaphaseAlignedVector);
LengthBinnedTbinnedProfile = max(LengthTbinnedVector);
NumAPbins = size(ltm.MeanProfiles{1}.AnaphaseAlignedCycleMeanTraces, 2);
%%
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleTimes = (1:LengthBinnedAnaphaseAlignedProfile)*deltaTbinWidth-deltaTbinWidth;
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleNumEmbryos = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleFractionOn = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleFractionOnStdErr = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleMeanTraces = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleStdErrors = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleNuclearMeanTraces = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleNuclearStdErrors = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleTotalNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleTotalOnNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleTotalOffNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleTotalQuiescentNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleTotalFinishedTranscribingNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleGradedTotalNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 10,6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleGradedTotalOnNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 10,6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleGradedTotalOffNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 10,6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedGradedCycleFractionOn = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 10,6, 5);

ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleTimes = (1:LengthBinnedAnaphaseAlignedProfile)*deltaTbinWidth-deltaTbinWidth;
ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleNumEmbryos = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleFractionOn = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleMeanTraces = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleStdErrors = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleNuclearMeanTraces = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleNuclearStdErrors = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleTotalNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleTotalOnNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleTotalOffNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleTotalQuiescentNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleTotalFinishedTranscribingNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleGradedTotalNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 10,6, 5);
ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleGradedTotalOnNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 10,6,  5);
ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleGradedTotalOffNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 10,6, 5);
ltm.BinnedMeanProfiles.AnaphaseAligned3DGradedCycleFractionOn = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 10,6, 5);



ltm.BinnedMeanProfiles.TbinnedCycleTimes = (1:LengthBinnedTbinnedProfile)*deltaTbinWidth-deltaTbinWidth;
ltm.BinnedMeanProfiles.TbinnedCycleNumEmbryos = NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.TbinnedCycleFractionOn = NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.TbinnedCycleMeanTraces = NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.TbinnedCycleStdErrors = NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.TbinnedCycleNuclearMeanTraces = NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.TbinnedCycleNuclearStdErrors = NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.TbinnedCycleTotalNuclei= NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.TbinnedCycleTotalOnNuclei= NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.TbinnedCycleTotalOffNuclei= NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.TbinnedCycleTotalQuiescentNuclei= NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.TbinnedCycleTotalFinishedTranscribingNuclei= NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.TbinnedCycleGradedTotalNuclei= NaN(LengthBinnedTbinnedProfile, NumAPbins, 10, 6,5);
ltm.BinnedMeanProfiles.TbinnedCycleGradedTotalOnNuclei= NaN(LengthBinnedTbinnedProfile, NumAPbins, 10, 6,5);
ltm.BinnedMeanProfiles.TbinnedCycleGradedTotalOffNuclei= NaN(LengthBinnedTbinnedProfile, NumAPbins, 10, 6,5);
ltm.BinnedMeanProfiles.TbinnedGradedCycleFractionOn = NaN(LengthBinnedTbinnedProfile, NumAPbins, 10, 6,5);
ltm.BinnedMeanProfiles.TbinnedCycleFractionOnStdErr = NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);

ltm.BinnedMeanProfiles.Tbinned3DCycleTimes = (1:LengthBinnedTbinnedProfile)*deltaTbinWidth-deltaTbinWidth;
ltm.BinnedMeanProfiles.Tbinned3DCycleNumEmbryos = NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.Tbinned3DCycleFractionOn = NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.Tbinned3DCycleMeanTraces = NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.Tbinned3DCycleStdErrors = NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.Tbinned3DCycleNuclearMeanTraces = NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.Tbinned3DCycleNuclearStdErrors = NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.Tbinned3DCycleTotalNuclei= NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.Tbinned3DCycleTotalOnNuclei= NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.Tbinned3DCycleTotalOffNuclei= NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.Tbinned3DCycleTotalQuiescentNuclei= NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.Tbinned3DCycleTotalFinishedTranscribingNuclei= NaN(LengthBinnedTbinnedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.Tbinned3DCycleGradedTotalNuclei= NaN(LengthBinnedTbinnedProfile, NumAPbins, 10,6, 5);
ltm.BinnedMeanProfiles.Tbinned3DCycleGradedTotalOnNuclei= NaN(LengthBinnedTbinnedProfile, NumAPbins,10,6,5);
ltm.BinnedMeanProfiles.Tbinned3DCycleGradedTotalOffNuclei= NaN(LengthBinnedTbinnedProfile, NumAPbins, 10, 6,5);
ltm.BinnedMeanProfiles.Tbinned3DGradedCycleFractionOn = NaN(LengthBinnedTbinnedProfile, NumAPbins, 10, 6,5);

%% First make average profiles binning anaphase aligned traces for each temperature 
for t_index = 1:NumTemperatures
    EmbryoIndices = find(EmbryoTemperatures == Temperatures(t_index) & IncludeExperimentTF);
    for nc_index = 1:6
        EmbryoProfiles = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfStds =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfsNumNuclei =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfsNumOnNuclei =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfsNumOffNuclei =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfsGradedNumNuclei =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins,10, length(EmbryoIndices));
        EmbryoProfsGradedNumOnNuclei =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 10,length(EmbryoIndices));
        EmbryoProfsGradedNumOffNuclei =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins,10, length(EmbryoIndices));
        EmbryoProfsNumQuiescentNuclei =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfsNumFinishedTranscribingNuclei =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoNuclearProfiles =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoFractionOns =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoGradedFractionOns =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins,10, length(EmbryoIndices));
        EmbryoTraceCounts =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, length(EmbryoIndices));
        for SetIndex = 1:length(EmbryoIndices)
            EmbryoIndex = EmbryoIndices(SetIndex);
            MatchingTimesForBinned = ismember(ltm.BinnedMeanProfiles.AnaphaseAlignedCycleTimes, ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleFrameTimes{nc_index});
            MatchingTimesForSet = ismember(ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleFrameTimes{nc_index}, ltm.BinnedMeanProfiles.AnaphaseAlignedCycleTimes);
            
            if isempty(MatchingTimesForBinned) | sum(MatchingTimesForBinned) == 0
                continue
            end
            EmbryoNCProfile = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins);
            EmbryoNCProfile(MatchingTimesForBinned,:) = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleMeanTraces(MatchingTimesForSet,:,nc_index);
            
            EmbryoNCNuclearProfile = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins);
            EmbryoNCNuclearProfile(MatchingTimesForBinned,:) = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleMeanPerNucleusTraces(MatchingTimesForSet,:,nc_index);
            
            EmbryoNCProfileStd = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins);
            EmbryoNCProfileStd(MatchingTimesForBinned,:) = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleTraceStdErrors(MatchingTimesForSet,:,nc_index);
            
            EmbryoNCProfileFractionOn = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins);
            EmbryoNCProfileFractionOn(MatchingTimesForBinned,:) = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleFractionOn(MatchingTimesForSet,:,nc_index);
            
            EmbryoNCProfileTraceCount = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins);
            EmbryoNCProfileTraceCount(MatchingTimesForBinned,:) = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleTraceCount(MatchingTimesForSet,:,nc_index);
%             EmbryoNCProfileTraceCount(isnan(EmbryoNCProfileTraceCount)) = 0;
            
            EmbryoNCProfileNumNuclei = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins);
            EmbryoNCProfileNumNuclei(MatchingTimesForBinned,:)  = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleNumNuclei(MatchingTimesForSet,:,nc_index);
            EmbryoNCProfileNumNuclei(isnan(EmbryoNCProfileNumNuclei)) = 0;
            EmbryoNCProfileNumNuclei(isnan(EmbryoNCProfileTraceCount)) = NaN;
            
            EmbryoNCProfileNumOnNuclei = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins);
            EmbryoNCProfileNumOnNuclei(MatchingTimesForBinned,:)  = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleNumOnNuclei(MatchingTimesForSet,:,nc_index);
            EmbryoNCProfileNumOnNuclei(isnan(EmbryoNCProfileNumOnNuclei)) = 0;
            EmbryoNCProfileNumOnNuclei(isnan(EmbryoNCProfileTraceCount)) = NaN;
            
            EmbryoNCProfileNumOffNuclei = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins);
            EmbryoNCProfileNumOffNuclei(MatchingTimesForBinned,:) = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleNumOffNuclei(MatchingTimesForSet,:,nc_index);
            EmbryoNCProfileNumOffNuclei(isnan(EmbryoNCProfileNumOffNuclei)) = 0;
            EmbryoNCProfileNumOffNuclei(isnan(EmbryoNCProfileTraceCount)) = NaN;
            
            
            EmbryoNCProfileNumQuiescentNuclei = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins);
            EmbryoNCProfileNumQuiescentNuclei(MatchingTimesForBinned,:)  = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleNumQuiescentNuclei(MatchingTimesForSet,:,nc_index);
            EmbryoNCProfileNumQuiescentNuclei(isnan(EmbryoNCProfileNumQuiescentNuclei)) = 0;
            EmbryoNCProfileNumQuiescentNuclei(isnan(EmbryoNCProfileTraceCount)) = NaN;
            
            EmbryoNCProfileNumFinishedTranscribingNuclei = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins);
            EmbryoNCProfileNumFinishedTranscribingNuclei(MatchingTimesForBinned,:) = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleNumFinishedTranscribingNuclei(MatchingTimesForSet,:,nc_index);
            EmbryoNCProfileNumFinishedTranscribingNuclei(isnan(EmbryoNCProfileNumFinishedTranscribingNuclei)) = 0; 
            EmbryoNCProfileNumFinishedTranscribingNuclei(isnan(EmbryoNCProfileTraceCount)) = NaN;
            
            EmbryoNCProfileFractionOn(EmbryoNCProfileNumNuclei == 0) = NaN;
            
            EmbryoNCProfileGradedNumNuclei = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 10);
            EmbryoNCProfileGradedNumNuclei(MatchingTimesForBinned,:,:)  = squeeze(ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleGradedTotalNuclei(MatchingTimesForSet,:,nc_index,:));
            EmbryoNCProfileGradedNumNuclei(isnan(EmbryoNCProfileGradedNumNuclei)) = 0;
            EmbryoNCProfileGradedNumNuclei(repmat(isnan(EmbryoNCProfileTraceCount), 1, 1, 10)) = NaN;
            
            EmbryoNCProfileGradedNumOnNuclei = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 10);
            EmbryoNCProfileGradedNumOnNuclei(MatchingTimesForBinned,:,:)  = squeeze(ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleGradedOnNuclei(MatchingTimesForSet,:,nc_index,:));
            EmbryoNCProfileGradedNumOnNuclei(isnan(EmbryoNCProfileGradedNumOnNuclei)) = 0;
            EmbryoNCProfileGradedNumOnNuclei(repmat(isnan(EmbryoNCProfileTraceCount), 1, 1, 10)) = NaN;
            
            
            EmbryoNCProfileGradedNumOffNuclei = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 10);
            EmbryoNCProfileGradedNumOffNuclei(MatchingTimesForBinned,:,:) = squeeze(ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleGradedOffNuclei(MatchingTimesForSet,:,nc_index,:));
            EmbryoNCProfileGradedNumOffNuclei(isnan(EmbryoNCProfileGradedNumOffNuclei)) = 0;
            EmbryoNCProfileGradedNumOffNuclei(repmat(isnan(EmbryoNCProfileTraceCount), 1, 1, 10)) = NaN;
            
            
            
            EmbryoNCProfileGradedFractionOn = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 10);
            EmbryoNCProfileGradedFractionOn(MatchingTimesForBinned,:,:) = squeeze(ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleGradedFractionOn(MatchingTimesForSet,:,nc_index,:));
            EmbryoNCProfileGradedFractionOn(EmbryoNCProfileGradedFractionOn == 0) = NaN;
            EmbryoNCProfileGradedFractionOn(repmat(isnan(EmbryoNCProfileTraceCount), 1, 1, 10)) = NaN;
            
            EmbryoProfiles(find(~all(isnan(EmbryoNCProfile), 2)),:,SetIndex) = ...
                EmbryoNCProfile(~all(isnan(EmbryoNCProfile), 2),:);
            EmbryoProfStds(find(~all(isnan(EmbryoNCProfile), 2)),:,SetIndex) = ...
                EmbryoNCProfileStd(~all(isnan(EmbryoNCProfile), 2),:);
            EmbryoNuclearProfiles(find(~all(isnan(EmbryoNCNuclearProfile), 2)),:,SetIndex) = ...
                EmbryoNCNuclearProfile(~all(isnan(EmbryoNCNuclearProfile), 2),:);
            EmbryoFractionOns(find(~all(isnan(EmbryoNCProfileFractionOn), 2)),:,SetIndex) = ...
                EmbryoNCProfileFractionOn(~all(isnan(EmbryoNCProfileFractionOn), 2),:);
            EmbryoTraceCounts(find(~all(isnan(EmbryoNCProfileTraceCount), 2)),:,SetIndex) = ...
                EmbryoNCProfileTraceCount(~all(isnan(EmbryoNCProfileTraceCount), 2),:);
            EmbryoProfsNumNuclei(find(~all(isnan(EmbryoNCProfileNumNuclei), 2)),:,SetIndex) = ...
                EmbryoNCProfileNumNuclei(~all(isnan(EmbryoNCProfileNumNuclei), 2),:);
            EmbryoProfsNumOnNuclei(find(~all(isnan(EmbryoNCProfileNumOnNuclei), 2)),:,SetIndex) = ...
                EmbryoNCProfileNumOnNuclei(~all(isnan(EmbryoNCProfileNumOnNuclei), 2),:);
            EmbryoProfsNumOffNuclei(find(~all(isnan(EmbryoNCProfileNumOffNuclei), 2)),:,SetIndex) = ...
                EmbryoNCProfileNumOffNuclei(~all(isnan(EmbryoNCProfileNumOffNuclei), 2),:);
            EmbryoProfsNumQuiescentNuclei(find(~all(isnan(EmbryoNCProfileNumQuiescentNuclei), 2)),:,SetIndex) = ...
                EmbryoNCProfileNumQuiescentNuclei(~all(isnan(EmbryoNCProfileNumQuiescentNuclei), 2),:);
            EmbryoProfsNumFinishedTranscribingNuclei(find(~all(isnan(EmbryoNCProfileNumFinishedTranscribingNuclei), 2)),:,SetIndex) = ...
                EmbryoNCProfileNumFinishedTranscribingNuclei(~all(isnan(EmbryoNCProfileNumFinishedTranscribingNuclei), 2),:);
            
            for graded_index = 1:10
                EmbryoProfsGradedNumNuclei(find(~all(isnan(EmbryoNCProfileGradedNumNuclei(:,:,graded_index)), 2)),:,graded_index,SetIndex) = ...
                    EmbryoNCProfileGradedNumNuclei(~all(isnan(EmbryoNCProfileGradedNumNuclei(:,:,graded_index)), 2),:,graded_index);
                EmbryoProfsGradedNumOnNuclei(find(~all(isnan(EmbryoNCProfileGradedNumOnNuclei(:,:,graded_index)), 2)),:,graded_index,SetIndex) = ...
                    EmbryoNCProfileGradedNumOnNuclei(~all(isnan(EmbryoNCProfileGradedNumOnNuclei(:,:,graded_index)), 2),:,graded_index);
                EmbryoProfsGradedNumOffNuclei(find(~all(isnan(EmbryoNCProfileGradedNumOffNuclei(:,:,graded_index)), 2)),:,graded_index,SetIndex) = ...
                    EmbryoNCProfileGradedNumOffNuclei(~all(isnan(EmbryoNCProfileGradedNumOffNuclei(:,:,graded_index)), 2),:,graded_index);
                EmbryoGradedFractionOns(find(~all(isnan(EmbryoNCProfileGradedFractionOn(:,:,graded_index)), 2)),:,graded_index,SetIndex) = ...
                    EmbryoNCProfileGradedFractionOn(~all(isnan(EmbryoNCProfileGradedFractionOn(:,:,graded_index)), 2),:,graded_index);
            end
        end
        
        BadEntries = EmbryoProfiles< ltm.MinimumSchnitzCount;
        EmbryoProfiles(EmbryoProfiles< ltm.MinimumSchnitzCount) = NaN; 
        TempMeanProfile = mean(EmbryoProfiles, 3, 'omitnan');
        TempStdError = std(EmbryoProfiles, 0, 3, 'omitnan');
        
        EmbryoNuclearProfiles(EmbryoNuclearProfiles< ltm.MinimumSchnitzCount) = NaN;
        TempNuclearMeanProfile = mean(EmbryoNuclearProfiles, 3, 'omitnan');
        TempNuclearStdError = std(EmbryoNuclearProfiles, 0, 3, 'omitnan');
        
        TempProfCount = sum(~isnan(EmbryoProfiles) & EmbryoProfiles >=  ltm.MinimumSchnitzCount , 3);
        
        EmbryoProfsNumNuclei(BadEntries) = NaN;
        TempTotalNuclei = sum(EmbryoProfsNumNuclei, 3, 'omitnan');
        
        EmbryoProfsNumOnNuclei(BadEntries) = NaN;
        TempTotalOnNuclei = sum(EmbryoProfsNumOnNuclei, 3, 'omitnan');
        
        %EmbryoProfsNumOffNuclei(BadEntries) = NaN;
        TempTotalOffNuclei = sum(EmbryoProfsNumOffNuclei, 3, 'omitnan');
        
        %EmbryoProfsNumQuiescentNuclei(BadEntries) = NaN;
        TempTotalQuiescentNuclei = sum(EmbryoProfsNumQuiescentNuclei, 3, 'omitnan');
        
        %EmbryoProfsNumFinishedTranscribingNuclei(BadEntries) = NaN;
        TempTotalFinishedTranscribingNuclei = sum(EmbryoProfsNumFinishedTranscribingNuclei, 3, 'omitnan');
        
        
        TempMeanFractionOnNuclei = mean(EmbryoFractionOns, 3, 'omitnan');
        TempStdErrFractionOnNuclei = std(EmbryoFractionOns, 0, 3, 'omitnan');
        
        TempMeanProfile(TempMeanFractionOnNuclei == 0) = 0;
        TempStdError(TempMeanFractionOnNuclei == 0) = 0;
        
        TempTotalGradedNuclei = sum(EmbryoProfsGradedNumNuclei, 4, 'omitnan');
        TempTotalGradedOnNuclei = sum(EmbryoProfsGradedNumOnNuclei, 4, 'omitnan');
        TempTotalGradedOffNuclei = sum(EmbryoProfsGradedNumOffNuclei, 4, 'omitnan');
        TempMeanGradedFractionOnNuclei = mean(EmbryoGradedFractionOns, 4, 'omitnan');
        
        ltm.BinnedMeanProfiles.AnaphaseAlignedCycleMeanTraces(:,:,nc_index,t_index) = TempMeanProfile;
        ltm.BinnedMeanProfiles.AnaphaseAlignedCycleNumEmbryos(:,:,nc_index,t_index) = TempProfCount;
        ltm.BinnedMeanProfiles.AnaphaseAlignedCycleStdErrors(:,:,nc_index,t_index) = TempStdError;
        ltm.BinnedMeanProfiles.AnaphaseAlignedCycleNuclearMeanTraces(:,:,nc_index,t_index) = TempNuclearMeanProfile;
        ltm.BinnedMeanProfiles.AnaphaseAlignedCycleNuclearStdErrors(:,:,nc_index,t_index) = TempNuclearStdError;
        ltm.BinnedMeanProfiles.AnaphaseAlignedCycleTotalNuclei(:,:,nc_index,t_index) = TempTotalNuclei;
        ltm.BinnedMeanProfiles.AnaphaseAlignedCycleTotalOnNuclei(:,:,nc_index,t_index) = TempTotalOnNuclei;
        ltm.BinnedMeanProfiles.AnaphaseAlignedCycleTotalOffNuclei(:,:,nc_index,t_index) = TempTotalOffNuclei;
        ltm.BinnedMeanProfiles.AnaphaseAlignedCycleTotalQuiescentNuclei(:,:,nc_index,t_index) = TempTotalQuiescentNuclei;
        ltm.BinnedMeanProfiles.AnaphaseAlignedCycleTotalFinishedTranscribingNuclei(:,:,nc_index,t_index) = TempTotalFinishedTranscribingNuclei;
        ltm.BinnedMeanProfiles.AnaphaseAlignedCycleFractionOn(:,:,nc_index,t_index) = TempMeanFractionOnNuclei;
        ltm.BinnedMeanProfiles.AnaphaseAlignedCycleFractionOnStdErr(:,:,nc_index,t_index) = TempStdErrFractionOnNuclei;
        
        ltm.BinnedMeanProfiles.AnaphaseAlignedCycleGradedTotalNuclei(:,:,:,nc_index,t_index) = TempTotalGradedNuclei;
        ltm.BinnedMeanProfiles.AnaphaseAlignedCycleGradedTotalOnNuclei(:,:,:,nc_index,t_index) = TempTotalGradedOnNuclei;
        ltm.BinnedMeanProfiles.AnaphaseAlignedCycleGradedTotalOffNuclei(:,:,:,nc_index,t_index) = TempTotalGradedOffNuclei;
        ltm.BinnedMeanProfiles.AnaphaseAlignedCycleGradedFractionOn(:,:,:,nc_index,t_index) = TempMeanGradedFractionOnNuclei;
     
    end
    
end

%% Then make average profiles binning anaphase aligned 3D traces for each temperature 
for t_index = 1:NumTemperatures
    EmbryoIndices = find(EmbryoTemperatures == Temperatures(t_index) & IncludeExperimentTF);
    for nc_index = 1:6
        EmbryoProfiles = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfStds =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfsNumNuclei =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfsNumOnNuclei =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfsNumOffNuclei =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfsNumQuiescentNuclei =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfsNumFinishedTranscribingNuclei =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoNuclearProfiles =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoFractionOns =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoTraceCounts =  NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, length(EmbryoIndices));
        for SetIndex = 1:length(EmbryoIndices)
            EmbryoIndex = EmbryoIndices(SetIndex);
            MatchingTimes = ismember(ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleTimes, ltm.MeanProfiles{EmbryoIndex}.AnaphaseAligned3DCycleFrameTimes{nc_index});
            EmbryoNCProfile = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAligned3DCycleMeanTraces(:,:,nc_index);
            EmbryoNCNuclearProfile = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAligned3DCycleMeanPerNucleusTraces(:,:,nc_index);
            EmbryoNCProfileStd = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAligned3DCycleTraceStdErrors(:,:,nc_index);
            EmbryoNCProfileFractionOn= ltm.MeanProfiles{EmbryoIndex}.AnaphaseAligned3DCycleFractionOn(:,:,nc_index);
            EmbryoNCProfileTraceCount = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAligned3DCycleTraceCount(:,:,nc_index);
            EmbryoNCProfileTraceCount(isnan(EmbryoNCProfileTraceCount)) = 0;
            EmbryoNCProfileNumNuclei = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAligned3DCycleNumNuclei(:,:,nc_index);
            EmbryoNCProfileNumNuclei(isnan(EmbryoNCProfileNumNuclei)) = 0;
            EmbryoNCProfileNumOnNuclei = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAligned3DCycleNumOnNuclei(:,:,nc_index);
            EmbryoNCProfileNumOnNuclei(isnan(EmbryoNCProfileNumOnNuclei)) = 0;
            EmbryoNCProfileNumOffNuclei = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAligned3DCycleNumOffNuclei(:,:,nc_index);
            EmbryoNCProfileNumOffNuclei(isnan(EmbryoNCProfileNumOffNuclei)) = 0;
            EmbryoNCProfileNumQuiescentNuclei = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAligned3DCycleNumQuiescentNuclei(:,:,nc_index);
            EmbryoNCProfileNumQuiescentNuclei(isnan(EmbryoNCProfileNumQuiescentNuclei)) = 0;
            EmbryoNCProfileNumFinishedTranscribingNuclei = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAligned3DCycleNumFinishedTranscribingNuclei(:,:,nc_index);
            EmbryoNCProfileNumFinishedTranscribingNuclei(isnan(EmbryoNCProfileNumFinishedTranscribingNuclei)) = 0; 
            EmbryoNCProfileFractionOn(EmbryoNCProfileNumNuclei == 0) = NaN;
            
            EmbryoProfiles(find(~all(isnan(EmbryoNCProfile), 2)),:,SetIndex) = ...
                EmbryoNCProfile(~all(isnan(EmbryoNCProfile), 2),:);
            EmbryoProfStds(find(~all(isnan(EmbryoNCProfile), 2)),:,SetIndex) = ...
                EmbryoNCProfileStd(~all(isnan(EmbryoNCProfile), 2),:);
            EmbryoNuclearProfiles(find(~all(isnan(EmbryoNCNuclearProfile), 2)),:,SetIndex) = ...
                EmbryoNCNuclearProfile(~all(isnan(EmbryoNCNuclearProfile), 2),:);
            EmbryoFractionOns(find(~all(isnan(EmbryoNCProfileFractionOn), 2)),:,SetIndex) = ...
                EmbryoNCProfileFractionOn(~all(isnan(EmbryoNCProfileFractionOn), 2),:);
            EmbryoTraceCounts(find(~all(isnan(EmbryoNCProfileTraceCount), 2)),:,SetIndex) = ...
                EmbryoNCProfileTraceCount(~all(isnan(EmbryoNCProfileTraceCount), 2),:);
            EmbryoProfsNumNuclei(find(~all(isnan(EmbryoNCProfileNumNuclei), 2)),:,SetIndex) = ...
                EmbryoNCProfileNumNuclei(~all(isnan(EmbryoNCProfileNumNuclei), 2),:);
            EmbryoProfsNumOnNuclei(find(~all(isnan(EmbryoNCProfileNumOnNuclei), 2)),:,SetIndex) = ...
                EmbryoNCProfileNumOnNuclei(~all(isnan(EmbryoNCProfileNumOnNuclei), 2),:);
            EmbryoProfsNumOffNuclei(find(~all(isnan(EmbryoNCProfileNumOffNuclei), 2)),:,SetIndex) = ...
                EmbryoNCProfileNumOffNuclei(~all(isnan(EmbryoNCProfileNumOffNuclei), 2),:);
            EmbryoProfsNumQuiescentNuclei(find(~all(isnan(EmbryoNCProfileNumQuiescentNuclei), 2)),:,SetIndex) = ...
                EmbryoNCProfileNumQuiescentNuclei(~all(isnan(EmbryoNCProfileNumQuiescentNuclei), 2),:);
            EmbryoProfsNumFinishedTranscribingNuclei(find(~all(isnan(EmbryoNCProfileNumFinishedTranscribingNuclei), 2)),:,SetIndex) = ...
                EmbryoNCProfileNumFinishedTranscribingNuclei(~all(isnan(EmbryoNCProfileNumFinishedTranscribingNuclei), 2),:);
        end
        
        BadEntries = EmbryoProfiles< ltm.MinimumSchnitzCount;
        EmbryoProfiles(EmbryoProfiles< ltm.MinimumSchnitzCount) = NaN; 
        TempMeanProfile = mean(EmbryoProfiles, 3, 'omitnan');
        TempStdError = std(EmbryoProfiles, 0, 3, 'omitnan');
        
        EmbryoNuclearProfiles(EmbryoNuclearProfiles< ltm.MinimumSchnitzCount) = NaN;
        TempNuclearMeanProfile = mean(EmbryoNuclearProfiles, 3, 'omitnan');
        TempNuclearStdError = std(EmbryoNuclearProfiles, 0, 3, 'omitnan');
        
        TempProfCount = sum(~isnan(EmbryoProfiles) & EmbryoProfiles >=  ltm.MinimumSchnitzCount , 3, 'omitnan');
        
        EmbryoProfsNumNuclei(BadEntries) = NaN;
        TempTotalNuclei = sum(EmbryoProfsNumNuclei, 3, 'omitnan');
        
        EmbryoProfsNumOnNuclei(BadEntries) = NaN;
        TempTotalOnNuclei = sum(EmbryoProfsNumOnNuclei, 3, 'omitnan');
        
        %EmbryoProfsNumOffNuclei(BadEntries) = NaN;
        TempTotalOffNuclei = sum(EmbryoProfsNumOffNuclei, 3, 'omitnan');
        
        %EmbryoProfsNumQuiescentNuclei(BadEntries) = NaN;
        TempTotalQuiescentNuclei = sum(EmbryoProfsNumQuiescentNuclei, 3, 'omitnan');
        
        %EmbryoProfsNumFinishedTranscribingNuclei(BadEntries) = NaN;
        TempTotalFinishedTranscribingNuclei = sum(EmbryoProfsNumFinishedTranscribingNuclei, 3, 'omitnan');
        
        
        TempMeanFractionOnNuclei = mean(EmbryoFractionOns, 3, 'omitnan');
        TempStdErrFractionOnNuclei = std(EmbryoFractionOns, 0, 3, 'omitnan');
        
        TempMeanProfile(TempMeanFractionOnNuclei == 0) = 0;
        TempStdError(TempMeanFractionOnNuclei == 0) = 0;
        
        ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleMeanTraces(:,:,nc_index,t_index) = TempMeanProfile;
        ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleNumEmbryos(:,:,nc_index,t_index) = TempProfCount;
        ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleStdErrors(:,:,nc_index,t_index) = TempStdError;
        ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleNuclearMeanTraces(:,:,nc_index,t_index) = TempNuclearMeanProfile;
        ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleNuclearStdErrors(:,:,nc_index,t_index) = TempNuclearStdError;
        ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleTotalNuclei(:,:,nc_index,t_index) = TempTotalNuclei;
        ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleTotalOnNuclei(:,:,nc_index,t_index) = TempTotalOnNuclei;
        ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleTotalOffNuclei(:,:,nc_index,t_index) = TempTotalOffNuclei;
        ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleTotalQuiescentNuclei(:,:,nc_index,t_index) = TempTotalQuiescentNuclei;
        ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleTotalFinishedTranscribingNuclei(:,:,nc_index,t_index) = TempTotalFinishedTranscribingNuclei;
        ltm.BinnedMeanProfiles.AnaphaseAligned3DCycleFractionOn(:,:,nc_index,t_index) = TempMeanFractionOnNuclei;
     
    end
    
end
   

%% First make average profiles binning anaphase aligned traces for each temperature 
for t_index = 1:NumTemperatures
    EmbryoIndices = find(EmbryoTemperatures == Temperatures(t_index) & IncludeExperimentTF);
    for nc_index = 1:6
        EmbryoProfiles = NaN(LengthBinnedTbinnedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfStds =  NaN(LengthBinnedTbinnedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfsNumNuclei =  NaN(LengthBinnedTbinnedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfsNumOnNuclei =  NaN(LengthBinnedTbinnedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfsNumOffNuclei =  NaN(LengthBinnedTbinnedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfsGradedNumNuclei =  NaN(LengthBinnedTbinnedProfile, NumAPbins,10, length(EmbryoIndices));
        EmbryoProfsGradedNumOnNuclei =  NaN(LengthBinnedTbinnedProfile, NumAPbins, 10,length(EmbryoIndices));
        EmbryoProfsGradedNumOffNuclei =  NaN(LengthBinnedTbinnedProfile, NumAPbins,10, length(EmbryoIndices));
        EmbryoProfsNumQuiescentNuclei =  NaN(LengthBinnedTbinnedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfsNumFinishedTranscribingNuclei =  NaN(LengthBinnedTbinnedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoNuclearProfiles =  NaN(LengthBinnedTbinnedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoFractionOns =  NaN(LengthBinnedTbinnedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoGradedFractionOns =  NaN(LengthBinnedTbinnedProfile, NumAPbins,10, length(EmbryoIndices));
        EmbryoTraceCounts =  NaN(LengthBinnedTbinnedProfile, NumAPbins, length(EmbryoIndices));
        for SetIndex = 1:length(EmbryoIndices)
            EmbryoIndex = EmbryoIndices(SetIndex);
            MatchingTimesForBinned = ismember(ltm.BinnedMeanProfiles.TbinnedCycleTimes, ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleFrameTimes{nc_index});
            MatchingTimesForSet = ismember(ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleFrameTimes{nc_index}, ltm.BinnedMeanProfiles.TbinnedCycleTimes);
            
            if isempty(MatchingTimesForBinned) | sum(MatchingTimesForBinned) == 0
                continue
            end
            EmbryoNCProfile = NaN(LengthBinnedTbinnedProfile, NumAPbins);
            EmbryoNCProfile(MatchingTimesForBinned,:) = ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleMeanTraces(MatchingTimesForSet,:,nc_index);
            
            EmbryoNCNuclearProfile = NaN(LengthBinnedTbinnedProfile, NumAPbins);
            EmbryoNCNuclearProfile(MatchingTimesForBinned,:) = ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleMeanPerNucleusTraces(MatchingTimesForSet,:,nc_index);
            
            EmbryoNCProfileStd = NaN(LengthBinnedTbinnedProfile, NumAPbins);
            EmbryoNCProfileStd(MatchingTimesForBinned,:) = ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleTraceStdErrors(MatchingTimesForSet,:,nc_index);
            
            EmbryoNCProfileFractionOn = NaN(LengthBinnedTbinnedProfile, NumAPbins);
            EmbryoNCProfileFractionOn(MatchingTimesForBinned,:) = ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleFractionOn(MatchingTimesForSet,:,nc_index);
            
            EmbryoNCProfileTraceCount = NaN(LengthBinnedTbinnedProfile, NumAPbins);
            EmbryoNCProfileTraceCount(MatchingTimesForBinned,:) = ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleTraceCount(MatchingTimesForSet,:,nc_index);
%             EmbryoNCProfileTraceCount(isnan(EmbryoNCProfileTraceCount)) = 0;
            
            EmbryoNCProfileNumNuclei = NaN(LengthBinnedTbinnedProfile, NumAPbins);
            EmbryoNCProfileNumNuclei(MatchingTimesForBinned,:)  = ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleNumNuclei(MatchingTimesForSet,:,nc_index);
            EmbryoNCProfileNumNuclei(isnan(EmbryoNCProfileNumNuclei)) = 0;
            EmbryoNCProfileNumNuclei(isnan(EmbryoNCProfileTraceCount)) = NaN;
            
            EmbryoNCProfileNumOnNuclei = NaN(LengthBinnedTbinnedProfile, NumAPbins);
            EmbryoNCProfileNumOnNuclei(MatchingTimesForBinned,:)  = ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleNumOnNuclei(MatchingTimesForSet,:,nc_index);
            EmbryoNCProfileNumOnNuclei(isnan(EmbryoNCProfileNumOnNuclei)) = 0;
            EmbryoNCProfileNumOnNuclei(isnan(EmbryoNCProfileTraceCount)) = NaN;
            
            EmbryoNCProfileNumOffNuclei = NaN(LengthBinnedTbinnedProfile, NumAPbins);
            EmbryoNCProfileNumOffNuclei(MatchingTimesForBinned,:) = ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleNumOffNuclei(MatchingTimesForSet,:,nc_index);
            EmbryoNCProfileNumOffNuclei(isnan(EmbryoNCProfileNumOffNuclei)) = 0;
            EmbryoNCProfileNumOffNuclei(isnan(EmbryoNCProfileTraceCount)) = NaN;
            
            
            EmbryoNCProfileNumQuiescentNuclei = NaN(LengthBinnedTbinnedProfile, NumAPbins);
            EmbryoNCProfileNumQuiescentNuclei(MatchingTimesForBinned,:)  = ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleNumQuiescentNuclei(MatchingTimesForSet,:,nc_index);
            EmbryoNCProfileNumQuiescentNuclei(isnan(EmbryoNCProfileNumQuiescentNuclei)) = 0;
            EmbryoNCProfileNumQuiescentNuclei(isnan(EmbryoNCProfileTraceCount)) = NaN;
            
            EmbryoNCProfileNumFinishedTranscribingNuclei = NaN(LengthBinnedTbinnedProfile, NumAPbins);
            EmbryoNCProfileNumFinishedTranscribingNuclei(MatchingTimesForBinned,:) = ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleNumFinishedTranscribingNuclei(MatchingTimesForSet,:,nc_index);
            EmbryoNCProfileNumFinishedTranscribingNuclei(isnan(EmbryoNCProfileNumFinishedTranscribingNuclei)) = 0; 
            EmbryoNCProfileNumFinishedTranscribingNuclei(isnan(EmbryoNCProfileTraceCount)) = NaN;
            
            EmbryoNCProfileFractionOn(EmbryoNCProfileNumNuclei == 0) = NaN;
            
            EmbryoNCProfileGradedNumNuclei = NaN(LengthBinnedTbinnedProfile, NumAPbins, 10);
            EmbryoNCProfileGradedNumNuclei(MatchingTimesForBinned,:,:)  = squeeze(ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleGradedTotalNuclei(MatchingTimesForSet,:,nc_index,:));
            EmbryoNCProfileGradedNumNuclei(isnan(EmbryoNCProfileGradedNumNuclei)) = 0;
            EmbryoNCProfileGradedNumNuclei(repmat(isnan(EmbryoNCProfileTraceCount), 1, 1, 10)) = NaN;
            
            EmbryoNCProfileGradedNumOnNuclei = NaN(LengthBinnedTbinnedProfile, NumAPbins, 10);
            EmbryoNCProfileGradedNumOnNuclei(MatchingTimesForBinned,:,:)  = squeeze(ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleGradedOnNuclei(MatchingTimesForSet,:,nc_index,:));
            EmbryoNCProfileGradedNumOnNuclei(isnan(EmbryoNCProfileGradedNumOnNuclei)) = 0;
            EmbryoNCProfileGradedNumOnNuclei(repmat(isnan(EmbryoNCProfileTraceCount), 1, 1, 10)) = NaN;
            
            
            EmbryoNCProfileGradedNumOffNuclei = NaN(LengthBinnedTbinnedProfile, NumAPbins, 10);
            EmbryoNCProfileGradedNumOffNuclei(MatchingTimesForBinned,:,:) = squeeze(ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleGradedOffNuclei(MatchingTimesForSet,:,nc_index,:));
            EmbryoNCProfileGradedNumOffNuclei(isnan(EmbryoNCProfileGradedNumOffNuclei)) = 0;
            EmbryoNCProfileGradedNumOffNuclei(repmat(isnan(EmbryoNCProfileTraceCount), 1, 1, 10)) = NaN;
            
            
            
            EmbryoNCProfileGradedFractionOn = NaN(LengthBinnedTbinnedProfile, NumAPbins, 10);
            EmbryoNCProfileGradedFractionOn(MatchingTimesForBinned,:,:) = squeeze(ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleGradedFractionOn(MatchingTimesForSet,:,nc_index,:));
            EmbryoNCProfileGradedFractionOn(EmbryoNCProfileGradedFractionOn == 0) = NaN;
            EmbryoNCProfileGradedFractionOn(repmat(isnan(EmbryoNCProfileTraceCount), 1, 1, 10)) = NaN;
            
            EmbryoProfiles(find(~all(isnan(EmbryoNCProfile), 2)),:,SetIndex) = ...
                EmbryoNCProfile(~all(isnan(EmbryoNCProfile), 2),:);
            EmbryoProfStds(find(~all(isnan(EmbryoNCProfile), 2)),:,SetIndex) = ...
                EmbryoNCProfileStd(~all(isnan(EmbryoNCProfile), 2),:);
            EmbryoNuclearProfiles(find(~all(isnan(EmbryoNCNuclearProfile), 2)),:,SetIndex) = ...
                EmbryoNCNuclearProfile(~all(isnan(EmbryoNCNuclearProfile), 2),:);
            EmbryoFractionOns(find(~all(isnan(EmbryoNCProfileFractionOn), 2)),:,SetIndex) = ...
                EmbryoNCProfileFractionOn(~all(isnan(EmbryoNCProfileFractionOn), 2),:);
            EmbryoTraceCounts(find(~all(isnan(EmbryoNCProfileTraceCount), 2)),:,SetIndex) = ...
                EmbryoNCProfileTraceCount(~all(isnan(EmbryoNCProfileTraceCount), 2),:);
            EmbryoProfsNumNuclei(find(~all(isnan(EmbryoNCProfileNumNuclei), 2)),:,SetIndex) = ...
                EmbryoNCProfileNumNuclei(~all(isnan(EmbryoNCProfileNumNuclei), 2),:);
            EmbryoProfsNumOnNuclei(find(~all(isnan(EmbryoNCProfileNumOnNuclei), 2)),:,SetIndex) = ...
                EmbryoNCProfileNumOnNuclei(~all(isnan(EmbryoNCProfileNumOnNuclei), 2),:);
            EmbryoProfsNumOffNuclei(find(~all(isnan(EmbryoNCProfileNumOffNuclei), 2)),:,SetIndex) = ...
                EmbryoNCProfileNumOffNuclei(~all(isnan(EmbryoNCProfileNumOffNuclei), 2),:);
            EmbryoProfsNumQuiescentNuclei(find(~all(isnan(EmbryoNCProfileNumQuiescentNuclei), 2)),:,SetIndex) = ...
                EmbryoNCProfileNumQuiescentNuclei(~all(isnan(EmbryoNCProfileNumQuiescentNuclei), 2),:);
            EmbryoProfsNumFinishedTranscribingNuclei(find(~all(isnan(EmbryoNCProfileNumFinishedTranscribingNuclei), 2)),:,SetIndex) = ...
                EmbryoNCProfileNumFinishedTranscribingNuclei(~all(isnan(EmbryoNCProfileNumFinishedTranscribingNuclei), 2),:);
            
            for graded_index = 1:10
                EmbryoProfsGradedNumNuclei(find(~all(isnan(EmbryoNCProfileGradedNumNuclei(:,:,graded_index)), 2)),:,graded_index,SetIndex) = ...
                    EmbryoNCProfileGradedNumNuclei(~all(isnan(EmbryoNCProfileGradedNumNuclei(:,:,graded_index)), 2),:,graded_index);
                EmbryoProfsGradedNumOnNuclei(find(~all(isnan(EmbryoNCProfileGradedNumOnNuclei(:,:,graded_index)), 2)),:,graded_index,SetIndex) = ...
                    EmbryoNCProfileGradedNumOnNuclei(~all(isnan(EmbryoNCProfileGradedNumOnNuclei(:,:,graded_index)), 2),:,graded_index);
                EmbryoProfsGradedNumOffNuclei(find(~all(isnan(EmbryoNCProfileGradedNumOffNuclei(:,:,graded_index)), 2)),:,graded_index,SetIndex) = ...
                    EmbryoNCProfileGradedNumOffNuclei(~all(isnan(EmbryoNCProfileGradedNumOffNuclei(:,:,graded_index)), 2),:,graded_index);
                EmbryoGradedFractionOns(find(~all(isnan(EmbryoNCProfileGradedFractionOn(:,:,graded_index)), 2)),:,graded_index,SetIndex) = ...
                    EmbryoNCProfileGradedFractionOn(~all(isnan(EmbryoNCProfileGradedFractionOn(:,:,graded_index)), 2),:,graded_index);
            end
        end
        
        BadEntries = EmbryoProfiles< ltm.MinimumSchnitzCount;
        EmbryoProfiles(EmbryoProfiles< ltm.MinimumSchnitzCount) = NaN; 
        TempMeanProfile = mean(EmbryoProfiles, 3, 'omitnan');
        TempStdError = std(EmbryoProfiles, 0, 3, 'omitnan');
        
        EmbryoNuclearProfiles(EmbryoNuclearProfiles< ltm.MinimumSchnitzCount) = NaN;
        TempNuclearMeanProfile = mean(EmbryoNuclearProfiles, 3, 'omitnan');
        TempNuclearStdError = std(EmbryoNuclearProfiles, 0, 3, 'omitnan');
        
        TempProfCount = sum(~isnan(EmbryoProfiles) & EmbryoProfiles >=  ltm.MinimumSchnitzCount , 3, 'omitnan');
        
        EmbryoProfsNumNuclei(BadEntries) = NaN;
        TempTotalNuclei = sum(EmbryoProfsNumNuclei, 3, 'omitnan');
        
        EmbryoProfsNumOnNuclei(BadEntries) = NaN;
        TempTotalOnNuclei = sum(EmbryoProfsNumOnNuclei, 3, 'omitnan');
        
        %EmbryoProfsNumOffNuclei(BadEntries) = NaN;
        TempTotalOffNuclei = sum(EmbryoProfsNumOffNuclei, 3, 'omitnan');
        
        %EmbryoProfsNumQuiescentNuclei(BadEntries) = NaN;
        TempTotalQuiescentNuclei = sum(EmbryoProfsNumQuiescentNuclei, 3, 'omitnan');
        
        %EmbryoProfsNumFinishedTranscribingNuclei(BadEntries) = NaN;
        TempTotalFinishedTranscribingNuclei = sum(EmbryoProfsNumFinishedTranscribingNuclei, 3, 'omitnan');
        
        
        TempMeanFractionOnNuclei = mean(EmbryoFractionOns, 3, 'omitnan');
        TempStdErrFractionOnNuclei = std(EmbryoFractionOns, 0, 3, 'omitnan');
        
        TempMeanProfile(TempMeanFractionOnNuclei == 0) = 0;
        TempStdError(TempMeanFractionOnNuclei == 0) = 0;
        
        TempTotalGradedNuclei = sum(EmbryoProfsGradedNumNuclei, 4, 'omitnan');
        TempTotalGradedOnNuclei = sum(EmbryoProfsGradedNumOnNuclei, 4, 'omitnan');
        TempTotalGradedOffNuclei = sum(EmbryoProfsGradedNumOffNuclei, 4, 'omitnan');
        TempMeanGradedFractionOnNuclei = mean(EmbryoGradedFractionOns, 4, 'omitnan');
        
        ltm.BinnedMeanProfiles.TbinnedCycleMeanTraces(:,:,nc_index,t_index) = TempMeanProfile;
        ltm.BinnedMeanProfiles.TbinnedCycleNumEmbryos(:,:,nc_index,t_index) = TempProfCount;
        ltm.BinnedMeanProfiles.TbinnedCycleStdErrors(:,:,nc_index,t_index) = TempStdError;
        ltm.BinnedMeanProfiles.TbinnedCycleNuclearMeanTraces(:,:,nc_index,t_index) = TempNuclearMeanProfile;
        ltm.BinnedMeanProfiles.TbinnedCycleNuclearStdErrors(:,:,nc_index,t_index) = TempNuclearStdError;
        ltm.BinnedMeanProfiles.TbinnedCycleTotalNuclei(:,:,nc_index,t_index) = TempTotalNuclei;
        ltm.BinnedMeanProfiles.TbinnedCycleTotalOnNuclei(:,:,nc_index,t_index) = TempTotalOnNuclei;
        ltm.BinnedMeanProfiles.TbinnedCycleTotalOffNuclei(:,:,nc_index,t_index) = TempTotalOffNuclei;
        ltm.BinnedMeanProfiles.TbinnedCycleTotalQuiescentNuclei(:,:,nc_index,t_index) = TempTotalQuiescentNuclei;
        ltm.BinnedMeanProfiles.TbinnedCycleTotalFinishedTranscribingNuclei(:,:,nc_index,t_index) = TempTotalFinishedTranscribingNuclei;
        ltm.BinnedMeanProfiles.TbinnedCycleFractionOn(:,:,nc_index,t_index) = TempMeanFractionOnNuclei;
        ltm.BinnedMeanProfiles.TbinnedCycleFractionOnStdErr(:,:,nc_index,t_index) = TempStdErrFractionOnNuclei;
        
        ltm.BinnedMeanProfiles.TbinnedCycleGradedTotalNuclei(:,:,:,nc_index,t_index) = TempTotalGradedNuclei;
        ltm.BinnedMeanProfiles.TbinnedCycleGradedTotalOnNuclei(:,:,:,nc_index,t_index) = TempTotalGradedOnNuclei;
        ltm.BinnedMeanProfiles.TbinnedCycleGradedTotalOffNuclei(:,:,:,nc_index,t_index) = TempTotalGradedOffNuclei;
        ltm.BinnedMeanProfiles.TbinnedCycleGradedFractionOn(:,:,:,nc_index,t_index) = TempMeanGradedFractionOnNuclei;
     
    end
    
end

%% Then make average profiles binning tbinned 3D traces for each temperature 
for t_index = 1:NumTemperatures
    EmbryoIndices = find(EmbryoTemperatures == Temperatures(t_index) & IncludeExperimentTF);
    for nc_index = 1:6
        EmbryoProfiles = NaN(LengthBinnedTbinnedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfStds =  NaN(LengthBinnedTbinnedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfsNumNuclei =  NaN(LengthBinnedTbinnedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfsNumOnNuclei =  NaN(LengthBinnedTbinnedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfsNumOffNuclei =  NaN(LengthBinnedTbinnedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfsNumQuiescentNuclei =  NaN(LengthBinnedTbinnedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoProfsNumFinishedTranscribingNuclei =  NaN(LengthBinnedTbinnedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoNuclearProfiles =  NaN(LengthBinnedTbinnedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoFractionOns =  NaN(LengthBinnedTbinnedProfile, NumAPbins, length(EmbryoIndices));
        EmbryoTraceCounts =  NaN(LengthBinnedTbinnedProfile, NumAPbins, length(EmbryoIndices));
        for SetIndex = 1:length(EmbryoIndices)
            EmbryoIndex = EmbryoIndices(SetIndex);
            MatchingTimes = ismember(ltm.BinnedMeanProfiles.Tbinned3DCycleTimes, ltm.MeanProfiles{EmbryoIndex}.Tbinned3DCycleFrameTimes{nc_index});
            EmbryoNCProfile = ltm.MeanProfiles{EmbryoIndex}.Tbinned3DCycleMeanTraces(:,:,nc_index);
            EmbryoNCNuclearProfile = ltm.MeanProfiles{EmbryoIndex}.Tbinned3DCycleMeanPerNucleusTraces(:,:,nc_index);
            EmbryoNCProfileStd = ltm.MeanProfiles{EmbryoIndex}.Tbinned3DCycleTraceStdErrors(:,:,nc_index);
            EmbryoNCProfileFractionOn= ltm.MeanProfiles{EmbryoIndex}.Tbinned3DCycleFractionOn(:,:,nc_index);
            EmbryoNCProfileTraceCount = ltm.MeanProfiles{EmbryoIndex}.Tbinned3DCycleTraceCount(:,:,nc_index);
            EmbryoNCProfileTraceCount(isnan(EmbryoNCProfileTraceCount)) = 0;
            EmbryoNCProfileNumNuclei = ltm.MeanProfiles{EmbryoIndex}.Tbinned3DCycleNumNuclei(:,:,nc_index);
            EmbryoNCProfileNumNuclei(isnan(EmbryoNCProfileNumNuclei)) = 0;
            EmbryoNCProfileNumOnNuclei = ltm.MeanProfiles{EmbryoIndex}.Tbinned3DCycleNumOnNuclei(:,:,nc_index);
            EmbryoNCProfileNumOnNuclei(isnan(EmbryoNCProfileNumOnNuclei)) = 0;
            EmbryoNCProfileNumOffNuclei = ltm.MeanProfiles{EmbryoIndex}.Tbinned3DCycleNumOffNuclei(:,:,nc_index);
            EmbryoNCProfileNumOffNuclei(isnan(EmbryoNCProfileNumOffNuclei)) = 0;
            EmbryoNCProfileNumQuiescentNuclei = ltm.MeanProfiles{EmbryoIndex}.Tbinned3DCycleNumQuiescentNuclei(:,:,nc_index);
            EmbryoNCProfileNumQuiescentNuclei(isnan(EmbryoNCProfileNumQuiescentNuclei)) = 0;
            EmbryoNCProfileNumFinishedTranscribingNuclei = ltm.MeanProfiles{EmbryoIndex}.Tbinned3DCycleNumFinishedTranscribingNuclei(:,:,nc_index);
            EmbryoNCProfileNumFinishedTranscribingNuclei(isnan(EmbryoNCProfileNumFinishedTranscribingNuclei)) = 0; 
            EmbryoNCProfileFractionOn(EmbryoNCProfileNumNuclei == 0) = NaN;
            
            EmbryoProfiles(find(~all(isnan(EmbryoNCProfile), 2)),:,SetIndex) = ...
                EmbryoNCProfile(~all(isnan(EmbryoNCProfile), 2),:);
            EmbryoProfStds(find(~all(isnan(EmbryoNCProfile), 2)),:,SetIndex) = ...
                EmbryoNCProfileStd(~all(isnan(EmbryoNCProfile), 2),:);
            EmbryoNuclearProfiles(find(~all(isnan(EmbryoNCNuclearProfile), 2)),:,SetIndex) = ...
                EmbryoNCNuclearProfile(~all(isnan(EmbryoNCNuclearProfile), 2),:);
            EmbryoFractionOns(find(~all(isnan(EmbryoNCProfileFractionOn), 2)),:,SetIndex) = ...
                EmbryoNCProfileFractionOn(~all(isnan(EmbryoNCProfileFractionOn), 2),:);
            EmbryoTraceCounts(find(~all(isnan(EmbryoNCProfileTraceCount), 2)),:,SetIndex) = ...
                EmbryoNCProfileTraceCount(~all(isnan(EmbryoNCProfileTraceCount), 2),:);
            EmbryoProfsNumNuclei(find(~all(isnan(EmbryoNCProfileNumNuclei), 2)),:,SetIndex) = ...
                EmbryoNCProfileNumNuclei(~all(isnan(EmbryoNCProfileNumNuclei), 2),:);
            EmbryoProfsNumOnNuclei(find(~all(isnan(EmbryoNCProfileNumOnNuclei), 2)),:,SetIndex) = ...
                EmbryoNCProfileNumOnNuclei(~all(isnan(EmbryoNCProfileNumOnNuclei), 2),:);
            EmbryoProfsNumOffNuclei(find(~all(isnan(EmbryoNCProfileNumOffNuclei), 2)),:,SetIndex) = ...
                EmbryoNCProfileNumOffNuclei(~all(isnan(EmbryoNCProfileNumOffNuclei), 2),:);
            EmbryoProfsNumQuiescentNuclei(find(~all(isnan(EmbryoNCProfileNumQuiescentNuclei), 2)),:,SetIndex) = ...
                EmbryoNCProfileNumQuiescentNuclei(~all(isnan(EmbryoNCProfileNumQuiescentNuclei), 2),:);
            EmbryoProfsNumFinishedTranscribingNuclei(find(~all(isnan(EmbryoNCProfileNumFinishedTranscribingNuclei), 2)),:,SetIndex) = ...
                EmbryoNCProfileNumFinishedTranscribingNuclei(~all(isnan(EmbryoNCProfileNumFinishedTranscribingNuclei), 2),:);
        end
        
        BadEntries = EmbryoProfiles< ltm.MinimumSchnitzCount;
        EmbryoProfiles(EmbryoProfiles< ltm.MinimumSchnitzCount) = NaN; 
        TempMeanProfile = mean(EmbryoProfiles, 3, 'omitnan');
        TempStdError = std(EmbryoProfiles, 0, 3, 'omitnan');
        
        EmbryoNuclearProfiles(EmbryoNuclearProfiles< ltm.MinimumSchnitzCount) = NaN;
        TempNuclearMeanProfile = mean(EmbryoNuclearProfiles, 3, 'omitnan');
        TempNuclearStdError = std(EmbryoNuclearProfiles, 0, 3, 'omitnan');
        
        TempProfCount = sum(~isnan(EmbryoProfiles) & EmbryoProfiles >=  ltm.MinimumSchnitzCount , 3, 'omitnan');
        
        EmbryoProfsNumNuclei(BadEntries) = NaN;
        TempTotalNuclei = sum(EmbryoProfsNumNuclei, 3, 'omitnan');
        
        EmbryoProfsNumOnNuclei(BadEntries) = NaN;
        TempTotalOnNuclei = sum(EmbryoProfsNumOnNuclei, 3, 'omitnan');
        
        %EmbryoProfsNumOffNuclei(BadEntries) = NaN;
        TempTotalOffNuclei = sum(EmbryoProfsNumOffNuclei, 3, 'omitnan');
        
        %EmbryoProfsNumQuiescentNuclei(BadEntries) = NaN;
        TempTotalQuiescentNuclei = sum(EmbryoProfsNumQuiescentNuclei, 3, 'omitnan');
        
        %EmbryoProfsNumFinishedTranscribingNuclei(BadEntries) = NaN;
        TempTotalFinishedTranscribingNuclei = sum(EmbryoProfsNumFinishedTranscribingNuclei, 3, 'omitnan');
        
        
        TempMeanFractionOnNuclei = mean(EmbryoFractionOns, 3, 'omitnan');
        TempStdErrFractionOnNuclei = std(EmbryoFractionOns, 0, 3, 'omitnan');
        
        TempMeanProfile(TempMeanFractionOnNuclei == 0) = 0;
        TempStdError(TempMeanFractionOnNuclei == 0) = 0;
        
        ltm.BinnedMeanProfiles.Tbinned3DCycleMeanTraces(:,:,nc_index,t_index) = TempMeanProfile;
        ltm.BinnedMeanProfiles.Tbinned3DCycleNumEmbryos(:,:,nc_index,t_index) = TempProfCount;
        ltm.BinnedMeanProfiles.Tbinned3DCycleStdErrors(:,:,nc_index,t_index) = TempStdError;
        ltm.BinnedMeanProfiles.Tbinned3DCycleNuclearMeanTraces(:,:,nc_index,t_index) = TempNuclearMeanProfile;
        ltm.BinnedMeanProfiles.Tbinned3DCycleNuclearStdErrors(:,:,nc_index,t_index) = TempNuclearStdError;
        ltm.BinnedMeanProfiles.Tbinned3DCycleTotalNuclei(:,:,nc_index,t_index) = TempTotalNuclei;
        ltm.BinnedMeanProfiles.Tbinned3DCycleTotalOnNuclei(:,:,nc_index,t_index) = TempTotalOnNuclei;
        ltm.BinnedMeanProfiles.Tbinned3DCycleTotalOffNuclei(:,:,nc_index,t_index) = TempTotalOffNuclei;
        ltm.BinnedMeanProfiles.Tbinned3DCycleTotalQuiescentNuclei(:,:,nc_index,t_index) = TempTotalQuiescentNuclei;
        ltm.BinnedMeanProfiles.Tbinned3DCycleTotalFinishedTranscribingNuclei(:,:,nc_index,t_index) = TempTotalFinishedTranscribingNuclei;
        ltm.BinnedMeanProfiles.Tbinned3DCycleFractionOn(:,:,nc_index,t_index) = TempMeanFractionOnNuclei;
     
    end
    
end