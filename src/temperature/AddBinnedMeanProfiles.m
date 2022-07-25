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
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleMeanTraces = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleStdErrors = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleNuclearMeanTraces = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleNuclearStdErrors = NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleTotalNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleTotalOnNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleTotalOffNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleTotalQuiescentNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);
ltm.BinnedMeanProfiles.AnaphaseAlignedCycleTotalFinishedTranscribingNuclei= NaN(LengthBinnedAnaphaseAlignedProfile, NumAPbins, 6, 5);

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

%% First make average profiles binning anaphase aligned traces for each temperature 
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
            MatchingTimes = ismember(ltm.BinnedMeanProfiles.AnaphaseAlignedCycleTimes, ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleFrameTimes{nc_index});
            EmbryoNCProfile = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleMeanTraces(:,:,nc_index);
            EmbryoNCNuclearProfile = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleMeanPerNucleusTraces(:,:,nc_index);
            EmbryoNCProfileStd = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleTraceStdErrors(:,:,nc_index);
            EmbryoNCProfileFractionOn= ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleFractionOn(:,:,nc_index);
            EmbryoNCProfileTraceCount = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleTraceCount(:,:,nc_index);
            EmbryoNCProfileTraceCount(isnan(EmbryoNCProfileTraceCount)) = 0;
            EmbryoNCProfileNumNuclei = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleNumNuclei(:,:,nc_index);
            EmbryoNCProfileNumNuclei(isnan(EmbryoNCProfileNumNuclei)) = 0;
            EmbryoNCProfileNumOnNuclei = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleNumOnNuclei(:,:,nc_index);
            EmbryoNCProfileNumOnNuclei(isnan(EmbryoNCProfileNumOnNuclei)) = 0;
            EmbryoNCProfileNumOffNuclei = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleNumOffNuclei(:,:,nc_index);
            EmbryoNCProfileNumOffNuclei(isnan(EmbryoNCProfileNumOffNuclei)) = 0;
            EmbryoNCProfileNumQuiescentNuclei = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleNumQuiescentNuclei(:,:,nc_index);
            EmbryoNCProfileNumQuiescentNuclei(isnan(EmbryoNCProfileNumQuiescentNuclei)) = 0;
            EmbryoNCProfileNumFinishedTranscribingNuclei = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleNumFinishedTranscribingNuclei(:,:,nc_index);
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
        
        TempMeanProfile = mean(EmbryoProfiles, 3, 'omitnan');
        TempStdError = std(EmbryoProfiles, 0, 3, 'omitnan');
        TempNuclearMeanProfile = mean(EmbryoNuclearProfiles, 3, 'omitnan');
        TempNuclearStdError = std(EmbryoNuclearProfiles, 0, 3, 'omitnan');
        TempProfCount = sum(~isnan(EmbryoProfiles) & EmbryoProfiles > 0 , 3);
        TempTotalNuclei = sum(EmbryoProfsNumNuclei, 3);
        TempTotalOnNuclei = sum(EmbryoProfsNumOnNuclei, 3);
        TempTotalOffNuclei = sum(EmbryoProfsNumOffNuclei, 3);
        TempTotalQuiescentNuclei = sum(EmbryoProfsNumQuiescentNuclei, 3);
        TempTotalFinishedTranscribingNuclei = sum(EmbryoProfsNumFinishedTranscribingNuclei, 3);
        TempMeanFractionOnNuclei = mean(EmbryoFractionOns, 3, 'omitnan');
        
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
        
        TempMeanProfile = mean(EmbryoProfiles, 3, 'omitnan');
        TempStdError = std(EmbryoProfiles, 0, 3, 'omitnan');
        TempNuclearMeanProfile = mean(EmbryoNuclearProfiles, 3, 'omitnan');
        TempNuclearStdError = std(EmbryoNuclearProfiles, 0, 3, 'omitnan');
        TempProfCount = sum(~isnan(EmbryoProfiles) & EmbryoProfiles > 0 , 3);
        TempTotalNuclei = sum(EmbryoProfsNumNuclei, 3);
        TempTotalOnNuclei = sum(EmbryoProfsNumOnNuclei, 3);
        TempTotalOffNuclei = sum(EmbryoProfsNumOffNuclei, 3);
        TempTotalQuiescentNuclei = sum(EmbryoProfsNumQuiescentNuclei, 3);
        TempTotalFinishedTranscribingNuclei = sum(EmbryoProfsNumFinishedTranscribingNuclei, 3);
        TempMeanFractionOnNuclei = mean(EmbryoFractionOns, 3, 'omitnan');
        
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
   
%% First make average profiles binning tbinned traces for each temperature 
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
            MatchingTimes = ismember(ltm.BinnedMeanProfiles.TbinnedCycleTimes, ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleFrameTimes{nc_index});
            EmbryoNCProfile = ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleMeanTraces(:,:,nc_index);
            EmbryoNCNuclearProfile = ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleMeanPerNucleusTraces(:,:,nc_index);
            EmbryoNCProfileStd = ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleTraceStdErrors(:,:,nc_index);
            EmbryoNCProfileFractionOn= ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleFractionOn(:,:,nc_index);
            EmbryoNCProfileTraceCount = ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleTraceCount(:,:,nc_index);
            EmbryoNCProfileTraceCount(isnan(EmbryoNCProfileTraceCount)) = 0;
            EmbryoNCProfileNumNuclei = ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleNumNuclei(:,:,nc_index);
            EmbryoNCProfileNumNuclei(isnan(EmbryoNCProfileNumNuclei)) = 0;
            EmbryoNCProfileNumOnNuclei = ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleNumOnNuclei(:,:,nc_index);
            EmbryoNCProfileNumOnNuclei(isnan(EmbryoNCProfileNumOnNuclei)) = 0;
            EmbryoNCProfileNumOffNuclei = ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleNumOffNuclei(:,:,nc_index);
            EmbryoNCProfileNumOffNuclei(isnan(EmbryoNCProfileNumOffNuclei)) = 0;
            EmbryoNCProfileNumQuiescentNuclei = ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleNumQuiescentNuclei(:,:,nc_index);
            EmbryoNCProfileNumQuiescentNuclei(isnan(EmbryoNCProfileNumQuiescentNuclei)) = 0;
            EmbryoNCProfileNumFinishedTranscribingNuclei = ltm.MeanProfiles{EmbryoIndex}.TbinnedCycleNumFinishedTranscribingNuclei(:,:,nc_index);
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
        
        TempMeanProfile = mean(EmbryoProfiles, 3, 'omitnan');
        TempStdError = std(EmbryoProfiles, 0, 3, 'omitnan');
        TempNuclearMeanProfile = mean(EmbryoNuclearProfiles, 3, 'omitnan');
        TempNuclearStdError = std(EmbryoNuclearProfiles, 0, 3, 'omitnan');
        TempProfCount = sum(~isnan(EmbryoProfiles) & EmbryoProfiles > 0 , 3);
        TempTotalNuclei = sum(EmbryoProfsNumNuclei, 3);
        TempTotalOnNuclei = sum(EmbryoProfsNumOnNuclei, 3);
        TempTotalOffNuclei = sum(EmbryoProfsNumOffNuclei, 3);
        TempTotalQuiescentNuclei = sum(EmbryoProfsNumQuiescentNuclei, 3);
        TempTotalFinishedTranscribingNuclei = sum(EmbryoProfsNumFinishedTranscribingNuclei, 3);
        TempMeanFractionOnNuclei = mean(EmbryoFractionOns, 3, 'omitnan');
        
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
        
        TempMeanProfile = mean(EmbryoProfiles, 3, 'omitnan');
        TempStdError = std(EmbryoProfiles, 0, 3, 'omitnan');
        TempNuclearMeanProfile = mean(EmbryoNuclearProfiles, 3, 'omitnan');
        TempNuclearStdError = std(EmbryoNuclearProfiles, 0, 3, 'omitnan');
        TempProfCount = sum(~isnan(EmbryoProfiles) & EmbryoProfiles > 0 , 3);
        TempTotalNuclei = sum(EmbryoProfsNumNuclei, 3);
        TempTotalOnNuclei = sum(EmbryoProfsNumOnNuclei, 3);
        TempTotalOffNuclei = sum(EmbryoProfsNumOffNuclei, 3);
        TempTotalQuiescentNuclei = sum(EmbryoProfsNumQuiescentNuclei, 3);
        TempTotalFinishedTranscribingNuclei = sum(EmbryoProfsNumFinishedTranscribingNuclei, 3);
        TempMeanFractionOnNuclei = mean(EmbryoFractionOns, 3, 'omitnan');
        
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