function CompiledEmbryos = AddNC13NormalizedProfiles(CompiledEmbryos, exp_index)
%%

AllSetInfo = GetFixedSetPrefixInfo;
AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';

SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];

APbins = 0:0.025:1;
NumAPbins = length(APbins);
NChannels = size(CompiledEmbryos.DorsalAvgAPProfiles, 3);

%%
ValidTestProfilesTF = CompiledEmbryos.IsNC13 & CompiledEmbryos.TestSetEmbryos;% &...


CompiledEmbryos.NC13 = {};
CompiledEmbryos.NC13.NormalizedProfiles = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Profiles = ...
    CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Profiles(ValidTestProfilesTF,:,:);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.MeanProfile = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.MeanProfile.mean = ...
    squeeze(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.MeanProfile.std = ...
    std(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.MeanProfile.count = ...
    size(CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Profiles(ValidTestProfilesTF,:,:), 1);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P50Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P50Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P50Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P50Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P75Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P75Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P75Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P75Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P80Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P80Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P80Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P80Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P90Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P90Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P90Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P90Profiles.count = zeros(1, NChannels);

NumProfs = size(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Profiles, 1);
ProfileIntegrals = NaN(NumProfs, NChannels);
for ch_index = 2:NChannels
    for pr =1:NumProfs

        NaNTF = isnan(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Profiles(pr,:,ch_index));
        prof_good = CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Profiles(pr,~NaNTF,ch_index);
        ap_good = APbins(~NaNTF);
        if ~isempty(prof_good)
            prof_sums = prof_good(1:end-1)+prof_good(2:end);
            ap_diffs = diff(ap_good);
            ProfileIntegrals(pr, ch_index) = sum(prof_sums.*ap_diffs)/2-min(prof_good)*sum(ap_diffs);
            
        end
    end
    
    P = prctile(ProfileIntegrals(:,ch_index), [50, 75, 80, 90]);
    
    TF50 = ProfileIntegrals(:,ch_index) >= P(1);
    if sum(TF50) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P50Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P50Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P50Profiles.count(ch_index) = sum(TF50);
    end
    
    TF75 = ProfileIntegrals(:,ch_index) >= P(2);
    if sum(TF75) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P75Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P75Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P75Profiles.count(ch_index) = sum(TF75);
    end
    
    TF80 = ProfileIntegrals(:,ch_index) >= P(3);
    if sum(TF80) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P80Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P80Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P80Profiles.count(ch_index) = sum(TF80);
    end
    
    TF90 = ProfileIntegrals(:,ch_index) >= P(4);
    if sum(TF90) > 0 
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P90Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P90Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.P90Profiles.count(ch_index) = sum(TF90);
    end
    
end
%%
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.Profiles = ...
    CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Profiles(ValidTestProfilesTF,:,:);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.MeanProfile = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.MeanProfile.mean = ...
    squeeze(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.MeanProfile.std = ...
    std(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.MeanProfile.count = ...
    size(CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P50Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P50Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P50Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P50Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P75Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P75Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P75Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P75Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P80Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P80Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P80Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P80Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P90Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P90Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P90Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P90Profiles.count = zeros(1, NChannels);

NumProfs = size(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.Profiles, 1);
ProfileIntegrals = NaN(NumProfs, NChannels);
for ch_index = 2:NChannels
    for pr =1:NumProfs

        NaNTF = isnan(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.Profiles(pr,:,ch_index));
        prof_good = CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.Profiles(pr,~NaNTF,ch_index);
        ap_good = APbins(~NaNTF);
        if ~isempty(prof_good)
            prof_sums = prof_good(1:end-1)+prof_good(2:end);
            ap_diffs = diff(ap_good);
            ProfileIntegrals(pr, ch_index) = sum(prof_sums.*ap_diffs)/2-min(prof_good)*sum(ap_diffs);
            
        end
    end
    
    P = prctile(ProfileIntegrals(:,ch_index), [50, 75, 80, 90]);
    
    TF50 = ProfileIntegrals(:,ch_index) >= P(1);
    if sum(TF50) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P50Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P50Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P50Profiles.count(ch_index) = sum(TF50);
    end
    
     TF75 = ProfileIntegrals(:,ch_index) >= P(2);
    if sum(TF75) > 0 
   
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P75Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P75Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P75Profiles.count(ch_index) = sum(TF75);
    end
    
TF80 = ProfileIntegrals(:,ch_index) >= P(3);
    if sum(TF80) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P80Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P80Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P80Profiles.count(ch_index) = sum(TF80);
    end
    
    TF90 = ProfileIntegrals(:,ch_index) >= P(4);
    if sum(TF90) > 0 
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P90Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P90Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.P90Profiles.count(ch_index) = sum(TF90);
    end
    
end

%% SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.Profiles = ...
    CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Profiles(ValidTestProfilesTF,:,:);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.MeanProfile = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.MeanProfile.mean = ...
    squeeze(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.MeanProfile.std = ...
    std(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.MeanProfile.count = ...
    size(CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P50Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P50Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P50Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P50Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P75Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P75Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P75Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P75Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P80Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P80Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P80Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P80Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P90Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P90Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P90Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P90Profiles.count = zeros(1, NChannels);

NumProfs = size(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.Profiles, 1);
ProfileIntegrals = NaN(NumProfs, NChannels);
for ch_index = 2:NChannels
    for pr =1:NumProfs

        NaNTF = isnan(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.Profiles(pr,:,ch_index));
        prof_good = CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.Profiles(pr,~NaNTF,ch_index);
        ap_good = APbins(~NaNTF);
        if ~isempty(prof_good)
            prof_sums = prof_good(1:end-1)+prof_good(2:end);
            ap_diffs = diff(ap_good);
            ProfileIntegrals(pr, ch_index) = sum(prof_sums.*ap_diffs)/2-min(prof_good)*sum(ap_diffs);
            
        end
    end
    
    P = prctile(ProfileIntegrals(:,ch_index), [50, 75, 80, 90]);
    
    TF50 = ProfileIntegrals(:,ch_index) >= P(1);
    if sum(TF50) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P50Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P50Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P50Profiles.count(ch_index) = sum(TF50);
    end
    
    TF75 = ProfileIntegrals(:,ch_index) >= P(2);
    if sum(TF75) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P75Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P75Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P75Profiles.count(ch_index) = sum(TF75);
    end
    
TF80 = ProfileIntegrals(:,ch_index) >= P(3);
    if sum(TF80) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P80Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P80Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P80Profiles.count(ch_index) = sum(TF80);
    end
    
    TF90 = ProfileIntegrals(:,ch_index) >= P(4);
    if sum(TF90) > 0 
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P90Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P90Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.P90Profiles.count(ch_index) = sum(TF90);
    end
end
%% SlideRescaledAvgAP.WindowedDeltaFCFactors
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.Profiles = ...
    CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Profiles(ValidTestProfilesTF,:,:);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.MeanProfile = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.MeanProfile.mean = ...
    squeeze(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.MeanProfile.std = ...
    std(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.MeanProfile.count = ...
    size(CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Profiles(ValidTestProfilesTF,:,:), 1);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P50Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P50Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P50Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P50Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P75Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P75Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P75Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P75Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P80Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P80Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P80Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P80Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P90Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P90Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P90Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P90Profiles.count = zeros(1, NChannels);

NumProfs = size(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.Profiles, 1);
ProfileIntegrals = NaN(NumProfs, NChannels);
for ch_index = 2:NChannels
    for pr =1:NumProfs

        NaNTF = isnan(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.Profiles(pr,:,ch_index));
        prof_good = CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.Profiles(pr,~NaNTF,ch_index);
        ap_good = APbins(~NaNTF);
        if ~isempty(prof_good)
            prof_sums = prof_good(1:end-1)+prof_good(2:end);
            ap_diffs = diff(ap_good);
            ProfileIntegrals(pr, ch_index) = sum(prof_sums.*ap_diffs)/2-min(prof_good)*sum(ap_diffs);
            
        end
    end
    
    P = prctile(ProfileIntegrals(:,ch_index), [50, 75, 80, 90]);
    
    TF50 = ProfileIntegrals(:,ch_index) >= P(1);
    if sum(TF50) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P50Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P50Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P50Profiles.count(ch_index) = sum(TF50);
    end
    
    TF75 = ProfileIntegrals(:,ch_index) >= P(2);
    if sum(TF75) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P75Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P75Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P75Profiles.count(ch_index) = sum(TF75);
end
        TF80 = ProfileIntegrals(:,ch_index) >= P(3);
    if sum(TF80) > 0 

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P80Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P80Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P80Profiles.count(ch_index) = sum(TF80);
    end
    
    TF90 = ProfileIntegrals(:,ch_index) >= P(4);
    if sum(TF90) > 0 
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P90Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P90Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.P90Profiles.count(ch_index) = sum(TF90);
    end
end

%% SlideRescaledAvgAP.WindowedDubuisTimesFactors
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.Profiles = ...
    CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Profiles(ValidTestProfilesTF,:,:);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.MeanProfile = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.MeanProfile.mean = ...
    squeeze(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.MeanProfile.std = ...
    std(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.MeanProfile.count = ...
    size(CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P50Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P50Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P50Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P50Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P75Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P75Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P75Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P75Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P80Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P80Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P80Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P80Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P90Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P90Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P90Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P90Profiles.count = zeros(1, NChannels);

NumProfs = size(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.Profiles, 1);
ProfileIntegrals = NaN(NumProfs, NChannels);
for ch_index = 2:NChannels
    for pr =1:NumProfs

        NaNTF = isnan(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.Profiles(pr,:,ch_index));
        prof_good = CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.Profiles(pr,~NaNTF,ch_index);
        ap_good = APbins(~NaNTF);
        if ~isempty(prof_good)
            prof_sums = prof_good(1:end-1)+prof_good(2:end);
            ap_diffs = diff(ap_good);
            ProfileIntegrals(pr, ch_index) = sum(prof_sums.*ap_diffs)/2-min(prof_good)*sum(ap_diffs);
            
        end
    end
    
    P = prctile(ProfileIntegrals(:,ch_index), [50, 75, 80, 90]);
    
    TF50 = ProfileIntegrals(:,ch_index) >= P(1);
    if sum(TF50) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P50Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P50Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P50Profiles.count(ch_index) = sum(TF50);
    end
    
     TF75 = ProfileIntegrals(:,ch_index) >= P(2);
    if sum(TF75) > 0 
   
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P75Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P75Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P75Profiles.count(ch_index) = sum(TF75);
end
        TF80 = ProfileIntegrals(:,ch_index) >= P(3);
    if sum(TF80) > 0 

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P80Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P80Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P80Profiles.count(ch_index) = sum(TF80);
    end
    
    TF90 = ProfileIntegrals(:,ch_index) >= P(4);
    if sum(TF90) > 0 
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P90Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P90Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.P90Profiles.count(ch_index) = sum(TF90);
    end
end
    
    
    %% SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.Profiles = ...
    CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Profiles(ValidTestProfilesTF,:,:);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.MeanProfile = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.MeanProfile.mean = ...
    squeeze(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.MeanProfile.std = ...
    std(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.MeanProfile.count = ...
    size(CompiledEmbryos.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P50Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P50Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P50Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P50Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P75Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P75Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P75Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P75Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P80Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P80Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P80Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P80Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P90Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P90Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P90Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P90Profiles.count = zeros(1, NChannels);

NumProfs = size(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.Profiles, 1);
ProfileIntegrals = NaN(NumProfs, NChannels);
for ch_index = 2:NChannels
    for pr =1:NumProfs

        NaNTF = isnan(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.Profiles(pr,:,ch_index));
        prof_good = CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.Profiles(pr,~NaNTF,ch_index);
        ap_good = APbins(~NaNTF);
        if ~isempty(prof_good)
            prof_sums = prof_good(1:end-1)+prof_good(2:end);
            ap_diffs = diff(ap_good);
            ProfileIntegrals(pr, ch_index) = sum(prof_sums.*ap_diffs)/2-min(prof_good)*sum(ap_diffs);
            
        end
    end
    
    P = prctile(ProfileIntegrals(:,ch_index), [50, 75, 80, 90]);
    
        TF50 = ProfileIntegrals(:,ch_index) >= P(1);
    if sum(TF50) > 0 

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P50Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P50Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P50Profiles.count(ch_index) = sum(TF50);
    end
    
    TF75 = ProfileIntegrals(:,ch_index) >= P(2);
    if sum(TF75) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P75Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P75Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P75Profiles.count(ch_index) = sum(TF75);
    end
    
TF80 = ProfileIntegrals(:,ch_index) >= P(3);
    if sum(TF80) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P80Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P80Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P80Profiles.count(ch_index) = sum(TF80);
    end
    
    TF90 = ProfileIntegrals(:,ch_index) >= P(4);
    if sum(TF90) > 0 
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P90Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P90Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.P90Profiles.count(ch_index) = sum(TF90);
    end
end
    
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.Profiles = ...
    CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Profiles(ValidTestProfilesTF,:,:);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.MeanProfile = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.MeanProfile.mean = ...
    squeeze(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.MeanProfile.std = ...
    std(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.MeanProfile.count = ...
    size(CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Profiles(ValidTestProfilesTF,:,:), 1);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P50Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P50Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P50Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P50Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P75Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P75Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P75Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P75Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P80Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P80Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P80Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P80Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P90Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P90Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P90Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P90Profiles.count = zeros(1, NChannels);

NumProfs = size(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.Profiles, 1);
ProfileIntegrals = NaN(NumProfs, NChannels);
for ch_index = 2:NChannels
    for pr =1:NumProfs

        NaNTF = isnan(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.Profiles(pr,:,ch_index));
        prof_good = CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.Profiles(pr,~NaNTF,ch_index);
        ap_good = APbins(~NaNTF);
        if ~isempty(prof_good)
            prof_sums = prof_good(1:end-1)+prof_good(2:end);
            ap_diffs = diff(ap_good);
            ProfileIntegrals(pr, ch_index) = sum(prof_sums.*ap_diffs)/2-min(prof_good)*sum(ap_diffs);
            
        end
    end
    
    P = prctile(ProfileIntegrals(:,ch_index), [50, 75, 80, 90]);
    
    TF50 = ProfileIntegrals(:,ch_index) >= P(1);
    if sum(TF50) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P50Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P50Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P50Profiles.count(ch_index) = sum(TF50);
    end
    
    TF75 = ProfileIntegrals(:,ch_index) >= P(2);
    if sum(TF75) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P75Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P75Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P75Profiles.count(ch_index) = sum(TF75);
    end
    
 TF80 = ProfileIntegrals(:,ch_index) >= P(3);
    if sum(TF80) > 0 
   
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P80Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P80Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P80Profiles.count(ch_index) = sum(TF80);
    end
    
    TF90 = ProfileIntegrals(:,ch_index) >= P(4);
    if sum(TF90) > 0 
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P90Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P90Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.P90Profiles.count(ch_index) = sum(TF90);
    end
    
end
%%
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.Profiles = ...
    CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Profiles(ValidTestProfilesTF,:,:);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.MeanProfile = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.MeanProfile.mean = ...
    squeeze(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.MeanProfile.std = ...
    std(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.MeanProfile.count = ...
    size(CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P50Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P50Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P50Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P50Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P75Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P75Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P75Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P75Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P80Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P80Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P80Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P80Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P90Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P90Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P90Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P90Profiles.count = zeros(1, NChannels);

NumProfs = size(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.Profiles, 1);
ProfileIntegrals = NaN(NumProfs, NChannels);
for ch_index = 2:NChannels
    for pr =1:NumProfs

        NaNTF = isnan(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.Profiles(pr,:,ch_index));
        prof_good = CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.Profiles(pr,~NaNTF,ch_index);
        ap_good = APbins(~NaNTF);
        if ~isempty(prof_good)
            prof_sums = prof_good(1:end-1)+prof_good(2:end);
            ap_diffs = diff(ap_good);
            ProfileIntegrals(pr, ch_index) = sum(prof_sums.*ap_diffs)/2-min(prof_good)*sum(ap_diffs);
            
        end
    end
    
    P = prctile(ProfileIntegrals(:,ch_index), [50, 75, 80, 90]);
    
    TF50 = ProfileIntegrals(:,ch_index) >= P(1);
    if sum(TF50) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P50Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P50Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P50Profiles.count(ch_index) = sum(TF50);
    end
    
        TF75 = ProfileIntegrals(:,ch_index) >= P(2);
    if sum(TF75) > 0 

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P75Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P75Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P75Profiles.count(ch_index) = sum(TF75);
    end
    
TF80 = ProfileIntegrals(:,ch_index) >= P(3);
    if sum(TF80) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P80Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P80Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P80Profiles.count(ch_index) = sum(TF80);
    end
    
    TF90 = ProfileIntegrals(:,ch_index) >= P(4);
    if sum(TF90) > 0 
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P90Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P90Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.P90Profiles.count(ch_index) = sum(TF90);
    end
    
end

%% SlideRescaledAP.SmoothedHisRFP25CTimesFactors
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.Profiles = ...
    CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Profiles(ValidTestProfilesTF,:,:);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.MeanProfile = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.MeanProfile.mean = ...
    squeeze(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.MeanProfile.std = ...
    std(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.MeanProfile.count = ...
    size(CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P50Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P50Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P50Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P50Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P75Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P75Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P75Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P75Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P80Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P80Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P80Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P80Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P90Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P90Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P90Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P90Profiles.count = zeros(1, NChannels);

NumProfs = size(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.Profiles, 1);
ProfileIntegrals = NaN(NumProfs, NChannels);
for ch_index = 2:NChannels
    for pr =1:NumProfs

        NaNTF = isnan(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.Profiles(pr,:,ch_index));
        prof_good = CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.Profiles(pr,~NaNTF,ch_index);
        ap_good = APbins(~NaNTF);
        if ~isempty(prof_good)
            prof_sums = prof_good(1:end-1)+prof_good(2:end);
            ap_diffs = diff(ap_good);
            ProfileIntegrals(pr, ch_index) = sum(prof_sums.*ap_diffs)/2-min(prof_good)*sum(ap_diffs);
            
        end
    end
    
    P = prctile(ProfileIntegrals(:,ch_index), [50, 75, 80, 90]);
    
    TF50 = ProfileIntegrals(:,ch_index) >= P(1);
    if sum(TF50) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P50Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P50Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P50Profiles.count(ch_index) = sum(TF50);
    end
    
     TF75 = ProfileIntegrals(:,ch_index) >= P(2);
    if sum(TF75) > 0 
   
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P75Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P75Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P75Profiles.count(ch_index) = sum(TF75);
    end
    
    TF80 = ProfileIntegrals(:,ch_index) >= P(3);
    if sum(TF80) > 0 

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P80Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P80Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P80Profiles.count(ch_index) = sum(TF80);
    end
    
    TF90 = ProfileIntegrals(:,ch_index) >= P(4);
    if sum(TF90) > 0 
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P90Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P90Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.P90Profiles.count(ch_index) = sum(TF90);
    end
end
%% SlideRescaledAP.WindowedDeltaFCFactors
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.Profiles = ...
    CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Profiles(ValidTestProfilesTF,:,:);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.MeanProfile = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.MeanProfile.mean = ...
    squeeze(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.MeanProfile.std = ...
    std(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.MeanProfile.count = ...
    size(CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Profiles(ValidTestProfilesTF,:,:), 1);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P50Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P50Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P50Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P50Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P75Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P75Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P75Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P75Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P80Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P80Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P80Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P80Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P90Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P90Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P90Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P90Profiles.count = zeros(1, NChannels);

NumProfs = size(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.Profiles, 1);
ProfileIntegrals = NaN(NumProfs, NChannels);
for ch_index = 2:NChannels
    for pr =1:NumProfs

        NaNTF = isnan(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.Profiles(pr,:,ch_index));
        prof_good = CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.Profiles(pr,~NaNTF,ch_index);
        ap_good = APbins(~NaNTF);
        if ~isempty(prof_good)
            prof_sums = prof_good(1:end-1)+prof_good(2:end);
            ap_diffs = diff(ap_good);
            ProfileIntegrals(pr, ch_index) = sum(prof_sums.*ap_diffs)/2-min(prof_good)*sum(ap_diffs);
            
        end
    end
    
    P = prctile(ProfileIntegrals(:,ch_index), [50, 75, 80, 90]);
    
    TF50 = ProfileIntegrals(:,ch_index) >= P(1);
    if sum(TF50) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P50Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P50Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P50Profiles.count(ch_index) = sum(TF50);
    end
    
       TF75 = ProfileIntegrals(:,ch_index) >= P(2);
    if sum(TF75) > 0 
 
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P75Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P75Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P75Profiles.count(ch_index) = sum(TF75);
    end
    
 TF80 = ProfileIntegrals(:,ch_index) >= P(3);
    if sum(TF80) > 0 
   
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P80Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P80Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P80Profiles.count(ch_index) = sum(TF80);
    end
    
    TF90 = ProfileIntegrals(:,ch_index) >= P(4);
    if sum(TF90) > 0 
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P90Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P90Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.P90Profiles.count(ch_index) = sum(TF90);
    end
end

%% SlideRescaledAP.WindowedDubuisTimesFactors
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.Profiles = ...
    CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Profiles(ValidTestProfilesTF,:,:);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.MeanProfile = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.MeanProfile.mean = ...
    squeeze(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.MeanProfile.std = ...
    std(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.MeanProfile.count = ...
    size(CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P50Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P50Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P50Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P50Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P75Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P75Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P75Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P75Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P80Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P80Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P80Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P80Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P90Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P90Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P90Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P90Profiles.count = zeros(1, NChannels);

NumProfs = size(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.Profiles, 1);
ProfileIntegrals = NaN(NumProfs, NChannels);
for ch_index = 2:NChannels
    for pr =1:NumProfs

        NaNTF = isnan(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.Profiles(pr,:,ch_index));
        prof_good = CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.Profiles(pr,~NaNTF,ch_index);
        ap_good = APbins(~NaNTF);
        if ~isempty(prof_good)
            prof_sums = prof_good(1:end-1)+prof_good(2:end);
            ap_diffs = diff(ap_good);
            ProfileIntegrals(pr, ch_index) = sum(prof_sums.*ap_diffs)/2-min(prof_good)*sum(ap_diffs);
            
        end
    end
    
    P = prctile(ProfileIntegrals(:,ch_index), [50, 75, 80, 90]);
    
    TF50 = ProfileIntegrals(:,ch_index) >= P(1);
    if sum(TF50) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P50Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P50Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P50Profiles.count(ch_index) = sum(TF50);
    end
      TF75 = ProfileIntegrals(:,ch_index) >= P(2);
    
    if sum(TF75) > 0 
  
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P75Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P75Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P75Profiles.count(ch_index) = sum(TF75);
    end
    
TF80 = ProfileIntegrals(:,ch_index) >= P(3);
    if sum(TF80) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P80Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P80Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P80Profiles.count(ch_index) = sum(TF80);
    end
    
    TF90 = ProfileIntegrals(:,ch_index) >= P(4);
    if sum(TF90) > 0 
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P90Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P90Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.P90Profiles.count(ch_index) = sum(TF90);
    end
end
    
    
    %% SlideRescaledAP.WindowedHisRFP25CTimesFactors
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.Profiles = ...
    CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Profiles(ValidTestProfilesTF,:,:);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.MeanProfile = {};
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.MeanProfile.mean = ...
    squeeze(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.MeanProfile.std = ...
    std(mean(CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1, 'omitnan'));
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.MeanProfile.count = ...
    size(CompiledEmbryos.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Profiles(ValidTestProfilesTF,:,:), 1);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P50Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P50Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P50Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P50Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P75Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P75Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P75Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P75Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P80Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P80Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P80Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P80Profiles.count = zeros(1, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P90Profiles = {}; 
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P90Profiles.mean = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P90Profiles.std = NaN(NumAPbins, NChannels);
CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P90Profiles.count = zeros(1, NChannels);

NumProfs = size(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.Profiles, 1);
ProfileIntegrals = NaN(NumProfs, NChannels);
for ch_index = 2:NChannels
    for pr =1:NumProfs

        NaNTF = isnan(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.Profiles(pr,:,ch_index));
        prof_good = CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.Profiles(pr,~NaNTF,ch_index);
        ap_good = APbins(~NaNTF);
        if ~isempty(prof_good)
            prof_sums = prof_good(1:end-1)+prof_good(2:end);
            ap_diffs = diff(ap_good);
            ProfileIntegrals(pr, ch_index) = sum(prof_sums.*ap_diffs)/2-min(prof_good)*sum(ap_diffs);
            
        end
    end
    
    P = prctile(ProfileIntegrals(:,ch_index), [50, 75, 80, 90]);
    
    TF50 = ProfileIntegrals(:,ch_index) >= P(1);
    if sum(TF50) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P50Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P50Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.Profiles(TF50,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P50Profiles.count(ch_index) = sum(TF50);
    end
    
     TF75 = ProfileIntegrals(:,ch_index) >= P(2);
    if sum(TF75) > 0 
   
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P75Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P75Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.Profiles(TF75,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P75Profiles.count(ch_index) = sum(TF75);
end
    TF80 = ProfileIntegrals(:,ch_index) >= P(3);
    if sum(TF80) > 0 
    
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P80Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P80Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.Profiles(TF80,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P80Profiles.count(ch_index) = sum(TF80);
    end
    
    TF90 = ProfileIntegrals(:,ch_index) >= P(4);
    if sum(TF90) > 0 
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P90Profiles.mean(:,ch_index) =...
        squeeze(mean(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';

    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P90Profiles.std(:,ch_index) =...
        squeeze(std(CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.Profiles(TF90,:,ch_index), 1, 'omitnan')).';
    CompiledEmbryos.NC13.NormalizedProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.P90Profiles.count(ch_index) = sum(TF90);
    end
end
    
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
save(CEoutpath, 'CompiledEmbryos');