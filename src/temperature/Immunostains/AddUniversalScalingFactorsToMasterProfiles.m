function CompiledEmbryos = AddUniversalScalingFactorsToMasterProfiles(CompiledEmbryos, exp_index)
%% List of Profiles to rescale


AllSetInfo = GetFixedSetPrefixInfo;
AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';

SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];


NChannels = size(CompiledEmbryos.DorsalAvgAPProfiles, 3);

%%
CompiledEmbryos.UnivBootstrappedScaledProfiles = {};

CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP = {};
%CompiledEmbryos.UnivBootstrappedScaledProfiles.ZeroCorrectedSlideRescaledAvgAP = {};
for i = 1:size(CompiledEmbryos.BootstrappedScaleFactors,1)
    SetString = ['Set', num2str(i)];
    CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).mean = NaN(size(CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles));
    %CompiledEmbryos.UnivBootstrappedScaledProfiles.ZeroCorrectedSlideRescaledAvgAP.(SetString).mean = NaN(size(CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles));
    for ch_index = 2:NChannels
        CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).mean(:,:,ch_index) = ...
            CompiledEmbryos.BootstrappedScaleFactors(i, ch_index)*CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + CompiledEmbryos.BootstrappedScaleIntercepts(i, ch_index);
%         CompiledEmbryos.UnivBootstrappedScaledProfiles.ZeroCorrectedSlideRescaledAvgAP.(SetString).mean (:,:,ch_index) = ...
%             CompiledEmbryos.BootstrappedScaleFactors(i, ch_index)*CompiledEmbryos.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + CompiledEmbryos.BootstrappedScaleIntercepts(i, ch_index);
    end


end




%%
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
save(CEoutpath, 'CompiledEmbryos');