function CompiledEmbryos = AddUniversalScaledTestMasterFitRescaledProfiles(CompiledEmbryos, exp_index)
%% List of Profiles to rescale
AllSetInfo = GetFixedSetPrefixInfo;

AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';


BootstrappedProfilePath = 'S:/Gabriella/Dropbox/BootstrappedTestData/25CBootstrappedGaussianSmoothedProfiles.mat';
load(BootstrappedProfilePath, 'MeanSmoothedProfiles', 'SmoothedProfileSEs', 'xfits');
NumMasterProfs = length(MeanSmoothedProfiles);

SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];


[NumEmbryos, NumAPbins, NChannels] = size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles);
 FitTypeStrings = {'Rep1Rep2ComboFitSlideRescaledDorsalAvgAPProfiles', 'Rep1FlippedComboFitSlideRescaledDorsalAvgAPProfiles',...
        'Rep2FlippedComboFitSlideRescaledDorsalAvgAPProfiles', 'AllCombinedFitSlideRescaledDorsalAvgAPProfiles',...
        'Rep1Rep2ComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles', 'Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles',...
        'Rep2FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles', 'AllCombinedZeroCorrectedSlideRescaledDorsalAvgAPProfiles',...
        'Rep1Rep2ComboSlideRescaledDorsalAvgAPProfiles', 'Rep1FlippedComboSlideRescaledDorsalAvgAPProfiles',...
        'Rep2FlippedComboSlideRescaledDorsalAvgAPProfiles', 'AllCombinedSlideRescaledDorsalAvgAPProfiles',...
        'Rep1Rep2ComboHybridRescaledDorsalAvgAPProfiles', 'Rep1FlippedComboHybridRescaledDorsalAvgAPProfiles',...
        'Rep2FlippedComboHybridRescaledDorsalAvgAPProfiles', 'AllCombinedHybridRescaledDorsalAvgAPProfiles'...
        };
SetStrings = {'Rep1', 'Rep2', 'Flipped'};
ProfileStrings = {'TestSetFitSlideRescaledDorsalAvgAPProfiles', 'TestSetFitSlideRescaledDorsalAvgAPProfiles',...
    'TestSetFitSlideRescaledDorsalAvgAPProfiles', 'TestSetFitSlideRescaledDorsalAvgAPProfiles',...
    'TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles', 'TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles',...
    'TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles', 'TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles',...
    'TestSetSlideRescaledDorsalAvgAPProfiles', 'TestSetSlideRescaledDorsalAvgAPProfiles',...
    'TestSetSlideRescaledDorsalAvgAPProfiles', 'TestSetSlideRescaledDorsalAvgAPProfiles',...
    'HybridTestSetRescaledDorsalAvgAPProfiles', 'HybridTestSetRescaledDorsalAvgAPProfiles',...
    'HybridTestSetRescaledDorsalAvgAPProfiles', 'HybridTestSetRescaledDorsalAvgAPProfiles'...
    };
%% 
if ~isfield(CompiledEmbryos, 'UnivScaledProfiles')
CompiledEmbryos.UnivScaledProfiles = {};
end
for i = 1:length(FitTypeStrings)
        CompiledEmbryos.UnivScaledProfiles.(FitTypeStrings{i}) = {};%..
        CompiledEmbryos.UnivScaledProfiles.(FitTypeStrings{i}).x = CompiledEmbryos.DubuisEmbryoTimes;
    for j = 1:length(SetStrings)
        CompiledEmbryos.UnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}) = NaN(NumEmbryos, NumAPbins, NChannels, NumMasterProfs);
        for ch_index = [3 5]
            for master_index = 1:NumMasterProfs
                slope_control = CompiledEmbryos.ScaleFits.(FitTypeStrings{i}).(SetStrings{j}).ScaleEstimate(master_index, ch_index);
                intercept_control = CompiledEmbryos.ScaleFits.(FitTypeStrings{i}).(SetStrings{j}).InterceptEstimate(master_index, ch_index);
                CompiledEmbryos.UnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j})(:,:,ch_index,master_index) = ...
                    slope_control*CompiledEmbryos.UnivScaledProfiles.(ProfileStrings{i}).(SetStrings{j})(:,:,ch_index)+intercept_control;
            end
        end
    end
end


%%
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
save(CEoutpath, 'CompiledEmbryos');