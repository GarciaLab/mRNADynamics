function CompiledEmbryos = AddUniversalScalingFactorsToMasterProfiles(CompiledEmbryos, exp_index)
%% List of Profiles to rescale
AllSetInfo = GetFixedSetPrefixInfo;
AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';

SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];


NChannels = size(CompiledEmbryos.DorsalAvgAPProfiles, 3);

NumMasterProfs = size(CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Control.ScaleEstimate, 1);


SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];


[NumEmbryos, NumAPbins, NChannels] = size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles);

FitTypeStrings = {'SlideRescaledDorsalAvgAPProfiles', 'ZeroCorrectedSlideRescaledDorsalAvgAPProfiles',...
    'FitSlideRescaledDorsalAvgAPProfiles','FitZeroedSlideRescaledDorsalAvgAPProfiles',...
    'HybridRescaledDorsalAvgAPProfiles'};
SetStrings = {'Control', 'Test'};
ProfileStrings = {'SlideRescaledDorsalAvgAPProfiles', 'ZeroCorrectedSlideRescaledDorsalAvgAPProfiles',...
    'FitSlideRescaledDorsalAvgAPProfiles','FitZeroedSlideRescaledDorsalAvgAPProfiles',...
    'ZeroCorrectedSlideRescaledDorsalAvgAPProfiles'};

chLists = cell(1, 5);
chLists{1} = [3 5];
chLists{2} = [3 5];
chLists{3} = [3];
chLists{4} = [3];
chLists{5} = [3 5];

UseTestScaling = false;
if AllSetInfo.Temperatures(exp_index) == 25
    UseTestScaling = true;
end
%%
if ~isfield(CompiledEmbryos, 'UnivScaledProfiles')
CompiledEmbryos.UnivScaledProfiles = {};
end

for i = 1:length(FitTypeStrings)
    CompiledEmbryos.UnivScaledProfiles.(FitTypeStrings{i}) = {};%..
    CompiledEmbryos.UnivScaledProfiles.(FitTypeStrings{i}).x = CompiledEmbryos.DubuisEmbryoTimes;
    for j = 1:length(SetStrings)
        CompiledEmbryos.UnivScaledProfiles.(FitTypeStrings{i}).([SetStrings{j}, 'Scaling'])=  NaN(NumEmbryos, NumAPbins, NChannels, NumMasterProfs);

        if strcmpi(lower(SetStrings{j}), 'test') & ~UseTestScaling
            continue
        end
        for ch_index = chLists{i}
            for master_index = 1:NumMasterProfs
                slope = CompiledEmbryos.ScaleFits.(ProfileStrings{i}).(SetStrings{j}).ScaleEstimate(master_index, ch_index);
                intercept= CompiledEmbryos.ScaleFits.(ProfileStrings{i}).(SetStrings{j}).InterceptEstimate(master_index, ch_index);
                CompiledEmbryos.UnivScaledProfiles.(FitTypeStrings{i}).([SetStrings{j}, 'Scaling'])(:,:,ch_index,master_index) = ...
                    slope*CompiledEmbryos.(ProfileStrings{i})(:,:,ch_index)+intercept;
            end
        end
        
    end
end



%%
% CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
% save(CEoutpath, 'CompiledEmbryos');