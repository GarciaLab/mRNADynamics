function CompiledEmbryos = AddUniversalScaleNC13Profiles(CompiledEmbryos, exp_index)
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
if ~isfield(CompiledEmbryos, 'NC13')
CompiledEmbryos.NC13 = {};
end
if ~isfield(CompiledEmbryos.NC13, 'UnivScaledProfiles')
CompiledEmbryos.NC13.UnivScaledProfiles = {};
end

for i = 1:length(FitTypeStrings)
    CompiledEmbryos.NC13.UnivScaledProfiles.(FitTypeStrings{i}) = {};%..
    for j = 1:length(SetStrings)
        if strcmpi(lower(SetStrings{j}), 'test') & ~UseTestScaling
            continue
        end
        TFTestNC13 = CompiledEmbryos.IsNC13 & CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.Approved;
        CompiledEmbryos.NC13.UnivScaledProfiles.(FitTypeStrings{i}).([SetStrings{j}, 'Scaling']).Test = CompiledEmbryos.UnivScaledProfiles.(FitTypeStrings{i}).([SetStrings{j}, 'Scaling'])(TFTestNC13,:,:,:) ;

    end
end



%%
% CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
% save(CEoutpath, 'CompiledEmbryos');