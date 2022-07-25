function CompiledEmbryos = AddTestMasterFitRescaledProfiles(CompiledEmbryos, exp_index)
%% List of Profiles to rescale
AllSetInfo = GetFixedSetPrefixInfo;

AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';

SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];


[NumEmbryos, NumAPbins, NChannels] = size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles);

%% 
if ~isfield(CompiledEmbryos, 'UnivScaledProfiles')
CompiledEmbryos.UnivScaledProfiles = {};
end
CompiledEmbryos.UnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.UnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.x = CompiledEmbryos.DubuisEmbryoTimes;
CompiledEmbryos.UnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1 = NaN(NumEmbryos, NumAPbins, NChannels);
CompiledEmbryos.UnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2 = NaN(NumEmbryos, NumAPbins, NChannels);
CompiledEmbryos.UnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped = NaN(NumEmbryos, NumAPbins, NChannels);

CompiledEmbryos.UnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.UnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.x = CompiledEmbryos.DubuisEmbryoTimes;
CompiledEmbryos.UnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1 = NaN(NumEmbryos, NumAPbins, NChannels);
CompiledEmbryos.UnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2 = NaN(NumEmbryos, NumAPbins, NChannels);
CompiledEmbryos.UnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped = NaN(NumEmbryos, NumAPbins, NChannels);

CompiledEmbryos.UnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.UnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.x = CompiledEmbryos.DubuisEmbryoTimes;
CompiledEmbryos.UnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1 = NaN(NumEmbryos, NumAPbins, NChannels);
CompiledEmbryos.UnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2 = NaN(NumEmbryos, NumAPbins, NChannels);
CompiledEmbryos.UnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped = NaN(NumEmbryos, NumAPbins, NChannels);

CompiledEmbryos.UnivScaledProfiles.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.UnivScaledProfiles.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.x = CompiledEmbryos.DubuisEmbryoTimes;
CompiledEmbryos.UnivScaledProfiles.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep1 = NaN(NumEmbryos, NumAPbins, NChannels);
CompiledEmbryos.UnivScaledProfiles.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep2 = NaN(NumEmbryos, NumAPbins, NChannels);
CompiledEmbryos.UnivScaledProfiles.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Flipped = NaN(NumEmbryos, NumAPbins, NChannels);

CompiledEmbryos.UnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.UnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.x = CompiledEmbryos.DubuisEmbryoTimes;
CompiledEmbryos.UnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep1 = NaN(NumEmbryos, NumAPbins, NChannels);
CompiledEmbryos.UnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep2 = NaN(NumEmbryos, NumAPbins, NChannels);
CompiledEmbryos.UnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Flipped = NaN(NumEmbryos, NumAPbins, NChannels);

%% First Add Slide Rescaling Fits of SlideRescaledDorsalAvgAPProfiles

for ch_index = [3 5]
    slope_control = CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1.ScaleEstimate(ch_index);
    intercept_control = CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1.InterceptEstimate(ch_index);
    CompiledEmbryos.UnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1(:,:,ch_index) = ...
        slope_control*CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_control;
    
    slope_control = CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2.ScaleEstimate(ch_index);
    intercept_control = CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2.InterceptEstimate(ch_index);
    CompiledEmbryos.UnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2(:,:,ch_index) = ...
        slope_control*CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_control;
    
    slope_control = CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped.ScaleEstimate(ch_index);
    intercept_control = CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped.InterceptEstimate(ch_index);
    CompiledEmbryos.UnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped(:,:,ch_index) = ...
        slope_control*CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_control;
end



%% Then Add Zero Corrected Fits of ZeroCorrectedSlideRescaledDorsalAvgAPProfiles
for ch_index = [3 5]
    slope_control = CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1.ScaleEstimate(ch_index);
    intercept_control = CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1.InterceptEstimate(ch_index);
    CompiledEmbryos.UnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1(:,:,ch_index) = ...
        slope_control*CompiledEmbryos.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_control;
    
    slope_control = CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.ScaleEstimate(ch_index);
    intercept_control = CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.InterceptEstimate(ch_index);
    CompiledEmbryos.UnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2(:,:,ch_index) = ...
        slope_control*CompiledEmbryos.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_control;
    
    slope_control = CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.ScaleEstimate(ch_index);
    intercept_control = CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.InterceptEstimate(ch_index);
    CompiledEmbryos.UnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped(:,:,ch_index) = ...
        slope_control*CompiledEmbryos.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_control;
end

%% Then add FitSlideRescaledDorsalAvgAPProfiles of FitSlideRescaledDorsalAvgAPProfiles
for ch_index = 3
    slope_control = CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1.ScaleEstimate(ch_index);
    intercept_control = CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1.InterceptEstimate(ch_index);
    CompiledEmbryos.UnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1(:,:,ch_index) = ...
        slope_control*CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_control;
    
    slope_control = CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2.ScaleEstimate(ch_index);
    intercept_control = CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2.InterceptEstimate(ch_index);
    CompiledEmbryos.UnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2(:,:,ch_index) = ...
        slope_control*CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_control;
    
    slope_control = CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped.ScaleEstimate(ch_index);
    intercept_control = CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped.InterceptEstimate(ch_index);
    CompiledEmbryos.UnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped(:,:,ch_index) = ...
        slope_control*CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_control;
end


%% Then add FitZeroedSlideRescaledDorsalAvgAPProfiles of FitZeroedSlideRescaledDorsalAvgAPProfiles
for ch_index = 3
    slope_control = CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep1.ScaleEstimate(ch_index);
    intercept_control = CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep1.InterceptEstimate(ch_index);
    CompiledEmbryos.UnivScaledProfiles.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep1(:,:,ch_index) = ...
        slope_control*CompiledEmbryos.FitZeroedSlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_control;
    
    slope_control = CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep2.ScaleEstimate(ch_index);
    intercept_control = CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep2.InterceptEstimate(ch_index);
    CompiledEmbryos.UnivScaledProfiles.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep2(:,:,ch_index) = ...
        slope_control*CompiledEmbryos.FitZeroedSlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_control;
    
    slope_control = CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Flipped.ScaleEstimate(ch_index);
    intercept_control = CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Flipped.InterceptEstimate(ch_index);
    CompiledEmbryos.UnivScaledProfiles.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Flipped(:,:,ch_index) = ...
        slope_control*CompiledEmbryos.FitZeroedSlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_control;
end
%% Then add FitZeroedSlideRescaledDorsalAvgAPProfiles of ZeroCorrectedSlideRescaledDorsalAvgAPProfiles
for ch_index = 3
    slope_control = CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep1.ScaleEstimate(ch_index);
    intercept_control = CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep1.InterceptEstimate(ch_index);
    CompiledEmbryos.UnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep1(:,:,ch_index) = ...
        slope_control*CompiledEmbryos.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_control;
    
    slope_control = CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep2.ScaleEstimate(ch_index);
    intercept_control = CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep2.InterceptEstimate(ch_index);
    CompiledEmbryos.UnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep2(:,:,ch_index) = ...
        slope_control*CompiledEmbryos.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_control;
    
    slope_control = CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Flipped.ScaleEstimate(ch_index);
    intercept_control = CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Flipped.InterceptEstimate(ch_index);
    CompiledEmbryos.UnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Flipped(:,:,ch_index) = ...
        slope_control*CompiledEmbryos.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_control;
end

%%
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
save(CEoutpath, 'CompiledEmbryos');