function CompiledEmbryos = AddMasterFitRescaledProfiles(CompiledEmbryos, exp_index)
%% List of Profiles to rescale
AllSetInfo = GetFixedSetPrefixInfo;

AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';

SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];


[NumEmbryos, NumAPbins, NChannels] = size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles);
NumMasterProfSets = size(CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Control.ScaleEstimate, 1);
%% 
CompiledEmbryos.UnivScaledProfiles = {};
CompiledEmbryos.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.x = CompiledEmbryos.DubuisEmbryoTimes;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling = NaN(NumEmbryos, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.TestScaling = NaN(NumEmbryos, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.UnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.UnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.x = CompiledEmbryos.DubuisEmbryoTimes;
CompiledEmbryos.UnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling = NaN(NumEmbryos, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.UnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling = NaN(NumEmbryos, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.UnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.UnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.x = CompiledEmbryos.DubuisEmbryoTimes;
CompiledEmbryos.UnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling = NaN(NumEmbryos, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.UnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling = NaN(NumEmbryos, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.UnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.UnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.x = CompiledEmbryos.DubuisEmbryoTimes;
CompiledEmbryos.UnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.ControlScaling = NaN(NumEmbryos, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.UnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.TestScaling = NaN(NumEmbryos, NumAPbins, NChannels, NumMasterProfSets);


%% First Add Slide Rescaling Fits of SlideRescaledDorsalAvgAPProfiles
for set_index = 1:NumMasterProfSets
    for ch_index = [3 5]
        slope_control = CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Control.ScaleEstimate(set_index, ch_index);
        intercept_control = CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Control.InterceptEstimate(set_index, ch_index);
        CompiledEmbryos.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(:,:,ch_index, set_index) = ...
            slope_control*CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_control;
        
        if AllSetInfo.Temperatures(exp_index) == 25
            slope_test = CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Test.ScaleEstimate(set_index, ch_index);
            intercept_test = CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Test.InterceptEstimate(set_index, ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.TestScaling(:,:,ch_index, set_index) = ...
                slope_test*CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_test;
        end
        
    end
end


%% Then Add Zero Corrected Fits of ZeroCorrectedSlideRescaledDorsalAvgAPProfiles
for set_index = 1:NumMasterProfSets
    for ch_index = [3 5]
        slope_control = CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Control.ScaleEstimate(set_index, ch_index);
        intercept_control = CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Control.InterceptEstimate(set_index, ch_index);
        CompiledEmbryos.UnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling(:,:,ch_index, set_index) = ...
            slope_control*CompiledEmbryos.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_control;
        
        if AllSetInfo.Temperatures(exp_index) == 25
            slope_test = CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Test.ScaleEstimate(set_index, ch_index);
            intercept_test = CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Test.InterceptEstimate(set_index, ch_index);
            CompiledEmbryos.UnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling(:,:,ch_index, set_index) = ...
                slope_test*CompiledEmbryos.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_test;
        end
        
    end
end

%% Then add FitSlideRescaledDorsalAvgAPProfiles of FitSlideRescaledDorsalAvgAPProfiles
for set_index = 1:NumMasterProfSets
    for ch_index = [3 5]
        slope_control = CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Control.ScaleEstimate(set_index, ch_index);
        intercept_control = CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Control.InterceptEstimate(set_index, ch_index);
        CompiledEmbryos.UnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling(:,:,ch_index, set_index) = ...
            slope_control*CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_control;
        
        if AllSetInfo.Temperatures(exp_index) == 25
            slope_test = CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Test.ScaleEstimate(set_index, ch_index);
            intercept_test = CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Test.InterceptEstimate(set_index, ch_index);
            CompiledEmbryos.UnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling(:,:,ch_index, set_index) = ...
                slope_test*CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_test;
        end
        
    end
end

%% Then add FitSlideRescaledDorsalAvgAPProfiles of SlideRescaledDorsalAvgAPProfiles
for set_index = 1:NumMasterProfSets
    for ch_index = [3 5]
        slope_control = CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Control.ScaleEstimate(set_index, ch_index);
        intercept_control = CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Control.InterceptEstimate(set_index, ch_index);
        CompiledEmbryos.UnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.ControlScaling(:,:,ch_index, set_index) = ...
            slope_control*CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_control;
        
        if AllSetInfo.Temperatures(exp_index) == 25
            slope_test = CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Test.ScaleEstimate(set_index, ch_index);
            intercept_test = CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Test.InterceptEstimate(set_index, ch_index);
            CompiledEmbryos.UnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.TestScaling(:,:,ch_index, set_index) = ...
                slope_test*CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + intercept_test;
        end
        
    end
end
%% Add Imin and Imax info  for SlideRescaledAvgAP.SmoothedDeltaFCFactors
% x = CompiledEmbryos.DubuisEmbryoTimes ;
% AllDeltaValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
%         CompiledEmbryos.TestSetEmbryos & ~isnan(x);% &...
% AllDeltaValidProfilesControlTF = CompiledEmbryos.IsNC14 &...
%         CompiledEmbryos.ControlSetEmbryos & ~isnan(x);% &...
% % 
% % CompiledEmbryos.UnivScaledProfiles.FitSlideRescaledAvgAP.(SetString).Test = {};
% x_test = x(AllDeltaValidProfilesTestTF);
% TestProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Profiles(AllDeltaValidProfilesTestTF,:,:);
% [x_test, test_order] = sort(x_test);
% TestProfiles = TestProfiles(test_order,:,:);
% CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.x = x_test;
% CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Profiles = TestProfiles;
% MTest = movmean(TestProfiles,10, 1, 'omitnan');
% MTest = MTest(6:end-5, :,:);
% CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Imin = NaN(1, NChannels); 
% CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Imax = NaN(1, NChannels); 
% for ch_index = 2:NChannels
%     MTestSub = MTest(:,:,ch_index);
%     CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Imin(ch_index) = min(MTestSub(:));
%     CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Imax(ch_index) = max(MTestSub(:));
% end
% 
% 
% 
% CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control = {};
% x_control = x(AllDeltaValidProfilesControlTF);
% ControlProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Profiles(AllDeltaValidProfilesControlTF,:,:);
% [x_control, control_order] = sort(x_control);
% 
% CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control.Imin = NaN(1, NChannels);
% CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control.Imax = NaN(1, NChannels);
% 
% ControlProfiles = ControlProfiles(control_order,:,:);
% if ~isempty(ControlProfiles)
%     CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control.x = x_control;
%     CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control.Profiles = ControlProfiles;
%     MControl = movmean(ControlProfiles,10, 1, 'omitnan');
%     
%     if size(MControl,1) >= 10
%         MControl = MControl(6:end-4, :,:);
%         
%         for ch_index = 2:NChannels
%             MControlSub = MControl(:,:,ch_index);
%             CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control.Imin(ch_index) = min(MControlSub(:));
%             CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control.Imax(ch_index) = max(MControlSub(:));
%         end
%     elseif mod(size(MControl,1),2) == 0
%         MControl = MControl((size(MControl,1)/2+1):(size(MControl,1)/2+1), :,:);
%         for ch_index = 2:NChannels
%             MControlSub = MControl(:,:,ch_index);
%             CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control.Imin(ch_index) = min(MControlSub(:));
%             CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control.Imax(ch_index) = max(MControlSub(:));
%         end
%     else
%         MControl = MControl(ceil(size(MControl,1)/2):ceil(size(MControl,1)/2), :,:);
%         for ch_index = 2:NChannels
%             MControlSub = MControl(:,:,ch_index);
%             CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control.Imin(ch_index) = min(MControlSub(:));
%             CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control.Imax(ch_index) = max(MControlSub(:));
%         end
%     end
% end




%%
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
save(CEoutpath, 'CompiledEmbryos');