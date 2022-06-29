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

%% Add Imin and Imax info  for SlideRescaledAvgAP.SmoothedDeltaFCFactors
% x = CompiledEmbryos.DubuisEmbryoTimes ;
% AllDeltaValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
%         CompiledEmbryos.TestSetEmbryos & ~isnan(x);% &...
% AllDeltaValidProfilesControlTF = CompiledEmbryos.IsNC14 &...
%         CompiledEmbryos.ControlSetEmbryos & ~isnan(x);% &...
% 
% CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).Test = {};
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