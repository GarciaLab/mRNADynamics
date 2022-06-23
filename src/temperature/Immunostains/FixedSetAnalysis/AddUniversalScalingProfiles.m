function CompiledEmbryos = AddUniversalScalingProfiles(CompiledEmbryos, exp_index)
%% List of Profiles to rescale
% SlideRescaledDorsalAvgAPProfiles
% SlideRescaledDorsalAvgNarrowAPProfiles
% SlideRescaledDorsalAPProfiles
% SlideRescaledDorsalNarrowAPProfiles
% --> The 4 above with FixCorrectedControlScalingFactors,
% DubuisTimesControlScalingFactors, HisRFP25CTimesControlScalingFactors,
% WindowedDubuisTimesControlScalingFactors,
% WindowedDeltaFCControlScalingFactors or WindowedHisRFP25CTimesControlScalingFactors

% FixCorrectedSmoothedAvgAPProfiles with FixCorrectedControlScalingFactors
% FixCorrectedSmoothedAvgNarrowAPProfiles
% DubuisTimeSmoothedAvgAPProfiles with DubuisTimesControlScalingFactors
% DubuisTimeSmoothedAvgNarrowAPProfiles
% yw25CTimeSmoothedAvgAPProfiles
% yw25CTimeSmoothedAvgNarrowAPProfiles
% HisRFP25CTimeSmoothedAvgAPProfiles with HisRFP25CTimesControlScalingFactors
% HisRFP25CTimeSmoothedAvgNarrowAPProfiles
% BinnedProfiles with DON'T HAVE THE RIGHT STUFF CALCULATED
% WindowedProfiles with WindowedDubuisTimesControlScalingFactors or
% WindowedDeltaFCControlScalingFactors or WindowedHisRFP25CTimesControlScalingFactors

%%


AllSetInfo = GetFixedSetPrefixInfo;
AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';

SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];


NChannels = size(CompiledEmbryos.DorsalAvgAPProfiles, 3);

%%
CompiledEmbryos.UnivScaledProfiles = {};

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP = {};
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors = {};
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.x = CompiledEmbryos.FixCorrectedDeltaFC_um.mean;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Profiles = ...
    NaN(size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles));
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors = {};
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.x = CompiledEmbryos.DubuisEmbryoTimes;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Profiles = ...
    NaN(size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles));
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors = {};
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.x = CompiledEmbryos.HisRFP25CEmbryoTimes;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Profiles = ...
    NaN(size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles));

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors = {};
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.x = CompiledEmbryos.FixCorrectedDeltaFC_um.mean;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Profiles = ...
    NaN(size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles));
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors = {};
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.x = CompiledEmbryos.DubuisEmbryoTimes;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Profiles = ...
    NaN(size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles));
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors = {};
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.x = CompiledEmbryos.HisRFP25CEmbryoTimes;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Profiles = ...
    NaN(size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles));

for ch_index = 2:NChannels
   CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Profiles(:,:,ch_index) = ...
       CompiledEmbryos.FixCorrectedControlScalingFactors(ch_index)*CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(:,:,ch_index);
   CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Profiles(:,:,ch_index) = ...
       CompiledEmbryos.DubuisTimesControlScalingFactors(ch_index)*CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(:,:,ch_index);
   CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Profiles(:,:,ch_index) = ...
       CompiledEmbryos.HisRFP25CTimesControlScalingFactors(ch_index)*CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(:,:,ch_index);
   CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Profiles(:,:,ch_index) = ...
       CompiledEmbryos.WindowedDeltaFCControlScalingFactors(ch_index)*CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Profiles(:,:,ch_index) = ...
       CompiledEmbryos.WindowedDubuisTimesControlScalingFactors(ch_index)*CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(:,:,ch_index);
   CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Profiles(:,:,ch_index) = ...
       CompiledEmbryos.WindowedHisRFP25CTimesControlScalingFactors(ch_index)*CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(:,:,ch_index);
end


%% Add Imin and Imax info  for SlideRescaledAvgAP.SmoothedDeltaFCFactors
x = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.x ;
AllDeltaValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan(x);% &...
AllDeltaValidProfilesControlTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.ControlSetEmbryos & ~isnan(x);% &...

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test = {};
x_test = x(AllDeltaValidProfilesTestTF);
TestProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Profiles(AllDeltaValidProfilesTestTF,:,:);
[x_test, test_order] = sort(x_test);
TestProfiles = TestProfiles(test_order,:,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.x = x_test;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Profiles = TestProfiles;
MTest = movmean(TestProfiles,10, 1, 'omitnan');
MTest = MTest(6:end-5, :,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Imin = NaN(1, NChannels); 
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Imax = NaN(1, NChannels); 
for ch_index = 2:NChannels
    MTestSub = MTest(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Imin(ch_index) = min(MTestSub(:));
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Test.Imax(ch_index) = max(MTestSub(:));
end



CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control = {};
x_control = x(AllDeltaValidProfilesControlTF);
ControlProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Profiles(AllDeltaValidProfilesControlTF,:,:);
[x_control, control_order] = sort(x_control);

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control.Imin = NaN(1, NChannels);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control.Imax = NaN(1, NChannels);

ControlProfiles = ControlProfiles(control_order,:,:);
if ~isempty(ControlProfiles)
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control.x = x_control;
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control.Profiles = ControlProfiles;
    MControl = movmean(ControlProfiles,10, 1, 'omitnan');
    
    if size(MControl,1) >= 10
        MControl = MControl(6:end-4, :,:);
        
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    elseif mod(size(MControl,1),2) == 0
        MControl = MControl((size(MControl,1)/2+1):(size(MControl,1)/2+1), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    else
        MControl = MControl(ceil(size(MControl,1)/2):ceil(size(MControl,1)/2), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDeltaFCFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    end
end

%% Add Imin and Imax info  for SlideRescaledAvgAP.SmoothedDubuisTimesFactors
x = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.x ;
AllDeltaValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan(x);% &...
AllDeltaValidProfilesControlTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.ControlSetEmbryos & ~isnan(x);% &...

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test = {};
x_test = x(AllDeltaValidProfilesTestTF);
TestProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Profiles(AllDeltaValidProfilesTestTF,:,:);
[x_test, test_order] = sort(x_test);
TestProfiles = TestProfiles(test_order,:,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.x = x_test;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.Profiles = TestProfiles;
MTest = movmean(TestProfiles,10, 1, 'omitnan');
MTest = MTest(6:end-5, :,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.Imin = NaN(1, NChannels); 
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.Imax = NaN(1, NChannels); 
for ch_index = 2:NChannels
    MTestSub = MTest(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.Imin(ch_index) = min(MTestSub(:));
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Test.Imax(ch_index) = max(MTestSub(:));
end



CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Control = {};
x_control = x(AllDeltaValidProfilesControlTF);
ControlProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Profiles(AllDeltaValidProfilesControlTF,:,:);
[x_control, control_order] = sort(x_control);

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Control.Imin = NaN(1, NChannels);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Control.Imax = NaN(1, NChannels);

ControlProfiles = ControlProfiles(control_order,:,:);
if ~isempty(ControlProfiles)
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Control.x = x_control;
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Control.Profiles = ControlProfiles;
    MControl = movmean(ControlProfiles,10, 1, 'omitnan');
    
    if size(MControl,1) >= 10
        MControl = MControl(6:end-4, :,:);
        
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    elseif mod(size(MControl,1),2) == 0
        MControl = MControl((size(MControl,1)/2+1):(size(MControl,1)/2+1), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    else
        MControl = MControl(ceil(size(MControl,1)/2):ceil(size(MControl,1)/2), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedDubuisTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    end
end
%% Add Imin and Imax info  for SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors
x = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.x ;
AllDeltaValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan(x);% &...
AllDeltaValidProfilesControlTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.ControlSetEmbryos & ~isnan(x);% &...

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test = {};
x_test = x(AllDeltaValidProfilesTestTF);
TestProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Profiles(AllDeltaValidProfilesTestTF,:,:);
[x_test, test_order] = sort(x_test);
TestProfiles = TestProfiles(test_order,:,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.x = x_test;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.Profiles = TestProfiles;
MTest = movmean(TestProfiles,10, 1, 'omitnan');
MTest = MTest(6:end-5, :,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.Imin = NaN(1, NChannels); 
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.Imax = NaN(1, NChannels); 
for ch_index = 2:NChannels
    MTestSub = MTest(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.Imin(ch_index) = min(MTestSub(:));
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Test.Imax(ch_index) = max(MTestSub(:));
end



CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Control = {};
x_control = x(AllDeltaValidProfilesControlTF);
ControlProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Profiles(AllDeltaValidProfilesControlTF,:,:);
[x_control, control_order] = sort(x_control);

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Control.Imin = NaN(1, NChannels);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Control.Imax = NaN(1, NChannels);

ControlProfiles = ControlProfiles(control_order,:,:);
if ~isempty(ControlProfiles)
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Control.x = x_control;
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Control.Profiles = ControlProfiles;
    MControl = movmean(ControlProfiles,10, 1, 'omitnan');
    
    if size(MControl,1) >= 10
        MControl = MControl(6:end-4, :,:);
        
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    elseif mod(size(MControl,1),2) == 0
        MControl = MControl((size(MControl,1)/2+1):(size(MControl,1)/2+1), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    else
        MControl = MControl(ceil(size(MControl,1)/2):ceil(size(MControl,1)/2), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.SmoothedHisRFP25CTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    end
end

%% Add Imin and Imax info  for SlideRescaledAvgAP.SmoothedDeltaFCFactors
x = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.x ;
AllDeltaValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan(x);% &...
AllDeltaValidProfilesControlTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.ControlSetEmbryos & ~isnan(x);% &...

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test = {};
x_test = x(AllDeltaValidProfilesTestTF);
TestProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Profiles(AllDeltaValidProfilesTestTF,:,:);
[x_test, test_order] = sort(x_test);
TestProfiles = TestProfiles(test_order,:,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.x = x_test;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.Profiles = TestProfiles;
MTest = movmean(TestProfiles,10, 1, 'omitnan');
MTest = MTest(6:end-5, :,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.Imin = NaN(1, NChannels); 
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.Imax = NaN(1, NChannels); 
for ch_index = 2:NChannels
    MTestSub = MTest(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.Imin(ch_index) = min(MTestSub(:));
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Test.Imax(ch_index) = max(MTestSub(:));
end



CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Control = {};
x_control = x(AllDeltaValidProfilesControlTF);
ControlProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Profiles(AllDeltaValidProfilesControlTF,:,:);
[x_control, control_order] = sort(x_control);

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Control.Imin = NaN(1, NChannels);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Control.Imax = NaN(1, NChannels);

ControlProfiles = ControlProfiles(control_order,:,:);
if ~isempty(ControlProfiles)
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Control.x = x_control;
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Control.Profiles = ControlProfiles;
    MControl = movmean(ControlProfiles,10, 1, 'omitnan');
    
    if size(MControl,1) >= 10
        MControl = MControl(6:end-4, :,:);
        
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    elseif mod(size(MControl,1),2) == 0
        MControl = MControl((size(MControl,1)/2+1):(size(MControl,1)/2+1), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    else
        MControl = MControl(ceil(size(MControl,1)/2):ceil(size(MControl,1)/2), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDeltaFCFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    end
end

%% Add Imin and Imax info  for SlideRescaledAvgAP.WindowedDubuisTimesFactors
x = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.x ;
AllDeltaValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan(x);% &...
AllDeltaValidProfilesControlTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.ControlSetEmbryos & ~isnan(x);% &...

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test = {};
x_test = x(AllDeltaValidProfilesTestTF);
TestProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Profiles(AllDeltaValidProfilesTestTF,:,:);
[x_test, test_order] = sort(x_test);
TestProfiles = TestProfiles(test_order,:,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.x = x_test;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.Profiles = TestProfiles;
MTest = movmean(TestProfiles,10, 1, 'omitnan');
MTest = MTest(6:end-5, :,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.Imin = NaN(1, NChannels); 
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.Imax = NaN(1, NChannels); 
for ch_index = 2:NChannels
    MTestSub = MTest(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.Imin(ch_index) = min(MTestSub(:));
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Test.Imax(ch_index) = max(MTestSub(:));
end



CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Control = {};
x_control = x(AllDeltaValidProfilesControlTF);
ControlProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Profiles(AllDeltaValidProfilesControlTF,:,:);
[x_control, control_order] = sort(x_control);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Control.Imin = NaN(1, NChannels);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Control.Imax = NaN(1, NChannels);

ControlProfiles = ControlProfiles(control_order,:,:);
if ~isempty(ControlProfiles)
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Control.x = x_control;
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Control.Profiles = ControlProfiles;
    MControl = movmean(ControlProfiles,10, 1, 'omitnan');
    
    if size(MControl,1) >= 10
        MControl = MControl(6:end-4, :,:);
        
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    elseif mod(size(MControl,1),2) == 0
        MControl = MControl((size(MControl,1)/2+1):(size(MControl,1)/2+1), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    else
        MControl = MControl(ceil(size(MControl,1)/2):ceil(size(MControl,1)/2), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedDubuisTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    end
end

%% Add Imin and Imax info  for SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors
x = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.x ;
AllDeltaValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan(x);% &...
AllDeltaValidProfilesControlTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.ControlSetEmbryos & ~isnan(x);% &...

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test = {};
x_test = x(AllDeltaValidProfilesTestTF);
TestProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Profiles(AllDeltaValidProfilesTestTF,:,:);
[x_test, test_order] = sort(x_test);
TestProfiles = TestProfiles(test_order,:,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.x = x_test;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.Profiles = TestProfiles;
MTest = movmean(TestProfiles,10, 1, 'omitnan');
MTest = MTest(6:end-5, :,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.Imin = NaN(1, NChannels); 
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.Imax = NaN(1, NChannels); 
for ch_index = 2:NChannels
    MTestSub = MTest(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.Imin(ch_index) = min(MTestSub(:));
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Test.Imax(ch_index) = max(MTestSub(:));
end



CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Control = {};
x_control = x(AllDeltaValidProfilesControlTF);
ControlProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Profiles(AllDeltaValidProfilesControlTF,:,:);
[x_control, control_order] = sort(x_control);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Control.Imin = NaN(1, NChannels);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Control.Imax = NaN(1, NChannels);

ControlProfiles = ControlProfiles(control_order,:,:);
if ~isempty(ControlProfiles)
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Control.x = x_control;
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Control.Profiles = ControlProfiles;
    MControl = movmean(ControlProfiles,10, 1, 'omitnan');
    
    if size(MControl,1) >= 10
        MControl = MControl(6:end-4, :,:);
        
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    elseif mod(size(MControl,1),2) == 0
        MControl = MControl((size(MControl,1)/2+1):(size(MControl,1)/2+1), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    else
        MControl = MControl(ceil(size(MControl,1)/2):ceil(size(MControl,1)/2), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAvgAP.WindowedHisRFP25CTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    end
end



%%
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP = {};
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors = {};
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.x = CompiledEmbryos.FixCorrectedDeltaFC_um.mean;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Profiles = ...
    NaN(size(CompiledEmbryos.SlideRescaledDorsalAPProfiles));
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors = {};
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.x = CompiledEmbryos.DubuisEmbryoTimes;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Profiles = ...
    NaN(size(CompiledEmbryos.SlideRescaledDorsalAPProfiles));

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors = {};
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.x = CompiledEmbryos.HisRFP25CEmbryoTimes;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Profiles = ...
    NaN(size(CompiledEmbryos.SlideRescaledDorsalAPProfiles));

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors = {};
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.x = CompiledEmbryos.FixCorrectedDeltaFC_um.mean;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Profiles = ...
    NaN(size(CompiledEmbryos.SlideRescaledDorsalAPProfiles));

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors = {};
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.x = CompiledEmbryos.DubuisEmbryoTimes;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Profiles = ...
    NaN(size(CompiledEmbryos.SlideRescaledDorsalAPProfiles));

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors = {};
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.x = CompiledEmbryos.HisRFP25CEmbryoTimes;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Profiles = ...
    NaN(size(CompiledEmbryos.SlideRescaledDorsalAPProfiles));

for ch_index = 2:NChannels
   CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Profiles(:,:,ch_index) = ...
       CompiledEmbryos.FixCorrectedControlScalingFactors(ch_index)*CompiledEmbryos.SlideRescaledDorsalAPProfiles(:,:,ch_index);
   CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Profiles(:,:,ch_index) = ...
       CompiledEmbryos.DubuisTimesControlScalingFactors(ch_index)*CompiledEmbryos.SlideRescaledDorsalAPProfiles(:,:,ch_index);
   CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Profiles(:,:,ch_index) = ...
       CompiledEmbryos.HisRFP25CTimesControlScalingFactors(ch_index)*CompiledEmbryos.SlideRescaledDorsalAPProfiles(:,:,ch_index);
   CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Profiles(:,:,ch_index) = ...
       CompiledEmbryos.WindowedDeltaFCControlScalingFactors(ch_index)*CompiledEmbryos.SlideRescaledDorsalAPProfiles(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Profiles(:,:,ch_index) = ...
       CompiledEmbryos.WindowedDubuisTimesControlScalingFactors(ch_index)*CompiledEmbryos.SlideRescaledDorsalAPProfiles(:,:,ch_index);
   CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Profiles(:,:,ch_index) = ...
       CompiledEmbryos.WindowedHisRFP25CTimesControlScalingFactors(ch_index)*CompiledEmbryos.SlideRescaledDorsalAPProfiles(:,:,ch_index);
end

%% Add Imin and Imax info  for SlideRescaledAP.SmoothedDeltaFCFactors
x = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.x ;
AllDeltaValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan(x);% &...
AllDeltaValidProfilesControlTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.ControlSetEmbryos & ~isnan(x);% &...

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test = {};
x_test = x(AllDeltaValidProfilesTestTF);
TestProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Profiles(AllDeltaValidProfilesTestTF,:,:);
[x_test, test_order] = sort(x_test);
TestProfiles = TestProfiles(test_order,:,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.x = x_test;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.Profiles = TestProfiles;
MTest = movmean(TestProfiles,10, 1, 'omitnan');
MTest = MTest(6:end-5, :,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.Imin = NaN(1, NChannels); 
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.Imax = NaN(1, NChannels); 
for ch_index = 2:NChannels
    MTestSub = MTest(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.Imin(ch_index) = min(MTestSub(:));
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Test.Imax(ch_index) = max(MTestSub(:));
end

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Control = {};
x_control = x(AllDeltaValidProfilesControlTF);
ControlProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Profiles(AllDeltaValidProfilesControlTF,:,:);
[x_control, control_order] = sort(x_control);

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Control.Imin = NaN(1, NChannels);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Control.Imax = NaN(1, NChannels);

ControlProfiles = ControlProfiles(control_order,:,:);
if ~isempty(ControlProfiles)
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Control.x = x_control;
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Control.Profiles = ControlProfiles;
    MControl = movmean(ControlProfiles,10, 1, 'omitnan');
    
    if size(MControl,1) >= 10
        MControl = MControl(6:end-4, :,:);
        
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    elseif mod(size(MControl,1),2) == 0
        MControl = MControl((size(MControl,1)/2+1):(size(MControl,1)/2+1), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    else
        MControl = MControl(ceil(size(MControl,1)/2):ceil(size(MControl,1)/2), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDeltaFCFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    end
end



%% Add Imin and Imax info  for SlideRescaledAP.SmoothedDubuisTimesFactors
x = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.x ;
AllDeltaValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan(x);% &...
AllDeltaValidProfilesControlTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.ControlSetEmbryos & ~isnan(x);% &...

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test = {};
x_test = x(AllDeltaValidProfilesTestTF);
TestProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Profiles(AllDeltaValidProfilesTestTF,:,:);
[x_test, test_order] = sort(x_test);
TestProfiles = TestProfiles(test_order,:,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.x = x_test;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.Profiles = TestProfiles;
MTest = movmean(TestProfiles,10, 1, 'omitnan');
MTest = MTest(6:end-5, :,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.Imin = NaN(1, NChannels); 
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.Imax = NaN(1, NChannels); 
for ch_index = 2:NChannels
    MTestSub = MTest(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.Imin(ch_index) = min(MTestSub(:));
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Test.Imax(ch_index) = max(MTestSub(:));
end



CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Control = {};
x_control = x(AllDeltaValidProfilesControlTF);
ControlProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Profiles(AllDeltaValidProfilesControlTF,:,:);
[x_control, control_order] = sort(x_control);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Control.Imin = NaN(1, NChannels);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Control.Imax = NaN(1, NChannels);

ControlProfiles = ControlProfiles(control_order,:,:);
if ~isempty(ControlProfiles)
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Control.x = x_control;
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Control.Profiles = ControlProfiles;
    MControl = movmean(ControlProfiles,10, 1, 'omitnan');
    
    if size(MControl,1) >= 10
        MControl = MControl(6:end-4, :,:);
        
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    elseif mod(size(MControl,1),2) == 0
        MControl = MControl((size(MControl,1)/2+1):(size(MControl,1)/2+1), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    else
        MControl = MControl(ceil(size(MControl,1)/2):ceil(size(MControl,1)/2), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedDubuisTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    end
end

%% Add Imin and Imax info  for SlideRescaledAP.SmoothedHisRFP25CTimesFactors
x = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.x ;
AllDeltaValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan(x);% &...
AllDeltaValidProfilesControlTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.ControlSetEmbryos & ~isnan(x);% &...

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test = {};
x_test = x(AllDeltaValidProfilesTestTF);
TestProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Profiles(AllDeltaValidProfilesTestTF,:,:);
[x_test, test_order] = sort(x_test);
TestProfiles = TestProfiles(test_order,:,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.x = x_test;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.Profiles = TestProfiles;
MTest = movmean(TestProfiles,10, 1, 'omitnan');
MTest = MTest(6:end-5, :,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.Imin = NaN(1, NChannels); 
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.Imax = NaN(1, NChannels); 
for ch_index = 2:NChannels
    MTestSub = MTest(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.Imin(ch_index) = min(MTestSub(:));
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Test.Imax(ch_index) = max(MTestSub(:));
end



CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Control = {};
x_control = x(AllDeltaValidProfilesControlTF);
ControlProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Profiles(AllDeltaValidProfilesControlTF,:,:);
[x_control, control_order] = sort(x_control);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Control.Imin = NaN(1, NChannels);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Control.Imax = NaN(1, NChannels);

ControlProfiles = ControlProfiles(control_order,:,:);
if ~isempty(ControlProfiles)
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Control.x = x_control;
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Control.Profiles = ControlProfiles;
    MControl = movmean(ControlProfiles,10, 1, 'omitnan');
    
    if size(MControl,1) >= 10
        MControl = MControl(6:end-4, :,:);
        
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    elseif mod(size(MControl,1),2) == 0
        MControl = MControl((size(MControl,1)/2+1):(size(MControl,1)/2+1), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    else
        MControl = MControl(ceil(size(MControl,1)/2):ceil(size(MControl,1)/2), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.SmoothedHisRFP25CTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    end
end

%% Add Imin and Imax info  for SlideRescaledAP.WindowedDeltaFCFactors
x = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.x ;
AllDeltaValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan(x);% &...
AllDeltaValidProfilesControlTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.ControlSetEmbryos & ~isnan(x);% &...

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test = {};
x_test = x(AllDeltaValidProfilesTestTF);
TestProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Profiles(AllDeltaValidProfilesTestTF,:,:);
[x_test, test_order] = sort(x_test);
TestProfiles = TestProfiles(test_order,:,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.x = x_test;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.Profiles = TestProfiles;
MTest = movmean(TestProfiles,10, 1, 'omitnan');
MTest = MTest(6:end-5, :,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.Imin = NaN(1, NChannels); 
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.Imax = NaN(1, NChannels); 
for ch_index = 2:NChannels
    MTestSub = MTest(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.Imin(ch_index) = min(MTestSub(:));
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Test.Imax(ch_index) = max(MTestSub(:));
end



CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Control = {};
x_control = x(AllDeltaValidProfilesControlTF);
ControlProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Profiles(AllDeltaValidProfilesControlTF,:,:);
[x_control, control_order] = sort(x_control);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Control.Imin = NaN(1, NChannels);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Control.Imax = NaN(1, NChannels);

ControlProfiles = ControlProfiles(control_order,:,:);
if ~isempty(ControlProfiles)
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Control.x = x_control;
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Control.Profiles = ControlProfiles;
    MControl = movmean(ControlProfiles,10, 1, 'omitnan');
    
    if size(MControl,1) >= 10
        MControl = MControl(6:end-4, :,:);
        
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    elseif mod(size(MControl,1),2) == 0
        MControl = MControl((size(MControl,1)/2+1):(size(MControl,1)/2+1), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    else
        MControl = MControl(ceil(size(MControl,1)/2):ceil(size(MControl,1)/2), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDeltaFCFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    end
end

%% Add Imin and Imax info  for SlideRescaledAP.WindowedDubuisTimesFactors
x = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.x ;
AllDeltaValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan(x);% &...
AllDeltaValidProfilesControlTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.ControlSetEmbryos & ~isnan(x);% &...

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test = {};
x_test = x(AllDeltaValidProfilesTestTF);
TestProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Profiles(AllDeltaValidProfilesTestTF,:,:);
[x_test, test_order] = sort(x_test);
TestProfiles = TestProfiles(test_order,:,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.x = x_test;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.Profiles = TestProfiles;
MTest = movmean(TestProfiles,10, 1, 'omitnan');
MTest = MTest(6:end-5, :,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.Imin = NaN(1, NChannels); 
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.Imax = NaN(1, NChannels); 
for ch_index = 2:NChannels
    MTestSub = MTest(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.Imin(ch_index) = min(MTestSub(:));
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Test.Imax(ch_index) = max(MTestSub(:));
end



CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Control = {};
x_control = x(AllDeltaValidProfilesControlTF);
ControlProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Profiles(AllDeltaValidProfilesControlTF,:,:);
[x_control, control_order] = sort(x_control);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Control.Imin = NaN(1, NChannels);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Control.Imax = NaN(1, NChannels);

ControlProfiles = ControlProfiles(control_order,:,:);
if ~isempty(ControlProfiles)
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Control.x = x_control;
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Control.Profiles = ControlProfiles;
    MControl = movmean(ControlProfiles,10, 1, 'omitnan');
    
    if size(MControl,1) >= 10
        MControl = MControl(6:end-4, :,:);
        
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    elseif mod(size(MControl,1),2) == 0
        MControl = MControl((size(MControl,1)/2+1):(size(MControl,1)/2+1), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    else
        MControl = MControl(ceil(size(MControl,1)/2):ceil(size(MControl,1)/2), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedDubuisTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    end
end

%% Add Imin and Imax info  for SlideRescaledAP.WindowedHisRFP25CTimesFactors
x = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.x ;
AllDeltaValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan(x);% &...
AllDeltaValidProfilesControlTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.ControlSetEmbryos & ~isnan(x);% &...

CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test = {};
x_test = x(AllDeltaValidProfilesTestTF);
TestProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Profiles(AllDeltaValidProfilesTestTF,:,:);
[x_test, test_order] = sort(x_test);
TestProfiles = TestProfiles(test_order,:,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.x = x_test;
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.Profiles = TestProfiles;
MTest = movmean(TestProfiles,10, 1, 'omitnan');
MTest = MTest(6:end-5, :,:);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.Imin = NaN(1, NChannels); 
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.Imax = NaN(1, NChannels); 
for ch_index = 2:NChannels
    MTestSub = MTest(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.Imin(ch_index) = min(MTestSub(:));
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Test.Imax(ch_index) = max(MTestSub(:));
end



CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Control = {};
x_control = x(AllDeltaValidProfilesControlTF);
ControlProfiles = CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Profiles(AllDeltaValidProfilesControlTF,:,:);
[x_control, control_order] = sort(x_control);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Control.Imin = NaN(1, NChannels);
CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Control.Imax = NaN(1, NChannels);

ControlProfiles = ControlProfiles(control_order,:,:);
if ~isempty(ControlProfiles)
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Control.x = x_control;
    CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Control.Profiles = ControlProfiles;
    MControl = movmean(ControlProfiles,10, 1, 'omitnan');
    
    if size(MControl,1) >= 10
        MControl = MControl(6:end-4, :,:);
        
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    elseif mod(size(MControl,1),2) == 0
        MControl = MControl((size(MControl,1)/2+1):(size(MControl,1)/2+1), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    else
        MControl = MControl(ceil(size(MControl,1)/2):ceil(size(MControl,1)/2), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.SlideRescaledAP.WindowedHisRFP25CTimesFactors.Control.Imax(ch_index) = max(MControlSub(:));
        end
    end
end

%%
CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP = {};
CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.x = CompiledEmbryos.FixCorrectedSmoothedDeltaFCs;
CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.Test = {};
CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.Test.Profiles = ...
    NaN(size(CompiledEmbryos.FixCorrectedSmoothedAvgAPProfiles.Test));
CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.Control = {};
CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.Control.Profiles  = ...
    NaN(size(CompiledEmbryos.FixCorrectedSmoothedAvgAPProfiles.Control));

CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP = {};
CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.x = CompiledEmbryos.DubuisSmoothedTimes;
CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.Test = {};
CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.Test.Profiles = ...
    NaN(size(CompiledEmbryos.DubuisTimeSmoothedAvgAPProfiles.Test));
CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.Control = {};
CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.Control.Profiles = ...
    NaN(size(CompiledEmbryos.DubuisTimeSmoothedAvgAPProfiles.Control));

CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP = {};
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.x = CompiledEmbryos.HisRFP25CSmoothedTimes;
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.Test = {};
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.Test.Profiles = ...
    NaN(size(CompiledEmbryos.HisRFP25CTimeSmoothedAvgAPProfiles.Test));
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.Control = {};
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.Control.Profiles = ...
    NaN(size(CompiledEmbryos.HisRFP25CTimeSmoothedAvgAPProfiles.Control));

for ch_index = 2:NChannels
    CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.Test.Profiles(:,:,ch_index) = ...
      CompiledEmbryos.FixCorrectedControlScalingFactors(ch_index)*CompiledEmbryos.FixCorrectedSmoothedAvgAPProfiles.Test(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.Control.Profiles(:,:,ch_index) = ...
      CompiledEmbryos.FixCorrectedControlScalingFactors(ch_index)*CompiledEmbryos.FixCorrectedSmoothedAvgAPProfiles.Control(:,:,ch_index);
    
    CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.Test.Profiles(:,:,ch_index) = ...
      CompiledEmbryos.DubuisTimesControlScalingFactors(ch_index)*CompiledEmbryos.DubuisTimeSmoothedAvgAPProfiles.Test(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.Control.Profiles(:,:,ch_index) = ...
      CompiledEmbryos.DubuisTimesControlScalingFactors(ch_index)*CompiledEmbryos.DubuisTimeSmoothedAvgAPProfiles.Control(:,:,ch_index);

     CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.Test.Profiles(:,:,ch_index) = ...
      CompiledEmbryos.HisRFP25CTimesControlScalingFactors(ch_index)*CompiledEmbryos.HisRFP25CTimeSmoothedAvgAPProfiles.Test(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.Control.Profiles(:,:,ch_index) = ...
      CompiledEmbryos.HisRFP25CTimesControlScalingFactors(ch_index)*CompiledEmbryos.HisRFP25CTimeSmoothedAvgAPProfiles.Control(:,:,ch_index);

end
%% Add Imin and Imax info  for DeltaFCSmoothedAvgAP
TestProfiles = CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.Test.Profiles;
MTest = movmean(TestProfiles,10, 1, 'omitnan');
MTest = MTest(6:end-5, :,:);
CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.Test.Imin = NaN(1, NChannels); 
CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.Test.Imax = NaN(1, NChannels); 
for ch_index = 2:NChannels
    MTestSub = MTest(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.Test.Imin(ch_index) = min(MTestSub(:));
    CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.Test.Imax(ch_index) = max(MTestSub(:));
end

ControlProfiles = CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.Control.Profiles;
CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.Control.Imin = NaN(1, NChannels);
CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.Control.Imax = NaN(1, NChannels);

if ~isempty(ControlProfiles)
    MControl = movmean(ControlProfiles,10, 1, 'omitnan');
    
    if size(MControl,1) >= 10
        MControl = MControl(6:end-4, :,:);
        
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.Control.Imax(ch_index) = max(MControlSub(:));
        end
    elseif mod(size(MControl,1),2) == 0
        MControl = MControl((size(MControl,1)/2+1):(size(MControl,1)/2+1), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.Control.Imax(ch_index) = max(MControlSub(:));
        end
    else
        MControl = MControl(ceil(size(MControl,1)/2):ceil(size(MControl,1)/2), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.DeltaFCSmoothedAvgAP.Control.Imax(ch_index) = max(MControlSub(:));
        end
    end
end

%% Add Imin and Imax info  for DubuisTimesSmoothedAvgAP
TestProfiles = CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.Test.Profiles;
MTest = movmean(TestProfiles,10, 1, 'omitnan');
MTest = MTest(6:end-5, :,:);
CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.Test.Imin = NaN(1, NChannels); 
CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.Test.Imax = NaN(1, NChannels); 
for ch_index = 2:NChannels
    MTestSub = MTest(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.Test.Imin(ch_index) = min(MTestSub(:));
    CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.Test.Imax(ch_index) = max(MTestSub(:));
end

ControlProfiles = CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.Control.Profiles;
CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.Control.Imin = NaN(1, NChannels);
CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.Control.Imax = NaN(1, NChannels);

if ~isempty(ControlProfiles)
    MControl = movmean(ControlProfiles,10, 1, 'omitnan');
    
    if size(MControl,1) >= 10
        MControl = MControl(6:end-4, :,:);
        
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.Control.Imax(ch_index) = max(MControlSub(:));
        end
    elseif mod(size(MControl,1),2) == 0
        MControl = MControl((size(MControl,1)/2+1):(size(MControl,1)/2+1), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.Control.Imax(ch_index) = max(MControlSub(:));
        end
    else
        MControl = MControl(ceil(size(MControl,1)/2):ceil(size(MControl,1)/2), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.DubuisTimesSmoothedAvgAP.Control.Imax(ch_index) = max(MControlSub(:));
        end
    end
end
%% Add Imin and Imax info  for HisRFP25CTimesSmoothedAvgAP
TestProfiles = CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.Test.Profiles;
MTest = movmean(TestProfiles,10, 1, 'omitnan');
MTest = MTest(6:end-5, :,:);
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.Test.Imin = NaN(1, NChannels); 
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.Test.Imax = NaN(1, NChannels); 
for ch_index = 2:NChannels
    MTestSub = MTest(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.Test.Imin(ch_index) = min(MTestSub(:));
    CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.Test.Imax(ch_index) = max(MTestSub(:));
end

ControlProfiles = CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.Control.Profiles;
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.Control.Imin = NaN(1, NChannels);
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.Control.Imax = NaN(1, NChannels);

if ~isempty(ControlProfiles)
    MControl = movmean(ControlProfiles,10, 1, 'omitnan');
    
    if size(MControl,1) >= 10
        MControl = MControl(6:end-4, :,:);
        
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.Control.Imax(ch_index) = max(MControlSub(:));
        end
    elseif mod(size(MControl,1),2) == 0
        MControl = MControl((size(MControl,1)/2+1):(size(MControl,1)/2+1), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.Control.Imax(ch_index) = max(MControlSub(:));
        end
    else
        MControl = MControl(ceil(size(MControl,1)/2):ceil(size(MControl,1)/2), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesSmoothedAvgAP.Control.Imax(ch_index) = max(MControlSub(:));
        end
    end
end

%%
CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP = {};
CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.x = CompiledEmbryos.WindowedProfiles.DeltaFC.x;
CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Test.mean = ...
    NaN(size(CompiledEmbryos.WindowedProfiles.DeltaFC.Test.mean));
CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Test.se = ...
    NaN(size(CompiledEmbryos.WindowedProfiles.DeltaFC.Test.se));
CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Test.count = ...
    NaN(size(CompiledEmbryos.WindowedProfiles.DeltaFC.Test.count));
CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Control.mean = ...
    NaN(size(CompiledEmbryos.WindowedProfiles.DeltaFC.Control.mean));
CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Control.se = ...
    NaN(size(CompiledEmbryos.WindowedProfiles.DeltaFC.Control.se));
CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Control.count = ...
    NaN(size(CompiledEmbryos.WindowedProfiles.DeltaFC.Control.count));
for ch_index = 2:NChannels
    CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Test.mean(:,:,ch_index) = ...
      CompiledEmbryos.WindowedDeltaFCControlScalingFactors(ch_index)*CompiledEmbryos.WindowedProfiles.DeltaFC.Test.mean(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Test.se(:,:,ch_index) = ...
      CompiledEmbryos.WindowedDeltaFCControlScalingFactors(ch_index)*CompiledEmbryos.WindowedProfiles.DeltaFC.Test.se(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Test.count(:,:,ch_index) = ...
      CompiledEmbryos.WindowedProfiles.DeltaFC.Test.count(:,:,ch_index);
  
    CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Control.mean(:,:,ch_index) = ...
      CompiledEmbryos.WindowedDeltaFCControlScalingFactors(ch_index)*CompiledEmbryos.WindowedProfiles.DeltaFC.Control.mean(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Control.se(:,:,ch_index) = ...
      CompiledEmbryos.WindowedDeltaFCControlScalingFactors(ch_index)*CompiledEmbryos.WindowedProfiles.DeltaFC.Control.se(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Control.count(:,:,ch_index) = ...
      CompiledEmbryos.WindowedProfiles.DeltaFC.Control.count(:,:,ch_index);
end

%% Add Imin and Imax info  for DeltaFCWindowedAvgAP
TestProfiles = CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Test.mean;
MTest = movmean(TestProfiles,10, 1, 'omitnan');
MTest = MTest(6:end-5, :,:);
CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Test.Imin = NaN(1, NChannels); 
CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Test.Imax = NaN(1, NChannels); 
for ch_index = 2:NChannels
    MTestSub = MTest(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Test.Imin(ch_index) = min(MTestSub(:));
    CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Test.Imax(ch_index) = max(MTestSub(:));
end


ControlProfiles = CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Control.mean;
CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Control.Imin = NaN(1, NChannels);
CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Control.Imax = NaN(1, NChannels);

if ~isempty(ControlProfiles)
    MControl = movmean(ControlProfiles,10, 1, 'omitnan');
    
    if size(MControl,1) >= 10
        MControl = MControl(6:end-4, :,:);
        
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Control.Imax(ch_index) = max(MControlSub(:));
        end
    elseif mod(size(MControl,1),2) == 0
        MControl = MControl((size(MControl,1)/2+1):(size(MControl,1)/2+1), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Control.Imax(ch_index) = max(MControlSub(:));
        end
    else
        MControl = MControl(ceil(size(MControl,1)/2):ceil(size(MControl,1)/2), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.DeltaFCWindowedAvgAP.Control.Imax(ch_index) = max(MControlSub(:));
        end
    end
end


%%
CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP = {};
CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.x = CompiledEmbryos.WindowedProfiles.DubuisTime.x;
CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Test.mean = ...
    NaN(size(CompiledEmbryos.WindowedProfiles.DubuisTime.Test.mean));
CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Test.se = ...
    NaN(size(CompiledEmbryos.WindowedProfiles.DubuisTime.Test.se));
CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Test.count = ...
    NaN(size(CompiledEmbryos.WindowedProfiles.DubuisTime.Test.count));
CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Control.mean = ...
    NaN(size(CompiledEmbryos.WindowedProfiles.DubuisTime.Control.mean));
CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Control.se = ...
    NaN(size(CompiledEmbryos.WindowedProfiles.DubuisTime.Control.se));
CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Control.count = ...
    NaN(size(CompiledEmbryos.WindowedProfiles.DubuisTime.Control.count));
for ch_index = 2:NChannels
    CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Test.mean(:,:,ch_index) = ...
      CompiledEmbryos.WindowedDubuisTimesControlScalingFactors(ch_index)*CompiledEmbryos.WindowedProfiles.DubuisTime.Test.mean(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Test.se(:,:,ch_index) = ...
      CompiledEmbryos.WindowedDubuisTimesControlScalingFactors(ch_index)*CompiledEmbryos.WindowedProfiles.DubuisTime.Test.se(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Test.count(:,:,ch_index) = ...
      CompiledEmbryos.WindowedProfiles.DubuisTime.Test.count(:,:,ch_index);
  
    CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Control.mean(:,:,ch_index) = ...
      CompiledEmbryos.WindowedDubuisTimesControlScalingFactors(ch_index)*CompiledEmbryos.WindowedProfiles.DubuisTime.Control.mean(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Control.se(:,:,ch_index) = ...
      CompiledEmbryos.WindowedDubuisTimesControlScalingFactors(ch_index)*CompiledEmbryos.WindowedProfiles.DubuisTime.Control.se(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Control.count(:,:,ch_index) = ...
      CompiledEmbryos.WindowedProfiles.DubuisTime.Control.count(:,:,ch_index);
end

%% Add Imin and Imax info  for DubuisTimesWindowedAvgAP
TestProfiles = CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Test.mean;
MTest = movmean(TestProfiles,10, 1, 'omitnan');
MTest = MTest(6:end-5, :,:);
CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Test.Imin = NaN(1, NChannels); 
CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Test.Imax = NaN(1, NChannels); 
for ch_index = 2:NChannels
    MTestSub = MTest(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Test.Imin(ch_index) = min(MTestSub(:));
    CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Test.Imax(ch_index) = max(MTestSub(:));
end

ControlProfiles = CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Control.mean;
CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Control.Imin = NaN(1, NChannels);
CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Control.Imax = NaN(1, NChannels);

if ~isempty(ControlProfiles)
    MControl = movmean(ControlProfiles,10, 1, 'omitnan');
    
    if size(MControl,1) >= 10
        MControl = MControl(6:end-4, :,:);
        
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Control.Imax(ch_index) = max(MControlSub(:));
        end
    elseif mod(size(MControl,1),2) == 0
        MControl = MControl((size(MControl,1)/2+1):(size(MControl,1)/2+1), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Control.Imax(ch_index) = max(MControlSub(:));
        end
    else
        MControl = MControl(ceil(size(MControl,1)/2):ceil(size(MControl,1)/2), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.DubuisTimesWindowedAvgAP.Control.Imax(ch_index) = max(MControlSub(:));
        end
    end
end
%%
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP = {};
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.x = CompiledEmbryos.WindowedProfiles.HisRFP25CTime.x;
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Test.mean = ...
    NaN(size(CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test.mean));
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Test.se = ...
    NaN(size(CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test.se));
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Test.count = ...
    NaN(size(CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test.count));
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Control.mean = ...
    NaN(size(CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control.mean));
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Test.se = ...
    NaN(size(CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control.se));
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Control.count = ...
    NaN(size(CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control.count));
for ch_index = 2:NChannels
    CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Test.mean(:,:,ch_index) = ...
      CompiledEmbryos.WindowedHisRFP25CTimesControlScalingFactors(ch_index)*CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test.mean(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Test.se(:,:,ch_index) = ...
      CompiledEmbryos.WindowedHisRFP25CTimesControlScalingFactors(ch_index)*CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test.se(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Test.count(:,:,ch_index) = ...
      CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test.count(:,:,ch_index);
  
    CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Control.mean(:,:,ch_index) = ...
      CompiledEmbryos.WindowedHisRFP25CTimesControlScalingFactors(ch_index)*CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control.mean(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Control.se(:,:,ch_index) = ...
      CompiledEmbryos.WindowedHisRFP25CTimesControlScalingFactors(ch_index)*CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control.se(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Control.count(:,:,ch_index) = ...
      CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control.count(:,:,ch_index);
end

%% Add Imin and Imax info  for HisRFP25CTimesWindowedAvgAP
TestProfiles = CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Test.mean;
MTest = movmean(TestProfiles,10, 1, 'omitnan');
MTest = MTest(6:end-5, :,:);
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Test.Imin = NaN(1, NChannels); 
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Test.Imax = NaN(1, NChannels); 
for ch_index = 2:NChannels
    MTestSub = MTest(:,:,ch_index);
    CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Test.Imin(ch_index) = min(MTestSub(:));
    CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Test.Imax(ch_index) = max(MTestSub(:));
end

ControlProfiles = CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Control.mean;
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Control.Imin = NaN(1, NChannels);
CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Control.Imax = NaN(1, NChannels);

if ~isempty(ControlProfiles)
    MControl = movmean(ControlProfiles,10, 1, 'omitnan');
    
    if size(MControl,1) >= 10
        MControl = MControl(6:end-4, :,:);
        
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Control.Imax(ch_index) = max(MControlSub(:));
        end
    elseif mod(size(MControl,1),2) == 0
        MControl = MControl((size(MControl,1)/2+1):(size(MControl,1)/2+1), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Control.Imax(ch_index) = max(MControlSub(:));
        end
    else
        MControl = MControl(ceil(size(MControl,1)/2):ceil(size(MControl,1)/2), :,:);
        for ch_index = 2:NChannels
            MControlSub = MControl(:,:,ch_index);
            CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Control.Imin(ch_index) = min(MControlSub(:));
            CompiledEmbryos.UnivScaledProfiles.HisRFP25CTimesWindowedAvgAP.Control.Imax(ch_index) = max(MControlSub(:));
        end
    end
end

CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
save(CEoutpath, 'CompiledEmbryos');