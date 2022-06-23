function AllSetInfo =  GetFixedSetPrefixInfo
%% author: G. Martini
% date created: 6/8/2022
AllSetInfo = {};
AllSetInfo.Temperatures = [];
AllSetInfo.Replicates = [];
AllSetInfo.Flipped = [];
AllSetInfo.SetLabels = {};
AllSetInfo.PlotLabels = {};
AllSetInfo.Prefixes = {};
AllSetInfo.DorsalAvgAPKnirpsBackgroundCutoff = {};
AllSetInfo.UsePeakRatioToIdentifyControls = {};
AllSetInfo.DorsalAvgAPKnirpsPeakValueCutoff = {};
AllSetInfo.UsePeakValueToIdentifyControls = {};
AllSetInfo.DorsalAvgAPKnirpsMeanValueCutoff = {};
AllSetInfo.UseMeanValueToIdentifyControls = {};
AllSetInfo.MinimumFixCorrectedDeltaFC = {};




% T=25C, Replicate 1

SetLabel = 'FixedT25CReplicate1';
PlotLabel = '25°C Replicate 1';
SetPrefixes = {'2022-04-03-yw104-180m25C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide3',...
    '2022-04-04-yw104-180m25C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide10',...
    '2022-04-07-yw104-180m25C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide9',...
    '2022-04-09-yw104-180m25C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide11'};
AllSetInfo.Temperatures(1, end+1) = 25;
AllSetInfo.Replicates(1, end+1) = 1;
AllSetInfo.Flipped(1, end+1) = 0;
AllSetInfo.SetLabels{end+1} = SetLabel;
AllSetInfo.PlotLabels{end+1} = PlotLabel;
AllSetInfo.Prefixes{end+1} = SetPrefixes;
AllSetInfo.DorsalAvgAPKnirpsBackgroundCutoff{end+1} = [1 NaN 6 6.5];
AllSetInfo.UsePeakRatioToIdentifyControls{end+1} = logical([1 0 1 1]);
AllSetInfo.DorsalAvgAPKnirpsPeakValueCutoff{end+1} = [1.5 15 10 10];
AllSetInfo.UsePeakValueToIdentifyControls{end+1} = logical([1 1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsMeanValueCutoff{end+1} = [1 NaN 8 7];
AllSetInfo.UseMeanValueToIdentifyControls{end+1} = logical([1 0 1 1]);
AllSetInfo.MinimumFixCorrectedDeltaFC{end+1} = [0 0 0 0];
% 
% T=25C, Replicate 2
SetLabel = 'FixedT25CReplicate2';
PlotLabel = '25°C Replicate 2';
SetPrefixes = {'2022-05-07-yw124-180m25C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep2Slide4',...
    '2022-05-11-yw124-180m25C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep2Slide6',...
    '2022-05-12-yw124-180m25C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep2Slide5'};
AllSetInfo.Temperatures(1, end+1) = 25;
AllSetInfo.Replicates(1, end+1) = 2;
AllSetInfo.Flipped(1, end+1) = 0;
AllSetInfo.SetLabels{end+1} = SetLabel;
AllSetInfo.PlotLabels{end+1} = PlotLabel;
AllSetInfo.Prefixes{end+1} = SetPrefixes;
AllSetInfo.DorsalAvgAPKnirpsBackgroundCutoff{end+1} = [6 5 5];
AllSetInfo.UsePeakRatioToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsPeakValueCutoff{end+1} = [15 NaN 10];
AllSetInfo.UsePeakValueToIdentifyControls{end+1} = logical([1 0 1]);
AllSetInfo.DorsalAvgAPKnirpsMeanValueCutoff{end+1} = [8 5 6];
AllSetInfo.UseMeanValueToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.MinimumFixCorrectedDeltaFC{end+1} = [0 0 0];

% % T=25C, Flipped Stain
SetLabel = 'FixedT25CFlipped';
PlotLabel = '25°C Flipped Stain';
SetPrefixes = {'2022-05-02-yw164-180m25C_Kni647_yw174-190m25C_Bcd488_Hb546_Hoechst_Flip1Slide1',...
    '2022-05-03-yw164-180m25C_Kni647_yw174-190m25C_Bcd488_Hb546_Hoechst_Flip1Slide2',...
    '2022-05-05-yw164-180m25C_Kni647_yw174-190m25C_Bcd488_Hb546_Hoechst_Flip1Slide4'};
AllSetInfo.Temperatures(1, end+1) = 25;
AllSetInfo.Replicates(1, end+1) = 0;
AllSetInfo.Flipped(1, end+1) = 1;
AllSetInfo.SetLabels{end+1} = SetLabel;
AllSetInfo.PlotLabels{end+1} = PlotLabel;
AllSetInfo.Prefixes{end+1} = SetPrefixes;
AllSetInfo.DorsalAvgAPKnirpsBackgroundCutoff{end+1} = [NaN 4 4.5];
AllSetInfo.UsePeakRatioToIdentifyControls{end+1} = logical([0 1 1]);
AllSetInfo.DorsalAvgAPKnirpsPeakValueCutoff{end+1} = [12 5 7.5];
AllSetInfo.UsePeakValueToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsMeanValueCutoff{end+1} = [9.5 4.5 5];
AllSetInfo.UseMeanValueToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.MinimumFixCorrectedDeltaFC{end+1} = [0 0 0];
%%
% % T=27.5C, Replicate 1
SetLabel = 'FixedT27_5CReplicate1';
PlotLabel = '27.5°C Replicate 1';
SetPrefixes = {'2022-04-02-yw84-160m27_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide10',...%'2022-04-05-yw84-160m27_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide9',...
    '2022-04-06-yw84-160m27_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide3',...
    '2022-04-08-yw84-160m27_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide7',...%'2022-04-10-yw84-160m27_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide11'...
    };

AllSetInfo.Temperatures(1, end+1) = 27.5;
AllSetInfo.Replicates(1, end+1) = 1;
AllSetInfo.Flipped(1, end+1) = 0;
AllSetInfo.SetLabels{end+1} = SetLabel;
AllSetInfo.PlotLabels{end+1} = PlotLabel;
AllSetInfo.Prefixes{end+1} = SetPrefixes;
AllSetInfo.DorsalAvgAPKnirpsBackgroundCutoff{end+1} = [0.8 7.5 6];
AllSetInfo.UsePeakRatioToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsPeakValueCutoff{end+1} = [1 10 10];
AllSetInfo.UsePeakValueToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsMeanValueCutoff{end+1} = [1 8 8];
AllSetInfo.UseMeanValueToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.MinimumFixCorrectedDeltaFC{end+1} = [0 0 0];
% % T=27.5C, Replicate 2 NEW
SetLabel = 'FixedT27_5CReplicate2';
PlotLabel = '27.5°C Replicate 2';
SetPrefixes = {'2022-05-08-yw104-160m27_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep2Slide1',...
'2022-05-10-yw104-160m27_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep2Slide7',...
'2022-05-10-yw104-160m27_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep2Slide7B'};
AllSetInfo.Temperatures(1, end+1) = 27.5;
AllSetInfo.Replicates(1, end+1) = 2;
AllSetInfo.Flipped(1, end+1) = 0;
AllSetInfo.SetLabels{end+1} = SetLabel;
AllSetInfo.PlotLabels{end+1} = PlotLabel;
AllSetInfo.Prefixes{end+1} = SetPrefixes;
AllSetInfo.DorsalAvgAPKnirpsBackgroundCutoff{end+1} = [6 5.5 5];
AllSetInfo.UsePeakRatioToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsPeakValueCutoff{end+1} = [10 10 NaN ];
AllSetInfo.UsePeakValueToIdentifyControls{end+1} = logical([1 1 0]);
AllSetInfo.DorsalAvgAPKnirpsMeanValueCutoff{end+1} = [8 6 6];
AllSetInfo.UseMeanValueToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.MinimumFixCorrectedDeltaFC{end+1} = [0 0 0];


% % T=27.5C, Flipped Stain
SetLabel = 'FixedT27_5CFlipped';
PlotLabel = '27.5°C Flipped Stain';
SetPrefixes = {'2022-05-01-yw144-160m27_5C_Kni647_yw174-190m25C_Bcd488_Hb546_Hoechst_Flip1Slide1',...
'2022-05-03-yw144-160m27_5C_Kni647_yw174-190m25C_Bcd488_Hb546_Hoechst_Flip1Slide3',...
'2022-05-06-yw144-160m27_5C_Kni647_yw174-190m25C_Bcd488HyD1_Hb546_Hoechst_Flip1Slide2'};
AllSetInfo.Temperatures(1, end+1) = 27.5;
AllSetInfo.Replicates(1, end+1) = 0;
AllSetInfo.Flipped(1, end+1) = 1;
AllSetInfo.SetLabels{end+1} = SetLabel;
AllSetInfo.PlotLabels{end+1} = PlotLabel;
AllSetInfo.Prefixes{end+1} = SetPrefixes;
AllSetInfo.DorsalAvgAPKnirpsBackgroundCutoff{end+1} = [8 7.5 4];
AllSetInfo.UsePeakRatioToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsPeakValueCutoff{end+1} = [12 10 6];
AllSetInfo.UsePeakValueToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsMeanValueCutoff{end+1} = [8.5 8 5];
AllSetInfo.UseMeanValueToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.MinimumFixCorrectedDeltaFC{end+1} = [4 0 0];
%%
% % T=22.5C, Replicate 1
SetLabel = 'FixedT22_5CReplicate1';
PlotLabel = '22.5°C Replicate 1';
SetPrefixes = {'2022-04-02-yw120-216m22_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide2',...
    '2022-04-06-yw120-216m22_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide3',...
    '2022-04-07-yw120-216m22_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide4',...
    '2022-04-08-yw120-216m22_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide11',...
    '2022-04-09-yw120-216m22_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide12'};
AllSetInfo.Temperatures(1, end+1) = 22.5;
AllSetInfo.Replicates(1, end+1) = 1;
AllSetInfo.Flipped(1, end+1) = 0;
AllSetInfo.SetLabels{end+1} = SetLabel;
AllSetInfo.PlotLabels{end+1} = PlotLabel;
AllSetInfo.Prefixes{end+1} = SetPrefixes;
AllSetInfo.DorsalAvgAPKnirpsBackgroundCutoff{end+1} = [0.8 7 6.5 5 7];
AllSetInfo.UsePeakRatioToIdentifyControls{end+1} = logical([1 1 1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsPeakValueCutoff{end+1} = [1 10 10 6 8];
AllSetInfo.UsePeakValueToIdentifyControls{end+1} = logical([1 1 1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsMeanValueCutoff{end+1} = [1 8 8 6 7.5];
AllSetInfo.UseMeanValueToIdentifyControls{end+1} = logical([1 1 1 1 1]);
AllSetInfo.MinimumFixCorrectedDeltaFC{end+1} = [0 0 0 0 0];

% % T=22.5C, Replicate 2 
SetLabel = 'FixedT22_5CReplicate2';
PlotLabel = '22.5°C Replicate 2';
SetPrefixes = {'2022-05-08-yw140-216m22_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep2Slide3',...
    '2022-05-10-yw140-216m22_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep2Slide6'};
AllSetInfo.Temperatures(1, end+1) = 22.5;
AllSetInfo.Replicates(1, end+1) = 2;
AllSetInfo.Flipped(1, end+1) = 0;
AllSetInfo.SetLabels{end+1} = SetLabel;
AllSetInfo.PlotLabels{end+1} = PlotLabel;
AllSetInfo.Prefixes{end+1} = SetPrefixes;
AllSetInfo.DorsalAvgAPKnirpsBackgroundCutoff{end+1} = [7 6];
AllSetInfo.UsePeakRatioToIdentifyControls{end+1} = logical([1 1]);
AllSetInfo.DorsalAvgAPKnirpsPeakValueCutoff{end+1} = [10 7];
AllSetInfo.UsePeakValueToIdentifyControls{end+1} = logical([1 1]);
AllSetInfo.DorsalAvgAPKnirpsMeanValueCutoff{end+1} = [7 7];
AllSetInfo.UseMeanValueToIdentifyControls{end+1} = logical([1 1]);
AllSetInfo.MinimumFixCorrectedDeltaFC{end+1} = [0 0];

% % T=22.5C, Flipped Stain 
SetLabel = 'FixedT22_5CFlipped';
PlotLabel = '22.5°C Flipped Stain';
SetPrefixes = {'2022-05-04-yw200-216m22_5C_Kni647_yw174-190m25C_Bcd488_Hb546_Hoechst_Flip1Slide2',...
    '2022-05-05-yw200-216m22_5C_Kni647_yw174-190m25C_Bcd488HyD1_Hb546_Hoechst_Flip1Slide1',...
    '2022-05-06-yw200-216m22_5C_Kni647_yw174-190m25C_Bcd488HyD1_Hb546_Hoechst_Flip1Slide3'};
AllSetInfo.Temperatures(1, end+1) = 22.5;
AllSetInfo.Replicates(1, end+1) = 0;
AllSetInfo.Flipped(1, end+1) = 1;
AllSetInfo.SetLabels{end+1} = SetLabel;
AllSetInfo.PlotLabels{end+1} = PlotLabel;
AllSetInfo.Prefixes{end+1} = SetPrefixes;
AllSetInfo.DorsalAvgAPKnirpsBackgroundCutoff{end+1} = [7 5 7.2];
AllSetInfo.UsePeakRatioToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsPeakValueCutoff{end+1} = [10 7.5 10];
AllSetInfo.UsePeakValueToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsMeanValueCutoff{end+1} = [8 6 8];
AllSetInfo.UseMeanValueToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.MinimumFixCorrectedDeltaFC{end+1} = [0 0 3];
%%
% % T=20C, Replicate 1
SetLabel = 'FixedT20CReplicate1';
PlotLabel = '20°C Replicate 1';
SetPrefixes = {'2022-04-02-yw144-280m20C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide6',...
    '2022-04-05-yw144-280m20C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide1',...
    '2022-04-06-yw144-280m20C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide9'};
AllSetInfo.Temperatures(1, end+1) = 20;
AllSetInfo.Replicates(1, end+1) = 1;
AllSetInfo.Flipped(1, end+1) = 0;
AllSetInfo.SetLabels{end+1} = SetLabel;
AllSetInfo.PlotLabels{end+1} = PlotLabel;
AllSetInfo.Prefixes{end+1} = SetPrefixes;
AllSetInfo.DorsalAvgAPKnirpsBackgroundCutoff{end+1} = [0.6 7.5 7];
AllSetInfo.UsePeakRatioToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsPeakValueCutoff{end+1} = [1 10 10];
AllSetInfo.UsePeakValueToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsMeanValueCutoff{end+1} = [1 9 8];
AllSetInfo.UseMeanValueToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.MinimumFixCorrectedDeltaFC{end+1} = [0 0 0];

% % T=20C, Replicate 2
SetLabel = 'FixedT20CReplicate2';
PlotLabel = '20°C Replicate 2';
SetPrefixes = {'2022-05-09-yw184-280m20C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep2Slide3',...
'2022-05-11-yw184-280m20C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep2Slide1',...
'2022-05-12-yw184-280m20C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep2Slide6'}; % Third set has no control embryos
AllSetInfo.Temperatures(1, end+1) = 20;
AllSetInfo.Replicates(1, end+1) = 2;
AllSetInfo.Flipped(1, end+1) = 0;
AllSetInfo.SetLabels{end+1} = SetLabel;
AllSetInfo.PlotLabels{end+1} = PlotLabel;
AllSetInfo.Prefixes{end+1} = SetPrefixes;
AllSetInfo.DorsalAvgAPKnirpsBackgroundCutoff{end+1} = [6 6 5];
AllSetInfo.UsePeakRatioToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsPeakValueCutoff{end+1} = [7 8 5];
AllSetInfo.UsePeakValueToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsMeanValueCutoff{end+1} = [6 6 5];
AllSetInfo.UseMeanValueToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.MinimumFixCorrectedDeltaFC{end+1} = [0 0 0];


% % T=20C, Flipped Stain
SetLabel = 'FixedT20CFlipped';
PlotLabel = '20°C Flipped Stain';
SetPrefixes = {'2022-04-30-yw244-280m20C_Kni647_yw174-190m25C_Bcd488_Hb546_Hoechst_Flip1Slide1',...
'2022-05-02-yw244-280m20C_Kni647_yw174-190m25C_Bcd488_Hb546_Hoechst_Flip1Slide4',...
'2022-05-02-yw244-280m20C_Kni647_yw174-190m25C_Bcd488_Hb546_Hoechst_Flip1Slide4B'};
AllSetInfo.Temperatures(1, end+1) = 20;
AllSetInfo.Replicates(1, end+1) = 0;
AllSetInfo.Flipped(1, end+1) = 1;
AllSetInfo.SetLabels{end+1} = SetLabel;
AllSetInfo.PlotLabels{end+1} = PlotLabel;
AllSetInfo.Prefixes{end+1} = SetPrefixes;
AllSetInfo.DorsalAvgAPKnirpsBackgroundCutoff{end+1} = [12 7 6.5];
AllSetInfo.UsePeakRatioToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsPeakValueCutoff{end+1} = [13 10 8];
AllSetInfo.UsePeakValueToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsMeanValueCutoff{end+1} = [12 8 8];
AllSetInfo.UseMeanValueToIdentifyControls{end+1} =logical([1 1 1]);
AllSetInfo.MinimumFixCorrectedDeltaFC{end+1} = [0 0 0];

%%
    
% % T=17.5C, Replicate 1
SetLabel = 'FixedT17_5CReplicate1';
PlotLabel = '17.5°C Replicate 1';
SetPrefixes ={'2022-04-03-yw168-344m17_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide6',...
    '2022-04-05-yw168-344m17_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide10',...
    '2022-04-08-yw168-344m17_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide7'};
AllSetInfo.Temperatures(1, end+1) = 17.5;
AllSetInfo.Replicates(1, end+1) = 1;
AllSetInfo.Flipped(1, end+1) = 0;
AllSetInfo.SetLabels{end+1} = SetLabel;
AllSetInfo.PlotLabels{end+1} = PlotLabel;
AllSetInfo.Prefixes{end+1} = SetPrefixes;
AllSetInfo.DorsalAvgAPKnirpsBackgroundCutoff{end+1} = [0.8 7.5 6];
AllSetInfo.UsePeakRatioToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsPeakValueCutoff{end+1} = [1 10 10];
AllSetInfo.UsePeakValueToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsMeanValueCutoff{end+1} = [1 8 7];
AllSetInfo.UseMeanValueToIdentifyControls{end+1} =logical([1 1 1]);
AllSetInfo.MinimumFixCorrectedDeltaFC{end+1} = [0 0 0];

% % T=17.5C, Replicate 2
SetLabel = 'FixedT17_5CReplicate2';
PlotLabel = '17.5°C Replicate 2';
SetPrefixes ={'2022-05-07-yw208-344m17_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep2Slide2',...
'2022-05-09-yw208-344m17_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep2Slide1',...
'2022-05-12-yw208-344m17_5C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep2Slide5'};
AllSetInfo.Temperatures(1, end+1) = 17.5;
AllSetInfo.Replicates(1, end+1) = 2;
AllSetInfo.Flipped(1, end+1) = 0;
AllSetInfo.SetLabels{end+1} = SetLabel;
AllSetInfo.PlotLabels{end+1} = PlotLabel;
AllSetInfo.Prefixes{end+1} = SetPrefixes;
AllSetInfo.DorsalAvgAPKnirpsBackgroundCutoff{end+1} = [5.5 6 4];
AllSetInfo.UsePeakRatioToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsPeakValueCutoff{end+1} = [7 7.5 5];
AllSetInfo.UsePeakValueToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsMeanValueCutoff{end+1} = [6 7 5];
AllSetInfo.UseMeanValueToIdentifyControls{end+1} =logical([1 1 1]);
AllSetInfo.MinimumFixCorrectedDeltaFC{end+1} = [0 0 0];

% % T=17.5C, Flipped Stain
SetLabel = 'FixedT17_5CFlipped';
PlotLabel = '17.5°C Flipped Stain';
SetPrefixes = {'2022-04-30-yw308-344m17_5C_Kni647_yw174-190m25C_Bcd488_Hb546_Hoechst_Flip1Slide1',...
'2022-05-01-yw308-344m17_5C_Kni647_yw174-190m25C_Bcd488_Hb546_Hoechst_Flip1Slide2',... % For second slide, use peak condition of 1.2 with FixCorrectedDeltaFC_um > 5
'2022-05-05-yw308-344m17_5C_Kni647_yw174-190m25C_Bcd488HyD1_Hb546_Hoechst_Flip1Slide4'};
AllSetInfo.Temperatures(1, end+1) = 17.5;
AllSetInfo.Replicates(1, end+1) = 0;
AllSetInfo.Flipped(1, end+1) = 1;
AllSetInfo.SetLabels{end+1} = SetLabel;
AllSetInfo.PlotLabels{end+1} = PlotLabel;
AllSetInfo.Prefixes{end+1} = SetPrefixes;
AllSetInfo.DorsalAvgAPKnirpsBackgroundCutoff{end+1} = [7 NaN 4.2];
AllSetInfo.UsePeakRatioToIdentifyControls{end+1} = logical([1 0 1]);
AllSetInfo.DorsalAvgAPKnirpsPeakValueCutoff{end+1} = [7.75 9.8 5];
AllSetInfo.UsePeakValueToIdentifyControls{end+1} = logical([1 1 1]);
AllSetInfo.DorsalAvgAPKnirpsMeanValueCutoff{end+1} = [7 8.9 4.5];
AllSetInfo.UseMeanValueToIdentifyControls{end+1} =logical([1 1 1]);
AllSetInfo.MinimumFixCorrectedDeltaFC{end+1} = [0 0 0];

AllSetInfo.Flipped = logical(AllSetInfo.Flipped);