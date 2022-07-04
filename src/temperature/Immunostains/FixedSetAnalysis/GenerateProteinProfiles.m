function  [AllCompiledEmbryos]   = GenerateProteinProfiles(version)
%%

if ~exist('version', 'var')
    version = 11;
end
SizeDataPath = 'S:/Gabriella/Dropbox/EmbryoSizeMeasurements/EmbryoSizeData.mat';
AllSetsProfFigPath = 'S:/Gabriella/Dropbox/ProteinProfiles/Figures/';
AllCompiledPath = ['S:/Gabriella/Dropbox/ProteinProfiles/V',num2str(version),'Profiles/'];
AllSetsVersionCombinedEmbryosPath = ['S:/Gabriella/Dropbox/ProteinProfiles/V',num2str(version),'CompiledEmbryos/'];
AllSetsCombinedEmbryosPath = ['S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/'];
load(SizeDataPath,'APLengths','DVLengths','VentralDistances','DorsalDistances','NGoodEmbryos',...
    'MeanAPLength','StdAPLength','MeanDVLength','StdDVLength','AspectRatios');
MeanAspectRatio = mean(AspectRatios);
StdAspectRatio = std(AspectRatios);
[DubuisTimes, DubuisMeanProfile, ~] = getMembraneFurrowProfiles( 'dubuis');
[yw25CTimes, yw25CProfile, yw25CSE] = getMembraneFurrowProfiles( 'yw25csquished_nopv');
[hisrfp25CTimes, hisrfp25CProfile, hisrfp25CSE] = getMembraneFurrowProfiles( 'hisrfp25c_nopv');
AllSetInfo = GetFixedSetPrefixInfo;
 
NumSets = length(AllSetInfo.Temperatures);

NChannels = 5;



APbins = 0:0.025:1;
NumAPbins = length(APbins);

NarrowAPbins = 0:0.0125:1;
NumNarrowAPbins = length(NarrowAPbins);
AllCompiledEmbryos = {};

%%
for exp_index = 1:length(AllSetInfo.Temperatures)

if AllSetInfo.Flipped(exp_index)
    FlipString = 'yes';
else
    FlipString = 'no';
end
disp(['T = ', num2str(AllSetInfo.Temperatures(exp_index)), ', Rep: ', num2str(AllSetInfo.Replicates(exp_index)), ', Flipped: ', FlipString])
SetLabel = AllSetInfo.SetLabels{exp_index};
PlotLabel = AllSetInfo.PlotLabels{exp_index};
SetPrefixes = AllSetInfo.Prefixes{exp_index};
SetIsFlipped = AllSetInfo.Flipped(exp_index);
ProfFigPath = [AllSetsProfFigPath, SetLabel];
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];
if ~isdir(ProfFigPath)
    mkdir(ProfFigPath)
end
if ~isdir(OutEmbryoPath)
    mkdir(OutEmbryoPath)
end
liveExperiments = cell(1, length(SetPrefixes));
for i = 1:length(SetPrefixes)
    liveExperiments{i} = LiveExperiment(SetPrefixes{i});
end
FixedPixelSize_um = liveExperiments{1}.pixelSize_um;
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
%load(CEoutpath, 'CompiledEmbryos');
CompiledEmbryos = CombineCompiledEmbryos(SetPrefixes);
CompiledEmbryos = AddEmbryoTimingInfo(CompiledEmbryos, exp_index);

NEmbryos = length(CompiledEmbryos.Approved);
AllEmbryos = 1:NEmbryos;
NC14Indices = find(CompiledEmbryos.IsNC14);
NC13Indices = find(CompiledEmbryos.IsNC13);
NC13NC14Indices = find(CompiledEmbryos.IsNC13orNC14);
NumEmbryosNC14 = CompiledEmbryos.NumEmbryosNC14;

KnirpsIndex = 4;
CompiledEmbryos = PartitionEmbryosTestControl(CompiledEmbryos, exp_index, KnirpsIndex);


CompiledEmbryos = RescaleSlideFluos(CompiledEmbryos, exp_index);
AllCompiledEmbryos{exp_index} = CompiledEmbryos;
end

AllCompiledPath = ['S:/Gabriella/Dropbox/ProteinProfiles/V',num2str(version),'Profiles/'];
if ~isdir(AllCompiledPath)
    mkdir(AllCompiledPath)
end
save([AllCompiledPath, 'AllCompiledEmbryos.Mat']);

if ~isdir(AllSetsVersionCombinedEmbryosPath)
    mkdir(AllSetsVersionCombinedEmbryosPath)
end
for exp_index=1:NumSets
CompiledEmbryos = AllCompiledEmbryos{exp_index};
SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsVersionCombinedEmbryosPath, filesep, SetLabel];
if ~isdir(OutEmbryoPath)
mkdir(OutEmbryoPath);
end
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
save(CEoutpath, 'CompiledEmbryos');
end
%%

for exp_index = 1:NumSets
    AllCompiledEmbryos{exp_index} = AddFitBicoidProfiles(AllCompiledEmbryos{exp_index} , exp_index);
    AllCompiledEmbryos{exp_index} = AddNonBicoidZeroCorrection(AllCompiledEmbryos{exp_index} , exp_index);
end
%%
for exp_index = 1:NumSets
    disp(['i = ', num2str(exp_index)])
    AllCompiledEmbryos{exp_index} = AddBootstrappedProfiles(AllCompiledEmbryos{exp_index} , exp_index);
end
AllCompiledEmbryos = FitScalingFactorsToTestProfiles(AllCompiledEmbryos);

for exp_index = 1:NumSets
    disp(['i = ', num2str(exp_index)])
    AllCompiledEmbryos{exp_index} = AddTestMasterFitRescaledProfiles(AllCompiledEmbryos{exp_index} , exp_index);
end


for exp_index = 1:NumSets
    disp(['i = ', num2str(exp_index)])
    AllCompiledEmbryos{exp_index} = AddBootstrappedTestFitProfiles(AllCompiledEmbryos{exp_index} , exp_index);
end




AllCompiledEmbryos = FitComboScalingFactorsToMasterProfiles(AllCompiledEmbryos);

for exp_index = 1:NumSets
    AllCompiledEmbryos{exp_index}  = AddUniversalScaledTestMasterFitRescaledProfiles(AllCompiledEmbryos{exp_index}, exp_index);
end
%%
parfor exp_index = 1:NumSets
    %disp(['Embryo Index = ', num2str(exp_index)])
    tic
    AllCompiledEmbryos{exp_index}  =  AddBootstrappedUnivScaledTestFitProfiles(AllCompiledEmbryos{exp_index} , exp_index);
    toc
end
%%
if ~isdir(AllCompiledPath)
    mkdir(AllCompiledPath)
end
save([AllCompiledPath, 'AllCompiledEmbryos.Mat'], 'AllCompiledEmbryos');%,


%%
% % 
% for exp_index = 1:NumSets
%     disp(['i = ', num2str(exp_index)])
%     SetLabel = AllSetInfo.SetLabels{exp_index};
%     PlotLabel = AllSetInfo.PlotLabels{exp_index};
%     SetPrefixes = AllSetInfo.Prefixes{exp_index};
%     SetIsFlipped = AllSetInfo.Flipped(exp_index);
%     ProfFigPath = [AllSetsProfFigPath, SetLabel];
%     OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];
%     CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
%     AllCompiledEmbryos{exp_index} = FitScalingFactorsToMasterProfiles(AllCompiledEmbryos{exp_index}, exp_index);
%     %AllCompiledEmbryos{exp_index} = AddUniversalScalingFactorsToMasterProfiles(AllCompiledEmbryos{exp_index}  , exp_index);
%     AllCompiledEmbryos{exp_index} = AddBootstrappedFitProfiles(AllCompiledEmbryos{exp_index}  , exp_index);
%     CompiledEmbryos = AllCompiledEmbryos{exp_index} ;
%     CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
%     save(CEoutpath, 'CompiledEmbryos');
% end

% TStarts = 0:10:50;
% TEnds = 10:10:60;
% NTimeBins = length(TStarts);
% CompiledEmbryos.Dubuis = {};
% CompiledEmbryos.Dubuis.TimeAveragedDorsalTestProfiles = {};
% CompiledEmbryos.Dubuis.StdTimeAveragedDorsalTestProfiles = {};
% CompiledEmbryos.Dubuis.CountTimeAveragedDorsalTestProfiles = {};
% CompiledEmbryos.Dubuis.TimeAveragedDorsalControlProfiles = {};
% CompiledEmbryos.Dubuis.StdTimeAveragedDorsalControlProfiles = {};
% CompiledEmbryos.Dubuis.CountTimeAveragedDorsalControlProfiles = {};
% CompiledEmbryos.Dubuis.TimeAveragedDorsalTestProfiles.NC14 = NaN(NTimeBins, NumAPbins, NChannels);
% CompiledEmbryos.Dubuis.TimeAveragedDorsalControlProfiles.NC14 = NaN(NTimeBins, NumAPbins, NChannels);
% CompiledEmbryos.Dubuis.StdTimeAveragedDorsalTestProfiles.NC14 = NaN(NTimeBins, NumAPbins, NChannels);
% CompiledEmbryos.Dubuis.StdTimeAveragedDorsalControlProfiles.NC14 = NaN(NTimeBins, NumAPbins, NChannels);
% CompiledEmbryos.Dubuis.CountTimeAveragedDorsalTestProfiles.NC14 = NaN(NTimeBins, NumAPbins, NChannels);
% CompiledEmbryos.Dubuis.CountTimeAveragedDorsalControlProfiles.NC14 = NaN(NTimeBins, NumAPbins, NChannels);
% CompiledEmbryos.Dubuis.TimeAveragedNarrowDorsalTestProfiles = {};
% CompiledEmbryos.Dubuis.StdTimeAveragedNarrowDorsalTestProfiles = {};
% CompiledEmbryos.Dubuis.CountTimeAveragedNarrowDorsalTestProfiles = {};
% CompiledEmbryos.Dubuis.TimeAveragedNarrowDorsalControlProfiles = {};
% CompiledEmbryos.Dubuis.StdTimeAveragedNarrowDorsalControlProfiles = {};
% CompiledEmbryos.Dubuis.CountTimeAveragedNarrowDorsalControlProfiles = {};
% CompiledEmbryos.Dubuis.TimeAveragedNarrowDorsalTestProfiles.NC14 = NaN(NTimeBins, NumNarrowAPbins, NChannels);
% CompiledEmbryos.Dubuis.TimeAveragedNarrowDorsalControlProfiles.NC14 = NaN(NTimeBins, NumNarrowAPbins, NChannels);
% CompiledEmbryos.Dubuis.StdTimeAveragedNarrowDorsalTestProfiles.NC14 = NaN(NTimeBins, NumNarrowAPbins, NChannels);
% CompiledEmbryos.Dubuis.StdTimeAveragedNarrowDorsalControlProfiles.NC14 = NaN(NTimeBins, NumNarrowAPbins, NChannels);
% 
% for i = 1:NTimeBins
%     TFtimeBin = CompiledEmbryos.DubuisEmbryoTimes >= TStarts(i) & CompiledEmbryos.DubuisEmbryoTimes < TEnds(i)  & CompiledEmbryos.IsNC14;
%     TFtestBin = TFtimeBin & CompiledEmbryos.TestSetEmbryos;
%     TFcontrolBin = TFtimeBin & CompiledEmbryos.ControlSetEmbryos;
%     for ch_index = 2:NChannels
%         if sum(TFtestBin) > 1
%             CompiledEmbryos.Dubuis.TimeAveragedDorsalTestProfiles.NC14(i,:,ch_index) = mean(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(TFtestBin,:,ch_index), 'omitnan');
%             CompiledEmbryos.Dubuis.StdTimeAveragedDorsalTestProfiles.NC14(i,:,ch_index) = std(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(TFtestBin,:,ch_index), 'omitnan');
%             CompiledEmbryos.Dubuis.TimeAveragedNarrowDorsalTestProfiles.NC14(i,:,ch_index) = mean(CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles(TFtestBin,:,ch_index), 'omitnan');
%             CompiledEmbryos.Dubuis.StdTimeAveragedNarrowDorsalTestProfiles.NC14(i,:,ch_index) = std(CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles(TFtestBin,:,ch_index), 'omitnan');
%         elseif sum(TFtestBin) == 1
%            CompiledEmbryos.Dubuis.TimeAveragedDorsalTestProfiles.NC14(i,:,ch_index) = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(TFtestBin,:,ch_index);
%             CompiledEmbryos.Dubuis.StdTimeAveragedDorsalTestProfiles.NC14(i,:,ch_index) = zeros(size(CompiledEmbryos.Dubuis.StdTimeAveragedDorsalTestProfiles.NC14(i,:,ch_index) )); 
%             CompiledEmbryos.Dubuis.TimeAveragedNarrowDorsalTestProfiles.NC14(i,:,ch_index) = CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles(TFtestBin,:,ch_index);
%             CompiledEmbryos.Dubuis.StdTimeAveragedNarrowDorsalTestProfiles.NC14(i,:,ch_index) = zeros(size(CompiledEmbryos.Dubuis.StdTimeAveragedNarrowDorsalTestProfiles.NC14(i,:,ch_index) )); 
%         end
%         CompiledEmbryos.Dubuis.CountTimeAveragedDorsalTestProfiles.NC14(i,:,ch_index)  = sum(TFtestBin);
%         
%         if sum(TFcontrolBin) > 1
%             CompiledEmbryos.Dubuis.TimeAveragedDorsalControlProfiles.NC14(i,:,ch_index) = mean(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(TFcontrolBin,:,ch_index), 'omitnan');
%             CompiledEmbryos.Dubuis.StdTimeAveragedDorsalControlProfiles.NC14(i,:,ch_index) = std(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(TFcontrolBin,:,ch_index), 'omitnan');
%             CompiledEmbryos.Dubuis.TimeAveragedNarrowDorsalControlProfiles.NC14(i,:,ch_index) = mean(CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles(TFcontrolBin,:,ch_index), 'omitnan');
%             CompiledEmbryos.Dubuis.StdTimeAveragedNarrowDorsalControlProfiles.NC14(i,:,ch_index) = std(CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles(TFcontrolBin,:,ch_index), 'omitnan');
%         elseif sum(TFcontrolBin) == 1
%            CompiledEmbryos.Dubuis.TimeAveragedDorsalControlProfiles.NC14(i,:,ch_index) = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(TFcontrolBin,:,ch_index);
%             CompiledEmbryos.Dubuis.StdTimeAveragedDorsalControlProfiles.NC14(i,:,ch_index) = zeros(size(CompiledEmbryos.Dubuis.StdTimeAveragedDorsalTestProfiles.NC14(i,:,ch_index) )); 
%             CompiledEmbryos.Dubuis.TimeAveragedNarrowDorsalControlProfiles.NC14(i,:,ch_index) = CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles(TFcontrolBin,:,ch_index);
%             CompiledEmbryos.Dubuis.StdTimeAveragedNarrowDorsalControlProfiles.NC14(i,:,ch_index) = zeros(size(CompiledEmbryos.Dubuis.StdTimeAveragedNarrowDorsalTestProfiles.NC14(i,:,ch_index) )); 
%         end
%         CompiledEmbryos.Dubuis.CountTimeAveragedDorsalControlProfiles.NC14(i,:,ch_index)  = sum(TFcontrolBin);
%     end
% end
% 
% CompiledEmbryos = AddSmoothedProfiles(CompiledEmbryos);
% CompiledEmbryos = AddBinnedProfiles(CompiledEmbryos);
% 
% CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
% save(CEoutpath, 'CompiledEmbryos');
% AllCompiledEmbryos{exp_index} = CompiledEmbryos;
% 
% end
% %%
% save([AllCompiledPath, 'AllCompiledEmbryos.Mat'], 'AllCompiledEmbryos');
% if ~isdir(AllSetsVersionCombinedEmbryosPath)
%     mkdir(AllSetsVersionCombinedEmbryosPath)
% end
% Construct25CMasterSet;
% 
% %%
% AllCompiledEmbryos = CalculateSetRescalingFactors(AllCompiledEmbryos);
% 
% 
% for exp_index = 1:NumSets
%     disp(['Embryo Index = ', num2str(exp_index)]);
%     AllCompiledEmbryos{exp_index} = AddUniversalScalingProfiles( AllCompiledEmbryos{exp_index}, exp_index);
% end
%  [AllCompiledEmbryos, Universal_Imins, Universal_Imaxs, Imins, Imaxs]  = AddNormalizedProfiles(AllCompiledEmbryos);
% 
% for exp_index = 1:NumSets
%     disp(['Embryo Index = ', num2str(exp_index)]);
%     AllCompiledEmbryos{exp_index} = AddNC13NormalizedProfiles( AllCompiledEmbryos{exp_index}, exp_index);
% end

% AllCompiledPath = ['S:/Gabriella/Dropbox/ProteinProfiles/V',num2str(version),'Profiles/'];
if ~isdir(AllCompiledPath)
    mkdir(AllCompiledPath)
end
save([AllCompiledPath, 'AllCompiledEmbryos.Mat'], 'AllCompiledEmbryos');%,'Universal_Imins','Universal_Imaxs','Imins','Imaxs');
% if ~isdir(AllSetsVersionCombinedEmbryosPath)
%     mkdir(AllSetsVersionCombinedEmbryosPath)
% end
% for exp_index=1:NumSets
% CompiledEmbryos = AllCompiledEmbryos{exp_index};
% SetLabel = AllSetInfo.SetLabels{exp_index};
% OutEmbryoPath = [AllSetsVersionCombinedEmbryosPath, filesep, SetLabel];
% if ~isdir(OutEmbryoPath)
% mkdir(OutEmbryoPath);
% end
% CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
% save(CEoutpath, 'CompiledEmbryos');
% end