clear all
SizeDataPath = 'S:/Gabriella/Dropbox/EmbryoSizeMeasurements/EmbryoSizeData.mat';
AllSetsProfFigPath = 'S:/Gabriella/Dropbox/ProteinProfiles/Figures/';
AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/V1CompiledEmbryos/';
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

ChannelNames = {'', 'Hoechst', 'Bicoid', 'Knirps', 'Hunchback'};

RefSets = [14, 3, 5, 14];
APbins = 0:0.025:1;
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
load(CEoutpath, 'CompiledEmbryos');
AllCompiledEmbryos{exp_index} = CompiledEmbryos;
end


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
for exp_index = 1:NumSets
    disp(['Embryo Index = ', num2str(exp_index)]);
    AllCompiledEmbryos{exp_index} = AddUniversalScalingProfiles( AllCompiledEmbryos{exp_index}, exp_index);
end
%%
 [AllCompiledEmbryos, Universal_Imins, Universal_Imaxs, Imins, Imaxs]  = AddNormalizedProfiles(AllCompiledEmbryos);
%%
for exp_index = 1:NumSets
    disp(['Embryo Index = ', num2str(exp_index)]);
    AllCompiledEmbryos{exp_index} = AddNC13NormalizedProfiles( AllCompiledEmbryos{exp_index}, exp_index);
end


%%
if ~isdir([ProfFigPath, filesep, 'SetComps'])
    mkdir([ProfFigPath, filesep, 'SetComps'])
end
close all
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
y_positions = [0, 0.5, 0.3, 0.6, 0.8];
unique_temperatures = unique(AllSetInfo.Temperatures);
for temp_index = 1:length(unique_temperatures)
    ReplicateIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 0);
    exp1 = ReplicateIndices(1); 
    exp2 = ReplicateIndices(2);
for ch_index = [3:5]
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.mean(:,:, ch_index);
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.mean(:,:, ch_index);
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 3 & ProfCounts2 >= 3);

    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);

    scatter(PlotProfSet1,PlotProfSet2,...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    hold on 


    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k')
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = [num2str(unique_temperatures(temp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp1))];
ylab = [num2str(unique_temperatures(temp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp2))];

xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, ChannelNames{ch_index}, 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [ProfFigPath, filesep, 'SetComps',filesep,'T', strrep(num2str( num2str(unique_temperatures(temp_index))),'.','_'),'Ctest_StandardReplicates_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end

%%
close all
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
y_positions = [0, 0.5, 0.3, 0.6, 0.8];
unique_temperatures = unique(AllSetInfo.Temperatures);
for temp_index = 1:length(unique_temperatures)
    ReplicateIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 0);
    exp1 = ReplicateIndices(1); 
    exp2 = ReplicateIndices(2);
for ch_index = 3:5
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.mean(:,:, ch_index);
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.mean(:,:, ch_index);
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    if isempty(PlotProfSet1)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        AddTitleStr = ' (Low Counts)';
    else
        AddTitleStr = '';
    end
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);

    scatter(PlotProfSet1,PlotProfSet2,...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    hold on 


    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k')
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = ['Control for ', num2str(unique_temperatures(temp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp1))];
ylab = ['Control for ', num2str(unique_temperatures(temp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp2))];

xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [ChannelNames{ch_index}, AddTitleStr], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [ProfFigPath, filesep, 'SetComps',filesep,'T', strrep(num2str( num2str(unique_temperatures(temp_index))),'.','_'),'Ccontrol_StandardReplicates_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end

%%
if ~isdir([ProfFigPath, filesep, 'FlippedSetComps'])
    mkdir([ProfFigPath, filesep, 'FlippedSetComps'])
end
close all
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
y_positions = [0, 0.5, 0.3, 0.6, 0.8];
unique_temperatures = unique(AllSetInfo.Temperatures);
for temp_index = 1:length(unique_temperatures)
    ReplicateIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 0);
    exp1 = ReplicateIndices(1); 
    FlippedIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 1);
    exp2 = FlippedIndices(1);
for ch_index = [3, 5]
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.mean(:,:, ch_index);
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.mean(:,:, ch_index);
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 3 & ProfCounts2 >= 3);
if isempty(PlotProfSet1)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        AddTitleStr = ' (Low Counts)';
    else
        AddTitleStr = '';
    end
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);

    scatter(PlotProfSet1,PlotProfSet2,...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    hold on 


    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k')
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = [num2str(unique_temperatures(temp_index)), 'ºC Replicate 1'];
ylab = [num2str(unique_temperatures(temp_index)), 'ºC Flipped Condition '];

xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [ChannelNames{ch_index}, AddTitleStr], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [ProfFigPath, filesep, 'FlippedSetComps',filesep,'T', strrep(num2str( num2str(unique_temperatures(temp_index))),'.','_'),'Ctest_FlippedRep1_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end

%%
for temp_index = 1:length(unique_temperatures)
    ReplicateIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 0);
    exp1 = ReplicateIndices(2); 
    FlippedIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 1);
    exp2 = FlippedIndices(1);
for ch_index = [3, 5]
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.mean(:,:, ch_index);
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.mean(:,:, ch_index);
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 3 & ProfCounts2 >= 3);
if isempty(PlotProfSet1)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        AddTitleStr = ' (Low Counts)';
    else
        AddTitleStr = '';
    end
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);

    scatter(PlotProfSet1,PlotProfSet2,...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    hold on 


    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k')
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = [num2str(unique_temperatures(temp_index)), 'ºC Replicate 2'];
ylab = [num2str(unique_temperatures(temp_index)), 'ºC Flipped Condition '];

xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [ChannelNames{ch_index}, AddTitleStr], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [ProfFigPath, filesep, 'FlippedSetComps',filesep,'T', strrep(num2str( num2str(unique_temperatures(temp_index))),'.','_'),'Ctest_FlippedRep2_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end

%%
close all
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
y_positions = [0, 0.5, 0.3, 0.6, 0.8];
unique_temperatures = unique(AllSetInfo.Temperatures);
for temp_index = 1:length(unique_temperatures)
    ReplicateIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 0);
    exp1 = ReplicateIndices(1); 
    FlippedIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 1);
    exp2 = FlippedIndices(1);
for ch_index = [3, 5]
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.mean(:,:, ch_index);
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.mean(:,:, ch_index);
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 3 & ProfCounts2 >= 3);
if isempty(PlotProfSet1)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        AddTitleStr = ' (Low Counts)';
    else
        AddTitleStr = '';
    end
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);

    scatter(PlotProfSet1,PlotProfSet2,...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    hold on 


    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k')
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = ['Control for ', num2str(unique_temperatures(temp_index)), 'ºC Replicate 1'];
ylab = ['Control for ',num2str(unique_temperatures(temp_index)), 'ºC Flipped Condition '];

xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [ChannelNames{ch_index}, AddTitleStr], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [ProfFigPath, filesep, 'FlippedSetComps',filesep,'T', strrep(num2str( num2str(unique_temperatures(temp_index))),'.','_'),'Ccontrol_FlippedRep1_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end
%%
close all
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
y_positions = [0, 0.5, 0.3, 0.6, 0.8];
unique_temperatures = unique(AllSetInfo.Temperatures);
for temp_index = 1:length(unique_temperatures)
    ReplicateIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 0);
    exp1 = ReplicateIndices(2); 
    FlippedIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 1);
    exp2 = FlippedIndices(1);
for ch_index = [3, 5]
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.mean(:,:, ch_index);
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.mean(:,:, ch_index);
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 3 & ProfCounts2 >= 3);
if isempty(PlotProfSet1)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        AddTitleStr = ' (Low Counts)';
    else
        AddTitleStr = '';
    end
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);

    scatter(PlotProfSet1,PlotProfSet2,...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    hold on 


    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k')
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = ['Control for ', num2str(unique_temperatures(temp_index)), 'ºC Replicate 2'];
ylab = ['Control for ',num2str(unique_temperatures(temp_index)), 'ºC Flipped Condition '];

xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [ChannelNames{ch_index}, AddTitleStr], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [ProfFigPath, filesep, 'FlippedSetComps',filesep,'T', strrep(num2str( num2str(unique_temperatures(temp_index))),'.','_'),'Ccontrol_FlippedRep2_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end