clear all
SizeDataPath = 'S:/Gabriella/Dropbox/EmbryoSizeMeasurements/EmbryoSizeData.mat';
AllSetsProfFigPath = 'S:/Gabriella/Dropbox/ProteinProfiles/Figures/';
AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';
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

exp_index = find(AllSetInfo.Temperatures == 20 & AllSetInfo.Replicates == 1);

ChannelNames = {'', 'Hoechst', 'Bicoid', 'Knirps', 'Hunchback'};

RefSets = [14, 14, 14, 14];
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
%%
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 177 32 ]/255;%0.9290 0.6940 0.1250]/255;
y_positions = [0, 0.5, 0.3, 0.6, 0.8];
for ch_index = 2:5 % 2:5
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    plot(APbins, CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDubuisValidProfilesTestTF,:,ch_index).',...
        'Color', [.7 .7 .7], 'LineWidth', 2.0);
    hold on 
    plot(APbins, CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDubuisValidProfilesTestTF & ...
        (CompiledEmbryos.DubuisEmbryoTimes >= 37) &...
        (CompiledEmbryos.DubuisEmbryoTimes <= 49),:,ch_index).',...
          'Color', gap_colors(ch_index,:), 'LineWidth', 2.0);
     MeanProfile = mean(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDubuisValidProfilesTestTF &...
         (CompiledEmbryos.DubuisEmbryoTimes >= 37) &...
         (CompiledEmbryos.DubuisEmbryoTimes <= 49),:,ch_index), 'omitnan');
    plot(APbins, MeanProfile, 'k', 'LineWidth', 2.0);
    MaxProfileValue = max(MeanProfile);
    MinProfileValue = min(MeanProfile);
    plot([0, 1], [MaxProfileValue, MaxProfileValue], 'k--','LineWidth', 2.0);
     plot([0, 1], [MinProfileValue, MinProfileValue], 'k--','LineWidth', 2.0);
     ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDubuisValidProfilesTestTF,:,ch_index)))*1.1;
    plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    xlabel('x/L)', 'FontSize', 16)
xlim([.1, .9])

ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, PlotLabel, 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [ProfFigPath, filesep, SetLabel, '_',ChannelNames{ch_index}, '_TestProfilesFig2C.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end

%%
for ch_index = 2:5 % 2:5
    RefBin = find(round(APbins,6) == y_positions(ch_index));
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    scatter(CompiledEmbryos.FixCorrectedDeltaFC_um.mean(CompiledEmbryos.AllDeltaValidProfilesTestTF),...
        CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,RefBin,ch_index).',...
        100, 'MarkerFaceColor', [.7 .7 .7],...
        'MarkerEdgeColor', [.7 .7 .7]);
    hold on 
    scatter(CompiledEmbryos.FixCorrectedDeltaFC_um.mean(CompiledEmbryos.AllDeltaValidProfilesTestTF& ...
        (CompiledEmbryos.DubuisEmbryoTimes >= 37) &...
        (CompiledEmbryos.DubuisEmbryoTimes <= 49)),...
        CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF & ...
        (CompiledEmbryos.DubuisEmbryoTimes >= 37) &...
        (CompiledEmbryos.DubuisEmbryoTimes <= 49),RefBin,ch_index).',...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    MeanProfile = CompiledEmbryos.FixCorrectedSmoothedAvgAPProfiles.Test(:,RefBin,ch_index).';
    plot(CompiledEmbryos.FixCorrectedSmoothedDeltaFCs, MeanProfile, 'k', 'LineWidth', 2.0);
    MaxProfileValue = max(MeanProfile);
    MinProfileValue = min(MeanProfile);
    ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(:,RefBin,ch_index)))*1.1;
    plot([10, 10], [0, ymax], 'k--','LineWidth', 2.0);
     plot([20, 20], [0, ymax], 'k--','LineWidth', 2.0);
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    xlabel('\delta_{FC} (\mum)', 'FontSize', 16)
%xlim([.1, .9])

ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 16)

ylim([0, ymax])
xlim([0 45])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [PlotLabel, ' (Test)'], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [ProfFigPath, filesep, SetLabel, '_',ChannelNames{ch_index}, '_TestProfilesFig2D.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end

%%
for ch_index = 2:5 % 2:5
    RefBin = find(round(APbins,6) == y_positions(ch_index));
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    scatter(CompiledEmbryos.FixCorrectedDeltaFC_um.mean(CompiledEmbryos.AllDeltaValidProfilesControlTF),...
        CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesControlTF,RefBin,ch_index).',...
        100, 'MarkerFaceColor', [.7 .7 .7],...
        'MarkerEdgeColor', [.7 .7 .7]);
    hold on 
    scatter(CompiledEmbryos.FixCorrectedDeltaFC_um.mean(CompiledEmbryos.AllDeltaValidProfilesControlTF& ...
        (CompiledEmbryos.DubuisEmbryoTimes >= 37) &...
        (CompiledEmbryos.DubuisEmbryoTimes <= 49)),...
        CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesControlTF & ...
        (CompiledEmbryos.DubuisEmbryoTimes >= 37) &...
        (CompiledEmbryos.DubuisEmbryoTimes <= 49),RefBin,ch_index).',...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    MeanProfile = CompiledEmbryos.FixCorrectedSmoothedAvgAPProfiles.Control(:,RefBin,ch_index).';
    plot(CompiledEmbryos.FixCorrectedSmoothedDeltaFCs, MeanProfile, 'k', 'LineWidth', 2.0);
    MaxProfileValue = max(MeanProfile);
    MinProfileValue = min(MeanProfile);
    ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(:,RefBin,ch_index)))*1.1;
    plot([10, 10], [0, ymax], 'k--','LineWidth', 2.0);
     plot([20, 20], [0, ymax], 'k--','LineWidth', 2.0);
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    xlabel('\delta_{FC} (\mum)', 'FontSize', 16)
%xlim([.1, .9])

ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 16)

ylim([0, ymax])
xlim([0 45])

DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [PlotLabel, ' (Control)'], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [ProfFigPath, filesep, SetLabel, '_',ChannelNames{ch_index}, '_ControlProfilesFig2D.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end

%%
for ch_index = 2:5 % 2:5
    RefBin = find(round(APbins,6) == y_positions(ch_index));
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    scatter(CompiledEmbryos.DubuisEmbryoTimes(CompiledEmbryos.AllDubuisValidProfilesTestTF),...
        CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDubuisValidProfilesTestTF,RefBin,ch_index).',...
        100, 'MarkerFaceColor', [.7 .7 .7],...
        'MarkerEdgeColor', [.7 .7 .7]);
    hold on 
    scatter(CompiledEmbryos.DubuisEmbryoTimes(CompiledEmbryos.AllDubuisValidProfilesTestTF& ...
        (CompiledEmbryos.FixCorrectedDeltaFC_um.mean >= 5) &...
        (CompiledEmbryos.FixCorrectedDeltaFC_um.mean <= 15)),...
        CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDubuisValidProfilesTestTF & ...
        (CompiledEmbryos.FixCorrectedDeltaFC_um.mean >= 5) &...
        (CompiledEmbryos.FixCorrectedDeltaFC_um.mean <= 15),RefBin,ch_index).',...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    MeanProfile = CompiledEmbryos.DubuisTimeSmoothedAvgAPProfiles.Test(:,RefBin,ch_index).';
    plot(CompiledEmbryos.DubuisSmoothedTimes, MeanProfile, 'k', 'LineWidth', 2.0);
    MaxProfileValue = max(MeanProfile);
    MinProfileValue = min(MeanProfile);
    MinTKeyRegion = min(CompiledEmbryos.DubuisEmbryoTimes(CompiledEmbryos.AllDubuisValidProfilesTestTF& ...
        (CompiledEmbryos.FixCorrectedDeltaFC_um.mean >= 5) &...
        (CompiledEmbryos.FixCorrectedDeltaFC_um.mean <= 15)));
    MaxTKeyRegion = max(CompiledEmbryos.DubuisEmbryoTimes(CompiledEmbryos.AllDubuisValidProfilesTestTF& ...
        (CompiledEmbryos.FixCorrectedDeltaFC_um.mean >= 5) &...
        (CompiledEmbryos.FixCorrectedDeltaFC_um.mean <= 15)));
    ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(:,RefBin,ch_index)))*1.1;
    plot([MinTKeyRegion, MinTKeyRegion], [0, ymax], 'k--','LineWidth', 2.0);
     plot([MaxTKeyRegion, MaxTKeyRegion], [0, ymax], 'k--','LineWidth', 2.0);
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
%xlim([.1, .9])

ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 16)

ylim([0, ymax])
xlim([0 60])

DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [PlotLabel, ' (Test)'], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [ProfFigPath, filesep, SetLabel, '_',ChannelNames{ch_index}, '_TestProfilesFig2D_TimeSmoothed.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end

%%
for ch_index = 2:5 % 2:5
    RefBin = find(round(APbins,6) == y_positions(ch_index));
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    scatter(CompiledEmbryos.DubuisEmbryoTimes(CompiledEmbryos.AllDubuisValidProfilesControlTF),...
        CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDubuisValidProfilesControlTF,RefBin,ch_index).',...
        100, 'MarkerFaceColor', [.7 .7 .7],...
        'MarkerEdgeColor', [.7 .7 .7]);
    hold on 
    scatter(CompiledEmbryos.DubuisEmbryoTimes(CompiledEmbryos.AllDubuisValidProfilesControlTF& ...
        (CompiledEmbryos.FixCorrectedDeltaFC_um.mean >= 5) &...
        (CompiledEmbryos.FixCorrectedDeltaFC_um.mean <= 15)),...
        CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDubuisValidProfilesControlTF & ...
        (CompiledEmbryos.FixCorrectedDeltaFC_um.mean >= 5) &...
        (CompiledEmbryos.FixCorrectedDeltaFC_um.mean <= 15),RefBin,ch_index).',...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    MeanProfile = CompiledEmbryos.DubuisTimeSmoothedAvgAPProfiles.Control(:,RefBin,ch_index).';
    plot(CompiledEmbryos.DubuisSmoothedTimes, MeanProfile, 'k', 'LineWidth', 2.0);
    MaxProfileValue = max(MeanProfile);
    MinProfileValue = min(MeanProfile);
    MinTKeyRegion = min(CompiledEmbryos.DubuisEmbryoTimes(CompiledEmbryos.AllDubuisValidProfilesControlTF& ...
        (CompiledEmbryos.FixCorrectedDeltaFC_um.mean >= 5) &...
        (CompiledEmbryos.FixCorrectedDeltaFC_um.mean <= 15)));
    MaxTKeyRegion = max(CompiledEmbryos.DubuisEmbryoTimes(CompiledEmbryos.AllDubuisValidProfilesTestTF& ...
        (CompiledEmbryos.FixCorrectedDeltaFC_um.mean >= 5) &...
        (CompiledEmbryos.FixCorrectedDeltaFC_um.mean <= 15)));
    ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(:,RefBin,ch_index)))*1.1;
    plot([MinTKeyRegion, MinTKeyRegion], [0, ymax], 'k--','LineWidth', 2.0);
     plot([MaxTKeyRegion, MaxTKeyRegion], [0, ymax], 'k--','LineWidth', 2.0);
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
%xlim([.1, .9])

ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 16)

ylim([0, ymax])
xlim([0 60])

DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [PlotLabel, ' (Control)'], 'FontSize', 18)



DeltaFCAx.FontSize = 16;

outpath = [ProfFigPath, filesep, SetLabel, '_',ChannelNames{ch_index}, '_ControlProfilesFig2D_TimeSmoothed.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end
%%
warning('off','stats:LinearModel:RankDefDesignMat')
exp_idx = 1:NumSets;
counts_1sigma = zeros(NumSets, size(AllCompiledEmbryos{1}.FixCorrectedSmoothedProfiles.ControlCounts.counts_within_1sigma, 2));
counts_2sigma = zeros(NumSets, size(AllCompiledEmbryos{1}.FixCorrectedSmoothedProfiles.ControlCounts.counts_within_2sigma, 2));
SmoothedDeltaFCs = AllCompiledEmbryos{1}.FixCorrectedSmoothedDeltaFCs;
ControlledProfiles = NaN([size(AllCompiledEmbryos{1}.FixCorrectedSmoothedAvgAPProfiles.Control), NumSets]);
for i = 1:NumSets
counts_1sigma(i,:) = AllCompiledEmbryos{exp_idx(i)}.FixCorrectedSmoothedProfiles.ControlCounts.counts_within_1sigma;
counts_2sigma(i,:) = AllCompiledEmbryos{exp_idx(i)}.FixCorrectedSmoothedProfiles.ControlCounts.counts_within_2sigma;
ControlledProfiles(:,:,:,i) = AllCompiledEmbryos{exp_idx(i)}.FixCorrectedSmoothedAvgAPProfiles.Control;
end
FitInfo = cell(NumSets, NumSets, NChannels);
ScalingFactors = NaN(NumSets, NumSets, NChannels);
ScalingSEs = NaN(NumSets, NumSets, NChannels);
GoodOverlapIndices = cell(NumSets, NumSets);

OkFitInfo = cell(NumSets, NumSets, NChannels);
OkScalingFactors = NaN(NumSets, NumSets, NChannels);
OkScalingSEs = NaN(NumSets, NumSets, NChannels);
OkOverlapIndices = cell(NumSets, NumSets);
for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        OkOverlapIndices{exp2, exp1} = (counts_1sigma(i,:)+counts_2sigma(i,:) >= 5) & (counts_1sigma(j,:)+counts_2sigma(j,:) >= 5);
        GoodOverlapIndices{exp2, exp1} = (counts_1sigma(i,:) >= 5) & (counts_1sigma(j,:) >= 5);
        if sum(GoodOverlapIndices{exp2, exp1}) >= 8
            for ch_index = 1:NChannels
                FitSet1 = ControlledProfiles(:,:,ch_index,i);
                FitSet1 = FitSet1( GoodOverlapIndices{exp2, exp1},:);
                FitSet1 = FitSet1(:);
                FitSet2 = ControlledProfiles(:,:,ch_index,j);
                FitSet2 = FitSet2( GoodOverlapIndices{exp2, exp1},:);
                FitSet2 = FitSet2(:);
                dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
                FitInfo{exp2, exp1, ch_index} = dlm;
                ScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate;
                ScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE;
            end
        else
            disp([num2str(exp1),', ', num2str(exp2)])
        end
        
        if sum(OkOverlapIndices{exp2, exp1}) >= 8
            for ch_index = 3:NChannels
                FitSet1 = ControlledProfiles(:,:,ch_index,i);
                FitSet1 = FitSet1( OkOverlapIndices{exp2, exp1},:);
                FitSet1 = FitSet1(:);
                FitSet2 = ControlledProfiles(:,:,ch_index,j);
                FitSet2 = FitSet2( OkOverlapIndices{exp2, exp1},:);
                FitSet2 = FitSet2(:);
                dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
                OkFitInfo{exp2, exp1, ch_index} = dlm;
                OkScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate;
                OkScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE;
            end
        end
    end
end

for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        for ch_index = 2:NChannels
            ScalingFactors(exp1, exp2, ch_index) = 1/ScalingFactors(exp2, exp1, ch_index);
            ScalingSEs(exp1, exp2, ch_index) = ScalingSEs(exp2, exp1, ch_index)/ScalingFactors(exp2, exp1, ch_index);
            OkScalingFactors(exp1, exp2, ch_index) = 1/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingSEs(exp1, exp2, ch_index) = OkScalingSEs(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
        end
        
    end
end

SetScalingFactors = NaN(NumSets, NChannels);
SetScalingSEs = NaN(NumSets, NChannels);
SetScalingFactors(RefSets(1), 2:5) = 1;
SetScalingSEs(RefSets(1), 2:5) = 1;
for ch_index = 2:NChannels 
    for i = AllSets(~ismember(AllSets, [RefSets(1)]))
        SetScalingFactors(i, ch_index) = ScalingFactors(i, RefSets(1), ch_index);
        SetScalingSEs(i, ch_index) = ScalingSEs(i, RefSets(1), ch_index);
        if isnan(ScalingFactors(i, RefSets(1), ch_index))
            SetScalingFactors(i, ch_index) = OkScalingFactors(i, RefSets(1), ch_index);
            SetScalingSEs(i, ch_index) = OkScalingSEs(i, RefSets(1), ch_index);
        end
    end
end
%%
% factor1 = [];
% factor2 = [];
% factor3 = [];
% channels = [];
% triple_multiplier = [];
% triple_multiplier_se = [];
% AllSets = 1:NumSets;
% for i = AllSets
%     for j = AllSets(~ismember(AllSets, [i]))
%         for k = AllSets(~ismember(AllSets, [i, j]))
%             for ch_index = 2:NChannels
%             factor1(end+1) = i;
%             factor2(end+1) = j;
%             factor3(end+1) = k;
%             channels(end+1) = ch_index;
%             triple_multiplier(end+1) = ScalingFactors(j, i, ch_index)*ScalingFactors(k, j, ch_index)*ScalingFactors(i, k, ch_index);
%             triple_multiplier_se(end+1) = sqrt(ScalingFactors(j, i, ch_index)^2*ScalingFactors(k, j, ch_index)^2*ScalingSEs(i, k, ch_index)^2 +...
%                 ScalingFactors(j, i, ch_index)^2*ScalingSEs(k, j, ch_index)^2*ScalingFactors(i, k, ch_index)^2+...
%                 ScalingSEs(j, i, ch_index)^2*ScalingFactors(k, j, ch_index)^2*ScalingFactors(i, k, ch_index)^2);
%             end
%         end
%     end
% end
% 
% % figure(2)
% % for i = AllSets
% %     scatter(i*ones(1, length(triple_multiplier(channels == 3 & (factor1 == i)))), triple_multiplier(channels == 3 & (factor1 == i)))
% %     hold on 
% % end
% % hold off
% % 
% % figure(3)
% % for i = AllSets
% %     scatter(i*ones(1, length(triple_multiplier(channels == 5 & (factor1 == i)))), triple_multiplier(channels == 5 & (factor1 == i)))
% %     hold on 
% % end
% % hold off
% 
% avg_triple_multiplier = NaN(NumSets,NChannels);
% bestsets = zeros(1, NChannels);
% for ch_index = 2:NChannels
%     for i = AllSets
%         
%         avg_triple_multiplier(i, ch_index) = mean(triple_multiplier(channels == ch_index & factor1 == i), 'omitnan');
%     end
%     [min_val, min_idx] = min(avg_triple_multiplier(:, ch_index));
%     bestsets(ch_index) = min_idx;
% end
% 
% SetScalingFactors = NaN(NumSets, NChannels);
% GoodFactors = zeros(NumSets, NumSets, NChannels, 'logical');
% % for ch_index = [3]%2:NChannels
% ch_index = 3;
% 
% BrokenCondition = false;
% IncludedSets = AllSets(any(~isnan(ScalingFactors(:,:,ch_index)), 1));
% SetsLeft = IncludedSets;
% SetsFound = [];
% 
% tm_notnan = ~isnan(triple_multiplier);
% factor1_sub = factor1(channels == ch_index & tm_notnan);
% factor2_sub = factor2(channels == ch_index & tm_notnan);
% factor3_sub = factor3(channels == ch_index & tm_notnan);
% 
% tm_sub = triple_multiplier(channels == ch_index & tm_notnan) ;
% tm_se_sub = triple_multiplier_se(channels == ch_index & tm_notnan) ;
% abs_tm_sub_diff = abs(tm_sub-1);
% [~, min_idx] = min(abs_tm_sub_diff);
% f1 = factor1_sub(min_idx);
% f2 = factor2_sub(min_idx);
% f3 = factor3_sub(min_idx);
% GoodFactors([f1 f1 f2 f2 f3 f3], [f2 f3 f3 f1 f1 f2], ch_index) = true;
% SetsFound = [SetsFound f1 f2 f3];
% SetsLeft = SetsLeft(~ismember(SetsLeft, SetsFound));
% RefSet = f1;
% SetScalingFactors(RefSet, ch_index) = 1;
% 
% counter = 0;
% while ~isempty(SetsLeft) & BrokenCondition == false
%     rm_cond = ((factor1_sub == f1)  & (factor2_sub == f2) & (factor3_sub == f3) ) | ...
%     ((factor1_sub == f1)  & (factor2_sub == f3) & (factor3_sub == f2) ) | ...
%     ((factor1_sub == f2)  & (factor2_sub == f3) & (factor3_sub == f1) ) | ...
%     ((factor1_sub == f2)  & (factor2_sub == f1) & (factor3_sub == f3) ) | ...
%     ((factor1_sub == f3)  & (factor2_sub == f1) & (factor3_sub == f2) ) | ...
%     ((factor1_sub == f3)  & (factor2_sub == f2) & (factor3_sub == f1) ) ;
% factor1_sub = factor1_sub(~rm_cond);
% factor2_sub = factor2_sub(~rm_cond);
% factor3_sub = factor3_sub(~rm_cond);
% 
% tm_sub = tm_sub(~rm_cond);
% tm_se_sub = tm_se_sub(~rm_cond);
% 
% is_good_tm = ismember(factor1_sub, SetsFound) | ismember(factor2_sub, SetsFound) | ismember(factor3_sub, SetsFound);
% if sum(is_good_tm) > 0
% factor1_temp = factor1_sub(is_good_tm);
% factor2_temp = factor2_sub(is_good_tm);
% factor3_temp = factor3_sub(is_good_tm);
% 
% tm_temp = tm_sub(is_good_tm);
% tm_se_temp = tm_se_sub(is_good_tm);
% 
% 
% abs_tm_temp_diff = abs(tm_temp-1);
% [min_val, min_idx] = min(abs_tm_temp_diff);
% f1 = factor1_temp(min_idx);
% f2 = factor2_temp(min_idx);
% f3 = factor3_temp(min_idx);
% GoodFactors(f1, f2, ch_index) = true;
% GoodFactors(f1, f3, ch_index) = true;
% GoodFactors(f2, f3, ch_index) = true;
% GoodFactors(f2, f1, ch_index) = true;
% GoodFactors(f3, f1, ch_index) = true;
% GoodFactors(f3, f2, ch_index) = true;
% NewFoundSets = [f1 f2 f3];
% SetsFound = [SetsFound NewFoundSets(~ismember(NewFoundSets, SetsFound))];
% SetsLeft = SetsLeft(~ismember(SetsLeft, SetsFound));
% else
%     BrokenCondition = true;
% end
% counter = counter + 1;
% end




%     
% end




