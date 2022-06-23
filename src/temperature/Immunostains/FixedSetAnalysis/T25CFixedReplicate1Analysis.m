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
 


NChannels = 5;

exp_index = find(AllSetInfo.Temperatures == 20 & AllSetInfo.Replicates == 1);


APbins = 0:0.025:1;
NumAPbins = length(APbins);

NarrowAPbins = 0:0.0125:1;
NumNarrowAPbins = length(NarrowAPbins);


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

%%

CompiledEmbryos = RescaleSlideFluos(CompiledEmbryos, exp_index);

TStarts = 0:10:50;
TEnds = 10:10:60;
NTimeBins = length(TStarts);
CompiledEmbryos.Dubuis = {};
CompiledEmbryos.Dubuis.TimeAveragedDorsalTestProfiles = {};
CompiledEmbryos.Dubuis.StdTimeAveragedDorsalTestProfiles = {};
CompiledEmbryos.Dubuis.CountTimeAveragedDorsalTestProfiles = {};
CompiledEmbryos.Dubuis.TimeAveragedDorsalControlProfiles = {};
CompiledEmbryos.Dubuis.StdTimeAveragedDorsalControlProfiles = {};
CompiledEmbryos.Dubuis.CountTimeAveragedDorsalControlProfiles = {};
CompiledEmbryos.Dubuis.TimeAveragedDorsalTestProfiles.NC14 = NaN(NTimeBins, NumAPbins, NChannels);
CompiledEmbryos.Dubuis.TimeAveragedDorsalControlProfiles.NC14 = NaN(NTimeBins, NumAPbins, NChannels);
CompiledEmbryos.Dubuis.StdTimeAveragedDorsalTestProfiles.NC14 = NaN(NTimeBins, NumAPbins, NChannels);
CompiledEmbryos.Dubuis.StdTimeAveragedDorsalControlProfiles.NC14 = NaN(NTimeBins, NumAPbins, NChannels);
CompiledEmbryos.Dubuis.CountTimeAveragedDorsalTestProfiles.NC14 = NaN(NTimeBins, NumAPbins, NChannels);
CompiledEmbryos.Dubuis.CountTimeAveragedDorsalControlProfiles.NC14 = NaN(NTimeBins, NumAPbins, NChannels);
CompiledEmbryos.Dubuis.TimeAveragedNarrowDorsalTestProfiles = {};
CompiledEmbryos.Dubuis.StdTimeAveragedNarrowDorsalTestProfiles = {};
CompiledEmbryos.Dubuis.CountTimeAveragedNarrowDorsalTestProfiles = {};
CompiledEmbryos.Dubuis.TimeAveragedNarrowDorsalControlProfiles = {};
CompiledEmbryos.Dubuis.StdTimeAveragedNarrowDorsalControlProfiles = {};
CompiledEmbryos.Dubuis.CountTimeAveragedNarrowDorsalControlProfiles = {};
CompiledEmbryos.Dubuis.TimeAveragedNarrowDorsalTestProfiles.NC14 = NaN(NTimeBins, NumNarrowAPbins, NChannels);
CompiledEmbryos.Dubuis.TimeAveragedNarrowDorsalControlProfiles.NC14 = NaN(NTimeBins, NumNarrowAPbins, NChannels);
CompiledEmbryos.Dubuis.StdTimeAveragedNarrowDorsalTestProfiles.NC14 = NaN(NTimeBins, NumNarrowAPbins, NChannels);
CompiledEmbryos.Dubuis.StdTimeAveragedNarrowDorsalControlProfiles.NC14 = NaN(NTimeBins, NumNarrowAPbins, NChannels);

for i = 1:NTimeBins
    TFtimeBin = CompiledEmbryos.DubuisEmbryoTimes >= TStarts(i) & CompiledEmbryos.DubuisEmbryoTimes < TEnds(i)  & CompiledEmbryos.IsNC14;
    TFtestBin = TFtimeBin & CompiledEmbryos.TestSetEmbryos;
    TFcontrolBin = TFtimeBin & CompiledEmbryos.ControlSetEmbryos;
    for ch_index = 2:NChannels
        if sum(TFtestBin) > 1
            CompiledEmbryos.Dubuis.TimeAveragedDorsalTestProfiles.NC14(i,:,ch_index) = mean(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(TFtestBin,:,ch_index), 'omitnan');
            CompiledEmbryos.Dubuis.StdTimeAveragedDorsalTestProfiles.NC14(i,:,ch_index) = std(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(TFtestBin,:,ch_index), 'omitnan');
            CompiledEmbryos.Dubuis.TimeAveragedNarrowDorsalTestProfiles.NC14(i,:,ch_index) = mean(CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles(TFtestBin,:,ch_index), 'omitnan');
            CompiledEmbryos.Dubuis.StdTimeAveragedNarrowDorsalTestProfiles.NC14(i,:,ch_index) = std(CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles(TFtestBin,:,ch_index), 'omitnan');
        elseif sum(TFtestBin) == 1
           CompiledEmbryos.Dubuis.TimeAveragedDorsalTestProfiles.NC14(i,:,ch_index) = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(TFtestBin,:,ch_index);
            CompiledEmbryos.Dubuis.StdTimeAveragedDorsalTestProfiles.NC14(i,:,ch_index) = zeros(size(CompiledEmbryos.Dubuis.StdTimeAveragedDorsalTestProfiles.NC14(i,:,ch_index) )); 
            CompiledEmbryos.Dubuis.TimeAveragedNarrowDorsalTestProfiles.NC14(i,:,ch_index) = CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles(TFtestBin,:,ch_index);
            CompiledEmbryos.Dubuis.StdTimeAveragedNarrowDorsalTestProfiles.NC14(i,:,ch_index) = zeros(size(CompiledEmbryos.Dubuis.StdTimeAveragedNarrowDorsalTestProfiles.NC14(i,:,ch_index) )); 
        end
        CompiledEmbryos.Dubuis.CountTimeAveragedDorsalTestProfiles.NC14(i,:,ch_index)  = sum(TFtestBin);
        
        if sum(TFcontrolBin) > 1
            CompiledEmbryos.Dubuis.TimeAveragedDorsalControlProfiles.NC14(i,:,ch_index) = mean(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(TFcontrolBin,:,ch_index), 'omitnan');
            CompiledEmbryos.Dubuis.StdTimeAveragedDorsalControlProfiles.NC14(i,:,ch_index) = std(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(TFcontrolBin,:,ch_index), 'omitnan');
            CompiledEmbryos.Dubuis.TimeAveragedNarrowDorsalControlProfiles.NC14(i,:,ch_index) = mean(CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles(TFcontrolBin,:,ch_index), 'omitnan');
            CompiledEmbryos.Dubuis.StdTimeAveragedNarrowDorsalControlProfiles.NC14(i,:,ch_index) = std(CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles(TFcontrolBin,:,ch_index), 'omitnan');
        elseif sum(TFcontrolBin) == 1
           CompiledEmbryos.Dubuis.TimeAveragedDorsalControlProfiles.NC14(i,:,ch_index) = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(TFcontrolBin,:,ch_index);
            CompiledEmbryos.Dubuis.StdTimeAveragedDorsalControlProfiles.NC14(i,:,ch_index) = zeros(size(CompiledEmbryos.Dubuis.StdTimeAveragedDorsalTestProfiles.NC14(i,:,ch_index) )); 
            CompiledEmbryos.Dubuis.TimeAveragedNarrowDorsalControlProfiles.NC14(i,:,ch_index) = CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles(TFcontrolBin,:,ch_index);
            CompiledEmbryos.Dubuis.StdTimeAveragedNarrowDorsalControlProfiles.NC14(i,:,ch_index) = zeros(size(CompiledEmbryos.Dubuis.StdTimeAveragedNarrowDorsalTestProfiles.NC14(i,:,ch_index) )); 
        end
        CompiledEmbryos.Dubuis.CountTimeAveragedDorsalControlProfiles.NC14(i,:,ch_index)  = sum(TFcontrolBin);
    end
end

CompiledEmbryos = AddSmoothedProfiles(CompiledEmbryos);
CompiledEmbryos = AddBinnedProfiles(CompiledEmbryos);
%%
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
save(CEoutpath, 'CompiledEmbryos');


end
AllCompiledEmbryos = CalculateSetRescalingFactors;
%%
for exp_index = 1:length(AllSetInfo.Temperatures)
    %%
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
CEpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
load(CEpath, 'CompiledEmbryos');
%%
NEmbryos = length(CompiledEmbryos.Approved);
AllEmbryos = 1:NEmbryos;
NC14Indices = find(CompiledEmbryos.IsNC14);
NC13Indices = find(CompiledEmbryos.IsNC13);
NC13NC14Indices = find(CompiledEmbryos.IsNC13orNC14);
NumEmbryosNC14 = CompiledEmbryos.NumEmbryosNC14;
%%
APbins = 0:0.025:1;
NumAPbins = length(APbins);

NarrowAPbins = 0:0.0125:1;
NumNarrowAPbins = length(NarrowAPbins);

ChannelNames = {'', 'Hoechst', 'Bicoid', 'Knirps', 'Hunchback'};
MinDeltaFC = min(CompiledEmbryos.FixMeanCorrectedDeltaFC_um.mean(NC14Indices));
MaxDeltaFC = max(CompiledEmbryos.FixMeanCorrectedDeltaFC_um.mean(NC14Indices));

DeltaFCRange = 0:0.5:45;
colors = hsv(length(DeltaFCRange)); % Colormap "jet" is another option
colors2 = brewermap(11,'Spectral');
colors2 = colors2([1 3 8 9 11], :);
MarkerStyles = {'o', 'd', 's', '>', '^','p', 'h', '*', 'x'};
FractionalDeltaFCRange = (DeltaFCRange)/45;

for ChannelIndex = [2,3,5,4]
close all
DeltaFCFixCorrectedFig = figure(1);
set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
DeltaFCAx = axes(DeltaFCFixCorrectedFig);
eb = cell(1,NumEmbryosNC14);
prof = cell(1, NumEmbryosNC14);

map = colormap(colors);
h = colorbar;
% %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
hold off
colorTitleHandle = get(h,'Title');
titleString = '\delta_{FC} (\mum)';
set(colorTitleHandle ,'String',titleString);
h.Ticks =  FractionalDeltaFCRange(1:10:91); %Create 8 ticks from zero to 1
h.TickLabels = {'0.0','5.0','10.0','15.0','20.0','25.0', '30.0', '35.0', '40.0', '45.0'} ;
hold on

for idx =1:NumEmbryosNC14
    if ~isnan(CompiledEmbryos.FixCorrectedDeltaFC_um.mean(NC14Indices(idx)))
    dfc_idx = find(abs(DeltaFCRange - CompiledEmbryos.FixCorrectedDeltaFC_um.mean(NC14Indices(idx))) == min(abs(DeltaFCRange - CompiledEmbryos.FixCorrectedDeltaFC_um.mean(NC14Indices(idx)))), 1);
    
    prof{idx} = plot(APbins,  CompiledEmbryos.DorsalAPProfiles(NC14Indices(idx),:,ChannelIndex).',...
        'LineStyle', '-', 'Color', colors(dfc_idx,:));
    
    hold on
    end
end


grid on

hold off

xlabel('Fraction Embryo Length', 'FontSize', 16)
xlim([0.1, 0.9])

ylabel([ChannelNames{ChannelIndex}, ' (AU)'], 'FontSize', 16)


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;
ylim([0, max(ceil(max(max(CompiledEmbryos.DorsalAPProfiles(NC14Indices,:,ChannelIndex).'))*1.1),1)])
title(DeltaFCAx, [PlotLabel, ' NC14'], 'FontSize', 18)
%


h.FontSize = 16;

outpath = [ProfFigPath, filesep, SetLabel, '_', ChannelNames{ChannelIndex}, 'NC14DeltaFCColor_UnbinnedProfiles.png'];
saveas(DeltaFCFixCorrectedFig,outpath);

DeltaFCFixCorrectedFig2 = figure(2);
set(DeltaFCFixCorrectedFig2,'units', 'normalized', 'position',[0.51, 0.05, .6, .6]);
set(gcf,'color','w');
DeltaFCAx2 = axes(DeltaFCFixCorrectedFig2);
prof2 = cell(1, NumEmbryosNC14);



for idx =1:NumEmbryosNC14
    

    prof2{idx} = plot(APbins,  CompiledEmbryos.DorsalAvgAPProfiles(NC14Indices(idx),:,ChannelIndex).',...
        'LineStyle', '-', 'Color', colors2(CompiledEmbryos.SlideIDs(NC14Indices(idx)),:));
    
    hold on
end


grid on

hold off

xlabel('Fraction Embryo Length', 'FontSize', 16)
xlim([0.1, 0.9])

ylabel([ChannelNames{ChannelIndex}, ' (AU)'], 'FontSize', 16)


DeltaFCAx2.YAxis.FontSize = 16;
DeltaFCAx2.XAxis.FontSize = 16;
ylim([0, max(ceil(max(max(CompiledEmbryos.DorsalAvgAPProfiles(NC14Indices,:,ChannelIndex).'))*1.1),1)])
title(DeltaFCAx2, [PlotLabel, ' NC14'], 'FontSize', 18)
%



outpath = [ProfFigPath, filesep, SetLabel, '_', ChannelNames{ChannelIndex}, 'NC14SlideIDColor_UnbinnedProfiles.png'];
saveas(DeltaFCFixCorrectedFig2,outpath);
end



%%
close all
NumEmbryosNC13 = CompiledEmbryos.NumEmbryosNC13;;
if NumEmbryosNC13 > 0
ChannelNames = {'', 'Hoechst', 'Bicoid', 'Knirps', 'Hunchback'};


for ChannelIndex = [2,3,5,4]
close all
DeltaFCFixCorrectedFig2 = figure(2);
set(DeltaFCFixCorrectedFig2,'units', 'normalized', 'position',[0.51, 0.05, .6, .6]);
set(gcf,'color','w');
DeltaFCAx2 = axes(DeltaFCFixCorrectedFig2);
prof2 = cell(1, NumEmbryosNC13);



for idx =1:NumEmbryosNC13
    

    prof2{idx} = plot(APbins,  CompiledEmbryos.DorsalAvgAPProfiles(NC13Indices(idx),:,ChannelIndex).',...
        'LineStyle', '-', 'Color', colors2(CompiledEmbryos.SlideIDs(NC13Indices(idx)),:));
    
    hold on
end


grid on

hold off

xlabel('Fraction Embryo Length', 'FontSize', 16)
xlim([0.1, 0.9])

ylabel([ChannelNames{ChannelIndex}, ' (AU)'], 'FontSize', 16)


DeltaFCAx2.YAxis.FontSize = 16;
DeltaFCAx2.XAxis.FontSize = 16;
ylim([0, max(ceil(max(max(CompiledEmbryos.DorsalAvgAPProfiles(NC13Indices,:,ChannelIndex).'))*1.1),1)])
title(DeltaFCAx2, [PlotLabel, ' NC13'], 'FontSize', 18)
%



outpath = [ProfFigPath, filesep, SetLabel, '_', ChannelNames{ChannelIndex}, 'NC13SlideIDColor_UnbinnedProfiles.png'];
saveas(DeltaFCFixCorrectedFig2,outpath);
end


close all
end

%% Add Partitioning
KnirpsIndex = 4;
%
for ChannelIndex = [2,3,5,4]
close all
DeltaFCFixCorrectedFig = figure(1);
set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.1, 0.05, .6, .6]);
set(gcf,'color','w');
DeltaFCAx = axes(DeltaFCFixCorrectedFig);
prof2 = cell(1, NumEmbryosNC14);

MaxYs = zeros(1, NumEmbryosNC14);

for idx =1:NumEmbryosNC14
    
    scale_factor = CompiledEmbryos.SlideRescalingFactors(CompiledEmbryos.SlideIDs(NC14Indices(idx)), ChannelIndex);
    prof2{idx} = plot(APbins,  scale_factor*CompiledEmbryos.DorsalAvgAPProfiles(NC14Indices(idx),:,ChannelIndex).',...
        'LineStyle', '-', 'Color', colors2(CompiledEmbryos.SlideIDs(NC14Indices(idx)),:));
    MaxYs(idx) = max(scale_factor*CompiledEmbryos.DorsalAvgAPProfiles(NC14Indices(idx),:,ChannelIndex).');
    hold on
end


grid on

hold off

xlabel('Fraction Embryo Length', 'FontSize', 16)
xlim([0.1, 0.9])

ylabel([ChannelNames{ChannelIndex}, ' (AU)'], 'FontSize', 16)


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;
ylim([0, max((max(MaxYs)*1.1),1)])
title(DeltaFCAx, [PlotLabel, ' NC14'], 'FontSize', 18)
%



outpath = [ProfFigPath, filesep, SetLabel, '_', ChannelNames{ChannelIndex}, 'SlideRescaledNC14SlideIDColor_UnbinnedProfiles.png'];
saveas(DeltaFCFixCorrectedFig,outpath);


end


%%


for ChannelIndex = 2:5
close all
DeltaFCFixCorrectedFig = figure(1);
set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
DeltaFCAx = axes(DeltaFCFixCorrectedFig);
eb = cell(1,NumEmbryosNC14);
prof = cell(1, NumEmbryosNC14);

map = colormap(colors);
h = colorbar;
% %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
hold off
colorTitleHandle = get(h,'Title');
titleString = '\delta_{FC} (\mum)';
set(colorTitleHandle ,'String',titleString);
h.Ticks =  FractionalDeltaFCRange(1:10:91); %Create 8 ticks from zero to 1
h.TickLabels = {'0.0','5.0','10.0','15.0','20.0','25.0', '30.0', '35.0', '40.0', '45.0'} ;
hold on

for idx =1:NumEmbryosNC14
    if ~isnan(CompiledEmbryos.FixCorrectedDeltaFC_um.mean(NC14Indices(idx))) & CompiledEmbryos.TestSetEmbryos(NC14Indices(idx))
    dfc_idx = find(abs(DeltaFCRange - CompiledEmbryos.FixCorrectedDeltaFC_um.mean(NC14Indices(idx))) == min(abs(DeltaFCRange - CompiledEmbryos.FixCorrectedDeltaFC_um.mean(NC14Indices(idx)))), 1);
    
    prof{idx} = plot(APbins,  CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(NC14Indices(idx),:,ChannelIndex).',...
        'LineStyle', '-', 'Color', colors(dfc_idx,:));
    
    hold on
    end
end


grid on

hold off

xlabel('Fraction Embryo Length', 'FontSize', 16)
xlim([0.1, 0.9])

ylabel([ChannelNames{ChannelIndex}, ' (AU)'], 'FontSize', 16)


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;
ylim([0, max(ceil(max(max(CompiledEmbryos.SlideRescaledDorsalAPProfiles(NC14Indices,:,ChannelIndex).'))*1.1),1)])
title(DeltaFCAx, [PlotLabel, ' NC14 (Test Set w/ Slide Rescaling)'], 'FontSize', 18)
%


h.FontSize = 16;

outpath = [ProfFigPath, filesep, SetLabel, '_', ChannelNames{ChannelIndex}, 'TestSetSlideRescaledNC14DeltaFCColor_UnbinnedProfiles.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end



%% Plot vs. Dubuis Times
MinDubuisTime = min(CompiledEmbryos.DubuisEmbryoTimes(NC14Indices).');
MaxDubuisTime = max(CompiledEmbryos.DubuisEmbryoTimes(NC14Indices).');

DubuisTimeRange = 0:1.5:60;
colors = hsv(length(DubuisTimeRange)); % Colormap "jet" is another option
colors2 = brewermap(11,'Spectral');
colors2 = colors2([1 3 8 9 11], :);
MarkerStyles = {'o', 'd', 's', '>', '^','p', 'h', '*', 'x'};
FractionalDubuisTimeRange = (DubuisTimeRange)/60;
for ChannelIndex = 2:5
close all
DeltaFCFixCorrectedFig = figure(1);
set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
DeltaFCAx = axes(DeltaFCFixCorrectedFig);
eb = cell(1,NumEmbryosNC14);
prof = cell(1, NumEmbryosNC14);

map = colormap(colors);
h = colorbar;
% %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
hold off
colorTitleHandle = get(h,'Title');
titleString = 'Time into NC14 (m)';
set(colorTitleHandle ,'String',titleString);
h.Ticks =  FractionalDubuisTimeRange(1:4:41); %Create 8 ticks from zero to 1
h.TickLabels = {'0.0','6.0','12.0','18.0','24.0','30.0', '36.0', '42.0', '48.0', '54.0', '60.0'} ;
hold on

for idx =1:NumEmbryosNC14
    if ~isnan(CompiledEmbryos.DubuisEmbryoTimes(NC14Indices(idx))) & CompiledEmbryos.TestSetEmbryos(NC14Indices(idx))
    dfc_idx = find(abs(DubuisTimeRange - CompiledEmbryos.DubuisEmbryoTimes(NC14Indices(idx))) == min(abs(DubuisTimeRange - CompiledEmbryos.DubuisEmbryoTimes(NC14Indices(idx)))), 1);
    
    prof{idx} = plot(APbins,  CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(NC14Indices(idx),:,ChannelIndex).',...
        'LineStyle', '-', 'Color', colors(dfc_idx,:));
    
    hold on
    end
end


grid on

hold off

xlabel('Fraction Embryo Length', 'FontSize', 16)
xlim([0.1, 0.9])

ylabel([ChannelNames{ChannelIndex}, ' (AU)'], 'FontSize', 16)


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;
ylim([0, max(ceil(max(max(CompiledEmbryos.SlideRescaledDorsalAPProfiles(NC14Indices,:,ChannelIndex).'))*1.1),1)])
title(DeltaFCAx, [PlotLabel, ' NC14 (Test Set)'], 'FontSize', 18)
%


h.FontSize = 16;

outpath = [ProfFigPath, filesep, SetLabel, '_', ChannelNames{ChannelIndex}, 'TestSetSlideRescaledNC14DubuisColor_UnbinnedProfiles.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end

%% Bin into 10-minute profiles - Note that this is very crude!
TStarts = 0:10:50;
TEnds = 10:10:60;
NTimeBins = length(TStarts);


%%
for ChannelIndex = 2:5
close all
DeltaFCFixCorrectedFig = figure(1);
set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .8, .6]);
set(gcf,'color','w');
DeltaFCAx = axes(DeltaFCFixCorrectedFig);
prof = cell(1,length(DubuisTimeRange));
prof2 = cell(1, length(DubuisTimeRange));
colors = hsv(NTimeBins); % Colormap "jet" is another option
map = colormap(colors);
h = colorbar;
% %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
hold off
colorTitleHandle = get(h,'Title');
titleString = 'Time into NC14 (m)';
set(colorTitleHandle ,'String',titleString);
h.Ticks =  (0.5:(NTimeBins+0.5))/NTimeBins; %Create 8 ticks from zero to 1
h.TickLabels = {'0-10', '10-20', '20-30', '30-40', '40-50', '50-60'} ;
hold on

for idx =1:length(TStarts)
    
    prof{idx} = plot(APbins,  CompiledEmbryos.Dubuis.TimeAveragedDorsalTestProfiles.NC14(idx,:,ChannelIndex).',...
        'LineStyle', '-', 'Color', colors(idx,:));
    
    hold on

end

for idx =1:length(TStarts)
    
    prof2{idx} = plot(APbins,  CompiledEmbryos.Dubuis.TimeAveragedDorsalControlProfiles.NC14(idx,:,ChannelIndex).',...
        'LineStyle', '--', 'Color', colors(idx,:));
    
    hold on

end


grid on

hold off

xlabel('Fraction Embryo Length', 'FontSize', 16)
xlim([0.1, 0.9])

ylabel([ChannelNames{ChannelIndex}, ' (AU)'], 'FontSize', 16)


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;
ylim([0, max([ceil(max(max(CompiledEmbryos.Dubuis.TimeAveragedDorsalTestProfiles.NC14(:,:,ChannelIndex).'))*1.1),...
    ceil(max(max(CompiledEmbryos.Dubuis.TimeAveragedDorsalControlProfiles.NC14(:,:,ChannelIndex).'))*1.1), 1])])
title(DeltaFCAx, [PlotLabel, ' NC14 (Test Set Solid, Control Set Dashed)'], 'FontSize', 18)
%


h.FontSize = 16;

outpath = [ProfFigPath, filesep, SetLabel, '_', ChannelNames{ChannelIndex}, '_TestControlSetSlideRescaledNC14DubuisColor_BinnedProfiles.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
% end
%% Use Bins from AddBinnedProfiles 
MinDeltaFC = min(CompiledEmbryos.BinnedProfiles.DeltaFC.x);
MaxDeltaFC = max(CompiledEmbryos.BinnedProfiles.DeltaFC.x);
DeltaFCBinRange = CompiledEmbryos.BinnedProfiles.DeltaFC.x;
colors = hsv(length(DeltaFCBinRange)); % Colormap "jet" is another option
FractionalDeltaFCBinRange = (DeltaFCBinRange)/(MaxDeltaFC+(CompiledEmbryos.BinnedProfiles.DeltaFC.x(2)-CompiledEmbryos.BinnedProfiles.DeltaFC.x(1)));
for ChannelIndex = 2:5
close all
DeltaFCFixCorrectedFig = figure(1);
set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .8, .6]);
set(gcf,'color','w');
DeltaFCAx = axes(DeltaFCFixCorrectedFig);
prof = cell(1,length(FractionalDeltaFCBinRange));
prof2 = cell(1, length(FractionalDeltaFCBinRange));
map = colormap(colors);
h = colorbar;
% %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
hold off
colorTitleHandle = get(h,'Title');
titleString = '\delta_{FC} (\mum)';
set(colorTitleHandle ,'String',titleString);
h.Ticks =  FractionalDeltaFCBinRange(1:2:end)+0.5*(FractionalDeltaFCBinRange(2)-FractionalDeltaFCBinRange(1)); %Create 8 ticks from zero to 1
h.TickLabels = {'0.0', '5.0', '10.0', '15.0', '20.0', '25.0', '30.0', '35.0', '40.0', '45.0'} ;
hold on

for idx =1:length(DeltaFCBinRange)
    if ~all(isnan(CompiledEmbryos.BinnedProfiles.DeltaFC.Test.mean(idx,:,ChannelIndex)))

    prof{idx} = plot(APbins,  CompiledEmbryos.BinnedProfiles.DeltaFC.Test.mean(idx,:,ChannelIndex).',...
        'LineStyle', '-', 'Color', colors(idx,:), 'LineWidth', 1.5);
    
    hold on
    end
end

for idx =1:length(DeltaFCBinRange)
    if ~all(isnan(CompiledEmbryos.BinnedProfiles.DeltaFC.Control.mean(idx,:,ChannelIndex)))

    prof{idx} = plot(APbins,  CompiledEmbryos.BinnedProfiles.DeltaFC.Control.mean(idx,:,ChannelIndex).',...
        'LineStyle', '--', 'Color', colors(idx,:), 'LineWidth', 1.5);
    
    hold on
    end
end


grid on

hold off

xlabel('Fraction Embryo Length', 'FontSize', 16)
xlim([0.1, 0.9])

ylabel([ChannelNames{ChannelIndex}, ' (AU)'], 'FontSize', 16)


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;
ylim([0, max([ceil(max(max(CompiledEmbryos.BinnedProfiles.DeltaFC.Test.mean(:,:,ChannelIndex).'))/5)*5,...
    ceil(max(max(CompiledEmbryos.BinnedProfiles.DeltaFC.Control.mean(:,:,ChannelIndex).'))/5)*5, 1])])
title(DeltaFCAx, [PlotLabel, ' NC14 (Test Set Solid, Control Set Dashed)'], 'FontSize', 18)
%


h.FontSize = 16;

% cbPos = h.Position;
% cbPos(1) = 0.94;
% cbPos(2) = 0.08;
% set(h,'Position',cbPos)

outpath = [ProfFigPath, filesep, SetLabel, '_', ChannelNames{ChannelIndex}, '_TestControlSetSlideRescaledNC14DeltaColor_FinerBinProfiles.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end


%% Use Bins from AddBinnedProfiles  - Dubuis Times
MinDubuisTimeBin = min(CompiledEmbryos.BinnedProfiles.DubuisTime.x);
MaxDubuisTimeBin = max(CompiledEmbryos.BinnedProfiles.DubuisTime.x);
DubuisTimeBinRange = CompiledEmbryos.BinnedProfiles.DubuisTime.x;
colors = hsv(length(DubuisTimeBinRange)); % Colormap "jet" is another option
FractionalDubuisTimeBinRange = (DubuisTimeBinRange)/(MaxDubuisTimeBin+(CompiledEmbryos.BinnedProfiles.DubuisTime.x(2)-CompiledEmbryos.BinnedProfiles.DubuisTime.x(1)));
for ChannelIndex = 2:5
close all
DeltaFCFixCorrectedFig = figure(1);
set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .8, .6]);
set(gcf,'color','w');
DeltaFCAx = axes(DeltaFCFixCorrectedFig);
prof = cell(1,length(FractionalDeltaFCBinRange));
prof2 = cell(1, length(FractionalDeltaFCBinRange));
map = colormap(colors);
h = colorbar;
% %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
hold off
colorTitleHandle = get(h,'Title');
titleString = {'cycle 14';' time (m)'};
set(colorTitleHandle ,'String',titleString);
h.Ticks =  FractionalDubuisTimeBinRange(1:end)+0.5*(FractionalDubuisTimeBinRange(2)-FractionalDubuisTimeBinRange(1)); %Create 8 ticks from zero to 1
h.TickLabels = {'0.0', '5.0', '10.0', '15.0', '20.0', '25.0', '30.0', '35.0', '40.0', '45.0', '50.0', '55.0', '60.0', '65.0'} ;
hold on

for idx =1:length(DubuisTimeBinRange)
    if ~all(isnan(CompiledEmbryos.BinnedProfiles.DubuisTime.Test.mean(idx,:,ChannelIndex)))

    prof{idx} = plot(APbins,  CompiledEmbryos.BinnedProfiles.DubuisTime.Test.mean(idx,:,ChannelIndex).',...
        'LineStyle', '-', 'Color', colors(idx,:), 'LineWidth', 1.5);
    
    hold on
    end
end

for idx =1:length(DubuisTimeBinRange)
    if ~all(isnan(CompiledEmbryos.BinnedProfiles.DubuisTime.Control.mean(idx,:,ChannelIndex)))

    prof{idx} = plot(APbins,  CompiledEmbryos.BinnedProfiles.DubuisTime.Control.mean(idx,:,ChannelIndex).',...
        'LineStyle', '--', 'Color', colors(idx,:), 'LineWidth', 1.5);
    
    hold on
    end
end


grid on

hold off

xlabel('Fraction Embryo Length', 'FontSize', 16)
xlim([0.1, 0.9])

ylabel([ChannelNames{ChannelIndex}, ' (AU)'], 'FontSize', 16)


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;
ylim([0, max([ceil(max(max(CompiledEmbryos.BinnedProfiles.DubuisTime.Test.mean(:,:,ChannelIndex).'))/5)*5,...
    ceil(max(max(CompiledEmbryos.BinnedProfiles.DubuisTime.Control.mean(:,:,ChannelIndex).'))/5)*5, 1])])
title(DeltaFCAx, [PlotLabel, ' NC14 (Test Set Solid, Control Set Dashed)'], 'FontSize', 18)
%


h.FontSize = 16;
cbPos = h.Position;
cbPos(1) = 0.94;
cbPos(2) = 0.08;
set(h,'Position',cbPos)

outpath = [ProfFigPath, filesep, SetLabel, '_', ChannelNames{ChannelIndex}, '_TestControlSetSlideRescaledNC14DubuisTimeColor_FinerBinProfiles.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end

%%
close all
DeltaFCFixCorrectedFig = figure(1);
set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
DeltaFCAx = axes(DeltaFCFixCorrectedFig);
h = histogram(CompiledEmbryos.FixCorrectedDeltaFC_um.mean(CompiledEmbryos.TestSetEmbryos), 0:5:50);
hold on 
h2 = histogram(CompiledEmbryos.FixCorrectedDeltaFC_um.mean(CompiledEmbryos.ControlSetEmbryos), 0:5:50);
grid on 
hold off

xlabel('\delta_{FC}', 'FontSize', 16)
xlim([0, 45])

ylabel('Counts', 'FontSize', 16)
ymax = ceil(max(h.Values)/5)*5+5;
ylim([0, ymax])

DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, PlotLabel, 'FontSize', 18)
%
legend({'Test Set', 'Control Set'});


DeltaFCAx.FontSize = 16;

outpath = [ProfFigPath, filesep, SetLabel, '_FixCorrectedDeltaFCHistogram.png'];
saveas(DeltaFCFixCorrectedFig,outpath);

%%
close all
DeltaFCFixCorrectedFig = figure(1);
set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
DeltaFCAx = axes(DeltaFCFixCorrectedFig);
h = histogram(CompiledEmbryos.DubuisEmbryoTimes(CompiledEmbryos.TestSetEmbryos), 0:5:60);
hold on 
h2 = histogram(CompiledEmbryos.DubuisEmbryoTimes(CompiledEmbryos.ControlSetEmbryos), 0:5:60);
grid on 
hold off

xlabel('Time into NC14 from Dubuis Profile (m)', 'FontSize', 16)
xlim([0, 60])

ylabel('Counts', 'FontSize', 16)
ymax = ceil(max(h.Values)/5)*5+5;
ylim([0, ymax])

legend({'Test Set', 'Control Set'});

DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, PlotLabel, 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [ProfFigPath, filesep, SetLabel, '_DubuisTimesHistogram.png'];
saveas(DeltaFCFixCorrectedFig,outpath);

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
    MeanProfile2 = CompiledEmbryos.WindowedProfiles.DeltaFC.Test.mean(:,RefBin,ch_index).';
    x2 = CompiledEmbryos.WindowedProfiles.DeltaFC.x;
    plot(x2(~isnan(MeanProfile2)), MeanProfile2(~isnan(MeanProfile2)), 'r', 'LineWidth', 2.0);
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
xlim([0, 45])

ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 16)

ylim([0, ymax])


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
    MeanProfile2 = CompiledEmbryos.WindowedProfiles.DeltaFC.Control.mean(:,RefBin,ch_index).';
    x2 = CompiledEmbryos.WindowedProfiles.DeltaFC.x;
    plot(x2(~isnan(MeanProfile2)), MeanProfile2(~isnan(MeanProfile2)), 'r', 'LineWidth', 2.0);
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
xlim([0, 45])

ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 16)

ylim([0, ymax])


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
    MeanProfile2 = CompiledEmbryos.WindowedProfiles.DubuisTime.Test.mean(:,RefBin,ch_index).';
    x2 = CompiledEmbryos.WindowedProfiles.DubuisTime.x;
    plot(x2(~isnan(MeanProfile2)), MeanProfile2(~isnan(MeanProfile2)), 'r', 'LineWidth', 2.0);
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
xlim([0, 65])

ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 16)

ylim([0, ymax])


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
    MeanProfile2 = CompiledEmbryos.WindowedProfiles.DubuisTime.Control.mean(:,RefBin,ch_index).';
    x2 = CompiledEmbryos.WindowedProfiles.DubuisTime.x;
    plot(x2(~isnan(MeanProfile2)), MeanProfile2(~isnan(MeanProfile2)), 'r', 'LineWidth', 2.0);
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
xlim([0, 65])

ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [PlotLabel, ' (Control)'], 'FontSize', 18)



DeltaFCAx.FontSize = 16;

outpath = [ProfFigPath, filesep, SetLabel, '_',ChannelNames{ch_index}, '_ControlProfilesFig2D_TimeSmoothed.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end