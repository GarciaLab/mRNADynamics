%clear all
SizeDataPath = 'S:/Gabriella/Dropbox/EmbryoSizeMeasurements/EmbryoSizeData.mat';
AllSetsProfFigPath = 'S:/Gabriella/Dropbox/ProteinProfiles/V5/Figures/';
AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/V5CompiledEmbryos/';
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

version = 3;
AllCompiledPath = ['S:/Gabriella/Dropbox/ProteinProfiles/V',num2str(version),'Profiles/'];
load([AllCompiledPath, 'AllCompiledEmbryos.Mat'], 'AllCompiledEmbryos','Universal_Imins','Universal_Imaxs','Imins','Imaxs');
AllSetInfo = GetFixedSetPrefixInfo;

NumSets = length(AllSetInfo.Temperatures);

NChannels = 5;

APbins = 0:0.025:1;
NumAPbins = length(APbins);

% 
% for i = 1:NumSets
%     AllCompiledEmbryos{i} = AddBinnedProfiles( AllCompiledEmbryos{i});
%     %AllCompiledEmbryos{i} = AddSmoothedProfiles( AllCompiledEmbryos{i});
% end

MasterSetPath = 'S:/Gabriella/Dropbox/ProteinProfiles/25CUnflippedMasterSets.mat';
load(MasterSetPath, 'CombinedMean', 'CombinedSE', 'CombinedCounts', 'Slopes', 'Intercepts', 'Fits', 'SubsetsIncluded',...
    'Ts', 'Reps', 'CTstrings', 'SubsetsIncluded');


%%
% AllCompiledEmbryos = {};
% for exp_index = 1:length(AllSetInfo.Temperatures)
% if AllSetInfo.Flipped(exp_index)
%     FlipString = 'yes';
% else
%     FlipString = 'no';
% end
% disp(['T = ', num2str(AllSetInfo.Temperatures(exp_index)), ', Rep: ', num2str(AllSetInfo.Replicates(exp_index)), ', Flipped: ', FlipString])
% SetLabel = AllSetInfo.SetLabels{exp_index};
% PlotLabel = AllSetInfo.PlotLabels{exp_index};
% SetPrefixes = AllSetInfo.Prefixes{exp_index};
% SetIsFlipped = AllSetInfo.Flipped(exp_index);
% ProfFigPath = [AllSetsProfFigPath, SetLabel];
% OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];
% if ~isdir(ProfFigPath)
%     mkdir(ProfFigPath)
% end
% if ~isdir(OutEmbryoPath)
%     mkdir(OutEmbryoPath)
% end
% liveExperiments = cell(1, length(SetPrefixes));
% for i = 1:length(SetPrefixes)
%     liveExperiments{i} = LiveExperiment(SetPrefixes{i});
% end
% FixedPixelSize_um = liveExperiments{1}.pixelSize_um;
% CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
% load(CEoutpath, 'CompiledEmbryos');
% AllCompiledEmbryos{exp_index} = CompiledEmbryos;
% end
% 
% 
% %% List of Profiles to rescale
% % SlideRescaledDorsalAvgAPProfiles
% % SlideRescaledDorsalAvgNarrowAPProfiles
% % SlideRescaledDorsalAPProfiles
% % SlideRescaledDorsalNarrowAPProfiles
% % --> The 4 above with FixCorrectedControlScalingFactors,
% % DubuisTimesControlScalingFactors, HisRFP25CTimesControlScalingFactors,
% % WindowedDubuisTimesControlScalingFactors,
% % WindowedDeltaFCControlScalingFactors or WindowedHisRFP25CTimesControlScalingFactors
% 
% % FixCorrectedSmoothedAvgAPProfiles with FixCorrectedControlScalingFactors
% % FixCorrectedSmoothedAvgNarrowAPProfiles
% % DubuisTimeSmoothedAvgAPProfiles with DubuisTimesControlScalingFactors
% % DubuisTimeSmoothedAvgNarrowAPProfiles
% % yw25CTimeSmoothedAvgAPProfiles
% % yw25CTimeSmoothedAvgNarrowAPProfiles
% % HisRFP25CTimeSmoothedAvgAPProfiles with HisRFP25CTimesControlScalingFactors
% % HisRFP25CTimeSmoothedAvgNarrowAPProfiles
% % BinnedProfiles with DON'T HAVE THE RIGHT STUFF CALCULATED
% % WindowedProfiles with WindowedDubuisTimesControlScalingFactors or
% % WindowedDeltaFCControlScalingFactors or WindowedHisRFP25CTimesControlScalingFactors
% for exp_index = 1:NumSets
%     disp(['Embryo Index = ', num2str(exp_index)]);
%     AllCompiledEmbryos{exp_index} = AddUniversalScalingProfiles( AllCompiledEmbryos{exp_index}, exp_index);
% end
% %%
%  [AllCompiledEmbryos, Universal_Imins, Universal_Imaxs, Imins, Imaxs]  = AddNormalizedProfiles(AllCompiledEmbryos);
% %%
% for exp_index = 1:NumSets
%     disp(['Embryo Index = ', num2str(exp_index)]);
%     AllCompiledEmbryos{exp_index} = AddNC13NormalizedProfiles( AllCompiledEmbryos{exp_index}, exp_index);
% end


%%
if ~isdir([AllSetsProfFigPath, 'SetComps'])
    mkdir([AllSetsProfFigPath, 'SetComps'])
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
    ProfSet1 = AllCompiledEmbryos{exp1}.Un.MasterDubuisTimesWindowedAvgAP.Test.mean(:,:, ch_index);
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.MasterDubuisTimesWindowedAvgAP.Test.mean(:,:, ch_index);
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.MasterDubuisTimesWindowedAvgAP.Test.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.MasterDubuisTimesWindowedAvgAP.Test.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    if isempty(PlotProfSet1) | isempty(PlotProfSet2)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
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

outpath = [AllSetsProfFigPath, filesep, 'SetComps',filesep,'T', strrep(num2str( num2str(unique_temperatures(temp_index))),'.','_'),'Ctest_StandardReplicates_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end

%%


SmoothedMasterSetPath = 'S:/Gabriella/Dropbox/ProteinProfiles/Smoothed25CUnflippedMasterSets.mat';
load(SmoothedMasterSetPath, 'CombinedMean', 'CombinedSE', 'CombinedCounts', 'Slopes', 'Intercepts', 'Fits', 'SubsetsIncluded',...
    'Ts', 'Reps', 'CTstrings', 'SubsetsIncluded');
close all
% SmoothMasterSetPath = 'S:/Gabriella/Dropbox/ProteinProfiles/25CMasterSets.mat';
% load(SmoothMasterSetPath, 'CombinedMean', 'CombinedSE', 'CombinedCounts', 'Slopes', 'Intercepts', 'Fits', 'SubsetsIncluded',...
%     'Ts', 'Reps', 'CTstrings', 'SubsetsIncluded');

unique_temperatures = unique(AllSetInfo.Temperatures);
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
DubuisTimeRange = 0:70;
colorsDubuisTimes = hsv(length(DubuisTimeRange)); % Colormap "jet" is another option
FractonalDubuisTimeRange = DubuisTimeRange/max(DubuisTimeRange);
BinsToFit = zeros(NChannels, NumAPbins, 'logical');
BinsToFit(2,9:33) = true;
BinsToFit(3,5:13) = true;
BinsToFit(4,22:37) = true;
BinsToFit(5,5:21) = true;
for exp_index = 1:15
for ch_index = [3 5]
    %ProfSet2 = AllCompiledEmbryos{exp_index}.NormalizedProfiles.MasterDubuisTimesSmoothedAvgAP.Control.Profiles(:,BinsToFit(ch_index,:), ch_index);
    ProfSet2 = AllCompiledEmbryos{exp_index}.UnivScaledProfiles.MasterDubuisTimesSmoothedAvgAP.Control.Profiles(:,BinsToFit(ch_index,:), ch_index);
    %ProfCounts2 = AllCompiledEmbryos{exp_index}.NormalizedProfiles.MasterDubuisTimesWindowedAvgAP.Control.count(:,BinsToFit(ch_index,:), ch_index);
    ProfSet2 = ProfSet2(:).';
    %ProfCounts2 = ProfCounts2(:).';
    
    
    ProfSet1 = CombinedMean(:,BinsToFit(ch_index,:),ch_index);
    ProfSet1 = ProfSet1(:).';
    ProfCounts1 = CombinedCounts(:,BinsToFit(ch_index,:),ch_index);
    ProfCounts1 = ProfCounts1(:).';
    TfNaN = isnan(ProfSet1) | isnan(ProfSet2);
    PlotProfSet1 = ProfSet1(~TfNaN);
    PlotProfSet2 = ProfSet2(~TfNaN);
%     PlotProfSet1 = ProfSet1(ProfCounts1 >= 5 & ProfCounts2 >= 3);
%     PlotProfSet2 = ProfSet2(ProfCounts1 >= 5 & ProfCounts2 >= 3);
%     
    
    DubuisTimesMat = repmat(AllCompiledEmbryos{exp_index}.DubuisSmoothedTimes.', 1, sum(BinsToFit(ch_index,:)));
    DubuisTimesFlat = DubuisTimesMat(:).';
    %PlotDubuisTimes = DubuisTimesFlat(ProfCounts1 >= 5 & ProfCounts2 >= 3);
    PlotDubuisTimes = DubuisTimesFlat(~TfNaN);
%     if isempty(PlotProfSet1) | isempty(PlotProfSet2)
%         PlotProfSet1 = ProfSet1(ProfCounts1 >= 3 & ProfCounts2 >= 1);
%         PlotProfSet2 = ProfSet2(ProfCounts1 >= 3 & ProfCounts2 >= 1);
%         PlotDubuisTimes = DubuisTimesFlat(ProfCounts1 >= 3 & ProfCounts2 >= 1);
%         AddTitleStr = ' (Low Counts)';
%         if isempty(PlotProfSet1)| isempty(PlotProfSet2)
%             PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
%             PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
%             PlotDubuisTimes = DubuisTimesFlat(ProfCounts1 >= 3 & ProfCounts2 >= 1);
%             AddTitleStr = ' (Very Low Counts)';
%             
%         end
%     else
%         AddTitleStr = '';
%     end
    
    if isempty(PlotProfSet1) | isempty(PlotProfSet2)
        continue
    end
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    
    map = colormap(colorsDubuisTimes);
    h = colorbar;
    % %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
    hold off
    colorTitleHandle = get(h,'Title');
    titleString = '\del (\mum)';
    titleString = 'time (min)';
    set(colorTitleHandle ,'String',titleString);
    h.Ticks =  FractonalDubuisTimeRange(1:5:71)+0.5*(FractonalDubuisTimeRange(2)-FractonalDubuisTimeRange(1)); %Create 8 ticks from zero to 1
    h.TickLabels = {'0.0', '5.0','10.0','15.0','20.0','25.0', '30.0', '35.0', '40.0','45.0','50.0','55.0', '60.0', '65.0', '70.0'};
    hold on
    
    for plot_index = 1:length(PlotProfSet1)
        scatter(PlotProfSet1(plot_index),PlotProfSet2(plot_index),...
            75, 'MarkerFaceColor', colorsDubuisTimes(PlotDubuisTimes(plot_index)+1,:),...
            'MarkerEdgeColor', colorsDubuisTimes(PlotDubuisTimes(plot_index)+1,:));
        hold on
    end
    
    ymax = ceil(max(PlotProfSet2)/5)*5;
    xmax = ceil(max(PlotProfSet1)/5)*5;
    ymax = max(PlotProfSet2);
    %             xmax = ceil(max([AllAPProfileFlatExp1 ; AllAPProfileFlatExp2])/.5)*.5;
    %             ymax = ceil(max(AllAPProfileFlatExp2)/5)*5;
    %             xmax = ceil(max(AllAPProfileFlatExp1)/5)*5;
    plot([0, xmax], [0, ymax], 'k');
    %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
    %     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
    xlim([0, xmax])
    xlab ='Master 25ºC Profile (AU)';
    ylab = ['Control for ', num2str(AllSetInfo.Temperatures(exp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp_index))];
    %ylab = [num2str(AllSetInfo.Temperatures(exp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp_index))];
    
    xlabel(xlab, 'FontSize', 16)
    ylabel(ylab, 'FontSize', 16)
    
    ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [ChannelNames{ch_index}], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;
if AllSetInfo.Replicates(exp_index) == 0
    outstring = ['Master25CProfile_T', strrep(num2str(AllSetInfo.Temperatures(exp_index)), '.', '_'),'Cflipped_control'];        
else
outstring = ['Master25CProfile_T', strrep(num2str(AllSetInfo.Temperatures(exp_index)), '.', '_'),'Crep', num2str(AllSetInfo.Replicates(exp_index)), '_test'];
end 
outpath = [AllSetsProfFigPath, filesep, 'SetComps',filesep,outstring, '_', ChannelNames{ch_index}, '_DubuisTimesColor.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end


%%


SmoothedMasterSetPath = 'S:/Gabriella/Dropbox/ProteinProfiles/Smoothed25CUnflippedMasterSets.mat';
load(SmoothedMasterSetPath, 'CombinedMean', 'CombinedSE', 'CombinedCounts', 'Slopes', 'Intercepts', 'Fits', 'SubsetsIncluded',...
    'Ts', 'Reps', 'CTstrings', 'SubsetsIncluded');
close all
% SmoothMasterSetPath = 'S:/Gabriella/Dropbox/ProteinProfiles/25CMasterSets.mat';
% load(SmoothMasterSetPath, 'CombinedMean', 'CombinedSE', 'CombinedCounts', 'Slopes', 'Intercepts', 'Fits', 'SubsetsIncluded',...
%     'Ts', 'Reps', 'CTstrings', 'SubsetsIncluded');

unique_temperatures = unique(AllSetInfo.Temperatures);
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
DubuisTimeRange = 0:70;
colorsDubuisTimes = hsv(length(DubuisTimeRange)); % Colormap "jet" is another option
FractonalDubuisTimeRange = DubuisTimeRange/max(DubuisTimeRange);
BinsToFit = zeros(NChannels, NumAPbins, 'logical');
BinsToFit(2,9:33) = true;
BinsToFit(3,5:13) = true;
BinsToFit(4,22:37) = true;
BinsToFit(5,5:21) = true;
for exp_index = 1:15
for ch_index = [3 5]
    %ProfSet2 = AllCompiledEmbryos{exp_index}.NormalizedProfiles.MasterDubuisTimesSmoothedAvgAP.Control.Profiles(:,BinsToFit(ch_index,:), ch_index);
    ProfSet2 = AllCompiledEmbryos{exp_index}.UnivScaledProfiles.MasterDubuisTimesSmoothedAvgAP.Control.Profiles(:,BinsToFit(ch_index,:), ch_index);
    %ProfCounts2 = AllCompiledEmbryos{exp_index}.NormalizedProfiles.MasterDubuisTimesWindowedAvgAP.Control.count(:,BinsToFit(ch_index,:), ch_index);
    ProfSet2 = ProfSet2(:).';
    %ProfCounts2 = ProfCounts2(:).';
    
    
    ProfSet1 = CombinedMean(:,BinsToFit(ch_index,:),ch_index);
    ProfSet1 = ProfSet1(:).';
    ProfCounts1 = CombinedCounts(:,BinsToFit(ch_index,:),ch_index);
    ProfCounts1 = ProfCounts1(:).';
    
    ProfSet1 = AllCompiledEmbryos{exp_index}.UnivScaledProfiles.MasterDubuisTimesSmoothedAvgAP.Test.Profiles(:,BinsToFit(ch_index,:), ch_index);
    %ProfCounts2 = AllCompiledEmbryos{exp_index}.NormalizedProfiles.MasterDubuisTimesWindowedAvgAP.Control.count(:,BinsToFit(ch_index,:), ch_index);
    ProfSet1 = ProfSet1(:).';
    TfNaN = isnan(ProfSet1) | isnan(ProfSet2);
    PlotProfSet1 = ProfSet1(~TfNaN);
    PlotProfSet2 = ProfSet2(~TfNaN);
%     PlotProfSet1 = ProfSet1(ProfCounts1 >= 5 & ProfCounts2 >= 3);
%     PlotProfSet2 = ProfSet2(ProfCounts1 >= 5 & ProfCounts2 >= 3);
%     
    
    DubuisTimesMat = repmat(AllCompiledEmbryos{exp_index}.DubuisSmoothedTimes.', 1, sum(BinsToFit(ch_index,:)));
    DubuisTimesFlat = DubuisTimesMat(:).';
    %PlotDubuisTimes = DubuisTimesFlat(ProfCounts1 >= 5 & ProfCounts2 >= 3);
    PlotDubuisTimes = DubuisTimesFlat(~TfNaN);
%     if isempty(PlotProfSet1) | isempty(PlotProfSet2)
%         PlotProfSet1 = ProfSet1(ProfCounts1 >= 3 & ProfCounts2 >= 1);
%         PlotProfSet2 = ProfSet2(ProfCounts1 >= 3 & ProfCounts2 >= 1);
%         PlotDubuisTimes = DubuisTimesFlat(ProfCounts1 >= 3 & ProfCounts2 >= 1);
%         AddTitleStr = ' (Low Counts)';
%         if isempty(PlotProfSet1)| isempty(PlotProfSet2)
%             PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
%             PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
%             PlotDubuisTimes = DubuisTimesFlat(ProfCounts1 >= 3 & ProfCounts2 >= 1);
%             AddTitleStr = ' (Very Low Counts)';
%             
%         end
%     else
%         AddTitleStr = '';
%     end
    
    if isempty(PlotProfSet1) | isempty(PlotProfSet2)
        continue
    end
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    
    map = colormap(colorsDubuisTimes);
    h = colorbar;
    % %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
    hold off
    colorTitleHandle = get(h,'Title');
    titleString = '\del (\mum)';
    titleString = 'time (min)';
    set(colorTitleHandle ,'String',titleString);
    h.Ticks =  FractonalDubuisTimeRange(1:5:71)+0.5*(FractonalDubuisTimeRange(2)-FractonalDubuisTimeRange(1)); %Create 8 ticks from zero to 1
    h.TickLabels = {'0.0', '5.0','10.0','15.0','20.0','25.0', '30.0', '35.0', '40.0','45.0','50.0','55.0', '60.0', '65.0', '70.0'};
    hold on
    
    for plot_index = 1:length(PlotProfSet1)
        scatter(PlotProfSet1(plot_index),PlotProfSet2(plot_index),...
            75, 'MarkerFaceColor', colorsDubuisTimes(PlotDubuisTimes(plot_index)+1,:),...
            'MarkerEdgeColor', colorsDubuisTimes(PlotDubuisTimes(plot_index)+1,:));
        hold on
    end
    
    ymax = ceil(max(PlotProfSet2)/5)*5;
    xmax = ceil(max(PlotProfSet1)/5)*5;
    ymax = max(PlotProfSet2);
    %             xmax = ceil(max([AllAPProfileFlatExp1 ; AllAPProfileFlatExp2])/.5)*.5;
    %             ymax = ceil(max(AllAPProfileFlatExp2)/5)*5;
    %             xmax = ceil(max(AllAPProfileFlatExp1)/5)*5;
    plot([0, xmax], [0, ymax], 'k');
    %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
    %     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
    xlim([0, xmax])
    xlab ='Master 25ºC Profile (AU)';
    xlab = [num2str(AllSetInfo.Temperatures(exp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp_index))];
    ylab = ['Control for ', num2str(AllSetInfo.Temperatures(exp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp_index))];
    %ylab = [num2str(AllSetInfo.Temperatures(exp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp_index))];
    
    xlabel(xlab, 'FontSize', 16)
    ylabel(ylab, 'FontSize', 16)
    
    ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [ChannelNames{ch_index}], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;
if AllSetInfo.Replicates(exp_index) == 0
    outstring = ['T', strrep(num2str(AllSetInfo.Temperatures(exp_index)), '.', '_'),'Cflipped_TestVsControl'];        
else
outstring = ['T', strrep(num2str(AllSetInfo.Temperatures(exp_index)), '.', '_'),'Crep', num2str(AllSetInfo.Replicates(exp_index)), '_TestVsControl'];
end 
outpath = [AllSetsProfFigPath, filesep, 'SetComps',filesep,outstring, '_', ChannelNames{ch_index}, '_DubuisTimesColor.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end

%%

SmoothedMasterSetPath = 'S:/Gabriella/Dropbox/ProteinProfiles/Smoothed25CUnflippedMasterSets.mat';
load(SmoothedMasterSetPath, 'CombinedMean', 'CombinedSE', 'CombinedCounts', 'Slopes', 'Intercepts', 'Fits', 'SubsetsIncluded',...
    'Ts', 'Reps', 'CTstrings', 'SubsetsIncluded');
close all
% SmoothMasterSetPath = 'S:/Gabriella/Dropbox/ProteinProfiles/25CMasterSets.mat';
% load(SmoothMasterSetPath, 'CombinedMean', 'CombinedSE', 'CombinedCounts', 'Slopes', 'Intercepts', 'Fits', 'SubsetsIncluded',...
%     'Ts', 'Reps', 'CTstrings', 'SubsetsIncluded');

unique_temperatures = unique(AllSetInfo.Temperatures);
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
DubuisTimeRange = 0:70;
colorsDubuisTimes = hsv(length(DubuisTimeRange)); % Colormap "jet" is another option
FractonalDubuisTimeRange = DubuisTimeRange/max(DubuisTimeRange);
BinsToFit = zeros(NChannels, NumAPbins, 'logical');
BinsToFit(2,9:33) = true;
BinsToFit(3,5:13) = true;
BinsToFit(4,22:37) = true;
BinsToFit(5,5:21) = true;
i_list = [1 1 2 4 4 5 7 7 8 10 10 11 13 13 14];
j_list =  [2 3 3 5 6 6 8 9 9 11 12 12 14 15 15];
for i = 1:15
    exp_index = i_list(i);
    exp_index2 =j_list(i);
for ch_index = [3 5]
    %ProfSet2 = AllCompiledEmbryos{exp_index}.NormalizedProfiles.MasterDubuisTimesSmoothedAvgAP.Control.Profiles(:,BinsToFit(ch_index,:), ch_index);
    ProfSet2 = AllCompiledEmbryos{exp_index2}.UnivScaledProfiles.MasterDubuisTimesSmoothedAvgAP.Test.Profiles(:,BinsToFit(ch_index,:), ch_index);
    %ProfCounts2 = AllCompiledEmbryos{exp_index}.NormalizedProfiles.MasterDubuisTimesWindowedAvgAP.Control.count(:,BinsToFit(ch_index,:), ch_index);
    ProfSet2 = ProfSet2(:).';
    %ProfCounts2 = ProfCounts2(:).';
    
    
    ProfSet1 = CombinedMean(:,BinsToFit(ch_index,:),ch_index);
    ProfSet1 = ProfSet1(:).';
    ProfCounts1 = CombinedCounts(:,BinsToFit(ch_index,:),ch_index);
    ProfCounts1 = ProfCounts1(:).';
    
    ProfSet1 = AllCompiledEmbryos{exp_index}.UnivScaledProfiles.MasterDubuisTimesSmoothedAvgAP.Test.Profiles(:,BinsToFit(ch_index,:), ch_index);
    %ProfCounts2 = AllCompiledEmbryos{exp_index}.NormalizedProfiles.MasterDubuisTimesWindowedAvgAP.Control.count(:,BinsToFit(ch_index,:), ch_index);
    ProfSet1 = ProfSet1(:).';
    TfNaN = isnan(ProfSet1) | isnan(ProfSet2);
    PlotProfSet1 = ProfSet1(~TfNaN);
    PlotProfSet2 = ProfSet2(~TfNaN);
%     PlotProfSet1 = ProfSet1(ProfCounts1 >= 5 & ProfCounts2 >= 3);
%     PlotProfSet2 = ProfSet2(ProfCounts1 >= 5 & ProfCounts2 >= 3);
%     
    
    DubuisTimesMat = repmat(AllCompiledEmbryos{exp_index}.DubuisSmoothedTimes.', 1, sum(BinsToFit(ch_index,:)));
    DubuisTimesFlat = DubuisTimesMat(:).';
    %PlotDubuisTimes = DubuisTimesFlat(ProfCounts1 >= 5 & ProfCounts2 >= 3);
    PlotDubuisTimes = DubuisTimesFlat(~TfNaN);
%     if isempty(PlotProfSet1) | isempty(PlotProfSet2)
%         PlotProfSet1 = ProfSet1(ProfCounts1 >= 3 & ProfCounts2 >= 1);
%         PlotProfSet2 = ProfSet2(ProfCounts1 >= 3 & ProfCounts2 >= 1);
%         PlotDubuisTimes = DubuisTimesFlat(ProfCounts1 >= 3 & ProfCounts2 >= 1);
%         AddTitleStr = ' (Low Counts)';
%         if isempty(PlotProfSet1)| isempty(PlotProfSet2)
%             PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
%             PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
%             PlotDubuisTimes = DubuisTimesFlat(ProfCounts1 >= 3 & ProfCounts2 >= 1);
%             AddTitleStr = ' (Very Low Counts)';
%             
%         end
%     else
%         AddTitleStr = '';
%     end
    
    if isempty(PlotProfSet1) | isempty(PlotProfSet2)
        continue
    end
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    
    map = colormap(colorsDubuisTimes);
    h = colorbar;
    % %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
    hold off
    colorTitleHandle = get(h,'Title');
    titleString = '\del (\mum)';
    titleString = 'time (min)';
    set(colorTitleHandle ,'String',titleString);
    h.Ticks =  FractonalDubuisTimeRange(1:5:71)+0.5*(FractonalDubuisTimeRange(2)-FractonalDubuisTimeRange(1)); %Create 8 ticks from zero to 1
    h.TickLabels = {'0.0', '5.0','10.0','15.0','20.0','25.0', '30.0', '35.0', '40.0','45.0','50.0','55.0', '60.0', '65.0', '70.0'};
    hold on
    
    for plot_index = 1:length(PlotProfSet1)
        scatter(PlotProfSet1(plot_index),PlotProfSet2(plot_index),...
            75, 'MarkerFaceColor', colorsDubuisTimes(PlotDubuisTimes(plot_index)+1,:),...
            'MarkerEdgeColor', colorsDubuisTimes(PlotDubuisTimes(plot_index)+1,:));
        hold on
    end
    
    ymax = ceil(max(PlotProfSet2)/5)*5;
    xmax = ceil(max(PlotProfSet1)/5)*5;
    ymax = max(PlotProfSet2);
    %             xmax = ceil(max([AllAPProfileFlatExp1 ; AllAPProfileFlatExp2])/.5)*.5;
    %             ymax = ceil(max(AllAPProfileFlatExp2)/5)*5;
    %             xmax = ceil(max(AllAPProfileFlatExp1)/5)*5;
    plot([0, max([xmax, ymax])], [0,  max([xmax, ymax])], 'k');
    %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
    %     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
    xlim([0, xmax])
    xlab ='Master 25ºC Profile (AU)';
    xlab = [num2str(AllSetInfo.Temperatures(exp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp_index))];
    ylab = [num2str(AllSetInfo.Temperatures(exp_index2)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp_index2))];
    %ylab = [num2str(AllSetInfo.Temperatures(exp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp_index))];
    
    xlabel(xlab, 'FontSize', 16)
    ylabel(ylab, 'FontSize', 16)
    
    ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [ChannelNames{ch_index}], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;
if AllSetInfo.Replicates(exp_index) == 0
    if AllSetInfo.Replicates(exp_index2) == 0
    outstring = ['T', strrep(num2str(AllSetInfo.Temperatures(exp_index)), '.', '_'),'Cflipped_T', strrep(num2str(AllSetInfo.Temperatures(exp_index2)), '.', '_'),'Cflipped'];     
    else
        outstring = ['T', strrep(num2str(AllSetInfo.Temperatures(exp_index)), '.', '_'),'Cflipped_T', strrep(num2str(AllSetInfo.Temperatures(exp_index2)), '.', '_'),'Crep', num2str(AllSetInfo.Replicates(exp_index2))];   
    end
else
    if AllSetInfo.Replicates(exp_index2) == 0
outstring = ['T', strrep(num2str(AllSetInfo.Temperatures(exp_index)), '.', '_'),'Crep', num2str(AllSetInfo.Replicates(exp_index)), '_T', strrep(num2str(AllSetInfo.Temperatures(exp_index2)), '.', '_'),'Cflipped'];
    else
       outstring = ['T', strrep(num2str(AllSetInfo.Temperatures(exp_index)), '.', '_'),'Crep', num2str(AllSetInfo.Replicates(exp_index)), '_T', strrep(num2str(AllSetInfo.Temperatures(exp_index2)), '.', '_'),'Crep', num2str(AllSetInfo.Replicates(exp_index2))];
    end
end 
outpath = [AllSetsProfFigPath, filesep, 'SetComps',filesep,outstring, '_', ChannelNames{ch_index}, '_DubuisTimesColor.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
   
end
