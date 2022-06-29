
version = 6;
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
NumAPbins = length(APbins);

ChannelNames = {'', 'Hoechst', 'Bicoid', 'Knirps', 'Hunchback'};

load([AllCompiledPath, 'AllCompiledEmbryos.Mat'], 'AllCompiledEmbryos')
%%
ch_index = 3;
for exp_index = 1:15
    
    CompiledEmbryos = AllCompiledEmbryos{exp_index};
    x_sample = CompiledEmbryos.DubuisEmbryoTimes;
    UseTF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample);
    x_sample =x_sample(UseTF);
    ys =  CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.Set1.AllProfiles.mean(UseTF, :, ch_index);
    BadTF =  sum(isnan(ys), 2).' >min( sum(isnan(ys), 2).');
    x_sample = x_sample(~BadTF);
    ys = ys(~BadTF,:);
    xfits = CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.Set1.Test.x;
    yfits = CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.Set1.Test.mean(:,:,ch_index);
    
    
    x_sample_control = CompiledEmbryos.DubuisEmbryoTimes;
    UseTF_control = CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample_control);
    x_sample_control =x_sample_control(UseTF_control);
    ys_control =  CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.Set1.AllProfiles.mean(UseTF_control, :, ch_index);
    BadTF_control=  sum(isnan(ys_control), 2).' >min( sum(isnan(ys_control), 2).');
    x_sample_control = x_sample_control(~BadTF_control);
    ys_control = ys_control(~BadTF_control,:);
    xfits_control = CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.Set1.Control.x;
    yfits_control = CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.Set1.Control.mean(:,:,ch_index);
    
    SetLabel = AllSetInfo.SetLabels{exp_index};
    
    for APindex =13:13%5:4:37
        
        close all
        DeltaFCFixCorrectedFig =figure(1);
        set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .4, .4]);
        set(gcf,'color','w');
        DeltaFCAx = axes(DeltaFCFixCorrectedFig);
        scatter(x_sample, ys(:,APindex),...
            100, 'MarkerFaceColor', [.7 .7 .7],...
            'MarkerEdgeColor', [.7 .7 .7]);
        hold on
        plot(xfits, yfits(:,APindex), 'k', 'LineWidth', 2.0);
        ymax = max([max(yfits(:, APindex)),max(ys(:,APindex))]) *1.1;
        
        grid on
        hold off
        xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
        %xlabel('\delta_{FC} (\mum)', 'FontSize', 16)
        xlim([0, 45])
        
        ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 16)
        
        %ylim([0, min([ymax, 100])])
        xlim([0 60])
        
        DeltaFCAx.YAxis.FontSize = 16;
        DeltaFCAx.XAxis.FontSize = 16;
        
        if AllSetInfo.Replicates(exp_index) ~= 0
            PlotLabel = [num2str(AllSetInfo.Temperatures(exp_index)), '°C Replicate ', num2str(AllSetInfo.Replicates(exp_index)),' Test Set'];
        else
            PlotLabel = [num2str(AllSetInfo.Temperatures(exp_index)), '°C Flipped Stain Test Set'];
        end
        title(DeltaFCAx, [PlotLabel,' AP: ',num2str(APbins(APindex)), ' (Test)'], 'FontSize', 18)
        %
        
        
        DeltaFCAx.FontSize = 16;
        %outstring = ['T', strrep(num2str(AllSetInfo.Temperatures(i)), '.', '_'),'CRep', num2str(AllSetInfo.Replicates(i)), 'test','_'];
        outpath = [AllSetsProfFigPath, filesep,'SetComps' , filesep, SetLabel, '_Test_',ChannelNames{ch_index}, '_DeltaFCWindowed_AP',num2str(APindex),'.png'];
        saveas(DeltaFCFixCorrectedFig,outpath);
        
        
        DeltaFCFixCorrectedFig2 =figure(2);
        set(DeltaFCFixCorrectedFig2,'units', 'normalized', 'position',[0.5, 0.05, .4, .4]);
        set(gcf,'color','w');
        DeltaFCAx2 = axes(DeltaFCFixCorrectedFig2);
        scatter(x_sample_control, ys_control(:,APindex),...
            100, 'MarkerFaceColor', [.7 .7 .7],...
            'MarkerEdgeColor', [.7 .7 .7]);
        hold on
        plot(xfits_control, yfits_control(:,APindex), 'k', 'LineWidth', 2.0);
        ymax = max([max(yfits_control(:, APindex)),max(ys_control(:,APindex))]) *1.1;
        
        grid on
        hold off
        xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
        %xlabel('\delta_{FC} (\mum)', 'FontSize', 16)
        %xlim([0, 45])
        
        ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 16)
        
        ylim([0, min([ymax, 100])])
        xlim([0 60])
        
        DeltaFCAx2.YAxis.FontSize = 16;
        DeltaFCAx2.XAxis.FontSize = 16;
        
        if AllSetInfo.Replicates(exp_index) ~= 0
            PlotLabel = [num2str(AllSetInfo.Temperatures(exp_index)), '°C Replicate ', num2str(AllSetInfo.Replicates(exp_index)),' Control Set'];
        else
            PlotLabel = [num2str(AllSetInfo.Temperatures(exp_index)), '°C Flipped Stain Control Set'];
        end
        title(DeltaFCAx2, [PlotLabel,' AP: ',num2str(APbins(APindex)), ' (Control)'], 'FontSize', 18)
        %
        
        
        DeltaFCAx2.FontSize = 16;
        %outstring = ['T', strrep(num2str(AllSetInfo.Temperatures(i)), '.', '_'),'CRep', num2str(AllSetInfo.Replicates(i)), 'test','_'];
        outpath2 = [AllSetsProfFigPath, filesep,'SetComps' , filesep, SetLabel, '_Control_',ChannelNames{ch_index}, '_DeltaFCWindowed_AP',num2str(APindex),'.png'];
        saveas(DeltaFCFixCorrectedFig2,outpath2);
    end
end

%%
i_list = [1 1 2 5 8 10 10 11 13 13 14];
j_list = [2 3 3 6 9 11 12 12 14 15 15];

ch_index = 3;
for list_index = 1:length(i_list)
    
    exp_index = i_list(list_index);
    CompiledEmbryos = AllCompiledEmbryos{exp_index};
    x_sample = CompiledEmbryos.DubuisEmbryoTimes;
    UseTF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample);
    x_sample =x_sample(UseTF);
    ys =  CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.Set1.AllProfiles.mean(UseTF, :, ch_index);
    BadTF =  sum(isnan(ys), 2).' >min( sum(isnan(ys), 2).');
    x_sample = x_sample(~BadTF);
    ys = ys(~BadTF,:);
    xfits = CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.Set1.Test.x;
    yfits = CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.Set1.Test.mean(:,:,ch_index);
    
    exp_index2 = j_list(list_index);
    CompiledEmbryos2 = AllCompiledEmbryos{exp_index2};
    x_sample_control = CompiledEmbryos2.DubuisEmbryoTimes;
    UseTF_control = CompiledEmbryos2.TestSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample_control);
    x_sample_control =x_sample_control(UseTF_control);
    ys_control =  CompiledEmbryos2.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.Set1.AllProfiles.mean(UseTF_control, :, ch_index);
    BadTF_control=  sum(isnan(ys_control), 2).' >min( sum(isnan(ys_control), 2).');
    x_sample_control = x_sample_control(~BadTF_control);
    ys_control = ys_control(~BadTF_control,:);
    xfits_control = CompiledEmbryos2.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.Set1.Test.x;
    yfits_control = CompiledEmbryos2.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.Set1.Test.mean(:,:,ch_index);
    APmat = repmat(APbins, length(xfits), 1);
    xmat = repmat(xfits.', 1, NumAPbins);
    xindex_mat = repmat((1:length(xfits)).', 1, NumAPbins);
    
    yfits = yfits(1:10:size(yfits, 1),:).';
    yfits = yfits(:).';
    yfits_control = yfits_control(1:10:size(yfits_control, 1),:).';
    yfits_control = yfits_control(:).';
    xflat = xmat(1:10:size(xmat, 1),:);
    xflat = xflat(:).';
    xindex_flat = xindex_mat(1:10:size(xindex_mat, 1),:);
    xindex_flat = xindex_flat(:).';

    TFkeep = ~isnan(yfits) & ~isnan(yfits_control);
    x_plot = xindex_flat(TFkeep);
    yfits = yfits(TFkeep);
    yfits2 = yfits_control(TFkeep);
    
    SetLabel = AllSetInfo.SetLabels{exp_index};
    SetLabel2 = AllSetInfo.SetLabels{exp_index2};
    DubuisTimeRange = xfits_control;
    colorsDubuisTimes = hsv(length(DubuisTimeRange)); % Colormap "jet" is another option
    FractonalDubuisTimeRange = DubuisTimeRange/max(DubuisTimeRange);
    BinsToFit = zeros(NChannels, NumAPbins, 'logical');
    BinsToFit(2,9:33) = true;
    BinsToFit(3,5:13) = true;
    BinsToFit(4,22:37) = true;
    BinsToFit(5,5:21) = true;
    
    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .4, .4]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    map = colormap(colorsDubuisTimes);
    h = colorbar;
    % %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])

    colorTitleHandle = get(h,'Title');
    titleString = '\del (\mum)';
    titleString = 'time (min)';
    set(colorTitleHandle ,'String',titleString);
    h.Ticks =  FractonalDubuisTimeRange(1:50:701)+0.5*(FractonalDubuisTimeRange(2)-FractonalDubuisTimeRange(1)); %Create 8 ticks from zero to 1
    h.TickLabels = {'0.0', '5.0','10.0','15.0','20.0','25.0', '30.0', '35.0', '40.0','45.0','50.0','55.0', '60.0','65.0', '70.0'};
    hold on 
    for plot_index = 1:length(x_plot)
        scatter(yfits(plot_index),yfits2(plot_index),...
            75, 'MarkerFaceColor', colorsDubuisTimes(x_plot(plot_index),:),...
            'MarkerEdgeColor', colorsDubuisTimes(x_plot(plot_index),:));
        hold on
    end
    
   ymax = max([max(yfits),max(yfits2)]) *1.1;
   plot([0 ymax], [0 ymax], 'k', 'LineWidth', 2.0);
 
    
    grid on
    hold off
    if AllSetInfo.Replicates(exp_index) ~= 0
        xLabel = [num2str(AllSetInfo.Temperatures(exp_index)), '°C Replicate ', num2str(AllSetInfo.Replicates(exp_index)),' Test Set'];
    else
        xLabel = [num2str(AllSetInfo.Temperatures(exp_index)), '°C Flipped Stain Test Set'];
    end
    xlabel(xLabel, 'FontSize', 16)
    
    if AllSetInfo.Replicates(exp_index2) ~= 0
        yLabel = [num2str(AllSetInfo.Temperatures(exp_index2)), '°C Replicate ', num2str(AllSetInfo.Replicates(exp_index2)),' Test Set'];
    else
        yLabel = [num2str(AllSetInfo.Temperatures(exp_index2)), '°C Flipped Stain Test Set'];
    end
    ylabel(yLabel, 'FontSize', 16)
    %xlabel('\delta_{FC} (\mum)', 'FontSize', 16)
    xlim([0, min([ymax, 150])])
    
    title(DeltaFCAx, [ChannelNames{ch_index}, ' (AU)'], 'FontSize', 18)
    
    ylim([0, min([ymax, 150])])

    DeltaFCAx.YAxis.FontSize = 16;
    DeltaFCAx.XAxis.FontSize = 16;
    
    %
    
    
    DeltaFCAx.FontSize = 16;
    %outstring = ['T', strrep(num2str(AllSetInfo.Temperatures(i)), '.', '_'),'CRep', num2str(AllSetInfo.Replicates(i)), 'test','_'];
    outpath = [AllSetsProfFigPath, filesep,'Set1Comps' , filesep, SetLabel, '_Test_',SetLabel2,'_Test_', ChannelNames{ch_index}, '_DubuisColored.png'];
    saveas(DeltaFCFixCorrectedFig,outpath);
    
end