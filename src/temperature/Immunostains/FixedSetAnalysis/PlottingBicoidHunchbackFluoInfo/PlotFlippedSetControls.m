function PlotFlippedSetControls(AllCompiledEmbryos)
%%
APbins = 0:0.025:1;
NumAPbins = length(APbins);

NarrowAPbins = 0:0.0125:1;
NumNarrowAPbins = length(NarrowAPbins);
NumAPbins = length(APbins);

ChannelNames = {'', 'Hoechst', 'Bicoid', 'Knirps', 'Hunchback'};
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32;
    76 0 153;
    255 0 255]/255;%0.9290 0.6940 0.1250]/255;

rainbow_colors = [212 84 66;
    255 128 0;
    237 176 32;
    100 149 93;
    66 87 168;
    76 0 153;
    255 0 255]/255;

SetLineStyles = {'-', '--', ':'};


temps = [25 27.5 22.5 20 17.5];
%% Make individual temperature plots with background subtracted normalized bootstrapped profiles

times_to_plot = 10:5:50;
colorsDubuisTimes = (hsv(round(length(times_to_plot)*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:length(times_to_plot),:);
FractonalDubuisTimeRange = (1:2:(2*length(times_to_plot)-1))/(2*length(times_to_plot));
for master_index = [1 2 7]
    for ch_index = [3 5]
        
        exp_indices = cell(1, length(temps));
        exp_index = zeros(1, length(temps));
        exp_index2 = zeros(1, length(temps));
        exp_index3 = zeros(1, length(temps));
        x_sample = cell(3, length(temps));
        UseTF = cell(3, length(temps));
        ys = cell(3, length(temps));
        BadTF = cell(3, length(temps));
        xfits = cell(3, length(temps));
        yfits = cell(3, length(temps));
        SetLabel = cell(3, length(temps));
        for temp_index = 1:length(temps)
            exp_indices{temp_index} = find(AllSetInfo.Temperatures == temps(temp_index));
            exp_index(temp_index) =  exp_indices{temp_index}(1);
            exp_index2(temp_index) =  exp_indices{temp_index}(2);
            exp_index3(temp_index) =  exp_indices{temp_index}(3);
            
            
            CompiledEmbryos = AllCompiledEmbryos{exp_index(temp_index)};
            x_sample{1, temp_index} = CompiledEmbryos.DubuisEmbryoTimes;
            UseTF{1, temp_index} = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample{1, temp_index});
            x_sample{1, temp_index} =x_sample{1, temp_index}(UseTF{1, temp_index});
            ys{1, temp_index} =  CompiledEmbryos.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{1, temp_index},:,ch_index, master_index);
            BadTF{1, temp_index} =  sum(isnan(ys{1, temp_index}), 2).' >NumAPbins;
            x_sample{1, temp_index} = x_sample{1, temp_index}(~BadTF{1, temp_index});
            ys{1, temp_index} = ys{1, temp_index}(~BadTF{1, temp_index},:);
            xfits{1, temp_index} = CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{1, temp_index} = CompiledEmbryos.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            SetLabel{1, temp_index} = AllSetInfo.SetLabels{exp_index(temp_index)};
            
            CompiledEmbryos2 = AllCompiledEmbryos{exp_index2(temp_index)};
            x_sample{2, temp_index} = CompiledEmbryos2.DubuisEmbryoTimes;
            UseTF{2, temp_index} = CompiledEmbryos2.TestSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample{2, temp_index});
            x_sample{2, temp_index} =x_sample{2, temp_index}(UseTF{2, temp_index});
            ys{2, temp_index} =  CompiledEmbryos2.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{2, temp_index},:,ch_index, master_index);
            BadTF{2, temp_index} =  sum(isnan(ys{2, temp_index}), 2).' >NumAPbins;
            x_sample{2, temp_index} = x_sample{2, temp_index}(~BadTF{2, temp_index});
            ys{2, temp_index} = ys{2, temp_index}(~BadTF{2, temp_index},:);
            xfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{2, temp_index} = CompiledEmbryos2.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            SetLabel{2, temp_index} = AllSetInfo.SetLabels{exp_index2(temp_index)};
            
            
            CompiledEmbryos3 = AllCompiledEmbryos{exp_index3(temp_index)};
            x_sample{3, temp_index} = CompiledEmbryos3.DubuisEmbryoTimes;
            UseTF{3, temp_index} = CompiledEmbryos3.TestSetEmbryos & CompiledEmbryos3.IsNC14 & ~isnan(x_sample{3, temp_index});
            x_sample{3, temp_index} =x_sample{3, temp_index}(UseTF{3, temp_index});
            ys{3, temp_index} =  CompiledEmbryos3.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{3, temp_index},:,ch_index, master_index);
            BadTF{3, temp_index} =  sum(isnan(ys{3, temp_index}), 2).' >NumAPbins;
            x_sample{3, temp_index} = x_sample{3, temp_index}(~BadTF{3, temp_index});
            ys{3, temp_index} = ys{3, temp_index}(~BadTF{3, temp_index},:);
            xfits{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{3, temp_index} = CompiledEmbryos3.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            SetLabel{3, temp_index} = AllSetInfo.SetLabels{exp_index3(temp_index)};
            
        end
        for temp_index = 1:length(temps)
            ymaxes(temp_index) = max([max( yfits{1, temp_index}(:)), max( yfits{2, temp_index}(:)), max( yfits{3, temp_index}(:))]);
            ymins(temp_index) = min([min( yfits{1, temp_index}(:)), min( yfits{2, temp_index}(:)), min( yfits{3, temp_index}(:))]);
        end
        plot_ymax = ceil(max(ymaxes)/0.2)*0.2;
        plot_ymin = floor(min(ymins)/0.2)*0.2;
        for temp_index = 1:length(temps)
            close all
            DeltaFCFixCorrectedFig =figure(1);
            set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .5, .6]);
            set(gcf,'color','w');
            DeltaFCAx = axes(DeltaFCFixCorrectedFig);
            map = colormap(colorsDubuisTimes);
            h = colorbar;
            % %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
            
            colorTitleHandle = get(h,'Title');
            titleString = '\del (\mum)';
            titleString = 'time (min)';
            set(colorTitleHandle ,'String',titleString);
            h.Ticks =  FractonalDubuisTimeRange; %Create 8 ticks from zero to 1
            h.TickLabels = {'10.0','15.0','20.0','25.0', '30.0', '35.0', '40.0','45.0','50.0'};
            hold on
            for time_index = 1:length(times_to_plot)
                
                for rep_index = 1:3
                    timebin = find(round(xfits{rep_index, temp_index}) == times_to_plot(time_index));
                    if ~isempty(timebin)
                        plot(APbins, yfits{rep_index, temp_index}(timebin,:),SetLineStyles{rep_index},...
                            'LineWidth', 2.0, 'Color', colorsDubuisTimes(time_index, :));
                        hold on
                    end
                end
                
            end
            grid on
            hold off
            xlabel('x/L', 'FontSize', 18)
            ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 18)
            ylim([plot_ymin,  plot_ymax])
            xlim([0.1 0.9])
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;
            title(DeltaFCAx, ['T = ',num2str(temps(temp_index)), '°C Replicates'], 'FontSize', 20)
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'FlippedStainReplicateComparisonAPProfiles', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'FlippedStainReplicateComparisonAPProfiles', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
            end
            outpath2 = [AllSetsProfFigPath, filesep, 'FlippedStainReplicateComparisonAPProfiles', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet_',num2str(master_index),'_T', strrep(num2str(temps(temp_index)), '.', '_'), 'C.png' ];
            saveas(DeltaFCFixCorrectedFig,outpath2);
        end
        
    end
end
%% Make individual temperature plots with background subtracted normalized bootstrapped profiles

times_to_plot = 10:5:50;
colorsDubuisTimes = (hsv(round(length(times_to_plot)*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:length(times_to_plot),:);
FractonalDubuisTimeRange = (1:2:(2*length(times_to_plot)-1))/(2*length(times_to_plot));
for master_index = [1 2 7]
    for ch_index = [3 5]
        
        exp_indices = cell(1, length(temps));
        exp_index = zeros(1, length(temps));
        exp_index2 = zeros(1, length(temps));
        exp_index3 = zeros(1, length(temps));
        x_sample = cell(3, length(temps));
        UseTF = cell(3, length(temps));
        ys = cell(3, length(temps));
        BadTF = cell(3, length(temps));
        xfits = cell(3, length(temps));
        yfits = cell(3, length(temps));
        SetLabel = cell(3, length(temps));
        for temp_index = 1:length(temps)
            exp_indices{temp_index} = find(AllSetInfo.Temperatures == temps(temp_index));
            exp_index(temp_index) =  exp_indices{temp_index}(1);
            exp_index2(temp_index) =  exp_indices{temp_index}(2);
            exp_index3(temp_index) =  exp_indices{temp_index}(3);
            
            
            CompiledEmbryos = AllCompiledEmbryos{exp_index(temp_index)};
            x_sample{1, temp_index} = CompiledEmbryos.DubuisEmbryoTimes;
            UseTF{1, temp_index} = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample{1, temp_index});
            x_sample{1, temp_index} =x_sample{1, temp_index}(UseTF{1, temp_index});
            ys{1, temp_index} =  CompiledEmbryos.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{1, temp_index},:,ch_index, master_index);
            BadTF{1, temp_index} =  sum(isnan(ys{1, temp_index}), 2).' >NumAPbins;
            x_sample{1, temp_index} = x_sample{1, temp_index}(~BadTF{1, temp_index});
            ys{1, temp_index} = ys{1, temp_index}(~BadTF{1, temp_index},:);
            xfits{1, temp_index} = CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{1, temp_index} = CompiledEmbryos.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            SetLabel{1, temp_index} = AllSetInfo.SetLabels{exp_index(temp_index)};
            
            CompiledEmbryos2 = AllCompiledEmbryos{exp_index2(temp_index)};
            x_sample{2, temp_index} = CompiledEmbryos2.DubuisEmbryoTimes;
            UseTF{2, temp_index} = CompiledEmbryos2.TestSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample{2, temp_index});
            x_sample{2, temp_index} =x_sample{2, temp_index}(UseTF{2, temp_index});
            ys{2, temp_index} =  CompiledEmbryos2.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{2, temp_index},:,ch_index, master_index);
            BadTF{2, temp_index} =  sum(isnan(ys{2, temp_index}), 2).' >NumAPbins;
            x_sample{2, temp_index} = x_sample{2, temp_index}(~BadTF{2, temp_index});
            ys{2, temp_index} = ys{2, temp_index}(~BadTF{2, temp_index},:);
            xfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{2, temp_index} = CompiledEmbryos2.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            SetLabel{2, temp_index} = AllSetInfo.SetLabels{exp_index2(temp_index)};
            
            
            CompiledEmbryos3 = AllCompiledEmbryos{exp_index3(temp_index)};
            x_sample{3, temp_index} = CompiledEmbryos3.DubuisEmbryoTimes;
            UseTF{3, temp_index} = CompiledEmbryos3.TestSetEmbryos & CompiledEmbryos3.IsNC14 & ~isnan(x_sample{3, temp_index});
            x_sample{3, temp_index} =x_sample{3, temp_index}(UseTF{3, temp_index});
            ys{3, temp_index} =  CompiledEmbryos3.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{3, temp_index},:,ch_index, master_index);
            BadTF{3, temp_index} =  sum(isnan(ys{3, temp_index}), 2).' >NumAPbins;
            x_sample{3, temp_index} = x_sample{3, temp_index}(~BadTF{3, temp_index});
            ys{3, temp_index} = ys{3, temp_index}(~BadTF{3, temp_index},:);
            xfits{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{3, temp_index} = CompiledEmbryos3.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            SetLabel{3, temp_index} = AllSetInfo.SetLabels{exp_index3(temp_index)};
            
        end
        for temp_index = 1:length(temps)
            ymaxes(temp_index) = max([max( yfits{1, temp_index}(:)), max( yfits{2, temp_index}(:)), max( yfits{3, temp_index}(:))]);
            ymins(temp_index) = min([min( yfits{1, temp_index}(:)), min( yfits{2, temp_index}(:)), min( yfits{3, temp_index}(:))]);
        end
        plot_ymax = ceil(max(ymaxes)/0.2)*0.2;
        plot_ymin = floor(min(ymins)/0.2)*0.2;
        for temp_index = 1:length(temps)
            close all
            DeltaFCFixCorrectedFig =figure(1);
            set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .5, .6]);
            set(gcf,'color','w');
            DeltaFCAx = axes(DeltaFCFixCorrectedFig);
            map = colormap(colorsDubuisTimes);
            h = colorbar;
            % %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
            
            colorTitleHandle = get(h,'Title');
            titleString = '\del (\mum)';
            titleString = 'time (min)';
            set(colorTitleHandle ,'String',titleString);
            h.Ticks =  FractonalDubuisTimeRange; %Create 8 ticks from zero to 1
            h.TickLabels = {'10.0','15.0','20.0','25.0', '30.0', '35.0', '40.0','45.0','50.0'};
            hold on
            for time_index = 1:length(times_to_plot)
                
                for rep_index = 1:3
                    timebin = find(round(xfits{rep_index, temp_index}) == times_to_plot(time_index));
                    if ~isempty(timebin)
                        plot(APbins, yfits{rep_index, temp_index}(timebin,:),SetLineStyles{rep_index},...
                            'LineWidth', 2.0, 'Color', colorsDubuisTimes(time_index, :));
                        hold on
                    end
                end
                
            end
            grid on
            hold off
            xlabel('x/L', 'FontSize', 18)
            ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 18)
            ylim([plot_ymin,  plot_ymax])
            xlim([0.1 0.9])
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;
            title(DeltaFCAx, ['T = ',num2str(temps(temp_index)), '°C Replicates'], 'FontSize', 20)
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'FlippedStainReplicateComparisonAPProfiles', filesep, 'NormalizedBootstrappedTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'FlippedStainReplicateComparisonAPProfiles', filesep, 'NormalizedBootstrappedTestProfiles'])
            end
            outpath2 = [AllSetsProfFigPath, filesep, 'FlippedStainReplicateComparisonAPProfiles', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet_',num2str(master_index),'_T', strrep(num2str(temps(temp_index)), '.', '_'), 'C.png' ];
            saveas(DeltaFCFixCorrectedFig,outpath2);
        end
        
    end
end
%% Make individual time plots with background subtracted normalized bootstrapped profiles

times_to_plot = 10:5:50;
colorsDubuisTimes = (hsv(round(length(times_to_plot)*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:length(times_to_plot),:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
TempTickLabels = {};
sorted_temps = [17.5, 20, 22.5, 25, 27.5];
for temp_index = 1:length(sorted_temps)
    TempTickLabels{temp_index} = num2str(round(sorted_temps(temp_index),1));
end
for master_index = [1 2 7]
    for ch_index = [3 5]
        
        exp_indices = cell(1, length(temps));
        exp_index = zeros(1, length(temps));
        exp_index2 = zeros(1, length(temps));
        x_sample = cell(3, length(temps));
        UseTF = cell(3, length(temps));
        ys = cell(3, length(temps));
        BadTF = cell(3, length(temps));
        xfits = cell(3, length(temps));
        yfits = cell(3, length(temps));
        x_sample_control = cell(3, length(temps));
        UseTF_control = cell(3, length(temps));
        ys_control = cell(3, length(temps));
        BadTF_control = cell(3, length(temps));
        xfits_control = cell(3, length(temps));
        yfits_control = cell(3, length(temps));
        SetLabel =cell(3, length(temps));
        for temp_index = 1:length(temps)
            exp_indices{temp_index} = find(AllSetInfo.Temperatures == sorted_temps(temp_index));
            exp_index(temp_index) =  exp_indices{temp_index}(1);
            exp_index2(temp_index) =  exp_indices{temp_index}(2);
            exp_index3(temp_index) =  exp_indices{temp_index}(3);
            
            
            CompiledEmbryos = AllCompiledEmbryos{exp_index(temp_index)};
            x_sample{1, temp_index} = CompiledEmbryos.DubuisEmbryoTimes;
            UseTF{1, temp_index} = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample{1, temp_index});
            x_sample{1, temp_index} =x_sample{1, temp_index}(UseTF{1, temp_index});
            ys{1, temp_index} =  CompiledEmbryos.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{1, temp_index},:,ch_index, master_index);
            BadTF{1, temp_index} =  sum(isnan(ys{1, temp_index}), 2).'> NumAPbins;
            x_sample{1, temp_index} = x_sample{1, temp_index}(~BadTF{1, temp_index});
            ys{1, temp_index} = ys{1, temp_index}(~BadTF{1, temp_index},:);
            xfits{1, temp_index} = CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{1, temp_index} = CompiledEmbryos.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            SetLabel{1, temp_index} = AllSetInfo.SetLabels{exp_index(temp_index)};
            
            CompiledEmbryos2 = AllCompiledEmbryos{exp_index2(temp_index)};
            x_sample{2, temp_index} = CompiledEmbryos2.DubuisEmbryoTimes;
            UseTF{2, temp_index} = CompiledEmbryos2.TestSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample{2, temp_index});
            x_sample{2, temp_index} =x_sample{2, temp_index}(UseTF{2, temp_index});
            ys{2, temp_index} =  CompiledEmbryos2.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{2, temp_index},:,ch_index, master_index);
            BadTF{2, temp_index} =  sum(isnan(ys{2, temp_index}), 2).'> NumAPbins;
            x_sample{2, temp_index} = x_sample{2, temp_index}(~BadTF{2, temp_index});
            ys{2, temp_index} = ys{2, temp_index}(~BadTF{2, temp_index},:);
            xfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{2, temp_index} = CompiledEmbryos2.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            SetLabel{2, temp_index} = AllSetInfo.SetLabels{exp_index2(temp_index)};
            
            CompiledEmbryos3 = AllCompiledEmbryos{exp_index3(temp_index)};
            x_sample{3, temp_index} = CompiledEmbryos3.DubuisEmbryoTimes;
            UseTF{3, temp_index} = CompiledEmbryos3.TestSetEmbryos & CompiledEmbryos3.IsNC14 & ~isnan(x_sample{3, temp_index});
            x_sample{3, temp_index} =x_sample{3, temp_index}(UseTF{3, temp_index});
            ys{3, temp_index} =  CompiledEmbryos3.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{3, temp_index},:,ch_index, master_index);
            BadTF{3, temp_index} =  sum(isnan(ys{3, temp_index}), 2).'> NumAPbins;
            x_sample{3, temp_index} = x_sample{3, temp_index}(~BadTF{3, temp_index});
            ys{3, temp_index} = ys{3, temp_index}(~BadTF{3, temp_index},:);
            xfits{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{3, temp_index} = CompiledEmbryos3.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            SetLabel{3, temp_index} = AllSetInfo.SetLabels{exp_index3(temp_index)};
            
        end
        for temp_index = 1:length(temps)
            ymaxes(temp_index) = max([max( yfits{1, temp_index}(:)), max( yfits{2, temp_index}(:)), max( yfits{3, temp_index}(:))]);
            ymins(temp_index) = min([min( yfits{1, temp_index}(:)), min( yfits{2, temp_index}(:)), min( yfits{2, temp_index}(:))]);
        end
        plot_ymax = ceil(max(ymaxes)/0.2)*0.2;
        plot_ymin = floor(min(ymins)/0.2)*0.2;
        for time_index = 1:length(times_to_plot)
            
            close all
            DeltaFCFixCorrectedFig =figure(1);
            set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .5, .6]);
            set(gcf,'color','w');
            DeltaFCAx = axes(DeltaFCFixCorrectedFig);
            map = colormap(colorsTemperatures);
            h = colorbar;
            % %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
            
            colorTitleHandle = get(h,'Title');
            titleString = '\del (\mum)';
            titleString = 'T (ºC)';
            set(colorTitleHandle ,'String',titleString);
            h.Ticks =  FractionalTemperatures; %Create 8 ticks from zero to 1
            
            h.TickLabels = TempTickLabels;
            hold on
            
            for temp_index = 1:length(sorted_temps)
                for rep_index = 1:3
                    timebin = find(round(xfits{rep_index, temp_index}) == times_to_plot(time_index));
                    if ~isempty(timebin)
                        plot(APbins, yfits{rep_index, temp_index}(timebin,:),SetLineStyles{rep_index},...
                            'LineWidth', 2.0, 'Color', colorsTemperatures(temp_index, :));
                        hold on
                    end
                end
                
            end
            grid on
            hold off
            xlabel('x/L', 'FontSize', 18)
            ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 18)
            ylim([plot_ymin,  plot_ymax])
            xlim([0.1 0.9])
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;
            title(DeltaFCAx, [num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 20)
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'FlippedTemperatureComparisonAPProfiles', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'FlippedTemperatureComparisonAPProfiles', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
            end
            outpath2 = [AllSetsProfFigPath, filesep, 'FlippedTemperatureComparisonAPProfiles', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
            saveas(DeltaFCFixCorrectedFig,outpath2);
        end
        
    end
end

%% Make individual time plots with normalized bootstrapped profiles

times_to_plot = 10:5:50;
colorsDubuisTimes = (hsv(round(length(times_to_plot)*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:length(times_to_plot),:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
TempTickLabels = {};
sorted_temps = [17.5, 20, 22.5, 25, 27.5];
for temp_index = 1:length(sorted_temps)
    TempTickLabels{temp_index} = num2str(round(sorted_temps(temp_index),1));
end
for master_index = [1 2 7]
    for ch_index = [3 5]
        
        exp_indices = cell(1, length(temps));
        exp_index = zeros(1, length(temps));
        exp_index2 = zeros(1, length(temps));
        x_sample = cell(3, length(temps));
        UseTF = cell(3, length(temps));
        ys = cell(3, length(temps));
        BadTF = cell(3, length(temps));
        xfits = cell(3, length(temps));
        yfits = cell(3, length(temps));
        x_sample_control = cell(3, length(temps));
        UseTF_control = cell(3, length(temps));
        ys_control = cell(3, length(temps));
        BadTF_control = cell(3, length(temps));
        xfits_control = cell(3, length(temps));
        yfits_control = cell(3, length(temps));
        SetLabel =cell(3, length(temps));
        for temp_index = 1:length(temps)
            exp_indices{temp_index} = find(AllSetInfo.Temperatures == sorted_temps(temp_index));
            exp_index(temp_index) =  exp_indices{temp_index}(1);
            exp_index2(temp_index) =  exp_indices{temp_index}(2);
            exp_index3(temp_index) =  exp_indices{temp_index}(3);
            
            
            CompiledEmbryos = AllCompiledEmbryos{exp_index(temp_index)};
            x_sample{1, temp_index} = CompiledEmbryos.DubuisEmbryoTimes;
            UseTF{1, temp_index} = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample{1, temp_index});
            x_sample{1, temp_index} =x_sample{1, temp_index}(UseTF{1, temp_index});
            ys{1, temp_index} =  CompiledEmbryos.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{1, temp_index},:,ch_index, master_index);
            BadTF{1, temp_index} =  sum(isnan(ys{1, temp_index}), 2).'> NumAPbins;
            x_sample{1, temp_index} = x_sample{1, temp_index}(~BadTF{1, temp_index});
            ys{1, temp_index} = ys{1, temp_index}(~BadTF{1, temp_index},:);
            xfits{1, temp_index} = CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{1, temp_index} = CompiledEmbryos.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            SetLabel{1, temp_index} = AllSetInfo.SetLabels{exp_index(temp_index)};
            
            CompiledEmbryos2 = AllCompiledEmbryos{exp_index2(temp_index)};
            x_sample{2, temp_index} = CompiledEmbryos2.DubuisEmbryoTimes;
            UseTF{2, temp_index} = CompiledEmbryos2.TestSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample{2, temp_index});
            x_sample{2, temp_index} =x_sample{2, temp_index}(UseTF{2, temp_index});
            ys{2, temp_index} =  CompiledEmbryos2.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{2, temp_index},:,ch_index, master_index);
            BadTF{2, temp_index} =  sum(isnan(ys{2, temp_index}), 2).'> NumAPbins;
            x_sample{2, temp_index} = x_sample{2, temp_index}(~BadTF{2, temp_index});
            ys{2, temp_index} = ys{2, temp_index}(~BadTF{2, temp_index},:);
            xfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{2, temp_index} = CompiledEmbryos2.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            SetLabel{2, temp_index} = AllSetInfo.SetLabels{exp_index2(temp_index)};
            
            CompiledEmbryos3 = AllCompiledEmbryos{exp_index3(temp_index)};
            x_sample{3, temp_index} = CompiledEmbryos3.DubuisEmbryoTimes;
            UseTF{3, temp_index} = CompiledEmbryos3.TestSetEmbryos & CompiledEmbryos3.IsNC14 & ~isnan(x_sample{3, temp_index});
            x_sample{3, temp_index} =x_sample{3, temp_index}(UseTF{3, temp_index});
            ys{3, temp_index} =  CompiledEmbryos3.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{3, temp_index},:,ch_index, master_index);
            BadTF{3, temp_index} =  sum(isnan(ys{3, temp_index}), 2).'> NumAPbins;
            x_sample{3, temp_index} = x_sample{3, temp_index}(~BadTF{3, temp_index});
            ys{3, temp_index} = ys{3, temp_index}(~BadTF{3, temp_index},:);
            xfits{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{3, temp_index} = CompiledEmbryos3.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            SetLabel{3, temp_index} = AllSetInfo.SetLabels{exp_index3(temp_index)};
            
        end
        for temp_index = 1:length(temps)
            ymaxes(temp_index) = max([max( yfits{1, temp_index}(:)), max( yfits{2, temp_index}(:)), max( yfits{3, temp_index}(:))]);
            ymins(temp_index) = min([min( yfits{1, temp_index}(:)), min( yfits{2, temp_index}(:)), min( yfits{2, temp_index}(:))]);
        end
        plot_ymax = ceil(max(ymaxes)/0.2)*0.2;
        plot_ymin = floor(min(ymins)/0.2)*0.2;
        for time_index = 1:length(times_to_plot)
            
            close all
            DeltaFCFixCorrectedFig =figure(1);
            set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .5, .6]);
            set(gcf,'color','w');
            DeltaFCAx = axes(DeltaFCFixCorrectedFig);
            map = colormap(colorsTemperatures);
            h = colorbar;
            % %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
            
            colorTitleHandle = get(h,'Title');
            titleString = '\del (\mum)';
            titleString = 'T (ºC)';
            set(colorTitleHandle ,'String',titleString);
            h.Ticks =  FractionalTemperatures; %Create 8 ticks from zero to 1
            
            h.TickLabels = TempTickLabels;
            hold on
            
            for temp_index = 1:length(sorted_temps)
                for rep_index = 1:3
                    timebin = find(round(xfits{rep_index, temp_index}) == times_to_plot(time_index));
                    if ~isempty(timebin)
                        plot(APbins, yfits{rep_index, temp_index}(timebin,:),SetLineStyles{rep_index},...
                            'LineWidth', 2.0, 'Color', colorsTemperatures(temp_index, :));
                        hold on
                    end
                end
                
            end
            grid on
            hold off
            xlabel('x/L', 'FontSize', 18)
            ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 18)
            ylim([plot_ymin,  plot_ymax])
            xlim([0.1 0.9])
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;
            title(DeltaFCAx, [num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 20)
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'FlippedTemperatureComparisonAPProfiles', filesep, 'NormalizedBootstrappedTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'FlippedTemperatureComparisonAPProfiles', filesep, 'NormalizedBootstrappedTestProfiles'])
            end
            outpath2 = [AllSetsProfFigPath, filesep, 'FlippedTemperatureComparisonAPProfiles', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
            saveas(DeltaFCFixCorrectedFig,outpath2);
        end
        
    end
end
%% Make individual temperature and time plots with background subtracted normalized bootstrapped profiles

times_to_plot = 10:5:50;
colorsDubuisTimes = (hsv(round(length(times_to_plot)*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:length(times_to_plot),:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
TempTickLabels = {};
sorted_temps = [17.5, 20, 22.5, 25, 27.5];
for temp_index = 1:length(sorted_temps)
    TempTickLabels{temp_index} = num2str(round(sorted_temps(temp_index),1));
end
for master_index = [1 2 7]
    for ch_index = [3 5]
        
        exp_indices = cell(1, length(temps));
        exp_index = zeros(1, length(temps));
        exp_index2 = zeros(1, length(temps));
        x_sample = cell(3, length(temps));
        UseTF = cell(3, length(temps));
        ys = cell(3, length(temps));
        BadTF = cell(3, length(temps));
        xfits = cell(3, length(temps));
        yfits = cell(3, length(temps));
        x_sample_control = cell(3, length(temps));
        UseTF_control = cell(3, length(temps));
        ys_control = cell(3, length(temps));
        BadTF_control = cell(3, length(temps));
        xfits_control = cell(3, length(temps));
        yfits_control = cell(3, length(temps));
        SetLabel =cell(3, length(temps));
        for temp_index = 1:length(temps)
            exp_indices{temp_index} = find(AllSetInfo.Temperatures == sorted_temps(temp_index));
            exp_index(temp_index) =  exp_indices{temp_index}(1);
            exp_index2(temp_index) =  exp_indices{temp_index}(2);
            exp_index3(temp_index) =  exp_indices{temp_index}(3);
            
            
            CompiledEmbryos = AllCompiledEmbryos{exp_index(temp_index)};
            x_sample{1, temp_index} = CompiledEmbryos.DubuisEmbryoTimes;
            UseTF{1, temp_index} = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample{1, temp_index});
            x_sample{1, temp_index} =x_sample{1, temp_index}(UseTF{1, temp_index});
            ys{1, temp_index} =  CompiledEmbryos.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{1, temp_index},:,ch_index, master_index);
            BadTF{1, temp_index} =  sum(isnan(ys{1, temp_index}), 2).'> NumAPbins;
            x_sample{1, temp_index} = x_sample{1, temp_index}(~BadTF{1, temp_index});
            ys{1, temp_index} = ys{1, temp_index}(~BadTF{1, temp_index},:);
            xfits{1, temp_index} = CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{1, temp_index} = CompiledEmbryos.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            SetLabel{1, temp_index} = AllSetInfo.SetLabels{exp_index(temp_index)};
            
            CompiledEmbryos2 = AllCompiledEmbryos{exp_index2(temp_index)};
            x_sample{2, temp_index} = CompiledEmbryos2.DubuisEmbryoTimes;
            UseTF{2, temp_index} = CompiledEmbryos2.TestSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample{2, temp_index});
            x_sample{2, temp_index} =x_sample{2, temp_index}(UseTF{2, temp_index});
            ys{2, temp_index} =  CompiledEmbryos2.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{2, temp_index},:,ch_index, master_index);
            BadTF{2, temp_index} =  sum(isnan(ys{2, temp_index}), 2).'> NumAPbins;
            x_sample{2, temp_index} = x_sample{2, temp_index}(~BadTF{2, temp_index});
            ys{2, temp_index} = ys{2, temp_index}(~BadTF{2, temp_index},:);
            xfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{2, temp_index} = CompiledEmbryos2.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            SetLabel{2, temp_index} = AllSetInfo.SetLabels{exp_index2(temp_index)};
            
            CompiledEmbryos3 = AllCompiledEmbryos{exp_index3(temp_index)};
            x_sample{3, temp_index} = CompiledEmbryos3.DubuisEmbryoTimes;
            UseTF{3, temp_index} = CompiledEmbryos3.TestSetEmbryos & CompiledEmbryos3.IsNC14 & ~isnan(x_sample{3, temp_index});
            x_sample{3, temp_index} =x_sample{3, temp_index}(UseTF{3, temp_index});
            ys{3, temp_index} =  CompiledEmbryos3.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{3, temp_index},:,ch_index, master_index);
            BadTF{3, temp_index} =  sum(isnan(ys{3, temp_index}), 2).'> NumAPbins;
            x_sample{3, temp_index} = x_sample{3, temp_index}(~BadTF{3, temp_index});
            ys{3, temp_index} = ys{3, temp_index}(~BadTF{3, temp_index},:);
            xfits{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{3, temp_index} = CompiledEmbryos3.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            SetLabel{3, temp_index} = AllSetInfo.SetLabels{exp_index3(temp_index)};
            
        end
        for temp_index = 1:length(temps)
            ymaxes(temp_index) = max([max( yfits{1, temp_index}(:)), max( yfits{2, temp_index}(:)), max( yfits{3, temp_index}(:))]);
            ymins(temp_index) = min([min( yfits{1, temp_index}(:)), min( yfits{2, temp_index}(:)), min( yfits{2, temp_index}(:))]);
        end
        plot_ymax = ceil(max(ymaxes)/0.2)*0.2;
        plot_ymin = floor(min(ymins)/0.2)*0.2;
        for time_index = 1:length(times_to_plot)
            
            
            
            for temp_index = 1:length(sorted_temps)
                close all
                DeltaFCFixCorrectedFig =figure(1);
                set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .5, .6]);
                set(gcf,'color','w');
                DeltaFCAx = axes(DeltaFCFixCorrectedFig);
                map = colormap(colorsTemperatures);
                hold on
                color_index  = [1, 4, 5];
                for rep_index = 1:3
                    
                    timebin = find(round(xfits{rep_index, temp_index}) == times_to_plot(time_index));
                    if ~isempty(timebin)
                        plot(APbins, yfits{rep_index, temp_index}(timebin,:),SetLineStyles{rep_index},...
                            'LineWidth', 3.0, 'Color', colorsTemperatures(color_index(rep_index), :));
                        hold on
                    end
                end
                grid on
                hold off
                xlabel('x/L', 'FontSize', 18)
                ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 18)
                legend({'Replicate 1', 'Replicate 2', 'Flipped Stain'}, 'FontSize', 18);
                ylim([plot_ymin,  plot_ymax])
                xlim([0.1 0.9])
                DeltaFCAx.YAxis.FontSize = 18;
                DeltaFCAx.XAxis.FontSize = 18;
                title(DeltaFCAx, ['T = ',num2str(sorted_temps(temp_index)),'ºC, ', num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 20)
                DeltaFCAx.FontSize = 18;
                if ~isdir([AllSetsProfFigPath, filesep, 'FlippedSingleTemperatureComparisonAPProfiles', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
                    mkdir([AllSetsProfFigPath, filesep, 'FlippedSingleTemperatureComparisonAPProfiles', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
                end
                outpath2 = [AllSetsProfFigPath, filesep, 'FlippedSingleTemperatureComparisonAPProfiles', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
                    ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_Time', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm_T',strrep(num2str(sorted_temps(temp_index)), '.', '_'),'C.png' ];
                saveas(DeltaFCFixCorrectedFig,outpath2);
                
            end
            
        end
        
    end
end

%% Make individual temperature and time plots with background subtracted normalized bootstrapped profiles

times_to_plot = 10:5:50;
colorsDubuisTimes = (hsv(round(length(times_to_plot)*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:length(times_to_plot),:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
TempTickLabels = {};
sorted_temps = [17.5, 20, 22.5, 25, 27.5];
for temp_index = 1:length(sorted_temps)
    TempTickLabels{temp_index} = num2str(round(sorted_temps(temp_index),1));
end
for master_index = [1 2 7]
    for ch_index = [3 5]
        
        exp_indices = cell(1, length(temps));
        exp_index = zeros(1, length(temps));
        exp_index2 = zeros(1, length(temps));
        x_sample = cell(3, length(temps));
        UseTF = cell(3, length(temps));
        ys = cell(3, length(temps));
        BadTF = cell(3, length(temps));
        xfits = cell(3, length(temps));
        yfits = cell(3, length(temps));
        x_sample_control = cell(3, length(temps));
        UseTF_control = cell(3, length(temps));
        ys_control = cell(3, length(temps));
        BadTF_control = cell(3, length(temps));
        xfits_control = cell(3, length(temps));
        yfits_control = cell(3, length(temps));
        SetLabel =cell(3, length(temps));
        for temp_index = 1:length(temps)
            exp_indices{temp_index} = find(AllSetInfo.Temperatures == sorted_temps(temp_index));
            exp_index(temp_index) =  exp_indices{temp_index}(1);
            exp_index2(temp_index) =  exp_indices{temp_index}(2);
            exp_index3(temp_index) =  exp_indices{temp_index}(3);
            
            
            CompiledEmbryos = AllCompiledEmbryos{exp_index(temp_index)};
            x_sample{1, temp_index} = CompiledEmbryos.DubuisEmbryoTimes;
            UseTF{1, temp_index} = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample{1, temp_index});
            x_sample{1, temp_index} =x_sample{1, temp_index}(UseTF{1, temp_index});
            ys{1, temp_index} =  CompiledEmbryos.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{1, temp_index},:,ch_index, master_index);
            BadTF{1, temp_index} =  sum(isnan(ys{1, temp_index}), 2).'> NumAPbins;
            x_sample{1, temp_index} = x_sample{1, temp_index}(~BadTF{1, temp_index});
            ys{1, temp_index} = ys{1, temp_index}(~BadTF{1, temp_index},:);
            xfits{1, temp_index} = CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{1, temp_index} = CompiledEmbryos.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            SetLabel{1, temp_index} = AllSetInfo.SetLabels{exp_index(temp_index)};
            
            CompiledEmbryos2 = AllCompiledEmbryos{exp_index2(temp_index)};
            x_sample{2, temp_index} = CompiledEmbryos2.DubuisEmbryoTimes;
            UseTF{2, temp_index} = CompiledEmbryos2.TestSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample{2, temp_index});
            x_sample{2, temp_index} =x_sample{2, temp_index}(UseTF{2, temp_index});
            ys{2, temp_index} =  CompiledEmbryos2.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{2, temp_index},:,ch_index, master_index);
            BadTF{2, temp_index} =  sum(isnan(ys{2, temp_index}), 2).'> NumAPbins;
            x_sample{2, temp_index} = x_sample{2, temp_index}(~BadTF{2, temp_index});
            ys{2, temp_index} = ys{2, temp_index}(~BadTF{2, temp_index},:);
            xfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{2, temp_index} = CompiledEmbryos2.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            SetLabel{2, temp_index} = AllSetInfo.SetLabels{exp_index2(temp_index)};
            
            CompiledEmbryos3 = AllCompiledEmbryos{exp_index3(temp_index)};
            x_sample{3, temp_index} = CompiledEmbryos3.DubuisEmbryoTimes;
            UseTF{3, temp_index} = CompiledEmbryos3.TestSetEmbryos & CompiledEmbryos3.IsNC14 & ~isnan(x_sample{3, temp_index});
            x_sample{3, temp_index} =x_sample{3, temp_index}(UseTF{3, temp_index});
            ys{3, temp_index} =  CompiledEmbryos3.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{3, temp_index},:,ch_index, master_index);
            BadTF{3, temp_index} =  sum(isnan(ys{3, temp_index}), 2).'> NumAPbins;
            x_sample{3, temp_index} = x_sample{3, temp_index}(~BadTF{3, temp_index});
            ys{3, temp_index} = ys{3, temp_index}(~BadTF{3, temp_index},:);
            xfits{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{3, temp_index} = CompiledEmbryos3.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            SetLabel{3, temp_index} = AllSetInfo.SetLabels{exp_index3(temp_index)};
            
        end
        for temp_index = 1:length(temps)
            ymaxes(temp_index) = max([max( yfits{1, temp_index}(:)), max( yfits{2, temp_index}(:)), max( yfits{3, temp_index}(:))]);
            ymins(temp_index) = min([min( yfits{1, temp_index}(:)), min( yfits{2, temp_index}(:)), min( yfits{2, temp_index}(:))]);
        end
        plot_ymax = ceil(max(ymaxes)/0.2)*0.2;
        plot_ymin = floor(min(ymins)/0.2)*0.2;
        for time_index = 1:length(times_to_plot)
            
            
            
            for temp_index = 1:length(sorted_temps)
                close all
                DeltaFCFixCorrectedFig =figure(1);
                set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .5, .6]);
                set(gcf,'color','w');
                DeltaFCAx = axes(DeltaFCFixCorrectedFig);
                map = colormap(colorsTemperatures);
                hold on
                color_index  = [1, 4, 5];
                for rep_index = 1:3
                    
                    timebin = find(round(xfits{rep_index, temp_index}) == times_to_plot(time_index));
                    if ~isempty(timebin)
                        plot(APbins, yfits{rep_index, temp_index}(timebin,:),SetLineStyles{rep_index},...
                            'LineWidth', 3.0, 'Color', colorsTemperatures(color_index(rep_index), :));
                        hold on
                    end
                end
                grid on
                hold off
                xlabel('x/L', 'FontSize', 18)
                ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 18)
                legend({'Replicate 1', 'Replicate 2', 'Flipped Stain'}, 'FontSize', 18);
                ylim([plot_ymin,  plot_ymax])
                xlim([0.1 0.9])
                DeltaFCAx.YAxis.FontSize = 18;
                DeltaFCAx.XAxis.FontSize = 18;
                title(DeltaFCAx, ['T = ',num2str(sorted_temps(temp_index)),'ºC, ', num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 20)
                DeltaFCAx.FontSize = 18;
                if ~isdir([AllSetsProfFigPath, filesep, 'FlippedSingleTemperatureComparisonAPProfiles', filesep, 'NormalizedBootstrappedTestProfiles'])
                    mkdir([AllSetsProfFigPath, filesep, 'FlippedSingleTemperatureComparisonAPProfiles', filesep, 'NormalizedBootstrappedTestProfiles'])
                end
                outpath2 = [AllSetsProfFigPath, filesep, 'FlippedSingleTemperatureComparisonAPProfiles', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
                    ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_Time', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm_T',strrep(num2str(sorted_temps(temp_index)), '.', '_'),'C.png' ];
                saveas(DeltaFCFixCorrectedFig,outpath2);
                
            end
            
        end
        
    end
end