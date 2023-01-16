function ReplicateComparisonPlots(AllCompiledEmbryos)
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
        x_sample = cell(2, length(temps));
        UseTF = cell(2, length(temps));
        ys = cell(2, length(temps));
        BadTF = cell(2, length(temps));
        xfits = cell(2, length(temps));
        yfits = cell(2, length(temps));
        x_sample_control = cell(2, length(temps));
        UseTF_control = cell(2, length(temps));
        ys_control = cell(2, length(temps));
        BadTF_control = cell(2, length(temps));
        xfits_control = cell(2, length(temps));
        yfits_control = cell(2, length(temps));
        SetLabel = cell(2, length(temps));
        for temp_index = 1:length(temps)
            exp_indices{temp_index} = find(AllSetInfo.Temperatures == temps(temp_index));
            exp_indices{temp_index}  = exp_indices{temp_index}(1:2);
            exp_index(temp_index) =  exp_indices{temp_index}(1);
            exp_index2(temp_index) =  exp_indices{temp_index}(2);
            
            
            CompiledEmbryos = AllCompiledEmbryos{exp_index(temp_index)};
            x_sample{1, temp_index} = CompiledEmbryos.DubuisEmbryoTimes;
            UseTF{1, temp_index} = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample{1, temp_index});
            x_sample{1, temp_index} =x_sample{1, temp_index}(UseTF{1, temp_index});
            ys{1, temp_index} =  CompiledEmbryos.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{1, temp_index},:,ch_index, master_index);
            BadTF{1, temp_index} =  sum(isnan(ys{1, temp_index}), 2).' > NumAPbins;
            x_sample{1, temp_index} = x_sample{1, temp_index}(~BadTF{1, temp_index});
            ys{1, temp_index} = ys{1, temp_index}(~BadTF{1, temp_index},:);
            xfits{1, temp_index} = CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{1, temp_index} = CompiledEmbryos.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            
            
            x_sample_control{1, temp_index} = CompiledEmbryos.DubuisEmbryoTimes;
            UseTF_control{1, temp_index} = CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample_control{1, temp_index});
            x_sample_control{1, temp_index} =x_sample_control{1, temp_index}(UseTF_control{1, temp_index});
            ys_control{1, temp_index} = CompiledEmbryos.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF_control{1, temp_index},:,ch_index, master_index);
            BadTF_control{1, temp_index} =  sum(isnan(ys_control{1, temp_index}), 2).' > NumAPbins;
            x_sample_control{1, temp_index} = x_sample_control{1, temp_index}(~BadTF_control{1, temp_index});
            ys_control{1, temp_index} = ys_control{1, temp_index}(~BadTF_control{1, temp_index},:);
            SetLabel{1, temp_index} = AllSetInfo.SetLabels{exp_index(temp_index)};
            
            CompiledEmbryos2 = AllCompiledEmbryos{exp_index2(temp_index)};
            x_sample{2, temp_index} = CompiledEmbryos2.DubuisEmbryoTimes;
            UseTF{2, temp_index} = CompiledEmbryos2.TestSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample{2, temp_index});
            x_sample{2, temp_index} =x_sample{2, temp_index}(UseTF{2, temp_index});
            %ys{2, temp_index} =  CompiledEmbryos2.UnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2(UseTF{2, temp_index},:,ch_index, master_index);
            ys{2, temp_index} =  CompiledEmbryos2.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{2, temp_index},:,ch_index, master_index);
            BadTF{2, temp_index} =  sum(isnan(ys{2, temp_index}), 2).' > NumAPbins;
            x_sample{2, temp_index} = x_sample{2, temp_index}(~BadTF{2, temp_index});
            ys{2, temp_index} = ys{2, temp_index}(~BadTF{2, temp_index},:);
            %xfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.x;
            %yfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.mean(:,:,ch_index, master_index);
            xfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{2, temp_index} = CompiledEmbryos2.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            
            
            x_sample_control{2, temp_index} = CompiledEmbryos2.DubuisEmbryoTimes;
            UseTF_control{2, temp_index} = CompiledEmbryos2.ControlSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample_control{2, temp_index});
            x_sample_control{2, temp_index} =x_sample_control{2, temp_index}(UseTF_control{2, temp_index});
            %ys_control{2, temp_index} =  CompiledEmbryos2.UnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2(UseTF_control{2, temp_index},:,ch_index, master_index);
            ys_control{2, temp_index} =  CompiledEmbryos2.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF_control{2, temp_index},:,ch_index, master_index);
            BadTF_control{2, temp_index}=  sum(isnan(ys_control{2, temp_index}), 2).' > NumAPbins;
            x_sample_control{2, temp_index} = x_sample_control{2, temp_index}(~BadTF_control{2, temp_index});
            ys_control{2, temp_index} = ys_control{2, temp_index}(~BadTF_control{2, temp_index},:);
            %xfits_control{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.x;
            %yfits_control{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.mean(:,:,ch_index, master_index);
            SetLabel{2, temp_index} = AllSetInfo.SetLabels{exp_index2(temp_index)};
            
        end
        for temp_index = 1:length(temps)
            ymaxes(temp_index) = max([max( yfits{1, temp_index}(:)), max( yfits{2, temp_index}(:))]);
            ymins(temp_index) = min([min( yfits{1, temp_index}(:)), min( yfits{2, temp_index}(:))]);
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
                
                for rep_index = 1:2
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
            xlabel('x/L', 'FontSize', 14)
            ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 14)
            ylim([plot_ymin,  plot_ymax])
            xlim([0.1 0.9])
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;
            title(DeltaFCAx, ['T = ',num2str(temps(temp_index)), '°C Replicates'], 'FontSize', 20)
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'ReplicateComparisonAPProfiles', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'ReplicateComparisonAPProfiles', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
            end
            outpath2 = [AllSetsProfFigPath, filesep, 'ReplicateComparisonAPProfiles', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet_',num2str(master_index),'_T', strrep(num2str(temps(temp_index)), '.', '_'), 'C.png' ];
            saveas(DeltaFCFixCorrectedFig,outpath2);
        end
        
    end
end
%% Make individual temperature plots with regular normalized bootstrapped profiles
for master_index = [1 2 7]
    for ch_index = [3 5]
        
        exp_indices = cell(1, length(temps));
        exp_index = zeros(1, length(temps));
        exp_index2 = zeros(1, length(temps));
        x_sample = cell(2, length(temps));
        UseTF = cell(2, length(temps));
        ys = cell(2, length(temps));
        BadTF = cell(2, length(temps));
        xfits = cell(2, length(temps));
        yfits = cell(2, length(temps));
        x_sample_control = cell(2, length(temps));
        UseTF_control = cell(2, length(temps));
        ys_control = cell(2, length(temps));
        BadTF_control = cell(2, length(temps));
        xfits_control = cell(2, length(temps));
        yfits_control = cell(2, length(temps));
        SetLabel = cell(2, length(temps));
        for temp_index = 1:length(temps)
            exp_indices{temp_index} = find(AllSetInfo.Temperatures == temps(temp_index));
            exp_indices{temp_index}  = exp_indices{temp_index}(1:2);
            exp_index(temp_index) =  exp_indices{temp_index}(1);
            exp_index2(temp_index) =  exp_indices{temp_index}(2);
            
            
            CompiledEmbryos = AllCompiledEmbryos{exp_index(temp_index)};
            x_sample{1, temp_index} = CompiledEmbryos.DubuisEmbryoTimes;
            UseTF{1, temp_index} = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample{1, temp_index});
            x_sample{1, temp_index} =x_sample{1, temp_index}(UseTF{1, temp_index});
            ys{1, temp_index} =  CompiledEmbryos.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{1, temp_index},:,ch_index, master_index);
            BadTF{1, temp_index} =  sum(isnan(ys{1, temp_index}), 2).' > NumAPbins;
            x_sample{1, temp_index} = x_sample{1, temp_index}(~BadTF{1, temp_index});
            ys{1, temp_index} = ys{1, temp_index}(~BadTF{1, temp_index},:);
            xfits{1, temp_index} = CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{1, temp_index} = CompiledEmbryos.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            
            
            x_sample_control{1, temp_index} = CompiledEmbryos.DubuisEmbryoTimes;
            UseTF_control{1, temp_index} = CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample_control{1, temp_index});
            x_sample_control{1, temp_index} =x_sample_control{1, temp_index}(UseTF_control{1, temp_index});
            ys_control{1, temp_index} = CompiledEmbryos.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF_control{1, temp_index},:,ch_index, master_index);
            BadTF_control{1, temp_index} =  sum(isnan(ys_control{1, temp_index}), 2).' > NumAPbins;
            x_sample_control{1, temp_index} = x_sample_control{1, temp_index}(~BadTF_control{1, temp_index});
            ys_control{1, temp_index} = ys_control{1, temp_index}(~BadTF_control{1, temp_index},:);
            SetLabel{1, temp_index} = AllSetInfo.SetLabels{exp_index(temp_index)};
            
            CompiledEmbryos2 = AllCompiledEmbryos{exp_index2(temp_index)};
            x_sample{2, temp_index} = CompiledEmbryos2.DubuisEmbryoTimes;
            UseTF{2, temp_index} = CompiledEmbryos2.TestSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample{2, temp_index});
            x_sample{2, temp_index} =x_sample{2, temp_index}(UseTF{2, temp_index});
            %ys{2, temp_index} =  CompiledEmbryos2.UnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2(UseTF{2, temp_index},:,ch_index, master_index);
            ys{2, temp_index} =  CompiledEmbryos2.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{2, temp_index},:,ch_index, master_index);
            BadTF{2, temp_index} =  sum(isnan(ys{2, temp_index}), 2).' > NumAPbins;
            x_sample{2, temp_index} = x_sample{2, temp_index}(~BadTF{2, temp_index});
            ys{2, temp_index} = ys{2, temp_index}(~BadTF{2, temp_index},:);
            %xfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.x;
            %yfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.mean(:,:,ch_index, master_index);
            xfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{2, temp_index} = CompiledEmbryos2.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            
            
            x_sample_control{2, temp_index} = CompiledEmbryos2.DubuisEmbryoTimes;
            UseTF_control{2, temp_index} = CompiledEmbryos2.ControlSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample_control{2, temp_index});
            x_sample_control{2, temp_index} =x_sample_control{2, temp_index}(UseTF_control{2, temp_index});
            %ys_control{2, temp_index} =  CompiledEmbryos2.UnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2(UseTF_control{2, temp_index},:,ch_index, master_index);
            ys_control{2, temp_index} =  CompiledEmbryos2.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF_control{2, temp_index},:,ch_index, master_index);
            BadTF_control{2, temp_index}=  sum(isnan(ys_control{2, temp_index}), 2).' > NumAPbins;
            x_sample_control{2, temp_index} = x_sample_control{2, temp_index}(~BadTF_control{2, temp_index});
            ys_control{2, temp_index} = ys_control{2, temp_index}(~BadTF_control{2, temp_index},:);
            %xfits_control{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.x;
            %yfits_control{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.mean(:,:,ch_index, master_index);
            SetLabel{2, temp_index} = AllSetInfo.SetLabels{exp_index2(temp_index)};
            
        end
        for temp_index = 1:length(temps)
            ymaxes(temp_index) = max([max( yfits{1, temp_index}(:)), max( yfits{2, temp_index}(:))]);
            ymins(temp_index) = min([min( yfits{1, temp_index}(:)), min( yfits{2, temp_index}(:))]);
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
                
                for rep_index = 1:2
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
            xlabel('x/L', 'FontSize', 14)
            ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 14)
            ylim([plot_ymin,  plot_ymax])
            xlim([0.1 0.9])
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;
            title(DeltaFCAx, ['T = ',num2str(temps(temp_index)), '°C Replicates'], 'FontSize', 20)
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'ReplicateComparisonAPProfiles', filesep, 'NormalizedBootstrappedTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'ReplicateComparisonAPProfiles', filesep, 'NormalizedBootstrappedTestProfiles'])
            end
            outpath2 = [AllSetsProfFigPath, filesep, 'ReplicateComparisonAPProfiles', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet_',num2str(master_index),'_T', strrep(num2str(temps(temp_index)), '.', '_'), 'C.png' ];
            saveas(DeltaFCFixCorrectedFig,outpath2);
        end
        
    end
end

%% Make individual temperature plots with background subtracted normalized single-embryo profiles

times_to_plot = 10:5:50;
colorsDubuisTimes = (hsv(round(length(times_to_plot)*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:length(times_to_plot),:);
FractonalDubuisTimeRange = (1:2:(2*length(times_to_plot)-1))/(2*length(times_to_plot));
RepMarkerStyles = {'o', 's'};
for master_index =  [1 2 7]
    for ch_index = [3 5]
        
        exp_indices = cell(1, length(temps));
        exp_index = zeros(1, length(temps));
        exp_index2 = zeros(1, length(temps));
        x_sample = cell(2, length(temps));
        UseTF = cell(2, length(temps));
        ys = cell(2, length(temps));
        BadTF = cell(2, length(temps));
        xfits = cell(2, length(temps));
        yfits = cell(2, length(temps));
        x_sample_control = cell(2, length(temps));
        UseTF_control = cell(2, length(temps));
        ys_control = cell(2, length(temps));
        BadTF_control = cell(2, length(temps));
        xfits_control = cell(2, length(temps));
        yfits_control = cell(2, length(temps));
        SetLabel = cell(2, length(temps));
        for temp_index = 1:length(temps)
            exp_indices{temp_index} = find(AllSetInfo.Temperatures == temps(temp_index));
            exp_indices{temp_index}  = exp_indices{temp_index}(1:2);
            exp_index(temp_index) =  exp_indices{temp_index}(1);
            exp_index2(temp_index) =  exp_indices{temp_index}(2);
            
            
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
            
            
            x_sample_control{1, temp_index} = CompiledEmbryos.DubuisEmbryoTimes;
            UseTF_control{1, temp_index} = CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample_control{1, temp_index});
            x_sample_control{1, temp_index} =x_sample_control{1, temp_index}(UseTF_control{1, temp_index});
            ys_control{1, temp_index} = CompiledEmbryos.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF_control{1, temp_index},:,ch_index, master_index);
            BadTF_control{1, temp_index} =  sum(isnan(ys_control{1, temp_index}), 2).' > NumAPbins;
            x_sample_control{1, temp_index} = x_sample_control{1, temp_index}(~BadTF_control{1, temp_index});
            ys_control{1, temp_index} = ys_control{1, temp_index}(~BadTF_control{1, temp_index},:);
            SetLabel{1, temp_index} = AllSetInfo.SetLabels{exp_index(temp_index)};
            
            CompiledEmbryos2 = AllCompiledEmbryos{exp_index2(temp_index)};
            x_sample{2, temp_index} = CompiledEmbryos2.DubuisEmbryoTimes;
            UseTF{2, temp_index} = CompiledEmbryos2.TestSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample{2, temp_index});
            x_sample{2, temp_index} =x_sample{2, temp_index}(UseTF{2, temp_index});
            %ys{2, temp_index} =  CompiledEmbryos2.UnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2(UseTF{2, temp_index},:,ch_index, master_index);
            ys{2, temp_index} =  CompiledEmbryos2.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{2, temp_index},:,ch_index, master_index);
            BadTF{2, temp_index} =  sum(isnan(ys{2, temp_index}), 2).' > NumAPbins;
            x_sample{2, temp_index} = x_sample{2, temp_index}(~BadTF{2, temp_index});
            ys{2, temp_index} = ys{2, temp_index}(~BadTF{2, temp_index},:);
            %xfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.x;
            %yfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.mean(:,:,ch_index, master_index);
            xfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{2, temp_index} = CompiledEmbryos2.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            
            
            x_sample_control{2, temp_index} = CompiledEmbryos2.DubuisEmbryoTimes;
            UseTF_control{2, temp_index} = CompiledEmbryos2.ControlSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample_control{2, temp_index});
            x_sample_control{2, temp_index} =x_sample_control{2, temp_index}(UseTF_control{2, temp_index});
            %ys_control{2, temp_index} =  CompiledEmbryos2.UnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2(UseTF_control{2, temp_index},:,ch_index, master_index);
            ys_control{2, temp_index} =  CompiledEmbryos2.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF_control{2, temp_index},:,ch_index, master_index);
            BadTF_control{2, temp_index}=  sum(isnan(ys_control{2, temp_index}), 2).' > NumAPbins;
            x_sample_control{2, temp_index} = x_sample_control{2, temp_index}(~BadTF_control{2, temp_index});
            ys_control{2, temp_index} = ys_control{2, temp_index}(~BadTF_control{2, temp_index},:);
            %xfits_control{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.x;
            %yfits_control{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.mean(:,:,ch_index, master_index);
            SetLabel{2, temp_index} = AllSetInfo.SetLabels{exp_index2(temp_index)};
            
        end
        for temp_index = 1:length(temps)
            ymaxes(temp_index) = max([max( ys{1, temp_index}(:)), max( ys{2, temp_index}(:))]);
            ymins(temp_index) = min([min( ys{1, temp_index}(:)), min( ys{2, temp_index}(:))]);
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
                
                for rep_index = 1:2
                    timebins = find((x_sample{rep_index, temp_index} <= times_to_plot(length(times_to_plot)+1-time_index) + 2.5)& ...
                        (x_sample{rep_index, temp_index} >= times_to_plot(length(times_to_plot)+1-time_index)-2.5));
                    if ~isempty(timebins)
                        for  timebin = timebins
                            if rep_index == 1
                            scatter(APbins-0.03+0.005*time_index, ys{rep_index, temp_index}(timebin,:),50,...
                                RepMarkerStyles{rep_index}, 'MarkerFaceColor', colorsDubuisTimes(length(times_to_plot)+1-time_index, :),...
                                'MarkerEdgeColor', colorsDubuisTimes(length(times_to_plot)+1-time_index, :),'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5 );
                            else
                                scatter(APbins-0.03+0.005*time_index, ys{rep_index, temp_index}(timebin,:),50,...
                                RepMarkerStyles{rep_index}, 'MarkerFaceColor', colorsDubuisTimes(length(times_to_plot)+1-time_index, :),...
                                'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
                            end
                            hold on
                        end
                    end
                end
                
            end
            grid on
            hold off
            xlabel('x/L', 'FontSize', 14)
            ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 14)
            ylim([plot_ymin,  plot_ymax])
            xlim([0.05 0.95])
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;
            title(DeltaFCAx, ['T = ',num2str(temps(temp_index)), '°C Replicates'], 'FontSize', 20)
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'ReplicateComparisonAPProfiles', filesep, 'BgdSubNormalizedSingleEmbryoTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'ReplicateComparisonAPProfiles', filesep, 'BgdSubNormalizedSingleEmbryoTestProfiles'])
            end
            outpath2 = [AllSetsProfFigPath, filesep, 'ReplicateComparisonAPProfiles', filesep, 'BgdSubNormalizedSingleEmbryoTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet_',num2str(master_index),'_T', strrep(num2str(temps(temp_index)), '.', '_'), 'C.png' ];
            saveas(DeltaFCFixCorrectedFig,outpath2);
        end
        
    end
end

%% Make individual temperature plots with standard normalized single-embryo profiles

times_to_plot = 10:5:50;
colorsDubuisTimes = (hsv(round(length(times_to_plot)*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:length(times_to_plot),:);
FractonalDubuisTimeRange = (1:2:(2*length(times_to_plot)-1))/(2*length(times_to_plot));
RepMarkerStyles = {'o', 's'};
for master_index =  [1 2 7]
    for ch_index = [3 5]
        
        exp_indices = cell(1, length(temps));
        exp_index = zeros(1, length(temps));
        exp_index2 = zeros(1, length(temps));
        x_sample = cell(2, length(temps));
        UseTF = cell(2, length(temps));
        ys = cell(2, length(temps));
        BadTF = cell(2, length(temps));
        xfits = cell(2, length(temps));
        yfits = cell(2, length(temps));
        x_sample_control = cell(2, length(temps));
        UseTF_control = cell(2, length(temps));
        ys_control = cell(2, length(temps));
        BadTF_control = cell(2, length(temps));
        xfits_control = cell(2, length(temps));
        yfits_control = cell(2, length(temps));
        SetLabel = cell(2, length(temps));
        for temp_index = 1:length(temps)
            exp_indices{temp_index} = find(AllSetInfo.Temperatures == temps(temp_index));
            exp_indices{temp_index}  = exp_indices{temp_index}(1:2);
            exp_index(temp_index) =  exp_indices{temp_index}(1);
            exp_index2(temp_index) =  exp_indices{temp_index}(2);
            
            
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
            
            
            x_sample_control{1, temp_index} = CompiledEmbryos.DubuisEmbryoTimes;
            UseTF_control{1, temp_index} = CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample_control{1, temp_index});
            x_sample_control{1, temp_index} =x_sample_control{1, temp_index}(UseTF_control{1, temp_index});
            ys_control{1, temp_index} = CompiledEmbryos.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF_control{1, temp_index},:,ch_index, master_index);
            BadTF_control{1, temp_index} =  sum(isnan(ys_control{1, temp_index}), 2).'> NumAPbins;
            x_sample_control{1, temp_index} = x_sample_control{1, temp_index}(~BadTF_control{1, temp_index});
            ys_control{1, temp_index} = ys_control{1, temp_index}(~BadTF_control{1, temp_index},:);
            SetLabel{1, temp_index} = AllSetInfo.SetLabels{exp_index(temp_index)};
            
            CompiledEmbryos2 = AllCompiledEmbryos{exp_index2(temp_index)};
            x_sample{2, temp_index} = CompiledEmbryos2.DubuisEmbryoTimes;
            UseTF{2, temp_index} = CompiledEmbryos2.TestSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample{2, temp_index});
            x_sample{2, temp_index} =x_sample{2, temp_index}(UseTF{2, temp_index});
            %ys{2, temp_index} =  CompiledEmbryos2.UnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2(UseTF{2, temp_index},:,ch_index, master_index);
            ys{2, temp_index} =  CompiledEmbryos2.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF{2, temp_index},:,ch_index, master_index);
            BadTF{2, temp_index} =  sum(isnan(ys{2, temp_index}), 2).' > NumAPbins;
            x_sample{2, temp_index} = x_sample{2, temp_index}(~BadTF{2, temp_index});
            ys{2, temp_index} = ys{2, temp_index}(~BadTF{2, temp_index},:);
            %xfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.x;
            %yfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.mean(:,:,ch_index, master_index);
            xfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{2, temp_index} = CompiledEmbryos2.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            
            
            x_sample_control{2, temp_index} = CompiledEmbryos2.DubuisEmbryoTimes;
            UseTF_control{2, temp_index} = CompiledEmbryos2.ControlSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample_control{2, temp_index});
            x_sample_control{2, temp_index} =x_sample_control{2, temp_index}(UseTF_control{2, temp_index});
            %ys_control{2, temp_index} =  CompiledEmbryos2.UnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2(UseTF_control{2, temp_index},:,ch_index, master_index);
            ys_control{2, temp_index} =  CompiledEmbryos2.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(UseTF_control{2, temp_index},:,ch_index, master_index);
            BadTF_control{2, temp_index}=  sum(isnan(ys_control{2, temp_index}), 2).' > NumAPbins;
            x_sample_control{2, temp_index} = x_sample_control{2, temp_index}(~BadTF_control{2, temp_index});
            ys_control{2, temp_index} = ys_control{2, temp_index}(~BadTF_control{2, temp_index},:);
            %xfits_control{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.x;
            %yfits_control{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.mean(:,:,ch_index, master_index);
            SetLabel{2, temp_index} = AllSetInfo.SetLabels{exp_index2(temp_index)};
            
        end
        for temp_index = 1:length(temps)
            ymaxes(temp_index) = max([max( ys{1, temp_index}(:)), max( ys{2, temp_index}(:))]);
            ymins(temp_index) = min([min( ys{1, temp_index}(:)), min( ys{2, temp_index}(:))]);
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
                
                for rep_index = 1:2
                    timebins = find((x_sample{rep_index, temp_index} <= times_to_plot(length(times_to_plot)+1-time_index) + 2.5)& ...
                        (x_sample{rep_index, temp_index} >= times_to_plot(length(times_to_plot)+1-time_index)-2.5));
                    if ~isempty(timebins)
                        for  timebin = timebins
                            if rep_index == 1
                            scatter(APbins-0.03+0.005*time_index, ys{rep_index, temp_index}(timebin,:),50,...
                                RepMarkerStyles{rep_index}, 'MarkerFaceColor', colorsDubuisTimes(length(times_to_plot)+1-time_index, :),...
                                'MarkerEdgeColor', colorsDubuisTimes(length(times_to_plot)+1-time_index, :),'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5 );
                            else
                                scatter(APbins-0.03+0.005*time_index, ys{rep_index, temp_index}(timebin,:),50,...
                                RepMarkerStyles{rep_index}, 'MarkerFaceColor', colorsDubuisTimes(length(times_to_plot)+1-time_index, :),...
                                'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
                            end
                            hold on
                        end
                    end
                end
                
            end
            grid on
            hold off
            xlabel('x/L', 'FontSize', 14)
            ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 14)
            ylim([plot_ymin,  plot_ymax])
            xlim([0.05 0.95])
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;
            title(DeltaFCAx, ['T = ',num2str(temps(temp_index)), '°C Replicates'], 'FontSize', 20)
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'ReplicateComparisonAPProfiles', filesep, 'NormalizedSingleEmbryoTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'ReplicateComparisonAPProfiles', filesep, 'NormalizedSingleEmbryoTestProfiles'])
            end
            outpath2 = [AllSetsProfFigPath, filesep, 'ReplicateComparisonAPProfiles', filesep, 'NormalizedSingleEmbryoTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet_',num2str(master_index),'_T', strrep(num2str(temps(temp_index)), '.', '_'), 'C.png' ];
            saveas(DeltaFCFixCorrectedFig,outpath2);
        end
        
    end
end