function CheungMaFigure1C(AllCompiledEmbryos)
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


temps = [17.5 20 22.5 25 27.5];

%% Make individual temperature plots with background subtracted normalized individual embryo profiles

times_to_plot = 10:5:50;
colorsDubuisTimes = (hsv(round(length(times_to_plot)*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:length(times_to_plot),:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
RepMarkerStyles = {'o', 's'};
TempTickLabels = {};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end
LegendLabels = {};
for temp_index = 1:length(temps)
    LegendLabels{temp_index} = [num2str(round(temps(temp_index),1)), 'ºC'];
end
times_to_plot = 0:15:60;
time_diff = times_to_plot(2)-times_to_plot(1);
for master_index = 7
    for ch_index = [3 5]
        exp_indices = cell(1, length(temps));
        exp_index = zeros(1, length(temps));
        exp_index2 = zeros(1, length(temps));
        exp_index3 = zeros(1, length(temps));
        x_sample = cell(3, length(temps));
        UseTF = cell(3, length(temps));
        ys = cell(3, length(temps));
        BadTF = cell(3, length(temps));
        xfits =cell(3, length(temps));
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
            BadTF{2, temp_index} =  sum(isnan(ys{2, temp_index}), 2).' > NumAPbins;
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
            BadTF{3, temp_index} =  sum(isnan(ys{3, temp_index}), 2).' > NumAPbins;
            x_sample{3, temp_index} = x_sample{3, temp_index}(~BadTF{3, temp_index});
            ys{3, temp_index} = ys{3, temp_index}(~BadTF{3, temp_index},:);
            xfits{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{3, temp_index} = CompiledEmbryos3.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            SetLabel{3, temp_index} = AllSetInfo.SetLabels{exp_index2(temp_index)};
            
        end
        for temp_index = 1:length(temps)
            ymaxes(temp_index) = max([max( ys{1, temp_index}(:)), max( ys{2, temp_index}(:)), max( ys{3, temp_index}(:))]);
            ymins(temp_index) = min([min( ys{1, temp_index}(:)), min( ys{2, temp_index}(:)), min( ys{3, temp_index}(:))]);
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
%             h = colorbar;
%             % %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
%             
%             colorTitleHandle = get(h,'Title');
%             titleString = '\del (\mum)';
%             titleString = 'T (ºC)';
%             set(colorTitleHandle ,'String',titleString);
%             h.Ticks =  FractionalTemperatures; %Create 8 ticks from zero to 1
%             
%             h.TickLabels = TempTickLabels;
%             hold on
%             
            plot_handles = zeros(1, 15);
            
            for temp_index = 1:length(temps)
                if temps(temp_index) == 17.5
                    rep_list = 1;
                else
                    rep_list = 1:2;
                end
                for rep_index = rep_list
                    timebins = find((x_sample{rep_index, temp_index} < times_to_plot(time_index) +time_diff/2)& ...
                        (x_sample{rep_index, temp_index} >= times_to_plot(time_index)-time_diff/2));
                    if ~isempty(timebins)
                        for  timebin = timebins
 
                            plot_handles((temp_index-1)*3+rep_index) = scatter(APbins-0.01+0.003*temp_index, ys{rep_index, temp_index}(timebin,:),30,...
                                'o', 'MarkerFaceColor', colorsTemperatures(temp_index, :),...
                                'MarkerEdgeColor', colorsTemperatures(temp_index, :),'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
                            
                            hold on
                        end
                    end
                end
                
            end
            grid on
            hold off
            xlabel('x/L', 'FontSize', 18)
            ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 18)
            ylim([plot_ymin,  plot_ymax])
            xlim([0.075 0.925])
            subplot_handles = plot_handles(1:3:13);
            SubLegendLabels = LegendLabels(subplot_handles ~= 0);
            subplot_handles = subplot_handles(subplot_handles ~= 0);
            legend(subplot_handles,SubLegendLabels, 'FontSize', 18);
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;
            title(DeltaFCAx, [num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 20)
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1C', filesep, 'BgdSubNormalizedSingleEmbryoTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1C', filesep, 'BgdSubNormalizedSingleEmbryoTestProfiles'])
            end
            outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFigure1C', filesep, 'BgdSubNormalizedSingleEmbryoTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_WidestTime', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
            saveas(DeltaFCFixCorrectedFig,outpath2);
        end
        
    end
end

%% Make individual temperature plots with normalized individual embryo profiles

times_to_plot = 10:5:50;
colorsDubuisTimes = (hsv(round(length(times_to_plot)*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:length(times_to_plot),:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
RepMarkerStyles = {'o', 's'};
TempTickLabels = {};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end
LegendLabels = {};
for temp_index = 1:length(temps)
    LegendLabels{temp_index} = [num2str(round(temps(temp_index),1)), 'ºC'];
end
times_to_plot = 0:15:60;
time_diff = times_to_plot(2)-times_to_plot(1);
for master_index = 7
    for ch_index = [3 5]
        exp_indices = cell(1, length(temps));
        exp_index = zeros(1, length(temps));
        exp_index2 = zeros(1, length(temps));
        exp_index3 = zeros(1, length(temps));
        x_sample = cell(3, length(temps));
        UseTF = cell(3, length(temps));
        ys = cell(3, length(temps));
        BadTF = cell(3, length(temps));
        xfits =cell(3, length(temps));
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
            BadTF{2, temp_index} =  sum(isnan(ys{2, temp_index}), 2).' > NumAPbins;
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
            BadTF{3, temp_index} =  sum(isnan(ys{3, temp_index}), 2).' > NumAPbins;
            x_sample{3, temp_index} = x_sample{3, temp_index}(~BadTF{3, temp_index});
            ys{3, temp_index} = ys{3, temp_index}(~BadTF{3, temp_index},:);
            xfits{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{3, temp_index} = CompiledEmbryos3.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
            SetLabel{3, temp_index} = AllSetInfo.SetLabels{exp_index2(temp_index)};
            
        end
        for temp_index = 1:length(temps)
            ymaxes(temp_index) = max([max( ys{1, temp_index}(:)), max( ys{2, temp_index}(:)), max( ys{3, temp_index}(:))]);
            ymins(temp_index) = min([min( ys{1, temp_index}(:)), min( ys{2, temp_index}(:)), min( ys{3, temp_index}(:))]);
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
%             h = colorbar;
%             % %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
%             
%             colorTitleHandle = get(h,'Title');
%             titleString = '\del (\mum)';
%             titleString = 'T (ºC)';
%             set(colorTitleHandle ,'String',titleString);
%             h.Ticks =  FractionalTemperatures; %Create 8 ticks from zero to 1
%             
%             h.TickLabels = TempTickLabels;
%             hold on
%             
            plot_handles = zeros(1, 15);
            
            for temp_index = 1:length(temps)
                if temps(temp_index) == 17.5
                    rep_list = 1;
                else
                    rep_list = 1:2;
                end
                for rep_index = rep_list
                    timebins = find((x_sample{rep_index, temp_index} < times_to_plot(time_index) + time_diff/2)& ...
                        (x_sample{rep_index, temp_index} >= times_to_plot(time_index)-time_diff/2));
                    if ~isempty(timebins)
                        for  timebin = timebins
 
                            plot_handles((temp_index-1)*3+rep_index) = scatter(APbins-0.01+0.003*temp_index, ys{rep_index, temp_index}(timebin,:),30,...
                                'o', 'MarkerFaceColor', colorsTemperatures(temp_index, :),...
                                'MarkerEdgeColor', colorsTemperatures(temp_index, :),'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
                            
                            hold on
                        end
                    end
                end
                
            end
            grid on
            hold off
            xlabel('x/L', 'FontSize', 18)
            ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 18)
            ylim([plot_ymin,  plot_ymax])
            xlim([0.075 0.925])
            subplot_handles = plot_handles(1:3:13);
            SubLegendLabels = LegendLabels(subplot_handles ~= 0);
            subplot_handles = subplot_handles(subplot_handles ~= 0);
            legend(subplot_handles,SubLegendLabels, 'FontSize', 18);
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;
            title(DeltaFCAx, [num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 20)
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1C', filesep, 'NormalizedSingleEmbryoTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1C', filesep, 'NormalizedSingleEmbryoTestProfiles'])
            end
            outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFigure1C', filesep, 'NormalizedSingleEmbryoTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_WidestTime', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
            saveas(DeltaFCFixCorrectedFig,outpath2);
        end
        
    end
end