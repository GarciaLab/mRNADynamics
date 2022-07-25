function CheungMaFigure1D_StandardNormalizedVersion(AllCompiledEmbryos)
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
RepMarkerStyles = {'o', 's', 'd'};
TempTickLabels = {};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end
LegendLabels = {};
for temp_index = 1:length(temps)
    LegendLabels{temp_index} = [num2str(round(temps(temp_index),1)), 'ºC'];
end

x_sample = cell(NumSets, length(times_to_plot), 2);
ys = cell(NumSets, length(times_to_plot), 2);
MeanProfiles = cell(NumSets, length(times_to_plot), 2);
SDProfiles = cell(NumSets, length(times_to_plot), 2);
CountProfiles = cell(NumSets, length(times_to_plot), 2);
yfits = cell(NumSets, length(times_to_plot), 2);
master_index = 2;
chList = [3 5];
for exp_index = 1:15
    CompiledEmbryos = AllCompiledEmbryos{exp_index};
    
    for idx = 1:length(chList)
        ch_index = chList(idx);
        for time_index = 1:length(times_to_plot)
            x = CompiledEmbryos.DubuisEmbryoTimes;
            TF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x) & (x >= times_to_plot(time_index)-2.5) & (x < times_to_plot(time_index)+2.5);
            x = x(TF);
            y = CompiledEmbryos.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(TF,:,ch_index, master_index);
            NonTF = sum(isnan(y), 2).' > NumAPbins;
            x = x(~NonTF);
            x_sample{exp_index, time_index, idx} = x;
            y = y(~NonTF, :);
            ys{exp_index, time_index, idx} = y;
            if ~isempty(y)
                MeanProfiles{exp_index, time_index, idx} = mean(y, 1, 'omitnan');
                if size(y, 1) > 1
                    SDProfiles{exp_index, time_index, idx} = std(y, 1, 'omitnan');
                else
                    SDProfiles{exp_index, time_index, idx} = NaN(size(y));
                end
            end
            CountProfiles{exp_index, time_index, idx} = size(ys, 1);
            xfit = CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            [minval, matchidx] = min(abs(xfit - times_to_plot(time_index)));
            yf = CompiledEmbryos.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(matchidx,:,ch_index, master_index);
            yfits{exp_index, time_index, idx} = yf;
        end
    end
end
       
CrossFitMeanProfiles =  cell(NumSets, NumSets, length(times_to_plot), 2);
AllSets = 1:NumSets;
for exp_index = 1:NumSets
    for exp_index2 = AllSets(~ismember(AllSets, exp_index))
        for idx = 1:length(chList)
            for time_index = 1:length(times_to_plot)
                MeanProf = MeanProfiles{exp_index, time_index, idx};
                MeanProf2 = MeanProfiles{exp_index2, time_index, idx};
                if isempty(MeanProf) |  isempty(MeanProf2) 
                    continue
                end
                dlm = fitlm(MeanProf, MeanProf2);
                CrossFitMeanProfiles{exp_index, exp_index2, time_index, idx} = dlm;
            end
        end
    end
end

CrossFitSmoothedProfiles =  cell(NumSets, NumSets, length(times_to_plot), 2);
AllSets = 1:NumSets;
for exp_index = 1:NumSets
    for exp_index2 = AllSets(~ismember(AllSets, exp_index))
        for idx = 1:length(chList)
            for time_index = 1:length(times_to_plot)
                MeanProf = yfits{exp_index, time_index, idx};
                MeanProf2 = yfits{exp_index2, time_index, idx};
                if isempty(MeanProf) |  isempty(MeanProf2) 
                    continue
                end
                dlm = fitlm(MeanProf, MeanProf2);
                CrossFitSmoothedProfiles{exp_index, exp_index2, time_index, idx} = dlm;
            end
        end
    end
end



%%
NC13_MeanProfiles = cell(NumSets, 5, 2);
NC13_SDProfiles = cell(NumSets,5, 2);
master_index = 2;
chList = [3 5];
for exp_index = 1:15
    CompiledEmbryos = AllCompiledEmbryos{exp_index};
    
    for idx = 1:length(chList)
        ch_index = chList(idx);
        for time_index = 1:5
            x = CompiledEmbryos.DubuisEmbryoTimes;
            if time_index == 1
                y = CompiledEmbryos.NC13.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.MeanProfile{ch_index, master_index};
                y_sd = CompiledEmbryos.NC13.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.ProfileSD{ch_index, master_index};
            elseif time_index == 2
                y = CompiledEmbryos.NC13.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.MeanProfileP50{ch_index, master_index};
                y_sd = CompiledEmbryos.NC13.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.ProfileSDP50{ch_index, master_index};
            elseif time_index == 3
                y = CompiledEmbryos.NC13.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.MeanProfileP75{ch_index, master_index};
                y_sd = CompiledEmbryos.NC13.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.ProfileSDP75{ch_index, master_index};
            
            elseif time_index == 4
                y = CompiledEmbryos.NC13.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.MeanProfileP90{ch_index, master_index};
                y_sd = CompiledEmbryos.NC13.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.ProfileSDP90{ch_index, master_index};
            
            elseif time_index == 5
                y = CompiledEmbryos.NC13.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.MeanProfileP25toP75{ch_index, master_index};
                y_sd = CompiledEmbryos.NC13.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.ProfileSDP25toP75{ch_index, master_index};
            
            end

            NC13_MeanProfiles{exp_index, time_index, idx} = y;
            NC13_SDProfiles{exp_index, time_index, idx} = y_sd;
   
        end
    end
end
       
NC13_CrossFitMeanProfiles =  cell(NumSets, NumSets, length(times_to_plot), 2);
AllSets = 1:NumSets;
for exp_index = 1:NumSets
    for exp_index2 = AllSets(~ismember(AllSets, exp_index))
        for idx = 1:length(chList)
            for time_index = 1:5
                MeanProf = NC13_MeanProfiles{exp_index, time_index, idx};
                MeanProf2 = NC13_MeanProfiles{exp_index2, time_index, idx};
                if isempty(MeanProf) |  isempty(MeanProf2) 
                    continue
                end
                dlm = fitlm(MeanProf, MeanProf2);
                NC13_CrossFitMeanProfiles{exp_index, exp_index2, time_index, idx} = dlm;
            end
        end
    end
end

%%

times_to_plot = 10:5:50;
colorsDubuisTimes = (hsv(round(length(times_to_plot)*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:length(times_to_plot),:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
RepMarkerStyles = {'o', 's', 'd'};
TempTickLabels = {};
RepList = {'Replicate 1', 'Replicate 2', 'FlippedStain'};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end
LegendLabels = {};
for temp_index = 1:length(temps)
    LegendLabels{temp_index} = [num2str(round(temps(temp_index),1)), 'ºC'];
end
LegendLabels = {'25.0ºC vs. 25.0ºC', '27.5ºC vs. 25.0ºC', '22.5ºC vs. 25.0ºC', '20.0ºC vs. 25.0ºC', '17.5ºC vs. 25.0ºC'};
temp_plot_order = [4 5 3 2 1];
for idx = 1:2
    ch_index = chList(idx);
    for time_index = 1:length(times_to_plot)
        for rep_type_index = 1:3
            set_indices = (0:3:12)+rep_type_index;
            PlotMeanProfs = cell(1, 5);
            PlotMeanSDs = cell(1, 5);
            PlotSlopes = NaN(1, 5);
            PlotSlopes(1) = 1;
            PlotIntercepts = NaN(1, 5);
            PlotIntercepts(1) = 0;
            ymaxes = NaN(1, 5);
            ymins = NaN(1, 5);
            xmaxes = NaN(1, 5);
            for idx2 = 1:5
                PlotMeanProfs{idx2} = MeanProfiles{set_indices(idx2), time_index, idx};
                PlotMeanSDs{idx2} = SDProfiles{set_indices(idx2), time_index, idx};
                if ~isempty(MeanProfiles{set_indices(idx2), time_index, idx}+SDProfiles{set_indices(idx2), time_index, idx}) 
                ymaxes(idx2) = max(MeanProfiles{set_indices(idx2), time_index, idx}+SDProfiles{set_indices(idx2), time_index, idx});
                xmaxes(idx2) = max(MeanProfiles{set_indices(idx2), time_index, idx});
                ymins(idx2) = min(MeanProfiles{set_indices(idx2), time_index, idx}+SDProfiles{set_indices(idx2), time_index, idx});
                end
                if idx2 ~= 1
                    dlm = CrossFitMeanProfiles{set_indices(1), set_indices(idx2), time_index, idx};
                    if ~isempty(dlm)
                    PlotSlopes(idx2) = dlm.Coefficients.Estimate(2);
                    PlotIntercepts(idx2) = dlm.Coefficients.Estimate(1);
                    end
                end
            end
            plot_ymax = (ceil(max(ymaxes)/0.2))*0.2;
            plot_xmax = (ceil(max(xmaxes)/0.2))*0.2;
            plot_ymin = floor(min(ymins)/0.05)*0.05;
            close all
            DeltaFCFixCorrectedFig =figure(1);
            set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
            set(gcf,'color','w');
            DeltaFCAx = axes(DeltaFCFixCorrectedFig);
            map = colormap(colorsTemperatures);
            for temp_index = 1:length(temps)
                plot([plot_ymin plot_ymax].',  ([plot_ymin plot_ymax]*PlotSlopes(temp_index)+PlotIntercepts(temp_index)).',...
                    'Color', colorsTemperatures(temp_plot_order(temp_index),:),...
                    'LineWidth', 2.0);
%                 plot(PlotMeanProfs{1}, PlotMeanProfs{temp_index}, colorsTemperatures(temp_plot_order(temp_index),:),...
%                     'LineWidth', 3.0);    
                hold on 
            end
            grid on 
            plot_handles = zeros(1, length(temps));
            for temp_index = 1:length(temps)
                if ~isempty(PlotMeanProfs{1}) & ~isempty( PlotMeanProfs{temp_index})
                 plot_handles(temp_index) = errorbar(PlotMeanProfs{1}, PlotMeanProfs{temp_index},PlotMeanSDs{temp_index},...
                     'o', 'MarkerSize', 10,'Color', colorsTemperatures(temp_plot_order(temp_index), :),...
                               'LineStyle', 'none', 'MarkerFaceColor', colorsTemperatures(temp_plot_order(temp_index), :),...
                                'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
                end
                            
            end
            
            hold off
            xlabel('Intensity 25.0ºC', 'FontSize', 18)
            ylabel('Intensity 27.5/25.0/22.5/20.0/17.5ºC', 'FontSize', 18);
            ylim([0,  plot_ymax])
            xlim([0 plot_xmax])
            subplot_handles = plot_handles([2 1 3 4 5]);
            SubLegendLabels = LegendLabels([2 1 3 4 5]);
            SubLegendLabels = SubLegendLabels(subplot_handles ~= 0);
            subplot_handles = subplot_handles(subplot_handles ~= 0);
            legend(subplot_handles,SubLegendLabels, 'FontSize', 18, 'location','northwest');
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;
            if rep_type_index < 3
            title(DeltaFCAx, [ChannelNames{ch_index},' Replicate ', num2str(rep_type_index),', ', num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 20)
            else
                title(DeltaFCAx, [ChannelNames{ch_index}, ' Flipped Stain, ', num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 20)
            end
            
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedMeanTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedMeanTestProfiles'])
            end
            if rep_type_index < 3
            outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedMeanTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_StainReplicate_', num2str(rep_type_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
            else
                 outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedMeanTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_FlippedStain_', num2str(rep_type_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
            end
            saveas(DeltaFCFixCorrectedFig,outpath2);
                
        end
    end
end

 %% Use Smoothed Profiles
 for idx = 1:2
    ch_index = chList(idx);
    for time_index = 1:length(times_to_plot)
        for rep_type_index = 1:3
            set_indices = (0:3:12)+rep_type_index;
            PlotMeanProfs = cell(1, 5);
            PlotMeanSDs = cell(1, 5);
            PlotSlopes = NaN(1, 5);
            PlotSlopes(1) = 1;
            PlotIntercepts = NaN(1, 5);
            PlotIntercepts(1) = 0;
            ymaxes = NaN(1, 5);
            ymins = NaN(1, 5);
            xmaxes = NaN(1, 5);
            for idx2 = 1:5
                PlotMeanProfs{idx2} = yfits{set_indices(idx2), time_index, idx};
                if ~isempty(yfits{set_indices(idx2), time_index, idx}) 
                ymaxes(idx2) = max(yfits{set_indices(idx2), time_index, idx});
                xmaxes(idx2) = max(yfits{set_indices(idx2), time_index, idx});
                ymins(idx2) = min(yfits{set_indices(idx2), time_index, idx});
                end
                if idx2 ~= 1
                    dlm = CrossFitSmoothedProfiles{set_indices(1), set_indices(idx2), time_index, idx};
                    if ~isempty(dlm)
                    PlotSlopes(idx2) = dlm.Coefficients.Estimate(2);
                    PlotIntercepts(idx2) = dlm.Coefficients.Estimate(1);
                    end
                end
            end
            if all(isnan(ymaxes))
                continue
            end
            plot_ymax = (ceil(max(ymaxes)/0.2))*0.2;
            plot_xmax = (ceil(max(xmaxes)/0.2))*0.2;
            plot_ymin = floor(min(ymins)/0.05)*0.05;
            close all
            DeltaFCFixCorrectedFig =figure(1);
            set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
            set(gcf,'color','w');
            DeltaFCAx = axes(DeltaFCFixCorrectedFig);
            map = colormap(colorsTemperatures);
            for temp_index = 1:length(temps)
                plot([plot_ymin plot_ymax].',  ([plot_ymin plot_ymax]*PlotSlopes(temp_index)+PlotIntercepts(temp_index)).',...
                    'Color', colorsTemperatures(temp_plot_order(temp_index),:),...
                    'LineWidth', 2.0);
%                 plot(PlotMeanProfs{1}, PlotMeanProfs{temp_index}, colorsTemperatures(temp_plot_order(temp_index),:),...
%                     'LineWidth', 3.0);    
                hold on 
            end
            grid on 
            plot_handles = zeros(1, length(temps));
            for temp_index = 1:length(temps)
                if ~isempty(PlotMeanProfs{1}) & ~isempty( PlotMeanProfs{temp_index})
                 plot_handles(temp_index) = scatter(PlotMeanProfs{1}, PlotMeanProfs{temp_index},...
                    100,'MarkerFaceColor', colorsTemperatures(temp_plot_order(temp_index), :),'MarkerEdgeColor', 'k');
                end
                            
            end
            
            hold off
            xlabel('Intensity 25.0ºC', 'FontSize', 18)
            ylabel('Intensity 27.5/25.0/22.5/20.0/17.5ºC', 'FontSize', 18);
            ylim([0,  plot_ymax])
            xlim([0 plot_xmax])
            subplot_handles = plot_handles([2 1 3 4 5]);
            SubLegendLabels = LegendLabels([2 1 3 4 5]);
            SubLegendLabels = SubLegendLabels(subplot_handles ~= 0);
            subplot_handles = subplot_handles(subplot_handles ~= 0);
            legend(subplot_handles,SubLegendLabels, 'FontSize', 18, 'location','northwest');
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;
            if rep_type_index < 3
            title(DeltaFCAx, [ChannelNames{ch_index},' Replicate ', num2str(rep_type_index),', ', num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 20)
            else
                title(DeltaFCAx, [ChannelNames{ch_index}, ' Flipped Stain, ', num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 20)
            end
            
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedSmoothedTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedSmoothedTestProfiles'])
            end
            if rep_type_index < 3
            outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedSmoothedTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_StainReplicate_', num2str(rep_type_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
            else
                 outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedSmoothedTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_FlippedStain_', num2str(rep_type_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
            end
            saveas(DeltaFCFixCorrectedFig,outpath2);
                
        end
    end
 end

 
%%
times_to_plot = 10:5:50;
colorsDubuisTimes = (hsv(round(length(times_to_plot)*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:length(times_to_plot),:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
RepMarkerStyles = {'o', 's', 'd'};
TempTickLabels = {};
RepList = {'Replicate 1', 'Replicate 2', 'FlippedStain'};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end
LegendLabels = {};
for temp_index = 1:length(temps)
    LegendLabels{temp_index} = [num2str(round(temps(temp_index),1)), 'ºC'];
end
LegendLabels = {'25.0ºC vs. 25.0ºC', '27.5ºC vs. 25.0ºC', '22.5ºC vs. 25.0ºC', '20.0ºC vs. 25.0ºC', '17.5ºC vs. 25.0ºC'};
temp_plot_order = [4 5 3 2 1];
for idx = 1:2
    ch_index = chList(idx);
    for time_index = 1:length(times_to_plot)
        for temp_index = 1:5
            LegendLabels = {};
            for i = 1:3
                LegendLabels{i} = [num2str(temps(temp_plot_order(temp_index))), 'ºC ', RepList{i}];
            end
            set_indices = (1:3)+(temp_index-1)*3;
            PlotMeanProfs = cell(1, 5);
            PlotMeanSDs = cell(1, 5);
            PlotSlopes = NaN(1, 5);
            PlotSlopes(1) = 1;
            PlotIntercepts = NaN(1, 5);
            PlotIntercepts(1) = 0;
            ymaxes = NaN(1, 5);
            ymins = NaN(1, 5);
            xmaxes = NaN(1, 5);
            for idx2 = 1:3
                PlotMeanProfs{idx2} = MeanProfiles{set_indices(idx2), time_index, idx};
                PlotMeanSDs{idx2} = SDProfiles{set_indices(idx2), time_index, idx};
                if ~isempty(MeanProfiles{set_indices(idx2), time_index, idx}+SDProfiles{set_indices(idx2), time_index, idx}) 
                ymaxes(idx2) = max(MeanProfiles{set_indices(idx2), time_index, idx}+SDProfiles{set_indices(idx2), time_index, idx});
                xmaxes(idx2) = max(MeanProfiles{set_indices(idx2), time_index, idx});
                ymins(idx2) = min(MeanProfiles{set_indices(idx2), time_index, idx}+SDProfiles{set_indices(idx2), time_index, idx});
                end
                if idx2 ~= 1
                    dlm = CrossFitMeanProfiles{set_indices(1), set_indices(idx2), time_index, idx};
                    if ~isempty(dlm)
                    PlotSlopes(idx2) = dlm.Coefficients.Estimate(2);
                    PlotIntercepts(idx2) = dlm.Coefficients.Estimate(1);
                    end
                end
            end
            plot_ymax = (ceil(max(ymaxes)/0.2))*0.2;
            plot_xmax = (ceil(max(xmaxes)/0.2))*0.2;
            plot_ymin = floor(min(ymins)/0.05)*0.05;
            close all
            DeltaFCFixCorrectedFig =figure(1);
            set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
            set(gcf,'color','w');
            DeltaFCAx = axes(DeltaFCFixCorrectedFig);
            map = colormap(colorsTemperatures);
            for rep_index = 1:3
                plot([plot_ymin plot_ymax].',  ([plot_ymin plot_ymax]*PlotSlopes(rep_index)+PlotIntercepts(rep_index)).',...
                    'Color', colorsTemperatures(rep_index,:),...
                    'LineWidth', 2.0);
%                 plot(PlotMeanProfs{1}, PlotMeanProfs{temp_index}, colorsTemperatures(temp_plot_order(temp_index),:),...
%                     'LineWidth', 3.0);    
                hold on 
            end
            grid on 
            plot_handles = zeros(1, 3);
            for rep_index = 1:3
                if ~isempty(PlotMeanProfs{1}) & ~isempty( PlotMeanProfs{rep_index})
                 plot_handles(rep_index) = errorbar(PlotMeanProfs{1}, PlotMeanProfs{rep_index},PlotMeanSDs{rep_index},...
                     'o', 'MarkerSize', 10,'Color', colorsTemperatures(rep_index, :),...
                               'LineStyle', 'none', 'MarkerFaceColor', colorsTemperatures(rep_index, :),...
                                'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
                end
                            
            end
            
            hold off
            xlabel(['Intensity ',num2str(temps(temp_plot_order(temp_index))),'ºC Replicate 1'], 'FontSize', 18)
            ylabel(['Intensity ',num2str(temps(temp_plot_order(temp_index))),'ºC Replicate 1/Replicate 2/Flipped'], 'FontSize', 18);
            ylim([0,  plot_ymax])
            xlim([0 plot_xmax])
            subplot_handles = plot_handles;
            SubLegendLabels = LegendLabels;
            SubLegendLabels = SubLegendLabels(subplot_handles ~= 0);
            subplot_handles = subplot_handles(subplot_handles ~= 0);
            legend(subplot_handles,SubLegendLabels, 'FontSize', 18, 'location','northwest');
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;

            title(DeltaFCAx, [ChannelNames{ch_index},' ', num2str(temps(temp_plot_order(temp_index))),'ºC ', num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 20)
           
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedMeanTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedMeanTestProfiles'])
            end

            outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedMeanTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_T',strrep(num2str(temps(temp_plot_order(temp_index))), '.', '_'),'C_Time', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
      
            saveas(DeltaFCFixCorrectedFig,outpath2);
                
        end
    end
end
   

%%
times_to_plot = 10:5:50;
colorsDubuisTimes = (hsv(round(length(times_to_plot)*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:length(times_to_plot),:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
RepMarkerStyles = {'o', 's', 'd'};
TempTickLabels = {};
RepList = {'Replicate 1', 'Replicate 2', 'FlippedStain'};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end
LegendLabels = {};
for temp_index = 1:length(temps)
    LegendLabels{temp_index} = [num2str(round(temps(temp_index),1)), 'ºC'];
end
LegendLabels = {'25.0ºC vs. 25.0ºC', '27.5ºC vs. 25.0ºC', '22.5ºC vs. 25.0ºC', '20.0ºC vs. 25.0ºC', '17.5ºC vs. 25.0ºC'};
temp_plot_order = [4 5 3 2 1];
for idx = 1:2
    ch_index = chList(idx);
    for time_index = 1:length(times_to_plot)
        for temp_index = 1:5
            LegendLabels = {};
            for i = 1:3
                LegendLabels{i} = [num2str(temps(temp_plot_order(temp_index))), 'ºC ', RepList{i}];
            end
            set_indices = (1:3)+(temp_index-1)*3;
            PlotMeanProfs = cell(1, 5);
            PlotSlopes = NaN(1, 5);
            PlotSlopes(1) = 1;
            PlotIntercepts = NaN(1, 5);
            PlotIntercepts(1) = 0;
            ymaxes = NaN(1, 5);
            ymins = NaN(1, 5);
            xmaxes = NaN(1, 5);
            for idx2 = 1:3
                PlotMeanProfs{idx2} = yfits{set_indices(idx2), time_index, idx};
                if ~isempty(yfits{set_indices(idx2), time_index, idx}) 
                ymaxes(idx2) = max(yfits{set_indices(idx2), time_index, idx});
                xmaxes(idx2) = max(yfits{set_indices(idx2), time_index, idx});
                ymins(idx2) = min(yfits{set_indices(idx2), time_index, idx});
                end
                if idx2 ~= 1
                    dlm = CrossFitSmoothedProfiles{set_indices(1), set_indices(idx2), time_index, idx};
                    if ~isempty(dlm)
                    PlotSlopes(idx2) = dlm.Coefficients.Estimate(2);
                    PlotIntercepts(idx2) = dlm.Coefficients.Estimate(1);
                    end
                end
            end
            plot_ymax = (ceil(max(ymaxes)/0.2))*0.2;
            plot_xmax = (ceil(max(xmaxes)/0.2))*0.2;
            plot_ymin = floor(min(ymins)/0.05)*0.05;
            close all
            DeltaFCFixCorrectedFig =figure(1);
            set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
            set(gcf,'color','w');
            DeltaFCAx = axes(DeltaFCFixCorrectedFig);
            map = colormap(colorsTemperatures);
            for rep_index = 1:3
                plot([plot_ymin plot_ymax].',  ([plot_ymin plot_ymax]*PlotSlopes(rep_index)+PlotIntercepts(rep_index)).',...
                    'Color', colorsTemperatures(rep_index,:),...
                    'LineWidth', 2.0);
%                 plot(PlotMeanProfs{1}, PlotMeanProfs{temp_index}, colorsTemperatures(temp_plot_order(temp_index),:),...
%                     'LineWidth', 3.0);    
                hold on 
            end
            grid on 
            plot_handles = zeros(1, 3);
            for rep_index = 1:3
                if ~isempty(PlotMeanProfs{1}) & ~isempty( PlotMeanProfs{rep_index})
                 plot_handles(rep_index) = scatter(PlotMeanProfs{1}, PlotMeanProfs{rep_index},100,...
                     'MarkerFaceColor', colorsTemperatures(rep_index, :),...
                                'MarkerEdgeColor', 'k');
                end
                            
            end
            
            hold off
            xlabel(['Intensity ',num2str(temps(temp_plot_order(temp_index))),'ºC Replicate 1'], 'FontSize', 18)
            ylabel(['Intensity ',num2str(temps(temp_plot_order(temp_index))),'ºC Replicate 1/Replicate 2/Flipped'], 'FontSize', 18);
            ylim([0,  plot_ymax])
            xlim([0 plot_xmax])
            subplot_handles = plot_handles;
            SubLegendLabels = LegendLabels;
            SubLegendLabels = SubLegendLabels(subplot_handles ~= 0);
            subplot_handles = subplot_handles(subplot_handles ~= 0);
            legend(subplot_handles,SubLegendLabels, 'FontSize', 18, 'location','northwest');
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;

            title(DeltaFCAx, [ChannelNames{ch_index},' ', num2str(temps(temp_plot_order(temp_index))),'ºC ', num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 20)
           
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedSmoothedTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedSmoothedTestProfiles'])
            end

            outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedSmoothedTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_T',strrep(num2str(temps(temp_plot_order(temp_index))), '.', '_'),'C_Time', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
      
            saveas(DeltaFCFixCorrectedFig,outpath2);
                
        end
    end
end

%%

AllSetInfo = GetFixedSetPrefixInfo;
CrossFitMeanSlopes = NaN(NumSets, NumSets, length(times_to_plot), 2);
CrossFitMeanIntercepts = NaN(NumSets, NumSets, length(times_to_plot), 2);
CrossFitMeanSlopeSEs = NaN(NumSets, NumSets, length(times_to_plot), 2);
CrossFitMeanInterceptSEs = NaN(NumSets, NumSets, length(times_to_plot), 2);
AllSets = 1:NumSets;

for exp_index = 1:NumSets
    for exp_index2 = 1:NumSets
        for idx = 1:length(chList)
            ch_index = chList(idx);
            for time_index = 1:length(times_to_plot)
                if exp_index == exp_index2 
                    CrossFitMeanSlopes(exp_index, exp_index2, time_index, idx) = 1;
                    CrossFitMeanIntercepts(exp_index, exp_index2, time_index, idx) = 0;
                    CrossFitMeanSlopeSEs(exp_index, exp_index2, time_index, idx) = 0;
                    CrossFitMeanInterceptSEs(exp_index, exp_index2, time_index, idx) = 0;
                    continue
                end
                dlm = CrossFitMeanProfiles{exp_index, exp_index2, time_index, idx};
                if ~isempty(dlm)
                   CrossFitMeanSlopes(exp_index, exp_index2, time_index, idx) = dlm.Coefficients.Estimate(2);
                   CrossFitMeanIntercepts(exp_index, exp_index2, time_index, idx) = dlm.Coefficients.Estimate(1);
                   CrossFitMeanSlopeSEs(exp_index, exp_index2, time_index, idx) = dlm.Coefficients.SE(2);
                   CrossFitMeanInterceptSEs(exp_index, exp_index2, time_index, idx) = dlm.Coefficients.SE(1);
                end
                
            end
        end
    end
end
%%
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
temps = [17 20 22.5 25 27.5];
MarkerStyles = {'o', 's', 'd'};
for idx =  1:length(chList)
for exp_index = 1:NumSets

    ch_index= chList(idx);
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
            ymaxes = zeros(1, NumSets);
            ymins = zeros(1, NumSets);
            for temp_index = 1:length(temps)
                MatchIndices = find(AllSetInfo.Temperatures == temps(temp_index));
                for rep_index = 1:length(MatchIndices)
                    errorbar(times_to_plot, squeeze(CrossFitMeanSlopes(exp_index, MatchIndices(rep_index), :,idx)).',...
                        squeeze(CrossFitMeanSlopeSEs(exp_index, MatchIndices(rep_index), :,idx)).', [MarkerStyles{rep_index}, '-'],...
                        'MarkerSize', 8,'LineWidth', 1.5, 'Color', colorsTemperatures(temp_index,:),...
                        'MarkerFaceColor', colorsTemperatures(temp_index,:),'MarkerEdgeColor', 'k');
                    ymaxes((temp_index-1)*3+rep_index) = max(squeeze(CrossFitMeanSlopes(exp_index, MatchIndices(rep_index), :,idx)).'+ squeeze(CrossFitMeanSlopeSEs(exp_index, MatchIndices(rep_index), :,idx)).');
                    ymins((temp_index-1)*3+rep_index) = min(squeeze(CrossFitMeanSlopes(exp_index, MatchIndices(rep_index), :,idx)).'- squeeze(CrossFitMeanSlopeSEs(exp_index, MatchIndices(rep_index), :,idx)).');
                    hold on 
                end
            end
            plot_ymax = (ceil(max(ymaxes)/0.5))*0.5;
            plot_xmax = (ceil(max(xmaxes)/0.5))*0.5;
            plot_ymin = floor(min(ymins)/0.1)*0.1;
            grid on 
            hold off
            xlabel(['Time into NC14 (min)'], 'FontSize', 18)
            ylabel(['Slope of Intensity Fit'], 'FontSize', 18);
            ylim([plot_ymin,  plot_ymax])
            xlim([5 max(times_to_plot)+5])
%             subplot_handles = plot_handles;
%             SubLegendLabels = LegendLabels;
%             SubLegendLabels = SubLegendLabels(subplot_handles ~= 0);
%             subplot_handles = subplot_handles(subplot_handles ~= 0);
%             legend(subplot_handles,SubLegendLabels, 'FontSize', 18, 'location','northwest');
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;
            
            if AllSetInfo.Replicates(exp_index) ~= 0
            title(DeltaFCAx, [ChannelNames{ch_index},' Fit Slopes to ', num2str(AllSetInfo.Temperatures(exp_index)), 'ºC Replicate ' num2str(AllSetInfo.Replicates(exp_index))], 'FontSize', 20)
            else
                title(DeltaFCAx, [ChannelNames{ch_index},' Fit Slopes to ', num2str(AllSetInfo.Temperatures(exp_index)), 'ºC Flipped Stain'], 'FontSize', 20);
            end
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedMeanTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedMeanTestProfiles'])
            end

            outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedMeanTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_Set', num2str(exp_index),'_AllTimesAndTemperatures.png' ];
      
            saveas(DeltaFCFixCorrectedFig,outpath2);
            
            
end
end
   

%% Re-do Slope v Time Plots with Standard Normalized 
CrossFitMeanSlopes = NaN(NumSets, NumSets, length(times_to_plot), 2);
CrossFitMeanIntercepts = NaN(NumSets, NumSets, length(times_to_plot), 2);
CrossFitMeanSlopeSEs = NaN(NumSets, NumSets, length(times_to_plot), 2);
CrossFitMeanInterceptSEs = NaN(NumSets, NumSets, length(times_to_plot), 2);
AllSets = 1:NumSets;

for exp_index = 1:NumSets
    for exp_index2 = 1:NumSets
        for idx = 1:length(chList)
            ch_index = chList(idx);
            for time_index = 1:length(times_to_plot)
                if exp_index == exp_index2 
                    CrossFitMeanSlopes(exp_index, exp_index2, time_index, idx) = 1;
                    CrossFitMeanIntercepts(exp_index, exp_index2, time_index, idx) = 0;
                    CrossFitMeanSlopeSEs(exp_index, exp_index2, time_index, idx) = 0;
                    CrossFitMeanInterceptSEs(exp_index, exp_index2, time_index, idx) = 0;
                    continue
                end
                dlm = CrossFitSmoothedProfiles{exp_index, exp_index2, time_index, idx};
                if ~isempty(dlm)
                   CrossFitMeanSlopes(exp_index, exp_index2, time_index, idx) = dlm.Coefficients.Estimate(2);
                   if CrossFitMeanSlopes(exp_index, exp_index2, time_index, idx) == 0
                       CrossFitMeanSlopes(exp_index, exp_index2, time_index, idx) = NaN;
                       CrossFitMeanSlopeSEs(exp_index, exp_index2, time_index, idx) = NaN;
                   end
                   CrossFitMeanIntercepts(exp_index, exp_index2, time_index, idx) = dlm.Coefficients.Estimate(1);
                   CrossFitMeanSlopeSEs(exp_index, exp_index2, time_index, idx) = dlm.Coefficients.SE(2);
                   CrossFitMeanInterceptSEs(exp_index, exp_index2, time_index, idx) = dlm.Coefficients.SE(1);
                end
                
            end
        end
    end
end
%%
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
temps = [17 20 22.5 25 27.5];
MarkerStyles = {'o', 's', 'd'};
for idx =  1:length(chList)
for exp_index = 1:NumSets

    ch_index= chList(idx);
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
            ymaxes = zeros(1, NumSets);
            ymins = zeros(1, NumSets);
            for temp_index = 1:length(temps)
                MatchIndices = find(AllSetInfo.Temperatures == temps(temp_index));
                for rep_index = 1:length(MatchIndices)
                    errorbar(times_to_plot, squeeze(CrossFitMeanSlopes(exp_index, MatchIndices(rep_index), :,idx)).',...
                        squeeze(CrossFitMeanSlopeSEs(exp_index, MatchIndices(rep_index), :,idx)).', [MarkerStyles{rep_index}, '-'],...
                        'MarkerSize', 8,'LineWidth', 1.5, 'Color', colorsTemperatures(temp_index,:),...
                        'MarkerFaceColor', colorsTemperatures(temp_index,:),'MarkerEdgeColor', 'k');
                    ymaxes((temp_index-1)*3+rep_index) = max(squeeze(CrossFitMeanSlopes(exp_index, MatchIndices(rep_index), :,idx)).'+ squeeze(CrossFitMeanSlopeSEs(exp_index, MatchIndices(rep_index), :,idx)).');
                    ymins((temp_index-1)*3+rep_index) = min(squeeze(CrossFitMeanSlopes(exp_index, MatchIndices(rep_index), :,idx)).'- squeeze(CrossFitMeanSlopeSEs(exp_index, MatchIndices(rep_index), :,idx)).');
                    hold on 
                end
            end
            plot_ymax = (ceil(max(ymaxes)/0.5))*0.5;
            plot_xmax = (ceil(max(xmaxes)/0.5))*0.5;
            plot_ymin = floor(min(ymins)/0.1)*0.1;
            grid on 
            hold off
            xlabel(['Time into NC14 (min)'], 'FontSize', 18)
            ylabel(['Slope of Intensity Fit'], 'FontSize', 18);
            ylim([plot_ymin,  plot_ymax])
            xlim([5 max(times_to_plot)+5])
%             subplot_handles = plot_handles;
%             SubLegendLabels = LegendLabels;
%             SubLegendLabels = SubLegendLabels(subplot_handles ~= 0);
%             subplot_handles = subplot_handles(subplot_handles ~= 0);
%             legend(subplot_handles,SubLegendLabels, 'FontSize', 18, 'location','northwest');
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;
            
            if AllSetInfo.Replicates(exp_index) ~= 0
            title(DeltaFCAx, [ChannelNames{ch_index},' Fit Slopes to ', num2str(AllSetInfo.Temperatures(exp_index)), 'ºC Replicate ' num2str(AllSetInfo.Replicates(exp_index))], 'FontSize', 20)
            else
                title(DeltaFCAx, [ChannelNames{ch_index},' Fit Slopes to ', num2str(AllSetInfo.Temperatures(exp_index)), 'ºC Flipped Stain'], 'FontSize', 20);
            end
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedSmoothedTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedSmoothedTestProfiles'])
            end

            outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedSmoothedTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_Set', num2str(exp_index),'_AllTimesAndTemperatures.png' ];
      
            saveas(DeltaFCFixCorrectedFig,outpath2);
            
            
end
end
   
%% REDO WITH BIG GROUPING

times_to_plot = [15, 45];
colorsDubuisTimes = (hsv(round(9*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:9,:);
colorsTemperatures = colorsDubuisTimes([9 5 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
RepMarkerStyles = {'o', 's', 'd'};
TempTickLabels = {};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end
LegendLabels = {};
for temp_index = 1:length(temps)
    LegendLabels{temp_index} = [num2str(round(temps(temp_index),1)), 'ºC'];
end

x_sample = cell(NumSets, length(times_to_plot), 2);
ys = cell(NumSets, length(times_to_plot), 2);
MeanProfiles = cell(NumSets, length(times_to_plot), 2);
SDProfiles = cell(NumSets, length(times_to_plot), 2);
CountProfiles = cell(NumSets, length(times_to_plot), 2);
yfits = cell(NumSets, length(times_to_plot), 2);
master_index = 7;
chList = [3 5];
for exp_index = 1:15
    CompiledEmbryos = AllCompiledEmbryos{exp_index};
    
    for idx = 1:length(chList)
        ch_index = chList(idx);
        for time_index = 1:length(times_to_plot)
            x = CompiledEmbryos.DubuisEmbryoTimes;
            TF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x) & (x >= times_to_plot(time_index)-15) & (x < times_to_plot(time_index)+15);
            x = x(TF);
            y = CompiledEmbryos.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(TF,:,ch_index, master_index);
            NonTF = sum(isnan(y), 2).' > NumAPbins;
            x = x(~NonTF);
            x_sample{exp_index, time_index, idx} = x;
            y = y(~NonTF, :);
            ys{exp_index, time_index, idx} = y;
            if ~isempty(y)
                MeanProfiles{exp_index, time_index, idx} = mean(y, 1, 'omitnan');
                if size(y, 1) > 1
                    SDProfiles{exp_index, time_index, idx} = std(y, 1, 'omitnan');
                else
                    SDProfiles{exp_index, time_index, idx} = NaN(size(y));
                end
            end
            CountProfiles{exp_index, time_index, idx} = size(ys, 1);
            xfit = CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            [minval, matchidx] = min(abs(xfit - times_to_plot(time_index)));
            yf = CompiledEmbryos.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(matchidx,:,ch_index, master_index);
            yfits{exp_index, time_index, idx} = yf;
        end
    end
end
       
CrossFitMeanProfiles =  cell(NumSets, NumSets, length(times_to_plot), 2);
AllSets = 1:NumSets;
for exp_index = 1:NumSets
    for exp_index2 = AllSets(~ismember(AllSets, exp_index))
        for idx = 1:length(chList)
            for time_index = 1:length(times_to_plot)
                MeanProf = MeanProfiles{exp_index, time_index, idx};
                MeanProf2 = MeanProfiles{exp_index2, time_index, idx};
                if isempty(MeanProf) |  isempty(MeanProf2) 
                    continue
                end
                dlm = fitlm(MeanProf, MeanProf2);
                CrossFitMeanProfiles{exp_index, exp_index2, time_index, idx} = dlm;
            end
        end
    end
end

CrossFitSmoothedProfiles =  cell(NumSets, NumSets, length(times_to_plot), 2);
AllSets = 1:NumSets;
for exp_index = 1:NumSets
    for exp_index2 = AllSets(~ismember(AllSets, exp_index))
        for idx = 1:length(chList)
            for time_index = 1:length(times_to_plot)
                MeanProf = yfits{exp_index, time_index, idx};
                MeanProf2 = yfits{exp_index2, time_index, idx};
                if isempty(MeanProf) |  isempty(MeanProf2) 
                    continue
                end
                dlm = fitlm(MeanProf, MeanProf2);
                CrossFitSmoothedProfiles{exp_index, exp_index2, time_index, idx} = dlm;
            end
        end
    end
end



%%

times_to_plot = [15, 45];
colorsDubuisTimes = (hsv(round(9*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:9,:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
RepMarkerStyles = {'o', 's', 'd'};
TempTickLabels = {};
RepList = {'Replicate 1', 'Replicate 2', 'FlippedStain'};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end
LegendLabels = {};
for temp_index = 1:length(temps)
    LegendLabels{temp_index} = [num2str(round(temps(temp_index),1)), 'ºC'];
end
LegendLabels = {'25.0ºC vs. 25.0ºC', '27.5ºC vs. 25.0ºC', '22.5ºC vs. 25.0ºC', '20.0ºC vs. 25.0ºC', '17.5ºC vs. 25.0ºC'};
temp_plot_order = [4 5 3 2 1];
for idx = 1:2
    ch_index = chList(idx);
    for time_index = 1:length(times_to_plot)
        for rep_type_index = 1:3
            set_indices = (0:3:12)+rep_type_index;
            PlotMeanProfs = cell(1, 5);
            PlotMeanSDs = cell(1, 5);
            PlotSlopes = NaN(1, 5);
            PlotSlopes(1) = 1;
            PlotIntercepts = NaN(1, 5);
            PlotIntercepts(1) = 0;
            ymaxes = NaN(1, 5);
            ymins = NaN(1, 5);
            xmaxes = NaN(1, 5);
            for idx2 = 1:5
                PlotMeanProfs{idx2} = MeanProfiles{set_indices(idx2), time_index, idx};
                PlotMeanSDs{idx2} = SDProfiles{set_indices(idx2), time_index, idx};
                if ~isempty(MeanProfiles{set_indices(idx2), time_index, idx}+SDProfiles{set_indices(idx2), time_index, idx}) 
                ymaxes(idx2) = max(MeanProfiles{set_indices(idx2), time_index, idx}+SDProfiles{set_indices(idx2), time_index, idx});
                xmaxes(idx2) = max(MeanProfiles{set_indices(idx2), time_index, idx});
                ymins(idx2) = min(MeanProfiles{set_indices(idx2), time_index, idx}+SDProfiles{set_indices(idx2), time_index, idx});
                end
                if idx2 ~= 1
                    dlm = CrossFitMeanProfiles{set_indices(1), set_indices(idx2), time_index, idx};
                    if ~isempty(dlm)
                    PlotSlopes(idx2) = dlm.Coefficients.Estimate(2);
                    PlotIntercepts(idx2) = dlm.Coefficients.Estimate(1);
                    end
                end
            end
            plot_ymax = (ceil(max(ymaxes)/0.2))*0.2;
            plot_xmax = (ceil(max(xmaxes)/0.2))*0.2;
            plot_ymin = floor(min(ymins)/0.05)*0.05;
            close all
            DeltaFCFixCorrectedFig =figure(1);
            set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
            set(gcf,'color','w');
            DeltaFCAx = axes(DeltaFCFixCorrectedFig);
            map = colormap(colorsTemperatures);
            for temp_index = 1:length(temps)
                plot([plot_ymin plot_ymax].',  ([plot_ymin plot_ymax]*PlotSlopes(temp_index)+PlotIntercepts(temp_index)).',...
                    'Color', colorsTemperatures(temp_plot_order(temp_index),:),...
                    'LineWidth', 2.0);
%                 plot(PlotMeanProfs{1}, PlotMeanProfs{temp_index}, colorsTemperatures(temp_plot_order(temp_index),:),...
%                     'LineWidth', 3.0);    
                hold on 
            end
            grid on 
            plot_handles = zeros(1, length(temps));
            for temp_index = 1:length(temps)
                if ~isempty(PlotMeanProfs{1}) & ~isempty( PlotMeanProfs{temp_index})
                 plot_handles(temp_index) = errorbar(PlotMeanProfs{1}, PlotMeanProfs{temp_index},PlotMeanSDs{temp_index},...
                     'o', 'MarkerSize', 10,'Color', colorsTemperatures(temp_plot_order(temp_index), :),...
                               'LineStyle', 'none', 'MarkerFaceColor', colorsTemperatures(temp_plot_order(temp_index), :),...
                                'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
                end
                            
            end
            
            hold off
            xlabel('Intensity 25.0ºC', 'FontSize', 18)
            ylabel('Intensity 27.5/25.0/22.5/20.0/17.5ºC', 'FontSize', 18);
            ylim([0,  plot_ymax])
            xlim([0 plot_xmax])
            subplot_handles = plot_handles([2 1 3 4 5]);
            SubLegendLabels = LegendLabels([2 1 3 4 5]);
            SubLegendLabels = SubLegendLabels(subplot_handles ~= 0);
            subplot_handles = subplot_handles(subplot_handles ~= 0);
            legend(subplot_handles,SubLegendLabels, 'FontSize', 18, 'location','northwest');
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;
            if rep_type_index < 3
            title(DeltaFCAx, [ChannelNames{ch_index},' Replicate ', num2str(rep_type_index),', ', num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 20)
            else
                title(DeltaFCAx, [ChannelNames{ch_index}, ' Flipped Stain, ', num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 20)
            end
            
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedMeanTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedMeanTestProfiles'])
            end
            if rep_type_index < 3
            outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedMeanTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_StainReplicate_', num2str(rep_type_index),'_WideTime', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
            else
                 outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedMeanTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_FlippedStain_', num2str(rep_type_index),'_WideTime', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
            end
            saveas(DeltaFCFixCorrectedFig,outpath2);
                
        end
    end
end

 %% Use Smoothed Profiles
 for idx = 1:2
    ch_index = chList(idx);
    for time_index = 1:length(times_to_plot)
        for rep_type_index = 1:3
            set_indices = (0:3:12)+rep_type_index;
            PlotMeanProfs = cell(1, 5);
            PlotMeanSDs = cell(1, 5);
            PlotSlopes = NaN(1, 5);
            PlotSlopes(1) = 1;
            PlotIntercepts = NaN(1, 5);
            PlotIntercepts(1) = 0;
            ymaxes = NaN(1, 5);
            ymins = NaN(1, 5);
            xmaxes = NaN(1, 5);
            for idx2 = 1:5
                PlotMeanProfs{idx2} = yfits{set_indices(idx2), time_index, idx};
                if ~isempty(yfits{set_indices(idx2), time_index, idx}) 
                ymaxes(idx2) = max(yfits{set_indices(idx2), time_index, idx});
                xmaxes(idx2) = max(yfits{set_indices(idx2), time_index, idx});
                ymins(idx2) = min(yfits{set_indices(idx2), time_index, idx});
                end
                if idx2 ~= 1
                    dlm = CrossFitSmoothedProfiles{set_indices(1), set_indices(idx2), time_index, idx};
                    if ~isempty(dlm)
                    PlotSlopes(idx2) = dlm.Coefficients.Estimate(2);
                    PlotIntercepts(idx2) = dlm.Coefficients.Estimate(1);
                    end
                end
            end
            if all(isnan(ymaxes))
                continue
            end
            plot_ymax = (ceil(max(ymaxes)/0.2))*0.2;
            plot_xmax = (ceil(max(xmaxes)/0.2))*0.2;
            plot_ymin = floor(min(ymins)/0.05)*0.05;
            close all
            DeltaFCFixCorrectedFig =figure(1);
            set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
            set(gcf,'color','w');
            DeltaFCAx = axes(DeltaFCFixCorrectedFig);
            map = colormap(colorsTemperatures);
            for temp_index = 1:length(temps)
                plot([plot_ymin plot_ymax].',  ([plot_ymin plot_ymax]*PlotSlopes(temp_index)+PlotIntercepts(temp_index)).',...
                    'Color', colorsTemperatures(temp_plot_order(temp_index),:),...
                    'LineWidth', 2.0);
%                 plot(PlotMeanProfs{1}, PlotMeanProfs{temp_index}, colorsTemperatures(temp_plot_order(temp_index),:),...
%                     'LineWidth', 3.0);    
                hold on 
            end
            grid on 
            plot_handles = zeros(1, length(temps));
            for temp_index = 1:length(temps)
                if ~isempty(PlotMeanProfs{1}) & ~isempty( PlotMeanProfs{temp_index})
                 plot_handles(temp_index) = scatter(PlotMeanProfs{1}, PlotMeanProfs{temp_index},...
                    100,'MarkerFaceColor', colorsTemperatures(temp_plot_order(temp_index), :),'MarkerEdgeColor', 'k');
                end
                            
            end
            
            hold off
            xlabel('Intensity 25.0ºC', 'FontSize', 18)
            ylabel('Intensity 27.5/25.0/22.5/20.0/17.5ºC', 'FontSize', 18);
            ylim([0,  plot_ymax])
            xlim([0 plot_xmax])
            subplot_handles = plot_handles([2 1 3 4 5]);
            SubLegendLabels = LegendLabels([2 1 3 4 5]);
            SubLegendLabels = SubLegendLabels(subplot_handles ~= 0);
            subplot_handles = subplot_handles(subplot_handles ~= 0);
            legend(subplot_handles,SubLegendLabels, 'FontSize', 18, 'location','northwest');
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;
            if rep_type_index < 3
            title(DeltaFCAx, [ChannelNames{ch_index},' Replicate ', num2str(rep_type_index),', ', num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 20)
            else
                title(DeltaFCAx, [ChannelNames{ch_index}, ' Flipped Stain, ', num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 20)
            end
            
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedSmoothedTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedSmoothedTestProfiles'])
            end
            if rep_type_index < 3
            outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedSmoothedTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_StainReplicate_', num2str(rep_type_index),'_WideTime', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
            else
                 outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedSmoothedTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_FlippedStain_', num2str(rep_type_index),'_WideTime', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
            end
            saveas(DeltaFCFixCorrectedFig,outpath2);
                
        end
    end
 end

 
%%
times_to_plot = [15 45];
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
RepMarkerStyles = {'o', 's', 'd'};
TempTickLabels = {};
RepList = {'Replicate 1', 'Replicate 2', 'FlippedStain'};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end
LegendLabels = {};
for temp_index = 1:length(temps)
    LegendLabels{temp_index} = [num2str(round(temps(temp_index),1)), 'ºC'];
end
LegendLabels = {'25.0ºC vs. 25.0ºC', '27.5ºC vs. 25.0ºC', '22.5ºC vs. 25.0ºC', '20.0ºC vs. 25.0ºC', '17.5ºC vs. 25.0ºC'};
temp_plot_order = [4 5 3 2 1];
for idx = 1:2
    ch_index = chList(idx);
    for time_index = 1:length(times_to_plot)
        for temp_index = 1:5
            LegendLabels = {};
            for i = 1:3
                LegendLabels{i} = [num2str(temps(temp_plot_order(temp_index))), 'ºC ', RepList{i}];
            end
            set_indices = (1:3)+(temp_index-1)*3;
            PlotMeanProfs = cell(1, 5);
            PlotMeanSDs = cell(1, 5);
            PlotSlopes = NaN(1, 5);
            PlotSlopes(1) = 1;
            PlotIntercepts = NaN(1, 5);
            PlotIntercepts(1) = 0;
            ymaxes = NaN(1, 5);
            ymins = NaN(1, 5);
            xmaxes = NaN(1, 5);
            for idx2 = 1:3
                PlotMeanProfs{idx2} = MeanProfiles{set_indices(idx2), time_index, idx};
                PlotMeanSDs{idx2} = SDProfiles{set_indices(idx2), time_index, idx};
                if ~isempty(MeanProfiles{set_indices(idx2), time_index, idx}+SDProfiles{set_indices(idx2), time_index, idx}) 
                ymaxes(idx2) = max(MeanProfiles{set_indices(idx2), time_index, idx}+SDProfiles{set_indices(idx2), time_index, idx});
                xmaxes(idx2) = max(MeanProfiles{set_indices(idx2), time_index, idx});
                ymins(idx2) = min(MeanProfiles{set_indices(idx2), time_index, idx}+SDProfiles{set_indices(idx2), time_index, idx});
                end
                if idx2 ~= 1
                    dlm = CrossFitMeanProfiles{set_indices(1), set_indices(idx2), time_index, idx};
                    if ~isempty(dlm)
                    PlotSlopes(idx2) = dlm.Coefficients.Estimate(2);
                    PlotIntercepts(idx2) = dlm.Coefficients.Estimate(1);
                    end
                end
            end
            plot_ymax = (ceil(max(ymaxes)/0.2))*0.2;
            plot_xmax = (ceil(max(xmaxes)/0.2))*0.2;
            plot_ymin = floor(min(ymins)/0.05)*0.05;
            close all
            DeltaFCFixCorrectedFig =figure(1);
            set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
            set(gcf,'color','w');
            DeltaFCAx = axes(DeltaFCFixCorrectedFig);
            map = colormap(colorsTemperatures);
            for rep_index = 1:3
                plot([plot_ymin plot_ymax].',  ([plot_ymin plot_ymax]*PlotSlopes(rep_index)+PlotIntercepts(rep_index)).',...
                    'Color', colorsTemperatures(rep_index,:),...
                    'LineWidth', 2.0);
%                 plot(PlotMeanProfs{1}, PlotMeanProfs{temp_index}, colorsTemperatures(temp_plot_order(temp_index),:),...
%                     'LineWidth', 3.0);    
                hold on 
            end
            grid on 
            plot_handles = zeros(1, 3);
            for rep_index = 1:3
                if ~isempty(PlotMeanProfs{1}) & ~isempty( PlotMeanProfs{rep_index})
                 plot_handles(rep_index) = errorbar(PlotMeanProfs{1}, PlotMeanProfs{rep_index},PlotMeanSDs{rep_index},...
                     'o', 'MarkerSize', 10,'Color', colorsTemperatures(rep_index, :),...
                               'LineStyle', 'none', 'MarkerFaceColor', colorsTemperatures(rep_index, :),...
                                'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
                end
                            
            end
            
            hold off
            xlabel(['Intensity ',num2str(temps(temp_plot_order(temp_index))),'ºC Replicate 1'], 'FontSize', 18)
            ylabel(['Intensity ',num2str(temps(temp_plot_order(temp_index))),'ºC Replicate 1/Replicate 2/Flipped'], 'FontSize', 18);
            ylim([0,  plot_ymax])
            xlim([0 plot_xmax])
            subplot_handles = plot_handles;
            SubLegendLabels = LegendLabels;
            SubLegendLabels = SubLegendLabels(subplot_handles ~= 0);
            subplot_handles = subplot_handles(subplot_handles ~= 0);
            legend(subplot_handles,SubLegendLabels, 'FontSize', 18, 'location','northwest');
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;

            title(DeltaFCAx, [ChannelNames{ch_index},' ', num2str(temps(temp_plot_order(temp_index))),'ºC ', num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 20)
           
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedMeanTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedMeanTestProfiles'])
            end

            outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedMeanTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_WideTime', strrep(num2str(times_to_plot(time_index)), '.', '_'),'_mT',strrep(num2str(temps(temp_plot_order(temp_index))), '.', '_'), 'C.png' ];
      
            saveas(DeltaFCFixCorrectedFig,outpath2);
                
        end
    end
end
   

%%

FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
RepMarkerStyles = {'o', 's', 'd'};
TempTickLabels = {};
RepList = {'Replicate 1', 'Replicate 2', 'FlippedStain'};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end
LegendLabels = {};
for temp_index = 1:length(temps)
    LegendLabels{temp_index} = [num2str(round(temps(temp_index),1)), 'ºC'];
end
LegendLabels = {'25.0ºC vs. 25.0ºC', '27.5ºC vs. 25.0ºC', '22.5ºC vs. 25.0ºC', '20.0ºC vs. 25.0ºC', '17.5ºC vs. 25.0ºC'};
temp_plot_order = [4 5 3 2 1];
for idx = 1:2
    ch_index = chList(idx);
    for time_index = 1:length(times_to_plot)
        for temp_index = 1:5
            LegendLabels = {};
            for i = 1:3
                LegendLabels{i} = [num2str(temps(temp_plot_order(temp_index))), 'ºC ', RepList{i}];
            end
            set_indices = (1:3)+(temp_index-1)*3;
            PlotMeanProfs = cell(1, 5);
            PlotSlopes = NaN(1, 5);
            PlotSlopes(1) = 1;
            PlotIntercepts = NaN(1, 5);
            PlotIntercepts(1) = 0;
            ymaxes = NaN(1, 5);
            ymins = NaN(1, 5);
            xmaxes = NaN(1, 5);
            for idx2 = 1:3
                PlotMeanProfs{idx2} = yfits{set_indices(idx2), time_index, idx};
                if ~isempty(yfits{set_indices(idx2), time_index, idx}) 
                ymaxes(idx2) = max(yfits{set_indices(idx2), time_index, idx});
                xmaxes(idx2) = max(yfits{set_indices(idx2), time_index, idx});
                ymins(idx2) = min(yfits{set_indices(idx2), time_index, idx});
                end
                if idx2 ~= 1
                    dlm = CrossFitSmoothedProfiles{set_indices(1), set_indices(idx2), time_index, idx};
                    if ~isempty(dlm)
                    PlotSlopes(idx2) = dlm.Coefficients.Estimate(2);
                    PlotIntercepts(idx2) = dlm.Coefficients.Estimate(1);
                    end
                end
            end
            plot_ymax = (ceil(max(ymaxes)/0.2))*0.2;
            plot_xmax = (ceil(max(xmaxes)/0.2))*0.2;
            plot_ymin = floor(min(ymins)/0.05)*0.05;
            close all
            DeltaFCFixCorrectedFig =figure(1);
            set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
            set(gcf,'color','w');
            DeltaFCAx = axes(DeltaFCFixCorrectedFig);
            map = colormap(colorsTemperatures);
            for rep_index = 1:3
                plot([plot_ymin plot_ymax].',  ([plot_ymin plot_ymax]*PlotSlopes(rep_index)+PlotIntercepts(rep_index)).',...
                    'Color', colorsTemperatures(rep_index,:),...
                    'LineWidth', 2.0);
%                 plot(PlotMeanProfs{1}, PlotMeanProfs{temp_index}, colorsTemperatures(temp_plot_order(temp_index),:),...
%                     'LineWidth', 3.0);    
                hold on 
            end
            grid on 
            plot_handles = zeros(1, 3);
            for rep_index = 1:3
                if ~isempty(PlotMeanProfs{1}) & ~isempty( PlotMeanProfs{rep_index})
                 plot_handles(rep_index) = scatter(PlotMeanProfs{1}, PlotMeanProfs{rep_index},100,...
                     'MarkerFaceColor', colorsTemperatures(rep_index, :),...
                                'MarkerEdgeColor', 'k');
                end
                            
            end
            
            hold off
            xlabel(['Intensity ',num2str(temps(temp_plot_order(temp_index))),'ºC Replicate 1'], 'FontSize', 18)
            ylabel(['Intensity ',num2str(temps(temp_plot_order(temp_index))),'ºC Replicate 1/Replicate 2/Flipped'], 'FontSize', 18);
            ylim([0,  plot_ymax])
            xlim([0 plot_xmax])
            subplot_handles = plot_handles;
            SubLegendLabels = LegendLabels;
            SubLegendLabels = SubLegendLabels(subplot_handles ~= 0);
            subplot_handles = subplot_handles(subplot_handles ~= 0);
            legend(subplot_handles,SubLegendLabels, 'FontSize', 18, 'location','northwest');
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;

            title(DeltaFCAx, [ChannelNames{ch_index},' ', num2str(temps(temp_plot_order(temp_index))),'ºC ', num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 20)
           
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedSmoothedTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedSmoothedTestProfiles'])
            end

            outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedSmoothedTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_T',strrep(num2str(temps(temp_plot_order(temp_index))), '.', '_'),'C_WideTime', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'C.png' ];
      
            saveas(DeltaFCFixCorrectedFig,outpath2);
                
        end
    end
end


%%
CrossFitMeanSlopes = NaN(NumSets, NumSets, length(times_to_plot), 2);
CrossFitMeanIntercepts = NaN(NumSets, NumSets, length(times_to_plot), 2);
CrossFitMeanSlopeSEs = NaN(NumSets, NumSets, length(times_to_plot), 2);
CrossFitMeanInterceptSEs = NaN(NumSets, NumSets, length(times_to_plot), 2);
AllSets = 1:NumSets;

for exp_index = 1:NumSets
    for exp_index2 = 1:NumSets
        for idx = 1:length(chList)
            ch_index = chList(idx);
            for time_index = 1:length(times_to_plot)
                if exp_index == exp_index2 
                    CrossFitMeanSlopes(exp_index, exp_index2, time_index, idx) = 1;
                    CrossFitMeanIntercepts(exp_index, exp_index2, time_index, idx) = 0;
                    CrossFitMeanSlopeSEs(exp_index, exp_index2, time_index, idx) = 0;
                    CrossFitMeanInterceptSEs(exp_index, exp_index2, time_index, idx) = 0;
                    continue
                end
                dlm = CrossFitMeanProfiles{exp_index, exp_index2, time_index, idx};
                if ~isempty(dlm)
                   CrossFitMeanSlopes(exp_index, exp_index2, time_index, idx) = dlm.Coefficients.Estimate(2);
                   CrossFitMeanIntercepts(exp_index, exp_index2, time_index, idx) = dlm.Coefficients.Estimate(1);
                   CrossFitMeanSlopeSEs(exp_index, exp_index2, time_index, idx) = dlm.Coefficients.SE(2);
                   CrossFitMeanInterceptSEs(exp_index, exp_index2, time_index, idx) = dlm.Coefficients.SE(1);
                end
                
            end
        end
    end
end
%%

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
temps = [17 20 22.5 25 27.5];
MarkerStyles = {'o', 's', 'd'};
for idx =  1:length(chList)
for exp_index = 1:NumSets

    ch_index= chList(idx);
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
            ymaxes = zeros(1, NumSets);
            ymins = zeros(1, NumSets);
            for temp_index = 1:length(temps)
                MatchIndices = find(AllSetInfo.Temperatures == temps(temp_index));
                for rep_index = 1:length(MatchIndices)
                    errorbar(times_to_plot, squeeze(CrossFitMeanSlopes(exp_index, MatchIndices(rep_index), :,idx)).',...
                        squeeze(CrossFitMeanSlopeSEs(exp_index, MatchIndices(rep_index), :,idx)).', [MarkerStyles{rep_index}, '-'],...
                        'MarkerSize', 8,'LineWidth', 1.5, 'Color', colorsTemperatures(temp_index,:),...
                        'MarkerFaceColor', colorsTemperatures(temp_index,:),'MarkerEdgeColor', 'k');
                    ymaxes((temp_index-1)*3+rep_index) = max(squeeze(CrossFitMeanSlopes(exp_index, MatchIndices(rep_index), :,idx)).'+ squeeze(CrossFitMeanSlopeSEs(exp_index, MatchIndices(rep_index), :,idx)).');
                    ymins((temp_index-1)*3+rep_index) = min(squeeze(CrossFitMeanSlopes(exp_index, MatchIndices(rep_index), :,idx)).'- squeeze(CrossFitMeanSlopeSEs(exp_index, MatchIndices(rep_index), :,idx)).');
                    hold on 
                end
            end
            plot_ymax = (ceil(max(ymaxes)/0.5))*0.5;
            plot_xmax = (ceil(max(xmaxes)/0.5))*0.5;
            plot_ymin = floor(min(ymins)/0.1)*0.1;
            grid on 
            hold off
            xlabel(['Time into NC14 (min)'], 'FontSize', 18)
            ylabel(['Slope of Intensity Fit'], 'FontSize', 18);
            ylim([plot_ymin,  plot_ymax])
            xlim([5 max(times_to_plot)+5])
%             subplot_handles = plot_handles;
%             SubLegendLabels = LegendLabels;
%             SubLegendLabels = SubLegendLabels(subplot_handles ~= 0);
%             subplot_handles = subplot_handles(subplot_handles ~= 0);
%             legend(subplot_handles,SubLegendLabels, 'FontSize', 18, 'location','northwest');
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;
            
            if AllSetInfo.Replicates(exp_index) ~= 0
            title(DeltaFCAx, [ChannelNames{ch_index},' Fit Slopes to ', num2str(AllSetInfo.Temperatures(exp_index)), 'ºC Replicate ' num2str(AllSetInfo.Replicates(exp_index))], 'FontSize', 20)
            else
                title(DeltaFCAx, [ChannelNames{ch_index},' Fit Slopes to ', num2str(AllSetInfo.Temperatures(exp_index)), 'ºC Flipped Stain'], 'FontSize', 20);
            end
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedMeanTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedMeanTestProfiles'])
            end

            outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedMeanTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_Set', num2str(exp_index),'_AllWideTimesAndTemperatures.png' ];
      
            saveas(DeltaFCFixCorrectedFig,outpath2);
            
            
end
end
   

%% Re-do Slope v Time Plots with Standard Normalized 
CrossFitMeanSlopes = NaN(NumSets, NumSets, length(times_to_plot), 2);
CrossFitMeanIntercepts = NaN(NumSets, NumSets, length(times_to_plot), 2);
CrossFitMeanSlopeSEs = NaN(NumSets, NumSets, length(times_to_plot), 2);
CrossFitMeanInterceptSEs = NaN(NumSets, NumSets, length(times_to_plot), 2);
AllSets = 1:NumSets;

for exp_index = 1:NumSets
    for exp_index2 = 1:NumSets
        for idx = 1:length(chList)
            ch_index = chList(idx);
            for time_index = 1:length(times_to_plot)
                if exp_index == exp_index2 
                    CrossFitMeanSlopes(exp_index, exp_index2, time_index, idx) = 1;
                    CrossFitMeanIntercepts(exp_index, exp_index2, time_index, idx) = 0;
                    CrossFitMeanSlopeSEs(exp_index, exp_index2, time_index, idx) = 0;
                    CrossFitMeanInterceptSEs(exp_index, exp_index2, time_index, idx) = 0;
                    continue
                end
                dlm = CrossFitSmoothedProfiles{exp_index, exp_index2, time_index, idx};
                if ~isempty(dlm)
                   CrossFitMeanSlopes(exp_index, exp_index2, time_index, idx) = dlm.Coefficients.Estimate(2);
                   if CrossFitMeanSlopes(exp_index, exp_index2, time_index, idx) == 0
                       CrossFitMeanSlopes(exp_index, exp_index2, time_index, idx) = NaN;
                       CrossFitMeanSlopeSEs(exp_index, exp_index2, time_index, idx) = NaN;
                   end
                   CrossFitMeanIntercepts(exp_index, exp_index2, time_index, idx) = dlm.Coefficients.Estimate(1);
                   CrossFitMeanSlopeSEs(exp_index, exp_index2, time_index, idx) = dlm.Coefficients.SE(2);
                   CrossFitMeanInterceptSEs(exp_index, exp_index2, time_index, idx) = dlm.Coefficients.SE(1);
                end
                
            end
        end
    end
end
%%
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
temps = [17 20 22.5 25 27.5];
MarkerStyles = {'o', 's', 'd'};
for idx =  1:length(chList)
for exp_index = 1:NumSets

    ch_index= chList(idx);
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
            ymaxes = zeros(1, NumSets);
            ymins = zeros(1, NumSets);
            for temp_index = 1:length(temps)
                MatchIndices = find(AllSetInfo.Temperatures == temps(temp_index));
                for rep_index = 1:length(MatchIndices)
                    errorbar(times_to_plot, squeeze(CrossFitMeanSlopes(exp_index, MatchIndices(rep_index), :,idx)).',...
                        squeeze(CrossFitMeanSlopeSEs(exp_index, MatchIndices(rep_index), :,idx)).', [MarkerStyles{rep_index}, '-'],...
                        'MarkerSize', 8,'LineWidth', 1.5, 'Color', colorsTemperatures(temp_index,:),...
                        'MarkerFaceColor', colorsTemperatures(temp_index,:),'MarkerEdgeColor', 'k');
                    ymaxes((temp_index-1)*3+rep_index) = max(squeeze(CrossFitMeanSlopes(exp_index, MatchIndices(rep_index), :,idx)).'+ squeeze(CrossFitMeanSlopeSEs(exp_index, MatchIndices(rep_index), :,idx)).');
                    ymins((temp_index-1)*3+rep_index) = min(squeeze(CrossFitMeanSlopes(exp_index, MatchIndices(rep_index), :,idx)).'- squeeze(CrossFitMeanSlopeSEs(exp_index, MatchIndices(rep_index), :,idx)).');
                    hold on 
                end
            end
            plot_ymax = (ceil(max(ymaxes)/0.5))*0.5;
            plot_xmax = (ceil(max(xmaxes)/0.5))*0.5;
            plot_ymin = floor(min(ymins)/0.1)*0.1;
            grid on 
            hold off
            xlabel(['Time into NC14 (min)'], 'FontSize', 18)
            ylabel(['Slope of Intensity Fit'], 'FontSize', 18);
            ylim([plot_ymin,  plot_ymax])
            xlim([5 max(times_to_plot)+5])
%             subplot_handles = plot_handles;
%             SubLegendLabels = LegendLabels;
%             SubLegendLabels = SubLegendLabels(subplot_handles ~= 0);
%             subplot_handles = subplot_handles(subplot_handles ~= 0);
%             legend(subplot_handles,SubLegendLabels, 'FontSize', 18, 'location','northwest');
            DeltaFCAx.YAxis.FontSize = 18;
            DeltaFCAx.XAxis.FontSize = 18;
            
            if AllSetInfo.Replicates(exp_index) ~= 0
            title(DeltaFCAx, [ChannelNames{ch_index},' Fit Slopes to ', num2str(AllSetInfo.Temperatures(exp_index)), 'ºC Replicate ' num2str(AllSetInfo.Replicates(exp_index))], 'FontSize', 20)
            else
                title(DeltaFCAx, [ChannelNames{ch_index},' Fit Slopes to ', num2str(AllSetInfo.Temperatures(exp_index)), 'ºC Flipped Stain'], 'FontSize', 20);
            end
            DeltaFCAx.FontSize = 18;
            if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep,'NormalizedSmoothedTestProfiles'])
                mkdir([AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedSmoothedTestProfiles'])
            end

            outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFigure1D', filesep, 'NormalizedSmoothedTestProfiles', filesep,...
                ChannelNames{ch_index}, '_MasterSet',num2str(master_index),'_Set', num2str(exp_index),'_AllWideTimesAndTemperatures.png' ];
      
            saveas(DeltaFCFixCorrectedFig,outpath2);
            
            
end
end
   
