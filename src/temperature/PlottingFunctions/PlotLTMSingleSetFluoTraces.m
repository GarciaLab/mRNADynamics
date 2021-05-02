function PlotLTMSingleSetFluoTraces(this, varargin)
%%

% PlotTitle, PlottingColors, UseDifferentColors,
% UseDiffProfiles, UsePhysicalAPLength
IncludeFits = true;
IncludedSets = this.ProcessedExperiments;

x = 1;
while x <= length(varargin)
    if strcmp(lower(varargin{x}), 'plottitle')
        PlotTitle = varargin{x+1};
        x = x+1;
    elseif strcmp(lower(varargin{x}), 'plottingcolors')
        PlottingColors = varargin{x+1};
        x = x+1;
    elseif strcmp(lower(varargin{x}), 'tracetype')
        TraceType = lower(varargin{x+1});
        x = x+1;
    elseif strcmp(lower(varargin{x}), 'excludefits')
        IncludeFits = false;
    elseif strcmp(lower(varargin{x}), 'includedsets')
        IncludedSets = lower(varargin{x+1});
        x = x+1;
    end
    x = x+1;
end

if ~exist('PlottingColors', 'var')
    PlottingColors = 'default';
elseif ~strcmp(lower(PlottingColors), 'default')  & ~strcmp(lower(PlottingColors), 'pboc')
    error('Invalid choice of plotting colors. Can use either "default", "pboc", or "gradient".') % change to error
end
if ~exist('TraceType', 'var')
    TraceType = 'AnaphaseAligned';
elseif  ~strcmp(lower(TraceType), 'tbinned')  & ~strcmp(lower(TraceType), 'anaphasealigned') &...
        ~strcmp(lower(TraceType), 'anaphasealigned3d') &  ~strcmp(lower(TraceType), 'tbinned3d')  & ...
        ~strcmp(lower(TraceType), 'fluo') &  ~strcmp(lower(TraceType), 'fluo3d')  & ...
        ~strcmp(lower(TraceType), 'unaligned') &  ~strcmp(lower(TraceType), 'unaligned3d')
    error('Invalid choice of trace type. Can use either "tbinned" or "anaphasealigned".') % change to error
end

if strcmpi(TraceType, 'anaphasealigned')
    TraceType = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'anaphasealigned3d')
    TraceType = 'AnaphaseAligned3D';
elseif strcmpi(TraceType, 'fluo') | strcmpi(TraceType, 'unaligned')
    TraceType = 'Unaligned';
elseif strcmpi(TraceType, 'fluo3d') | strcmpi(TraceType, 'unaligned3d')
    TraceType = 'Unaligned3D';
elseif strcmpi(TraceType, 'tbinned')
    TraceType = 'Tbinned';
elseif strcmpi(TraceType, 'tbinned3d')
    TraceType = 'Tbinned3D';
end

%%

if IncludeFits
    FitString = 'Fitted';
else
    FitString = '';
end
Temp_obs = this.Temp_obs;
Temp_sp = this.Temp_sps;


if strcmp(lower(PlottingColors), 'default')
    [~, colors] = getColorPalettes();
elseif strcmp(lower(PlottingColors), 'pboc')
    [colors, ~] = getColorPalettes();
end
%%
timeSubDir = 'MeanTraces';
Nsigfigs = 3;
%%
NumSets = length(this.ExperimentPrefixes);
temperatures = flip(unique(this.Temp_sps(this.ProcessedExperiments)));
NumTemperatures = length(temperatures);


APResolution = this.Experiments{1}.APResolution;
APbins = 0:APResolution:1;
NumAPbins = length(APbins);


legend_labels = this.LegendLabels;
MinimumTraceCount = this.MinimumTraceCount;

for SetIndex=1:NumSets%length(temperatures)
    if ~ismember(SetIndex, this.ProcessedExperiments)
        continue
    end
    resultsFolder = this.Experiments{SetIndex}.resultsFolder;
    
    outdir2 = [resultsFolder, filesep, timeSubDir];
    
    if ~exist(outdir2, 'dir')
        mkdir(outdir2)
    end
    
    
    
    
    CurrentTemperature = this.Temp_sps(SetIndex);
    
    for NC = this.IncludedNCs

        % Prepare Traces for plotting
        
        
        
        NCTimes = this.MeanProfiles{SetIndex}.([TraceType, 'CycleFrameTimes']){NC-8};
        IncludedRows = 1:length(NCTimes);
        MeanFluoMat = squeeze(this.MeanProfiles{SetIndex}.([TraceType, 'CycleMeanTraces'])(IncludedRows,:,NC-8));
        StdFluoMat = squeeze(this.MeanProfiles{SetIndex}.([TraceType, 'CycleTraceStdErrors'])(IncludedRows,:,NC-8));
        NumNucMat = squeeze(this.MeanProfiles{SetIndex}.([TraceType, 'CycleNumOnNuclei'])(IncludedRows,:,NC-8));
        
        IncludedRows = find(sum(~isnan(MeanFluoMat),2).' > 0);
        if isempty(IncludedRows)
            continue
        end
        
        IncludedColumns = find(sum(~isnan(MeanFluoMat),1).' > 0);
        if ~isempty(IncludedColumns)
            MinAPbin = min(IncludedColumns);
            MaxAPbin = max(IncludedColumns);
        end
        MeanFluoMat = MeanFluoMat(IncludedRows,:);
        StdFluoMat = StdFluoMat(IncludedRows,:);
        NumNucMat = NumNucMat(IncludedRows,:);
        
        NCTimes = NCTimes(IncludedRows)/60;
        
        MaximumNCTime = max(NCTimes);
        MaxFluo = max(max(MeanFluoMat+StdFluoMat));
        MinFluo = min(min(MeanFluoMat+StdFluoMat));
        NumFrames = length(NCTimes);
        
        
        
        
        
       
        
        outdir3 = [outdir2, filesep, 'NC', num2str(NC)];
        
        if ~exist(outdir3, 'dir')
            mkdir(outdir3)
        end
        
        
        outdir4 = [outdir3, filesep, TraceType];
        
        
        if ~exist(outdir4, 'dir')
            mkdir(outdir4)
        end
        
        
        outdir5 = [outdir4, filesep, datestr(now, 'yyyymmdd')];
        
        if ~exist(outdir5, 'dir')
            mkdir(outdir5)
        end
        close all
        
        
        FrameProfFig = figure(1);
        set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
        set(gcf,'color','w');
        FrameProfAx = axes(FrameProfFig);
        
        
        
        temp_idx = find(temperatures == CurrentTemperature);
        if IncludeFits
            ci_plotline = fill([0, 1, 1, 0], [0, 0, max(MaxFluo*1.2,1), max(MaxFluo*1.2,1)], colors(temp_idx,:));
            ci_plotline.FaceAlpha = 0.2;
            set(ci_plotline, 'EdgeColor', 'none');
            hold on
            t_on_plotline = xline(0,'-', 'color', [0,0,0]);
            t_off_plotline = xline(0, '--',  'color', [0,0,0]);
            t_peak_plotline = xline(0, '-.',  'color', [0,0,0]);
            %t_peak_plotline = xline(0, '--',  'color', [0,0,0]+0.5);
            
            
            fitted_prof = plot(APbins, ones(1, length(APbins)), '-', 'Color', colors(temp_idx,:));
            set(fitted_prof,'Visible','off'); %'off' or 'on'
            set(t_on_plotline,'Visible','off'); %'off' or 'on'
            set(t_off_plotline,'Visible','off'); %'off' or 'on'
            set(t_peak_plotline,'Visible','off'); %'off' or 'on'
            set(ci_plotline,'Visible','off'); %'off' or 'on'
        end
        
       
        
        eb = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
        hold on
        set(eb, 'color', colors(temp_idx,:), 'LineWidth', 1);
        set(get(get(eb, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        %set(get(get(eb, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        
        
        prof = plot(APbins, ones(1, length(APbins)),'o',...
                'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colors(temp_idx,:),...
                'MarkerSize', 6, 'LineStyle', '-', 'Color', colors(temp_idx,:));
        
        set(eb,'Visible','off'); %'off' or 'on'
        set(prof,'Visible','off'); %'off' or 'on'
        
        
        hold off
        if strcmp(lower(TraceType), 'anaphasealigned') | strcmp(lower(TraceType), 'anaphasealigned3d')
            xlabel('Time since anaphase (min)')
        else
            xlabel('Time since NC start (min)')
        end
        
       xlim([0-0.01* max(MaximumNCTime), max(MaximumNCTime)*1.1]);
        
        if strcmp(lower(TraceType), 'unaligned') | strcmp(lower(TraceType), 'anaphasealigned')  | strcmp(lower(TraceType), 'tbinned')
            ylabel('Fluo (AU)')
        else
            ylabel('3D Fluo (AU)')
        end
        ylim([max(0, min(MinFluo*0.95)), max(MaxFluo*1.05,1)])
        %
        
        title(FrameProfAx, {'',...
            ['Nuclear Cycle ',num2str(NC),', Fraction Embryo Length: ',num2str(-1) ]})
        
        
        
        if (NC == 14) & ~IncludeFits
            legend(FrameProfAx, legend_labels(SetIndex), 'Location', 'northeast')
        elseif (NC < 14) & ~IncludeFits
            legend(FrameProfAx, legend_labels(SetIndex), 'Location', 'northwest')
        end
        for i = MinAPbin:MaxAPbin
            APBinHasData = false;
            set(FrameProfAx.Children(2),'Visible','off'); %'off' or 'on'
            set(FrameProfAx.Children(1),'Visible','off'); %'off' or 'on'
            set(get(get(FrameProfAx.Children(1), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            
            
            use_idx = NumNucMat(:,i) >= this.MinimumTraceCount;
            
            
            if sum(use_idx) == 0 %| sum(DiffMeanFluoMat(i, use_idx, j) == 0)
                FrameProfAx.Children(1).XData = NCTimes(use_idx);
                FrameProfAx.Children(1).YData = zeros(1, length(NCTimes(use_idx)));
                FrameProfAx.Children(2).XData = NCTimes(use_idx);
                FrameProfAx.Children(2).YData = zeros(1, length(NCTimes(use_idx)));
                FrameProfAx.Children(2).YPositiveDelta = zeros(1, length(NCTimes(use_idx)));
                FrameProfAx.Children(2).YNegativeDelta = zeros(1, length(NCTimes(use_idx)));
                set(FrameProfAx.Children(1),'Visible','off'); %'off' or 'on'
                set(FrameProfAx.Children(2),'Visible','off'); %'off' or 'on'
                set(get(get(FrameProfAx.Children(1), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            else
                APBinHasData = true;
                set(get(get(FrameProfAx.Children(1), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
                set(FrameProfAx.Children(1),'Visible','on'); %'off' or 'on'
                set(FrameProfAx.Children(2),'Visible','on'); %'off' or 'on'
                if IncludeFits
                    plot_legend_labels = cell(1, 5);
                    
                    plot_legend_labels{1} = legend_labels{SetIndex};
                    
                    [t_vector, fit_solution, ci, pos_slope, se_pos_slope, neg_slope,...
                        se_neg_slope, time_on, se_time_on, time_off,...
                        se_time_off, time_peak, se_time_peak, R2] =...
                        getFittedTrace(this, SetIndex, i, NC, TraceType, MaximumNCTime);
                    if ~isempty(t_vector)
                        
                        FrameProfAx.Children(3).YData = fit_solution;
                        FrameProfAx.Children(3).XData =t_vector;
                        set(FrameProfAx.Children(3),'Visible','on'); %'off' or 'on'
                        if isnan(neg_slope) | (NC == 14)
                            lab2 = MeanSE_num2str(pos_slope, se_pos_slope, Nsigfigs);
                            plot_legend_labels{2} = ['Loading Rate: ', lab2.m,...
                                ' \pm ', lab2.se, ' AU/min, R^2: ', num2str(R2, 3)];
                    
                            
                        else
                            lab2a = MeanSE_num2str(pos_slope, se_pos_slope, Nsigfigs);
                            lab2b = MeanSE_num2str(neg_slope, se_neg_slope, Nsigfigs);
                            plot_legend_labels{2} = ['Loading Rate: ', lab2a.m,...
                                ' \pm ', lab2a.se, ' AU/min, R^2: ', num2str(R2, 3),'\newlineUnloading Rate: ',...
                                lab2b.m,' \pm ', lab2b.se, ' AU/min'];
                        end
                        if ~isempty(ci)
                            curve1 = ci(:,1).';
                            curve2 = ci(:,2).';
                            curve1(curve1 > max(MaxFluo*1.05,1)) = max(MaxFluo*1.05,1);
                            curve1(curve1 < 0) = 0;
                            curve2(curve2 > max(MaxFluo*1.05,1)) = max(MaxFluo*1.05,1);
                            curve2(curve2 < 0) = 0;
                            t_vector2 = [t_vector, fliplr(t_vector)];
                            inBetween = [curve1, fliplr(curve2)];
                            FrameProfAx.Children(7).Faces= 1:length(t_vector2);
                            FrameProfAx.Children(7).Vertices = [t_vector2; inBetween].';
                            set(FrameProfAx.Children(7),'Visible','on');
                            plot_legend_labels{3} = 'Prediction 95% Confidence Interval';
                        else
                            set(FrameProfAx.Children(7),'Visible','off');
                            plot_legend_labels{3} = '';
                        end
                        
               
                        
                       
                    else
                        FrameProfAx.Children(3).XData = NCTimes(use_idx);
                        FrameProfAx.Children(3).YData = zeros(1, length(NCTimes(use_idx)));
                        set(FrameProfAx.Children(3),'Visible','off'); %'off' or 'on'
                        plot_legend_labels{2} = 'No Fitting Info';
                        set(FrameProfAx.Children(7),'Visible','off');
                        plot_legend_labels{3} = '';
                    end
                    if ~isnan(time_peak) & ~isnan(time_on)
%                         FrameProfAx.Children(4).Vertices(:, 1) = [time_on; time_peak; time_peak; time_on];
%                         set(FrameProfAx.Children(4),'Visible','on');
%                         se_time_elongation = sqrt(se_time_peak^2+se_time_on^2);
%                         lab3 = MeanSE_num2str(time_peak-time_on, se_time_elongation, Nsigfigs);
%                         plot_legend_labels{5} = ['t_{elongation}: ', lab3.m, ' \pm ', lab3.se,...
%                             ' min'];%, t_{elongation}: ', num2str(time_peak-time_on, 3), ' min'];
                        FrameProfAx.Children(4).Value = time_peak;
                        set(FrameProfAx.Children(4),'Visible','on');
                        lab3 = MeanSE_num2str(time_peak, se_time_peak, Nsigfigs);
                        time_elongation = time_peak -time_on;
                        se_time_elongation = sqrt(se_time_peak^2 + se_time_on^2);
                        lab8 = MeanSE_num2str(time_elongation, se_time_elongation, Nsigfigs);
                        plot_legend_labels{6} = ['t_{peak}: ', lab3.m,' \pm ', lab3.se,...
                            ' min\newlinet_{elongation}: ', lab8.m,' \pm ', lab8.se, ' min'];
                       
                    else
                        FrameProfAx.Children(4).Value = 0;
                        set(FrameProfAx.Children(4),'Visible','off');
                        plot_legend_labels{6} = 'No t_{peak} info';
                    end
                    
                    if ~isnan(time_off)
                        FrameProfAx.Children(5).Value = time_off;
                        set(FrameProfAx.Children(5),'Visible','on');
                        lab4 = MeanSE_num2str(time_off, se_time_off, Nsigfigs);
                        plot_legend_labels{5} = ['t_{off}: ', lab4.m,' \pm ', lab4.se,  ' min'];
                    else
                        FrameProfAx.Children(5).Value = 0;
                        set(FrameProfAx.Children(5),'Visible','off');
                        plot_legend_labels{5} = 'No t_{off} info';
                    end
                    
                    if ~isnan(time_on)
                        FrameProfAx.Children(6).Value = time_on;
                        set(FrameProfAx.Children(6),'Visible','on');
                        lab5 = MeanSE_num2str(time_on, se_time_on, Nsigfigs);
                        plot_legend_labels{4} =  ['t_{on}: ', lab5.m,...
                            ' \pm ', lab5.se, ' min'];
                    else
                        FrameProfAx.Children(6).Value = 0;
                        set(FrameProfAx.Children(6),'Visible','off');
                        plot_legend_labels{4} = 'No t_{on} info';
                    end
                    
                    
                    
                    
                end
                
                FrameProfAx.Children(1).YData = MeanFluoMat(use_idx, i);
                FrameProfAx.Children(1).XData = NCTimes(use_idx);
                FrameProfAx.Children(2).YData = MeanFluoMat(use_idx, i);
                FrameProfAx.Children(2).XData = NCTimes(use_idx);
                FrameProfAx.Children(2).YPositiveDelta = StdFluoMat(use_idx, i);
                FrameProfAx.Children(2).YNegativeDelta  = StdFluoMat(use_idx, i);
                
                set(FrameProfAx.Children(1),'Visible','on'); %'off' or 'on'
                set(FrameProfAx.Children(2),'Visible','on'); %'off' or 'on'
                
            end
            
            %try
            if exist('PlotTitle', 'var')
                title(FrameProfAx, {PlotTitle,...
                    ['Nuclear Cycle ',num2str(NC),', Fraction Embryo Length: ',num2str(APbins(i)) ]})
                
                
            else
                title(FrameProfAx, ['Nuclear Cycle ',num2str(NC),', Fraction Embryo Length: ',num2str(APbins(i)) ])
                
                
            end
            %end
            if ~APBinHasData
                continue
            end
            if (NC == 14) & IncludeFits
                plot_handles ={prof, fitted_prof, ci_plotline, t_on_plotline, t_off_plotline,  t_peak_plotline};
                hlegend = legend([plot_handles{:}], plot_legend_labels, 'Location', 'eastoutside',...
                    'FontSize', 10);
                
            elseif (NC < 14) & IncludeFits
                plot_handles ={prof, fitted_prof, ci_plotline, t_on_plotline, t_off_plotline,  t_peak_plotline};
                hlegend = legend([plot_handles{:}], plot_legend_labels, 'Location', 'eastoutside',...
                    'FontSize', 10);
            end
            
            saveas(FrameProfFig,[outdir5, filesep,...
                FitString,'FluoTrace_NC',num2str(NC),'_Bin', num2str(i),'.png']);
            
        end
    end
end
close all
end