function GeneralizedScatterPlotLTMEmbryoStats(this, parameter, outdir, varargin)
%%

% PlotTitle, PlottingColors, UseDifferentColors,
% UseDiffProfiles, UsePhysicalAPLength
SkipParamsVsTemp = false;
SkipBinnedParamsVsTemp = false;
SkipAllSubplots = false;
SkipAllBinnedSubplots = false;
UsePerNucleusTraces = false;
UseBinnedTraces = false;
UseBinnedPerNucleusTraces = false;


x = 1;
while x <= length(varargin)
    if strcmp(lower(varargin{x}), 'plottitle')
        PlotTitle = varargin{x+1};
        x = x+1;
    elseif strcmp(lower(varargin{x}), 'plottingcolors')
        PlottingColors = varargin{x+1};
        x = x+1;
    elseif strcmpi(varargin{x}, 'SkipParamsVsTemp')
        SkipParamsVsTemp = true;
    elseif strcmpi(varargin{x}, 'SkipBinnedParamsVsTemp')
        SkipBinnedParamsVsTemp = true;
    elseif strcmpi(varargin{x}, 'SkipAllSubplots')
        SkipAllSubplots = true;
    elseif strcmpi(varargin{x}, 'SkipAllBinnedSubplots')
        SkipAllBinnedSubplots = true;
    end
    x = x+1;
end


if ~exist('PlottingColors', 'var')
    PlottingColors = 'default';
elseif ~strcmpi(PlottingColors, 'gradient') &~strcmp(lower(PlottingColors), 'default')  & ~strcmp(lower(PlottingColors), 'pboc')
    error('Invalid choice of plotting colors. Can use either "default", "pboc", or "gradient".') % change to error
end

if ~exist(outdir, 'dir')
    mkdir(outdir)
end

Temp_obs = this.Temp_obs;
Temp_sp = this.Temp_sps;

if strcmp(lower(PlottingColors), 'default')
    [~, colors] = getColorPalettes();
    GradString = '';
elseif strcmp(lower(PlottingColors), 'pboc')
    [colors, ~] = getColorPalettes();
    GradString = '';
else
    
    Temp_range = 15:0.1:max(Temp_obs);
    colors = jet(length(Temp_range));
    FractionalTempRange = (Temp_range-min(Temp_range))/(max(Temp_range)-min(Temp_range));
    GradString = 'Gradient';
end

MarkerStyles = {'o', 'd', 's', '>', '^', '*', 'x', 'p', 'o', 'd', 's', '>', '^', '*', 'x', 'p'};




%%

R = 8.314*10^(-3); % kJ * K^(-1)*mol^(-1)
NumSets = length(this.ExperimentPrefixes);

temperatures = flip(unique(this.Temp_sps(this.ProcessedExperiments)));
NumTemperatures = length(temperatures);
TempMatches = cell(1, NumTemperatures);
UseSet = ismember(1:NumSets, this.ProcessedExperiments);
for t_index = 1:NumTemperatures
    TempMatches{t_index} = find((this.Temp_sps == temperatures(t_index)) & UseSet);
end


legend_labels = this.LegendLabels;



%% Load relevant parameters into memory
[PlottedParams, PlottedParamSEs,ylab,OutputString,GlobalPlotYmax,GlobalPlotYmin,LogPlotYmin, IncludeLogPlots] = ...
    getEmbryoStatsPlottingVariables(this, parameter);
[BinnedParams, BinnedSEParams, Counts, ParamTemperatures, ParamSETemperatures] = ...
    getBinnedPlottingVariables(this, PlottedParams, PlottedParamSEs);

PlottedHasSEs = ~all(all(isnan(PlottedParams)));

%% Calculate x  limits for 2nd subplot


TemperatureVector = 1./(R*(this.Temp_obs + 273));
Subplot2DataXmin = min(TemperatureVector);
Subplot2DataXmax = max(TemperatureVector);
Subplot2Xspan = Subplot2DataXmax-Subplot2DataXmin;
Plot2Xmin = Subplot2DataXmin - Subplot2Xspan*.05;
Plot2Xmax = Subplot2DataXmax + Subplot2Xspan*.05;

if ~SkipParamsVsTemp
    outdir2 = [outdir,filesep,OutputString];
    if ~exist(outdir2, 'dir')
        mkdir(outdir2)
    end
    
    outdir3 = [outdir2,filesep, datestr(now, 'yyyymmdd')];
    if ~exist(outdir3, 'dir')
        mkdir(outdir3)
    end
    close all
    
    eb = cell(1, NumSets);
    prof = cell(1, NumSets);
    FrameProfFig = figure(1);
    set(gcf,'color','w');
    if NumSets > 25
        set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .95, .6]);
    else
        set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .8, .6]);
    end
    if IncludeLogPlots
        FrameProfAx = subplot(1, 3, 1, gca);
        if NumSets > 25
            pos1 = get(FrameProfAx, 'Position');
            posnew = pos1;
            posnew(1) = posnew(1) -0.05;
            set(FrameProfAx, 'Position', posnew)
        end
    else
        FrameProfAx = axes(FrameProfFig);
    end
    
    for SetIndex =1:NumSets
        temp_idx = find(temperatures == Temp_sp(SetIndex));
        marker_idx = find(TempMatches{temp_idx} == SetIndex);
        eb{SetIndex} = errorbar([0], [0], [0], 'vertical', 'LineStyle', 'none');
        hold on
        set(eb{SetIndex}, 'color', colors(temp_idx,:), 'LineWidth', 1);
        set(get(get(eb{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        set(eb{SetIndex},'Visible','off'); %'off' or 'on'
        
        
        if marker_idx <= length(MarkerStyles)/2
            prof{SetIndex} = plot([0], [0], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
                'linestyle', 'none');
        else
            prof{SetIndex} = plot([0], [0], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                'MarkerEdgeColor', 'k','MarkerFaceColor', colors(temp_idx,:),...
                'linestyle', 'none');
        end
        set(prof{SetIndex},'Visible','off'); %'off' or 'on'
    end
    
    hold off
    
    xlabel('Temperature (°C)')
    xlim([15, 30])
    
    ylabel(ylab)
    
    ylim([GlobalPlotYmin, GlobalPlotYmax])
    
    if IncludeLogPlots
        eb2 = cell(1, NumSets);
        prof2 = cell(1, NumSets);
        FrameProfAx2 = subplot(1, 3, 2);
        if NumSets > 25
            pos2 = get(FrameProfAx2, 'Position');
            posnew = pos2;
            posnew(1) = posnew(1) -0.05;
            set(FrameProfAx2, 'Position', posnew)
        end
        for SetIndex =1:NumSets
            temp_idx = find(temperatures == Temp_sp(SetIndex));
            marker_idx = find(TempMatches{temp_idx} == SetIndex);
            
            eb2{SetIndex} = errorbar([0], [0], [0], 'vertical', 'LineStyle', 'none');
            hold on
            set(eb2{SetIndex}, 'color', colors(temp_idx,:), 'LineWidth', 1);
            set(get(get(eb2{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            set(eb2{SetIndex},'Visible','off'); %'off' or 'on'
            if marker_idx <= length(MarkerStyles)/2
                prof2{SetIndex} = plot([0], [0], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                    'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
                    'linestyle', 'none');
            else
                prof2{SetIndex} = plot([0], [0], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                    'MarkerEdgeColor', 'k','MarkerFaceColor', colors(temp_idx,:),...
                    'linestyle', 'none');
            end
            set(prof2{SetIndex},'Visible','off'); %'off' or 'on'
        end
        
        hold off
        
        xlabel('1/(RT) (mol/kJ)')
        xlim([Plot2Xmin, Plot2Xmax])
        
        
        
        ylabel(ylab)
        ylim([LogPlotYmin, GlobalPlotYmax])
        
        set(FrameProfAx2, 'YScale', 'log')
    end
    
    for nc_idx=1:length(this.IncludedNCs)
        NC = this.IncludedNCs(nc_idx);
        AllNCParams = PlottedParams(:, NC-8).';
        if all(isnan(AllNCParams))
            continue
        end
        
        if PlottedHasSEs
            AllNCParamSEs = PlottedParamSEs(:, NC-8).';
            TempSEs = AllNCParamSEs;
            TempSEs(isnan(TempSEs)) = 0;
            MaxYvalue = max(AllNCParams+TempSEs);
        else
            MaxYvalue = max(AllNCParams);
        end
        
        PlottedSets = zeros(1, NumSets, 'logical');
        if strcmpi(parameter, 'NCDivisionTimes') | strcmpi(parameter, 'NCDivisionDurations') | strcmpi(parameter, 'NCDivisions')
            set(FrameProfAx, 'ylim', [GlobalPlotYmin, MaxYvalue*1.05]);
            if IncludeLogPlots
                set(FrameProfAx2, 'ylim', [LogPlotYmin, MaxYvalue*1.05]);
            end
        end
        for SetIndex = 1:NumSets
            if ~ismember(SetIndex, this.ProcessedExperiments)
                continue
            end
            if isnan(AllNCParams(SetIndex))
                
                FrameProfAx.Children(end-(2*(SetIndex-1)+1)).XData = 25;
                FrameProfAx.Children(end-(2*(SetIndex-1)+1)).YData = 1;
                FrameProfAx.Children(end-(2*(SetIndex-1))).XData = 25;
                FrameProfAx.Children(end-(2*(SetIndex-1))).YData = 1;
                FrameProfAx.Children(end-(2*(SetIndex-1))).YPositiveDelta = 1;
                FrameProfAx.Children(end-(2*(SetIndex-1))).YNegativeDelta = 1;
                set(FrameProfAx.Children(end-(2*(SetIndex-1)+1)),'Visible','off'); %'off' or 'on'
                set(FrameProfAx.Children(end-(2*(SetIndex-1))),'Visible','off'); %'off' or 'on'
                set(get(get(prof{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                if IncludeLogPlots
                    FrameProfAx2.Children(end-(2*(SetIndex-1)+1)).XData = 1/(R*(25+273));
                    FrameProfAx2.Children(end-(2*(SetIndex-1)+1)).YData = 1;
                    FrameProfAx2.Children(end-(2*(SetIndex-1))).XData = 1/(R*(25+273));
                    FrameProfAx2.Children(end-(2*(SetIndex-1))).YData = 1;
                    FrameProfAx2.Children(end-(2*(SetIndex-1))).YPositiveDelta = 1;
                    FrameProfAx2.Children(end-(2*(SetIndex-1))).YNegativeDelta = 1;
                    set(FrameProfAx2.Children(end-(2*(SetIndex-1)+1)),'Visible','off'); %'off' or 'on'
                    set(FrameProfAx2.Children(end-(2*(SetIndex-1))),'Visible','off'); %'off' or 'on'
                    set(get(get(prof2{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                end
            else
                PlottedSets(SetIndex) = 1;
                FrameProfAx.Children(end-(2*(SetIndex-1)+1)).YData = AllNCParams(SetIndex);
                FrameProfAx.Children(end-(2*(SetIndex-1)+1)).XData =this.Temp_obs(SetIndex);
                set(FrameProfAx.Children(end-(2*(SetIndex-1)+1)),'Visible','on'); %'off' or 'on'
                set(get(get(prof{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
                if PlottedHasSEs
                    FrameProfAx.Children(end-(2*(SetIndex-1))).YData =  AllNCParams(SetIndex);
                    FrameProfAx.Children(end-(2*(SetIndex-1))).XData = this.Temp_obs(SetIndex);
                    FrameProfAx.Children(end-(2*(SetIndex-1))).YPositiveDelta = AllNCParamSEs(SetIndex);
                    FrameProfAx.Children(end-(2*(SetIndex-1))).YNegativeDelta  = AllNCParamSEs(SetIndex);
                    set(FrameProfAx.Children(end-(2*(SetIndex-1))),'Visible','on'); %'off' or 'on'
                    
                else
                    set(FrameProfAx.Children(end-(2*(SetIndex-1))),'Visible','off'); %'off' or 'on'
                end
                if IncludeLogPlots
                    FrameProfAx2.Children(end-(2*(SetIndex-1)+1)).YData = AllNCParams(SetIndex);
                    FrameProfAx2.Children(end-(2*(SetIndex-1)+1)).XData =1/(R*(this.Temp_obs(SetIndex)+273));
                    set(FrameProfAx2.Children(end-(2*(SetIndex-1)+1)),'Visible','on'); %'off' or 'on'
                    set(get(get(prof2{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
                    if PlottedHasSEs
                        FrameProfAx2.Children(end-(2*(SetIndex-1))).YData =  AllNCParams(SetIndex);
                        FrameProfAx2.Children(end-(2*(SetIndex-1))).XData = 1/(R*(this.Temp_obs(SetIndex)+273));
                        FrameProfAx2.Children(end-(2*(SetIndex-1))).YPositiveDelta = AllNCParamSEs(SetIndex);
                        FrameProfAx2.Children(end-(2*(SetIndex-1))).YNegativeDelta  = AllNCParamSEs(SetIndex);
                        set(FrameProfAx2.Children(end-(2*(SetIndex-1))),'Visible','on'); %'off' or 'on'
                    else
                        set(FrameProfAx2.Children(end-(2*(SetIndex-1))),'Visible','off'); %'off' or 'on'
                    end
                end
            end
        end
        
        
        if all(~PlottedSets)
            continue
        end
        
        if IncludeLogPlots
            LegendAx = subplot(1, 3, 3);
            legend_profs = cell(1, NumSets);
            legend_labels2 = {};
            for SetIndex = 1:NumSets
                temp_idx = find(temperatures == Temp_sp(SetIndex));
                marker_idx = find(TempMatches{temp_idx} == SetIndex);
                if marker_idx <= length(MarkerStyles)/2
                    legend_profs{SetIndex} = plot([0], [0], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                        'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
                        'linestyle', 'none');
                else
                    legend_profs{SetIndex} = plot([0], [0], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                        'MarkerEdgeColor', 'k','MarkerFaceColor', colors(temp_idx,:),...
                        'linestyle', 'none');
                end
                hold on
                if ~PlottedSets(SetIndex)
                    set(get(get(legend_profs{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                else
                    legend_labels2{1, length(legend_labels2)+1} = legend_labels{SetIndex};
                end
            end
            hold off
            axis off
            
            if NumSets > 25
                hlegend = legend(LegendAx, legend_labels2,...
                    'FontSize', 10, 'NumColumns', 2);%,
            else
                hlegend = legend(LegendAx, legend_labels2,...
                    'FontSize', 10);%, 'Location', );
            end
            hlegend.Position(1) = LegendAx.Position(1)+ LegendAx.Position(3)/2 - hlegend.Position(3)/2;
            hlegend.Position(2) = LegendAx.Position(2)+ LegendAx.Position(4)/2 - hlegend.Position(4)/2;
            
            %try
        else
            legend_labels2 = {};
            for SetIndex = 1:NumSets
                if ~PlottedSets(SetIndex)
                    set(get(get(prof{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                else
                    legend_labels2{1, length(legend_labels2)+1} = legend_labels{SetIndex};
                end
            end
            if NumSets > 25
                hlegend = legend(legend_labels2,...
                    'FontSize', 10, 'Location', 'eastoutside', 'NumColumns', 2);
            else
                hlegend = legend(legend_labels2,...
                    'FontSize', 10, 'Location', 'eastoutside');
            end
            
        end
        
        if exist('PlotTitle', 'var')
            sgtitle([PlotTitle,'Nuclear Cycle ',num2str(NC)])
        else
            sgtitle(['Nuclear Cycle ',num2str(NC)])
        end
        
        
        saveas(FrameProfFig,[outdir3, filesep,OutputString,'_NC',num2str(NC), '.png']);
        
        
    end
    close all
end


%% Single Temperature version


TemperatureVector = 1./(R*(this.Temp_obs + 273));
Subplot2DataXmin = min(TemperatureVector);
Subplot2DataXmax = max(TemperatureVector);
Subplot2Xspan = Subplot2DataXmax-Subplot2DataXmin;
Plot2Xmin = Subplot2DataXmin - Subplot2Xspan*.05;
Plot2Xmax = Subplot2DataXmax + Subplot2Xspan*.05;

if ~SkipParamsVsTemp
    outdir2 = [outdir,filesep,OutputString];
    if ~exist(outdir2, 'dir')
        mkdir(outdir2)
    end
    
    outdir3 = [outdir2,filesep, datestr(now, 'yyyymmdd')];
    if ~exist(outdir3, 'dir')
        mkdir(outdir3)
    end
    
    
    for temp_idx = 1:NumTemperatures
        close all
        CurrentTemp = temperatures(temp_idx);
        NumTempSets = sum(this.Temp_sps == CurrentTemp);
        TempSetList = find(this.Temp_sps == CurrentTemp);
        eb = cell(1, NumTempSets);
        prof = cell(1, NumTempSets);
        FrameProfFig = figure(1);
        set(gcf,'color','w');
        if NumTempSets > 25
            set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .95, .6]);
        else
            set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .8, .6]);
        end
        if IncludeLogPlots
            FrameProfAx = subplot(1, 3, 1, gca);
            if NumTempSets > 25
                pos1 = get(FrameProfAx, 'Position');
                posnew = pos1;
                posnew(1) = posnew(1) -0.05;
                set(FrameProfAx, 'Position', posnew)
            end
        else
            FrameProfAx = axes(FrameProfFig);
        end
        
        for SetIndex =1:NumTempSets
            marker_idx =  SetIndex;
            eb{SetIndex} = errorbar([0], [0], [0], 'vertical', 'LineStyle', 'none');
            hold on
            set(eb{SetIndex}, 'color', colors(temp_idx,:), 'LineWidth', 1);
            set(get(get(eb{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            set(eb{SetIndex},'Visible','off'); %'off' or 'on'
            
            
            if marker_idx <= length(MarkerStyles)/2
                prof{SetIndex} = plot([0], [0], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                    'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
                    'linestyle', 'none');
            else
                prof{SetIndex} = plot([0], [0], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                    'MarkerEdgeColor', 'k','MarkerFaceColor', colors(temp_idx,:),...
                    'linestyle', 'none');
            end
            set(prof{SetIndex},'Visible','off'); %'off' or 'on'
        end
        grid on
        hold off
        
        xlabel('Temperature (°C)')
        xlim([CurrentTemp-1, CurrentTemp+1])
        
        ylabel(ylab)
        
        %ylim([GlobalPlotYmin, GlobalPlotYmax])
        
        if IncludeLogPlots
            eb2 = cell(1, NumTempSets);
            prof2 = cell(1, NumTempSets);
            FrameProfAx2 = subplot(1, 3, 2);
            if NumTempSets > 25
                pos2 = get(FrameProfAx2, 'Position');
                posnew = pos2;
                posnew(1) = posnew(1) -0.05;
                set(FrameProfAx2, 'Position', posnew)
            end
            for SetIndex =1:NumTempSets
                marker_idx = SetIndex;
                
                eb2{SetIndex} = errorbar([0], [0], [0], 'vertical', 'LineStyle', 'none');
                hold on
                set(eb2{SetIndex}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                set(get(get(eb2{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                set(eb2{SetIndex},'Visible','off'); %'off' or 'on'
                
                if marker_idx <= length(MarkerStyles)/2
                    prof2{SetIndex} = plot([0], [0], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                        'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
                        'linestyle', 'none');
                else
                    prof2{SetIndex} = plot([0], [0], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                        'MarkerEdgeColor', 'k','MarkerFaceColor', colors(temp_idx,:),...
                        'linestyle', 'none');
                end
                set(prof2{SetIndex},'Visible','off'); %'off' or 'on'
            end
            grid on
            hold off
            
            xlabel('1/(RT) (mol/kJ)')
            %xlim([Plot2Xmin, Plot2Xmax])
            
            
            
            ylabel(ylab)
            %ylim([LogPlotYmin, GlobalPlotYmax])
            
            set(FrameProfAx2, 'YScale', 'log')
        end
        
        for nc_idx=1:length(this.IncludedNCs)
            NC = this.IncludedNCs(nc_idx);
            AllNCParams = PlottedParams(TempSetList, NC-8).';
            if all(isnan(AllNCParams))
                continue
            end
            
            if PlottedHasSEs
                AllNCParamSEs = PlottedParamSEs(TempSetList, NC-8).';
                TempSEs = AllNCParamSEs;
                TempSEs(isnan(TempSEs)) = 0;
                MaxYvalue = max(AllNCParams+TempSEs);
                MinYValue = min(AllNCParams-TempSEs);
            else
                MaxYvalue = max(AllNCParams);
                MinYValue = min(AllNCParams);
            end
            
            PlottedSets = zeros(1, NumTempSets, 'logical');
            if strcmpi(parameter, 'NCDivisionTimes') | strcmpi(parameter, 'NCDivisionDurations') | strcmpi(parameter, 'NCDivisions')
                set(FrameProfAx, 'ylim', [floor(MinYValue)-1,ceil(MaxYvalue)+1]);
                if IncludeLogPlots
                    set(FrameProfAx2, 'ylim', [floor(MinYValue)-1,ceil(MaxYvalue)+1]);
                end
            end
            for SetIndex = 1:NumTempSets
                if ~ismember(TempSetList(SetIndex), this.ProcessedExperiments)
                    continue
                end
                if isnan(AllNCParams(SetIndex))
                    
                    FrameProfAx.Children(end-(2*(SetIndex-1)+1)).XData = 25;
                    FrameProfAx.Children(end-(2*(SetIndex-1)+1)).YData = 1;
                    FrameProfAx.Children(end-(2*(SetIndex-1))).XData = 25;
                    FrameProfAx.Children(end-(2*(SetIndex-1))).YData = 1;
                    FrameProfAx.Children(end-(2*(SetIndex-1))).YPositiveDelta = 1;
                    FrameProfAx.Children(end-(2*(SetIndex-1))).YNegativeDelta = 1;
                    set(FrameProfAx.Children(end-(2*(SetIndex-1)+1)),'Visible','off'); %'off' or 'on'
                    set(FrameProfAx.Children(end-(2*(SetIndex-1))),'Visible','off'); %'off' or 'on'
                    set(get(get(prof{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                    if IncludeLogPlots
                        FrameProfAx2.Children(end-(2*(SetIndex-1)+1)).XData = 1/(R*(25+273));
                        FrameProfAx2.Children(end-(2*(SetIndex-1)+1)).YData = 1;
                        FrameProfAx2.Children(end-(2*(SetIndex-1))).XData = 1/(R*(25+273));
                        FrameProfAx2.Children(end-(2*(SetIndex-1))).YData = 1;
                        FrameProfAx2.Children(end-(2*(SetIndex-1))).YPositiveDelta = 1;
                        FrameProfAx2.Children(end-(2*(SetIndex-1))).YNegativeDelta = 1;
                        set(FrameProfAx2.Children(end-(2*(SetIndex-1)+1)),'Visible','off'); %'off' or 'on'
                        set(FrameProfAx2.Children(end-(2*(SetIndex-1))),'Visible','off'); %'off' or 'on'
                        set(get(get(prof2{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                    end
                else
                    PlottedSets(SetIndex) = 1;
                    FrameProfAx.Children(end-(2*(SetIndex-1)+1)).YData = AllNCParams(SetIndex);
                    FrameProfAx.Children(end-(2*(SetIndex-1)+1)).XData =this.Temp_obs(TempSetList(SetIndex));
                    set(FrameProfAx.Children(end-(2*(SetIndex-1)+1)),'Visible','on'); %'off' or 'on'
                    set(get(get(prof{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
                    if PlottedHasSEs
                        FrameProfAx.Children(end-(2*(SetIndex-1))).YData =  AllNCParams(SetIndex);
                        FrameProfAx.Children(end-(2*(SetIndex-1))).XData = this.Temp_obs(TempSetList(SetIndex));
                        FrameProfAx.Children(end-(2*(SetIndex-1))).YPositiveDelta = AllNCParamSEs(SetIndex);
                        FrameProfAx.Children(end-(2*(SetIndex-1))).YNegativeDelta  = AllNCParamSEs(SetIndex);
                        set(FrameProfAx.Children(end-(2*(SetIndex-1))),'Visible','on'); %'off' or 'on'
                        
                    else
                        set(FrameProfAx.Children(end-(2*(SetIndex-1))),'Visible','off'); %'off' or 'on'
                    end
                    if IncludeLogPlots
                        FrameProfAx2.Children(end-(2*(SetIndex-1)+1)).YData = AllNCParams(SetIndex);
                        FrameProfAx2.Children(end-(2*(SetIndex-1)+1)).XData =1/(R*(this.Temp_obs(TempSetList(SetIndex))+273));
                        set(FrameProfAx2.Children(end-(2*(SetIndex-1)+1)),'Visible','on'); %'off' or 'on'
                        set(get(get(prof2{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
                        if PlottedHasSEs
                            FrameProfAx2.Children(end-(2*(SetIndex-1))).YData =  AllNCParams(SetIndex);
                            FrameProfAx2.Children(end-(2*(SetIndex-1))).XData = 1/(R*(this.Temp_obs(TempSetList(SetIndex))+273));
                            FrameProfAx2.Children(end-(2*(SetIndex-1))).YPositiveDelta = AllNCParamSEs(SetIndex);
                            FrameProfAx2.Children(end-(2*(SetIndex-1))).YNegativeDelta  = AllNCParamSEs(SetIndex);
                            set(FrameProfAx2.Children(end-(2*(SetIndex-1))),'Visible','on'); %'off' or 'on'
                        else
                            set(FrameProfAx2.Children(end-(2*(SetIndex-1))),'Visible','off'); %'off' or 'on'
                        end
                    end
                end
            end
            
            
            if all(~PlottedSets)
                continue
            end
            
            if IncludeLogPlots
                LegendAx = subplot(1, 3, 3);
                legend_profs = cell(1, NumSets);
                legend_labels2 = {};
                for SetIndex = 1:NumTempSets
                    marker_idx = SetIndex;
                    if marker_idx <= length(MarkerStyles)/2
                        legend_profs{SetIndex} = plot([0], [0], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                            'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
                            'linestyle', 'none');
                    else
                        legend_profs{SetIndex} = plot([0], [0], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                            'MarkerEdgeColor', 'k','MarkerFaceColor', colors(temp_idx,:),...
                            'linestyle', 'none');
                    end
                    hold on
                    if ~PlottedSets(SetIndex)
                        set(get(get(legend_profs{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                    else
                        legend_labels2{1, length(legend_labels2)+1} = legend_labels{TempSetList(SetIndex)};
                    end
                end
                hold off
                axis off
                
                if NumTempSets > 25
                    hlegend = legend(LegendAx, legend_labels2,...
                        'FontSize', 10, 'NumColumns', 2);%,
                else
                    hlegend = legend(LegendAx, legend_labels2,...
                        'FontSize', 10);%, 'Location', );
                end
                hlegend.Position(1) = LegendAx.Position(1)+ LegendAx.Position(3)/2 - hlegend.Position(3)/2;
                hlegend.Position(2) = LegendAx.Position(2)+ LegendAx.Position(4)/2 - hlegend.Position(4)/2;
                
                %try
            else
                legend_labels2 = {};
                for SetIndex = 1:NumTempSets
                    if ~PlottedSets(SetIndex)
                        set(get(get(prof{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                    else
                        legend_labels2{1, length(legend_labels2)+1} = legend_labels{TempSetList(SetIndex)};
                    end
                end
                if NumTempSets > 25
                    hlegend = legend(legend_labels2,...
                        'FontSize', 10, 'Location', 'eastoutside', 'NumColumns', 2);
                else
                    hlegend = legend(legend_labels2,...
                        'FontSize', 10, 'Location', 'eastoutside');
                end
                
            end
            
            if exist('PlotTitle', 'var')
                sgtitle([PlotTitle,'T = ',num2str(CurrentTemp),'ºC, Nuclear Cycle ',num2str(NC)])
            else
                sgtitle(['T = ',num2str(CurrentTemp),'ºC, Nuclear Cycle ',num2str(NC)])
            end
            
            
            saveas(FrameProfFig,[outdir3, filesep,OutputString,'T',strrep(num2str(CurrentTemp), '.', '_'),'C_NC',num2str(NC), '.png']);
            
        end
    end
    close all
end


%% Calculate x  limits for 2nd subplot

AllTpoints = reshape(ParamTemperatures, 1,...
    size(ParamTemperatures, 1)*size(ParamTemperatures, 2)*size(ParamTemperatures, 3));
AllSETpoints = reshape(ParamSETemperatures, 1,...
    size(ParamTemperatures, 1)*size(ParamTemperatures, 2)*size(ParamTemperatures, 3));
AllCombinedTpoints = reshape(ParamTemperatures+ParamSETemperatures, 1,...
    size(ParamTemperatures, 1)*size(ParamTemperatures, 2)*size(ParamTemperatures, 3));
AllCombinedTpoints = AllCombinedTpoints(~isnan(AllCombinedTpoints));
AllSETpoints = AllSETpoints(~isnan(AllTpoints));
AllTpoints = AllTpoints(~isnan(AllTpoints));
AllSETpoints(isnan(AllSETpoints)) = 0;
TemperatureVector = 1./(R*(AllTpoints + 273))+sqrt(((1./(R*(AllTpoints + 273).^2)).^2).*AllSETpoints.^2);
Subplot1DataXmin = min(AllCombinedTpoints)*0.95;
Subplot1DataXmax = max(AllCombinedTpoints)*1.05;
Subplot2DataXmin = min(TemperatureVector);
Subplot2DataXmax = max(TemperatureVector);
Subplot2Xspan = Subplot2DataXmax-Subplot2DataXmin;
Plot2Xmin = Subplot2DataXmin - Subplot2Xspan*.05;
Plot2Xmax = Subplot2DataXmax + Subplot2Xspan*.05;

if ~SkipBinnedParamsVsTemp
    outdir2 = [outdir,filesep,OutputString];
    if ~exist(outdir2, 'dir')
        mkdir(outdir2)
    end
    
    outdir3 = [outdir2,filesep, datestr(now, 'yyyymmdd')];
    if ~exist(outdir3, 'dir')
        mkdir(outdir3)
    end
    close all
    
    eb = cell(1, NumTemperatures);
    prof = cell(1, NumTemperatures);
    FrameProfFig = figure(1);
    set(gcf,'color','w');
    set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .8, .5]);
    if IncludeLogPlots
        FrameProfAx = subplot(1, 3, 1, gca);
    else
        FrameProfAx = axes(FrameProfFig);
    end
    for SetIndex =1:NumTemperatures
        temp_idx = SetIndex;
        marker_idx = 1;
        eb{SetIndex} = errorbar([0], [0], [0], [0], [0], [0], 'LineStyle', 'none');
        hold on
        set(eb{SetIndex}, 'color', colors(temp_idx,:), 'LineWidth', 1);
        set(get(get(eb{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        set(eb{SetIndex},'Visible','off'); %'off' or 'on'
        
        prof{SetIndex} = plot([0], [0],MarkerStyles{1},...
            'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
            'linestyle', 'none');
        set(prof{SetIndex},'Visible','off'); %'off' or 'on'
    end
    
    hold off
    
    xlabel('Temperature (°C)')
    xlim([min(15, Subplot1DataXmin), max(30, Subplot1DataXmax)])
    
    ylabel(ylab)
    ylim([GlobalPlotYmin, GlobalPlotYmax])
    
    if IncludeLogPlots
        eb2 = cell(1, NumTemperatures);
        prof2 = cell(1, NumTemperatures);
        FrameProfAx2 = subplot(1, 3, 2);
        
        for SetIndex =1:NumTemperatures
            temp_idx = SetIndex;
            marker_idx = 1;
            
            eb2{SetIndex} = errorbar([0], [0], [0], [0], [0], [0],  'LineStyle', 'none');
            hold on
            set(eb2{SetIndex}, 'color', colors(temp_idx,:), 'LineWidth', 1);
            set(get(get(eb2{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            set(eb2{SetIndex},'Visible','off'); %'off' or 'on'
            prof2{SetIndex} = plot([0], [0], MarkerStyles{1},...
                'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
                'linestyle', 'none');
            set(prof2{SetIndex},'Visible','off'); %'off' or 'on'
        end
        
        hold off
        
        xlabel('1/(RT) (mol/kJ)')
        xlim([Plot2Xmin, Plot2Xmax])
        
        
        
        ylabel(ylab)
        ylim([LogPlotYmin, GlobalPlotYmax])
        
        set(FrameProfAx2, 'YScale', 'log')
    end
    
    for nc_idx=1:length(this.IncludedNCs)
        NC = this.IncludedNCs(nc_idx);
        AllNCParams = BinnedParams(:, NC-8).';
        AllNCParamSEs = BinnedSEParams(:, NC-8).';
        TemperatureCounts  = Counts(:, NC-8).';
        NCParamTemperatures  = ParamTemperatures(:, NC-8).';
        NCParamTemperatureSEs  = ParamSETemperatures(:, NC-8).';
        
        AllNCParams(TemperatureCounts < this.MinimumBinCount) = NaN;
        AllNCParamSEs(TemperatureCounts < this.MinimumBinCount) = NaN;
        NCParamTemperatures(TemperatureCounts < this.MinimumBinCount) = NaN;
        NCParamTemperatureSEs(TemperatureCounts < this.MinimumBinCount) = NaN;
        
        
        TemporarySEs = AllNCParamSEs;
        TemporarySEs(isnan(TemporarySEs)) = 0;
        MaxYvalue = max(AllNCParams+TemporarySEs);
        
        if all(isnan(AllNCParams))
            continue
        end
        
        if strcmpi(parameter, 'NCDivisionTimes') | strcmpi(parameter, 'NCDivisionDurations') | strcmpi(parameter, 'NCDivisions')
            set(FrameProfAx, 'ylim', [GlobalPlotYmin, MaxYvalue*1.05]);
            if IncludeLogPlots
                set(FrameProfAx2, 'ylim', [LogPlotYmin, MaxYvalue*1.05]);
            end
        end
        
        PlottedSets = zeros(1, NumTemperatures, 'logical');
        for SetIndex = 1:NumTemperatures
            if isnan(AllNCParams(SetIndex))
                
                FrameProfAx.Children(end-(2*(SetIndex-1)+1)).XData = 25;
                FrameProfAx.Children(end-(2*(SetIndex-1)+1)).YData = 1;
                FrameProfAx.Children(end-(2*(SetIndex-1))).XData = 25;
                FrameProfAx.Children(end-(2*(SetIndex-1))).YData = 1;
                FrameProfAx.Children(end-(2*(SetIndex-1))).YPositiveDelta = 1;
                FrameProfAx.Children(end-(2*(SetIndex-1))).YNegativeDelta = 1;
                FrameProfAx.Children(end-(2*(SetIndex-1))).XPositiveDelta = 1;
                FrameProfAx.Children(end-(2*(SetIndex-1))).XNegativeDelta = 1;
                set(FrameProfAx.Children(end-(2*(SetIndex-1)+1)),'Visible','off'); %'off' or 'on'
                set(FrameProfAx.Children(end-(2*(SetIndex-1))),'Visible','off'); %'off' or 'on'
                set(get(get(prof{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                if IncludeLogPlots
                    FrameProfAx2.Children(end-(2*(SetIndex-1)+1)).XData = 1/(R*(25+273));
                    FrameProfAx2.Children(end-(2*(SetIndex-1)+1)).YData = 1;
                    FrameProfAx2.Children(end-(2*(SetIndex-1))).XData = 1/(R*(25+273));
                    FrameProfAx2.Children(end-(2*(SetIndex-1))).YData = 1;
                    FrameProfAx2.Children(end-(2*(SetIndex-1))).YPositiveDelta = 1;
                    FrameProfAx2.Children(end-(2*(SetIndex-1))).YNegativeDelta = 1;
                    FrameProfAx2.Children(end-(2*(SetIndex-1))).XPositiveDelta = 1;
                    FrameProfAx2.Children(end-(2*(SetIndex-1))).XNegativeDelta = 1;
                    set(FrameProfAx2.Children(end-(2*(SetIndex-1)+1)),'Visible','off'); %'off' or 'on'
                    set(FrameProfAx2.Children(end-(2*(SetIndex-1))),'Visible','off'); %'off' or 'on'
                    set(get(get(prof2{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                end
            else
                PlottedSets(SetIndex) = 1;
                FrameProfAx.Children(end-(2*(SetIndex-1)+1)).YData = AllNCParams(SetIndex);
                FrameProfAx.Children(end-(2*(SetIndex-1)+1)).XData = NCParamTemperatures(SetIndex);
                set(FrameProfAx.Children(end-(2*(SetIndex-1)+1)),'Visible','on'); %'off' or 'on'
                set(get(get(prof{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
                if PlottedHasSEs
                    FrameProfAx.Children(end-(2*(SetIndex-1))).YData =  AllNCParams(SetIndex);
                    FrameProfAx.Children(end-(2*(SetIndex-1))).XData = NCParamTemperatures(SetIndex);
                    FrameProfAx.Children(end-(2*(SetIndex-1))).YPositiveDelta = AllNCParamSEs(SetIndex);
                    FrameProfAx.Children(end-(2*(SetIndex-1))).YNegativeDelta  = AllNCParamSEs(SetIndex);
                    FrameProfAx.Children(end-(2*(SetIndex-1))).XPositiveDelta = NCParamTemperatureSEs(SetIndex);
                    FrameProfAx.Children(end-(2*(SetIndex-1))).XNegativeDelta  = NCParamTemperatureSEs(SetIndex);
                    set(FrameProfAx.Children(end-(2*(SetIndex-1))),'Visible','on'); %'off' or 'on'
                    
                else
                    set(FrameProfAx.Children(end-(2*(SetIndex-1))),'Visible','off'); %'off' or 'on'
                end
                if IncludeLogPlots
                    FrameProfAx2.Children(end-(2*(SetIndex-1)+1)).YData = AllNCParams(SetIndex);
                    FrameProfAx2.Children(end-(2*(SetIndex-1)+1)).XData =1/(R*(NCParamTemperatures(SetIndex)+273));
                    set(FrameProfAx2.Children(end-(2*(SetIndex-1)+1)),'Visible','on'); %'off' or 'on'
                    set(get(get(prof2{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
                    if PlottedHasSEs
                        FrameProfAx2.Children(end-(2*(SetIndex-1))).YData =  AllNCParams(SetIndex);
                        FrameProfAx2.Children(end-(2*(SetIndex-1))).XData = 1/(R*(NCParamTemperatures(SetIndex)+273));
                        FrameProfAx2.Children(end-(2*(SetIndex-1))).YPositiveDelta = AllNCParamSEs(SetIndex);
                        FrameProfAx2.Children(end-(2*(SetIndex-1))).YNegativeDelta  = AllNCParamSEs(SetIndex);
                        FrameProfAx2.Children(end-(2*(SetIndex-1))).XPositiveDelta =  sqrt(NCParamTemperatureSEs(SetIndex)^2*(1/(R*(NCParamTemperatures(SetIndex)+273)^2))^2);
                        FrameProfAx2.Children(end-(2*(SetIndex-1))).XNegativeDelta  = sqrt(NCParamTemperatureSEs(SetIndex)^2*(1/(R*(NCParamTemperatures(SetIndex)+273)^2))^2);
                        set(FrameProfAx2.Children(end-(2*(SetIndex-1))),'Visible','on'); %'off' or 'on'
                    else
                        set(FrameProfAx2.Children(end-(2*(SetIndex-1))),'Visible','off'); %'off' or 'on'
                    end
                end
            end
        end
        
        
        if all(~PlottedSets)
            continue
        end
        
        if IncludeLogPlots
            LegendAx = subplot(1, 3, 3);
            legend_profs = cell(1, NumSets);
            legend_labels2 = {};
            for SetIndex = 1:NumTemperatures
                temp_idx = SetIndex;
                marker_idx = 1;
                if ~PlottedSets(SetIndex)
                    set(get(get(prof{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                    %                 else
                    %                     %legend_labels2{1, length(legend_labels2)+1} = legend_labels{SetIndex};
                    %                     legend_labels2{1, length(legend_labels2)+1} = [num2str(temperatures(SetIndex)), ' °C'];
                end
                legend_profs{SetIndex} = plot([0], [0], MarkerStyles{1},...
                    'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
                    'linestyle', 'none');
                hold on
                legend_labels2{1, length(legend_labels2)+1} = [num2str(temperatures(SetIndex)), ' °C'];
                
                
                
            end
            hold off
            axis off
            
            
            hlegend = legend(LegendAx, legend_labels2,...
                'FontSize', 10);%, 'Location', );
            
            hlegend.Position(1) = LegendAx.Position(1)+ LegendAx.Position(3)/2 - hlegend.Position(3)/2;
            hlegend.Position(2) = LegendAx.Position(2)+ LegendAx.Position(4)/2 - hlegend.Position(4)/2;
            
            %try
        else
            legend_labels2 = {};
            for SetIndex = 1:NumTemperatures
                legend_labels2{1, length(legend_labels2)+1} = [num2str(temperatures(SetIndex)), ' °C'];
            end
            
            hlegend = legend(legend_labels2,...
                'FontSize', 10, 'Location', 'eastoutside');
            
        end
        
        if exist('PlotTitle', 'var')
            sgtitle([PlotTitle,'Nuclear Cycle ',num2str(NC)])
        else
            sgtitle(['Nuclear Cycle ',num2str(NC)])
        end
        
        
        saveas(FrameProfFig,[outdir3, filesep,'Binned', OutputString,'_NC',num2str(NC), '.png']);
        
        
    end
    close all
end

%%
if ~SkipAllSubplots
    outdir2 = [outdir,filesep, OutputString];
    if ~exist(outdir2, 'dir')
        mkdir(outdir2)
    end
    
    outdir3 = [outdir2,filesep, datestr(now, 'yyyymmdd')];
    if ~exist(outdir3, 'dir')
        mkdir(outdir3)
    end
    if  strcmp(lower(PlottingColors), 'gradient')
        [~, colors] = getColorPalettes();
        GradString = '';
    end
    
    
    PlottedParamSEs((PlottedParams > GlobalPlotYmax) | (PlottedParams < GlobalPlotYmin)) = NaN;
    PlottedParams((PlottedParams > GlobalPlotYmax) | (PlottedParams < GlobalPlotYmin)) = NaN;
    
    WhereValidNC = squeeze(sum(~isnan(PlottedParams), 1));
    ValidNCIndices = find(WhereValidNC);
    
    
    SubplotDims = [1, length(ValidNCIndices)+1];
    
    SubFigDims = [0.95, 0.95*3072/1920*SubplotDims(1)/SubplotDims(2)*1.2];
    SubFigDims(2) = min([SubFigDims(2), 0.95]);
    SubFigDims(2) = max([SubFigDims(2), 0.5]);
    SubFigDims = round(SubFigDims, 2);
    SubFigBuffer = 0.02;
    SubFigWidth = (1-3*SubFigBuffer)/SubplotDims(2);
    
    SubFigPositions = 0:(SubplotDims(2)-1);
    SubFigPositions = SubFigBuffer*2+SubFigPositions*SubFigWidth;
    eb = cell(length(ValidNCIndices), NumSets);
    prof = cell(length(ValidNCIndices), NumSets);
    FrameProfAx = cell(1,  length(ValidNCIndices));
    close all
    FrameProfFig = figure(1);
    set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.01, SubFigDims(1), SubFigDims(2)]);
    set(gcf,'color','w');
    
    for nc_idx=1:length(ValidNCIndices)
        NC = ValidNCIndices(nc_idx)+8;
        % Prepare Traces for plotting
        
        NCMaxParams = NaN(1,NumSets);
        AllNCParams = NaN(1,NumSets);
        AllNCParamSEs = NaN(1,NumSets);
        AllNCTemperatures = NaN(1,NumSets);
        for SetIndex=1:NumSets
            AllNCParams(SetIndex) = PlottedParams(SetIndex,NC-8);
            AllNCParamSEs(SetIndex) = PlottedParamSEs(SetIndex,NC-8);
            AllNCTemperatures(SetIndex) = this.Temp_obs(SetIndex);
            TemporarySEs = AllNCParamSEs(SetIndex);
            TemporarySEs(isnan(TemporarySEs)) = 0;
            NCMaxParams(SetIndex) = AllNCParams(SetIndex) + TemporarySEs;
        end
        
        if all(isnan(NCMaxParams))
            continue
        end
        
        if nc_idx == 1
            FrameProfAx{nc_idx} = subplot(SubplotDims(1), SubplotDims(2), nc_idx, gca);
        else
            FrameProfAx{nc_idx} = subplot(SubplotDims(1), SubplotDims(2), nc_idx);
        end
        
        PlottedSets = zeros(1, NumSets, 'logical');
        for SetIndex = 1:NumSets
            temp_idx = find(temperatures == this.Temp_sps(SetIndex));
            marker_idx = find(TempMatches{temp_idx} == SetIndex);
            if isnan(AllNCParams(SetIndex))
                eb{nc_idx, SetIndex} = errorbar([0], [0], [0], 'vertical',  'LineStyle', 'none');
                hold on
                set(eb{nc_idx,SetIndex}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                set(get(get(eb{nc_idx,SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                set(eb{nc_idx,SetIndex},'Visible','off'); %'off' or 'on'
                
                prof{nc_idx,SetIndex} = plot([0], [0], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                    'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
                    'linestyle', 'none');
                set(prof{nc_idx,SetIndex},'Visible','off'); %'off' or 'on'
                if ~ismember(SetIndex, this.ProcessedExperiments)
                    set(get(get(prof{nc_idx,SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                end
                
            else
                PlottedSets(SetIndex) = 1;
                eb{nc_idx,SetIndex} = errorbar([AllNCTemperatures(SetIndex)],...
                    [AllNCParams(SetIndex)], [AllNCParamSEs(SetIndex)], 'vertical',  'LineStyle', 'none');
                hold on
                set(eb{nc_idx,SetIndex}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                set(get(get(eb{nc_idx,SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                set(eb{nc_idx,SetIndex},'Visible','on'); %'off' or 'on'
                prof{nc_idx,SetIndex} = plot([AllNCTemperatures(SetIndex)],...
                    [AllNCParams(SetIndex)], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                    'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
                    'linestyle', 'none');
                set(prof{nc_idx,SetIndex},'Visible','on'); %'off' or 'on'
            end
        end
        
        hold off
        
        xlabel('Temperature (°C)')
        %xlim([Subplot1DataXmin, Subplot1DataXmax])
        xlim([15, 30])
        
        ylabel(ylab)
        ylim([GlobalPlotYmin, GlobalPlotYmax])
        title(['NC: ', num2str(NC)])
        
        pos = get(gca, 'Position');
        pos(1) = SubFigPositions(nc_idx);
        pos(3) = SubFigWidth-2*SubFigBuffer;
        set(gca, 'Position', pos)
        
    end
    %         if exist('PlotTitle', 'var')
    %             sgtitle({PlotTitle, ['Nuclear Cycle ',num2str(NC)]})
    %         else
    %             sgtitle(['Nuclear Cycle ',num2str(NC)])
    %         end
    LegendAx = subplot(SubplotDims(1), SubplotDims(2), SubplotDims(1)*SubplotDims(2));
    legend_profs = cell(1, NumSets);
    legend_labels2 = {};
    for SetIndex = 1:NumSets
        temp_idx = find(temperatures == Temp_sp(SetIndex));
        marker_idx = find(TempMatches{temp_idx} == SetIndex);
        
        legend_profs{SetIndex} = plot([0], [0], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
            'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
            'linestyle', 'none');
        hold on
        if ~PlottedSets(SetIndex)
            set(get(get(legend_profs{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        else
            legend_labels2{1, length(legend_labels2)+1} = legend_labels{SetIndex};
        end
    end
    xlim([-5, -1])
    hold off
    axis off
    if length(legend_labels2)
        hlegend = legend(LegendAx, legend_labels2,...
            'FontSize', 6, 'Location', 'eastoutside', 'NumColumns', 2);%, 'Location', );
        legpos = get(hlegend, 'Position');
        legposnew = legpos;
        legposnew(1) = legposnew(1)-0.035;
        set(hlegend, 'Position', legposnew);
    else
        hlegend = legend(LegendAx, legend_labels2,...
            'FontSize', 10, 'Location', 'eastoutside');%, 'Location', );
    end
    
    saveas(FrameProfFig,[outdir3, filesep,OutputString,'_Subplots.png']);
    
    
    %%
    if IncludeLogPlots
        outdir2 = [outdir,filesep, OutputString];
        if ~exist(outdir2, 'dir')
            mkdir(outdir2)
        end
        
        outdir3 = [outdir2,filesep, datestr(now, 'yyyymmdd')];
        if ~exist(outdir3, 'dir')
            mkdir(outdir3)
        end
        if  strcmp(lower(PlottingColors), 'gradient')
            [~, colors] = getColorPalettes();
            GradString = '';
        end
        
        
        PlottedParamSEs((PlottedParams > GlobalPlotYmax) | (PlottedParams < GlobalPlotYmin)) = NaN;
        PlottedParams((PlottedParams > GlobalPlotYmax) | (PlottedParams < GlobalPlotYmin)) = NaN;
        
        WhereValidNC = squeeze(sum(~isnan(PlottedParams), 1));
        ValidNCIndices = find(WhereValidNC);
        
        
        SubplotDims = [1, length(ValidNCIndices)+1];
        
        SubFigDims = [0.95, 0.95*3072/1920*SubplotDims(1)/SubplotDims(2)*1.2];
        SubFigDims(2) = min([SubFigDims(2), 0.95]);
        SubFigDims(2) = max([SubFigDims(2), 0.5]);
        SubFigDims = round(SubFigDims, 2);
        SubFigBuffer = 0.02;
        SubFigWidth = (1-3*SubFigBuffer)/SubplotDims(2);
        
        SubFigPositions = 0:(SubplotDims(2)-1);
        SubFigPositions = SubFigBuffer*2+SubFigPositions*SubFigWidth;
        eb = cell(length(ValidNCIndices), NumSets);
        prof = cell(length(ValidNCIndices), NumSets);
        FrameProfAx = cell(1,  length(ValidNCIndices));
        close all
        FrameProfFig = figure(1);
        set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.01, SubFigDims(1), SubFigDims(2)]);
        set(gcf,'color','w');
        
        for nc_idx=1:length(ValidNCIndices)
            NC = ValidNCIndices(nc_idx)+8;
            % Prepare Traces for plotting
            
            NCMaxParams = NaN(1,NumSets);
            AllNCParams = NaN(1,NumSets);
            AllNCParamSEs = NaN(1,NumSets);
            AllNCTemperatures = NaN(1,NumSets);
            for SetIndex=1:NumSets
                AllNCParams(SetIndex) = PlottedParams(SetIndex,NC-8);
                AllNCParamSEs(SetIndex) = PlottedParamSEs(SetIndex,NC-8);
                AllNCTemperatures(SetIndex) = this.Temp_obs(SetIndex);
                TemporarySEs = AllNCParamSEs(SetIndex);
                TemporarySEs(isnan(TemporarySEs)) = 0;
                NCMaxParams(SetIndex) = AllNCParams(SetIndex) + TemporarySEs;
            end
            
            if all(isnan(NCMaxParams))
                continue
            end
            
            if nc_idx == 1
                FrameProfAx{nc_idx} = subplot(SubplotDims(1), SubplotDims(2), nc_idx, gca);
            else
                FrameProfAx{nc_idx} = subplot(SubplotDims(1), SubplotDims(2), nc_idx);
            end
            
            PlottedSets = zeros(1, NumSets, 'logical');
            for SetIndex = 1:NumSets
                temp_idx = find(temperatures == this.Temp_sps(SetIndex));
                marker_idx = find(TempMatches{temp_idx} == SetIndex);
                if isnan(AllNCParams(SetIndex))
                    eb{nc_idx, SetIndex} = errorbar([0], [0], [0], 'vertical', 'LineStyle', 'none');
                    hold on
                    set(eb{nc_idx,SetIndex}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                    set(get(get(eb{nc_idx,SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                    set(eb{nc_idx,SetIndex},'Visible','off'); %'off' or 'on'
                    
                    prof{nc_idx,SetIndex} = plot([0], [0], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                        'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
                        'linestyle', 'none');
                    set(prof{nc_idx,SetIndex},'Visible','off'); %'off' or 'on'
                    if ~ismember(SetIndex, this.ProcessedExperiments)
                        set(get(get(prof{nc_idx,SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                    end
                    
                else
                    PlottedSets(SetIndex) = 1;
                    eb{nc_idx,SetIndex} = errorbar([1/(R*(AllNCTemperatures(SetIndex)+273))],...
                        [AllNCParams(SetIndex)], [AllNCParamSEs(SetIndex)],'vertical',  'LineStyle', 'none');
                    hold on
                    set(eb{nc_idx,SetIndex}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                    set(get(get(eb{nc_idx,SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                    set(eb{nc_idx,SetIndex},'Visible','on'); %'off' or 'on'
                    prof{nc_idx,SetIndex} = plot([1/(R*(AllNCTemperatures(SetIndex)+273))],...
                        [AllNCParams(SetIndex)], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                        'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
                        'linestyle', 'none');
                    set(prof{nc_idx,SetIndex},'Visible','on'); %'off' or 'on'
                end
            end
            
            hold off
            
            
            title(['NC: ', num2str(NC)])
            
            xlabel('1/(RT) (mol/kJ)')
            xlim([Plot2Xmin, Plot2Xmax])
            
            
            ylabel(ylab)
            ylim([LogPlotYmin, GlobalPlotYmax])
            
            set(FrameProfAx{nc_idx} , 'YScale', 'log')
            
            pos = get(gca, 'Position');
            pos(1) = SubFigPositions(nc_idx);
            pos(3) = SubFigWidth-2*SubFigBuffer;
            set(gca, 'Position', pos)
        end
        %         if exist('PlotTitle', 'var')
        %             sgtitle({PlotTitle, ['Nuclear Cycle ',num2str(NC)]})
        %         else
        %             sgtitle(['Nuclear Cycle ',num2str(NC)])
        %         end
        LegendAx = subplot(SubplotDims(1), SubplotDims(2), SubplotDims(1)*SubplotDims(2));
        legend_profs = cell(1, NumSets);
        legend_labels2 = {};
        for SetIndex = 1:NumSets
            temp_idx = find(temperatures == Temp_sp(SetIndex));
            marker_idx = find(TempMatches{temp_idx} == SetIndex);
            
            legend_profs{SetIndex} = plot([0], [0], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
                'linestyle', 'none');
            hold on
            if ~PlottedSets(SetIndex)
                set(get(get(legend_profs{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            else
                legend_labels2{1, length(legend_labels2)+1} = legend_labels{SetIndex};
            end
        end
        xlim([-5, -1])
        hold off
        axis off
        
        if length(legend_labels2)
            hlegend = legend(LegendAx, legend_labels2,...
                'FontSize', 6, 'Location', 'eastoutside', 'NumColumns', 2);%, 'Location', );
            legpos = get(hlegend, 'Position');
            legposnew = legpos;
            legposnew(1) = legposnew(1)-0.035;
            set(hlegend, 'Position', legposnew);
        else
            hlegend = legend(LegendAx, legend_labels2,...
                'FontSize', 10, 'Location', 'eastoutside');%, 'Location', );
        end
        
        
        saveas(FrameProfFig,[outdir3, filesep,OutputString,'_LogSubplots.png']);
    end
end


%%
if ~SkipAllBinnedSubplots
    outdir2 = [outdir,filesep, OutputString];
    if ~exist(outdir2, 'dir')
        mkdir(outdir2)
    end
    
    outdir3 = [outdir2,filesep, datestr(now, 'yyyymmdd')];
    if ~exist(outdir3, 'dir')
        mkdir(outdir3)
    end
    if  strcmp(lower(PlottingColors), 'gradient')
        [~, colors] = getColorPalettes();
        GradString = '';
    end
    
    
    BinnedSEParams((Counts < this.MinimumBinCount) | (BinnedParams > GlobalPlotYmax) | (BinnedParams < GlobalPlotYmin)) = NaN;
    BinnedParams((Counts < this.MinimumBinCount) |(BinnedParams > GlobalPlotYmax) | (BinnedParams < GlobalPlotYmin)) = NaN;
    
    WhereValidNC = squeeze(sum(~isnan(BinnedParams), 1));
    ValidNCIndices = find(WhereValidNC);
    
    
    SubplotDims = [1, length(ValidNCIndices)+1];
    
    SubFigDims = [0.95, 0.95*3072/1920*SubplotDims(1)/SubplotDims(2)*1.2];
    SubFigDims(2) = min([SubFigDims(2), 0.95]);
    SubFigDims(2) = max([SubFigDims(2), 0.5]);
    SubFigDims = round(SubFigDims, 2);
    SubFigBuffer = 0.02;
    SubFigWidth = (1-3*SubFigBuffer)/SubplotDims(2);
    
    SubFigPositions = 0:(SubplotDims(2)-1);
    SubFigPositions = SubFigBuffer*2+SubFigPositions*SubFigWidth;
    eb = cell(length(ValidNCIndices), NumTemperatures);
    prof = cell(length(ValidNCIndices), NumTemperatures);
    FrameProfAx = cell(1,  length(ValidNCIndices));
    close all
    FrameProfFig = figure(1);
    set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.01, SubFigDims(1), SubFigDims(2)]);
    set(gcf,'color','w');
    
    for nc_idx=1:length(ValidNCIndices)
        NC = ValidNCIndices(nc_idx)+8;
        % Prepare Traces for plotting
        NCMaxParams = NaN(1,NumTemperatures);
        AllNCParams = NaN(1,NumTemperatures);
        AllNCParamSEs = NaN(1,NumTemperatures);
        AllNCTemperatures = NaN(1,NumTemperatures);
        AllNCTemperatureSEs = NaN(1, NumTemperatures);
        for SetIndex=1:NumTemperatures
            AllNCParams(SetIndex) = BinnedParams(SetIndex,NC-8);
            AllNCParamSEs(SetIndex) = BinnedSEParams(SetIndex,NC-8);
            AllNCTemperatures(SetIndex) = ParamTemperatures(SetIndex, NC-8);
            AllNCTemperatureSEs(SetIndex) = ParamSETemperatures(SetIndex, NC-8);
            TemporarySEs = AllNCParamSEs(SetIndex);
            TemporarySEs(isnan(TemporarySEs)) = 0;
            NCMaxParams(SetIndex) = AllNCParams(SetIndex) + TemporarySEs;
        end
        
        if all(isnan(NCMaxParams))
            continue
        end
        
        if nc_idx == 1
            FrameProfAx{nc_idx} = subplot(SubplotDims(1), SubplotDims(2), nc_idx, gca);
        else
            FrameProfAx{nc_idx} = subplot(SubplotDims(1), SubplotDims(2), nc_idx);
        end
        
        PlottedSets = zeros(1, NumTemperatures, 'logical');
        for SetIndex = 1:NumTemperatures
            temp_idx = SetIndex;
            marker_idx = 1;
            if isnan(AllNCParams(SetIndex))
                eb{nc_idx, SetIndex} = errorbar([0],[0],[0],[0],[0],[0], 'LineStyle', 'none');
                hold on
                set(eb{nc_idx,SetIndex}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                set(get(get(eb{nc_idx,SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                set(eb{nc_idx,SetIndex},'Visible','off'); %'off' or 'on'
                
                prof{nc_idx,SetIndex} = plot([0], [0], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                    'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
                    'linestyle', 'none');
                set(prof{nc_idx,SetIndex},'Visible','off'); %'off' or 'on'
                if ~ismember(SetIndex, this.ProcessedExperiments)
                    set(get(get(prof{nc_idx,SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                end
                
            else
                PlottedSets(SetIndex) = 1;
                eb{nc_idx,SetIndex} = errorbar([AllNCTemperatures(SetIndex)],...
                    [AllNCParams(SetIndex)], [AllNCParamSEs(SetIndex)],...
                    [AllNCParamSEs(SetIndex)],[AllNCTemperatureSEs(SetIndex)],...
                    [AllNCTemperatureSEs(SetIndex)], 'LineStyle', 'none');
                hold on
                set(eb{nc_idx,SetIndex}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                set(get(get(eb{nc_idx,SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                set(eb{nc_idx,SetIndex},'Visible','on'); %'off' or 'on'
                prof{nc_idx,SetIndex} = plot([AllNCTemperatures(SetIndex)],...
                    [AllNCParams(SetIndex)], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                    'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
                    'linestyle', 'none');
                set(prof{nc_idx,SetIndex},'Visible','on'); %'off' or 'on'
            end
        end
        
        hold off
        
        xlabel('Temperature (°C)')
        %xlim([Subplot1DataXmin, Subplot1DataXmax])
        xlim([min(15, Subplot1DataXmin), max(Subplot1DataXmax, 30)])
        
        ylabel(ylab)
        ylim([GlobalPlotYmin, GlobalPlotYmax])
        title(['NC: ', num2str(NC)])
        
        pos = get(gca, 'Position');
        pos(1) = SubFigPositions(nc_idx);
        pos(3) = SubFigWidth-2*SubFigBuffer;
        set(gca, 'Position', pos)
        %disp('pause')
    end
    %         if exist('PlotTitle', 'var')
    %             sgtitle({PlotTitle, ['Nuclear Cycle ',num2str(NC)]})
    %         else
    %             sgtitle(['Nuclear Cycle ',num2str(NC)])
    %         end
    LegendAx = subplot(SubplotDims(1), SubplotDims(2), SubplotDims(1)*SubplotDims(2));
    legend_profs = cell(1, NumTemperatures);
    legend_labels2 = {};
    for SetIndex = 1:NumTemperatures
        temp_idx = SetIndex;
        marker_idx = 1;
        
        legend_profs{SetIndex} = plot([0], [0], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
            'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
            'linestyle', 'none');
        hold on
        if ~PlottedSets(SetIndex)
            set(get(get(legend_profs{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        else
            legend_labels2{1, length(legend_labels2)+1} = [num2str(temperatures(SetIndex)), ' °C'];
        end
    end
    xlim([-5, -1])
    hold off
    axis off
    
    if length(legend_labels2)
        hlegend = legend(LegendAx, legend_labels2,...
            'FontSize', 10, 'Location', 'eastoutside', 'NumColumns', 2);%, 'Location', );
    else
        hlegend = legend(LegendAx, legend_labels2,...
            'FontSize', 10, 'Location', 'eastoutside');
    end
    
    saveas(FrameProfFig,[outdir3, filesep,OutputString,'_BinnedSubplots.png']);
    
    
    %%
    if IncludeLogPlots
        outdir2 = [outdir,filesep, OutputString];
        if ~exist(outdir2, 'dir')
            mkdir(outdir2)
        end
        
        outdir3 = [outdir2,filesep, datestr(now, 'yyyymmdd')];
        if ~exist(outdir3, 'dir')
            mkdir(outdir3)
        end
        if  strcmp(lower(PlottingColors), 'gradient')
            [~, colors] = getColorPalettes();
            GradString = '';
        end
        
        
        BinnedSEParams((Counts < this.MinimumBinCount) | (BinnedParams > GlobalPlotYmax) | (BinnedParams < GlobalPlotYmin)) = NaN;
        BinnedParams((Counts < this.MinimumBinCount) | (BinnedParams > GlobalPlotYmax) | (BinnedParams < GlobalPlotYmin)) = NaN;
        
        WhereValidNC = squeeze(sum(~isnan(BinnedParams), 1));
        ValidNCIndices = find(WhereValidNC);
        
        
        SubplotDims = [1, length(ValidNCIndices)+1];
        
        SubFigDims = [0.95, 0.95*3072/1920*SubplotDims(1)/SubplotDims(2)*1.2];
        SubFigDims(2) = min([SubFigDims(2), 0.95]);
        SubFigDims(2) = max([SubFigDims(2), 0.5]);
        SubFigDims = round(SubFigDims, 2);
        SubFigBuffer = 0.02;
        SubFigWidth = (1-3*SubFigBuffer)/SubplotDims(2);
        
        SubFigPositions = 0:(SubplotDims(2)-1);
        SubFigPositions = SubFigBuffer*2+SubFigPositions*SubFigWidth;
        eb = cell(length(ValidNCIndices), NumTemperatures);
        prof = cell(length(ValidNCIndices), NumTemperatures);
        FrameProfAx = cell(1,  length(ValidNCIndices));
        close all
        FrameProfFig = figure(1);
        set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.01, SubFigDims(1), SubFigDims(2)]);
        set(gcf,'color','w');
        
        for nc_idx=1:length(ValidNCIndices)
            NC = ValidNCIndices(nc_idx)+8;
            % Prepare Traces for plotting
            
            NCMaxParams = NaN(1,NumTemperatures);
            AllNCParams = NaN(1,NumTemperatures);
            AllNCParamSEs = NaN(1,NumTemperatures);
            AllNCTemperatures = NaN(1,NumTemperatures);
            AllNCTemperatureSEs = NaN(1,NumTemperatures);
            for SetIndex=1:NumTemperatures
                AllNCParams(SetIndex) = BinnedParams(SetIndex,NC-8);
                AllNCParamSEs(SetIndex) = BinnedSEParams(SetIndex,NC-8);
                AllNCTemperatures(SetIndex) = ParamTemperatures(SetIndex, NC-8);
                AllNCTemperatureSEs(SetIndex) = ParamSETemperatures(SetIndex, NC-8);
                TemporarySEs = AllNCParamSEs(SetIndex);
                TemporarySEs(isnan(TemporarySEs)) = 0;
                NCMaxParams(SetIndex) = AllNCParams(SetIndex) + TemporarySEs;
            end
            
            if all(isnan(NCMaxParams))
                continue
            end
            
            if nc_idx == 1
                FrameProfAx{nc_idx} = subplot(SubplotDims(1), SubplotDims(2), nc_idx, gca);
            else
                FrameProfAx{nc_idx} = subplot(SubplotDims(1), SubplotDims(2), nc_idx);
            end
            
            PlottedSets = zeros(1, NumTemperatures, 'logical');
            for SetIndex = 1:NumTemperatures
                temp_idx = SetIndex;
                marker_idx = 1;
                if isnan(AllNCParams(SetIndex))
                    eb{nc_idx, SetIndex} = errorbar([0], [0], [0],...
                        [0], [0], [0], 'LineStyle', 'none');
                    hold on
                    set(eb{nc_idx,SetIndex}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                    set(get(get(eb{nc_idx,SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                    set(eb{nc_idx,SetIndex},'Visible','off'); %'off' or 'on'
                    
                    prof{nc_idx,SetIndex} = plot([0], [0], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                        'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
                        'linestyle', 'none');
                    set(prof{nc_idx,SetIndex},'Visible','off'); %'off' or 'on'
                    if ~ismember(SetIndex, this.ProcessedExperiments)
                        set(get(get(prof{nc_idx,SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                    end
                    
                else
                    PlottedSets(SetIndex) = 1;
                    Xdelta = sqrt(AllNCTemperatureSEs(SetIndex)^2*(1/(R*(AllNCTemperatures(SetIndex)+273)^2))^2);
                    eb{nc_idx,SetIndex} = errorbar([1/(R*(AllNCTemperatures(SetIndex)+273))],...
                        [AllNCParams(SetIndex)], [AllNCParamSEs(SetIndex)],[AllNCParamSEs(SetIndex)],...
                        [Xdelta], [Xdelta], 'LineStyle', 'none');
                    hold on
                    set(eb{nc_idx,SetIndex}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                    set(get(get(eb{nc_idx,SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                    set(eb{nc_idx,SetIndex},'Visible','on'); %'off' or 'on'
                    prof{nc_idx,SetIndex} = plot([1/(R*(AllNCTemperatures(SetIndex)+273))],...
                        [AllNCParams(SetIndex)], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                        'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
                        'linestyle', 'none');
                    set(prof{nc_idx,SetIndex},'Visible','on'); %'off' or 'on'
                end
            end
            
            hold off
            
            
            title(['NC: ', num2str(NC)])
            
            xlabel('1/(RT) (mol/kJ)')
            xlim([Plot2Xmin, Plot2Xmax])
            
            
            ylabel(ylab)
            ylim([LogPlotYmin, GlobalPlotYmax])
            
            set(FrameProfAx{nc_idx} , 'YScale', 'log')
            
            pos = get(gca, 'Position');
            pos(1) = SubFigPositions(nc_idx);
            pos(3) = SubFigWidth-2*SubFigBuffer;
            set(gca, 'Position', pos)
        end
        %         if exist('PlotTitle', 'var')
        %             sgtitle({PlotTitle, ['Nuclear Cycle ',num2str(NC)]})
        %         else
        %             sgtitle(['Nuclear Cycle ',num2str(NC)])
        %         end
        LegendAx = subplot(SubplotDims(1), SubplotDims(2), SubplotDims(1)*SubplotDims(2));
        legend_profs = cell(1, NumSets);
        legend_labels2 = {};
        for SetIndex = 1:NumTemperatures
            temp_idx = SetIndex;
            marker_idx = 1;
            
            legend_profs{SetIndex} = plot([0], [0], MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                'MarkerEdgeColor', colors(temp_idx,:),'MarkerFaceColor', colors(temp_idx,:),...
                'linestyle', 'none');
            hold on
            if ~PlottedSets(SetIndex)
                set(get(get(legend_profs{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            else
                legend_labels2{1, length(legend_labels2)+1} = [num2str(temperatures(SetIndex)), ' °C'];
            end
        end
        xlim([-5, -1])
        hold off
        axis off
        
        if length(legend_labels2) > 25
            hlegend = legend(LegendAx, legend_labels2,...
                'FontSize', 10, 'Location', 'eastoutside','NumColumns', 2);%, 'Location', );
        else
            hlegend = legend(LegendAx, legend_labels2,...
                'FontSize', 10, 'Location', 'eastoutside');
        end
        
        
        saveas(FrameProfFig,[outdir3, filesep,OutputString,'_BinnedLogSubplots.png']);
    end
end


end








