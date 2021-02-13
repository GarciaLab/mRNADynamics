function PlotLTMSingleTempTrapParamsVsAP(this, parameter, outdir, varargin)
%%

% PlotTitle, PlottingColors, UseDifferentColors,
% UseDiffProfiles, UsePhysicalAPLength
UsePhysicalAPLength = false;
UseLines = true;
UseDifferentColors = true;

x = 1;
while x <= length(varargin)
    if strcmp(lower(varargin{x}), 'plottitle')
        PlotTitle = varargin{x+1};
        x = x+1;
    elseif strcmp(lower(varargin{x}), 'plottingcolors')
        PlottingColors = varargin{x+1};
        x = x+1;
    elseif strcmp(lower(varargin{x}), 'usephysicalaplength')
        UsePhysicalAPLength = true;
    elseif strcmp(lower(varargin{x}), 'noline')
        UseLines = false;
    elseif strcmp(lower(varargin{x}), 'tracetype')
        TraceType = lower(varargin{x+1});
        x = x+1;
    elseif strcmp(lower(varargin{x}), 'uniformcolor')
        UseDifferentColors = false;
    end
    x = x+1;
end


if ~exist('PlottingColors', 'var')
    PlottingColors = 'default';
elseif ~strcmpi(PlottingColors, 'gradient') &~strcmp(lower(PlottingColors), 'default')  & ~strcmp(lower(PlottingColors), 'pboc')
    error('Invalid choice of plotting colors. Can use either "default", "pboc", or "gradient".') % change to error
end
if ~exist('TraceType', 'var')
    TraceType = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'Fluo3D') | strcmpi(TraceType, 'Unaligned3D')
    TraceType = 'Unaligned3D';
elseif strcmpi(TraceType, 'Fluo') | strcmpi(TraceType, 'Unaligned')
    TraceType = 'Unaligned';
elseif strcmpi(TraceType, 'AnaphaseAligned')
    TraceType = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'AnaphaseAligned3D')
    TraceType = 'AnaphaseAligned3D';
elseif strcmpi(TraceType, 'Tbinned')
    TraceType = 'Tbinned';
elseif strcmpi(TraceType, 'Tbinned3D')
    TraceType = 'Tbinned3D';
else
    error('Invalid choice of trace type. Can use either "fluo", "fluo3d", "anaphasealigned", or "anaphasealigned3d".') % change to error
end

%%

if ~exist(outdir, 'dir')
    mkdir(outdir)
end
if UsePhysicalAPLength
    PhysicalAPString = 'PhysAP';
else
    PhysicalAPString = '';
end
Temp_obs = this.Temp_obs;
Temp_sp = this.Temp_sps;


if strcmp(lower(PlottingColors), 'pboc')
    [colors, ~] = getColorPalettes();
else
    [~, colors] = getColorPalettes();
end

%%

R = this.R;
R2bound = this.R2bound;

temperatures = flip(unique(this.Temp_sps(this.ProcessedExperiments)));
NumTemperatures = length(temperatures);


APResolution = this.Experiments{1}.APResolution;
APbins = 0:APResolution:1;
NumAPbins = length(APbins);

NumSets = length(this.ExperimentPrefixes);
UseSet = ismember(1:NumSets, this.ProcessedExperiments);
Nsigfigs = 3;
legend_labels = this.LegendLabels;


MarkerStyles = {'o', 'd', 's', '>', '^','p', 'h', '*', 'x'};

TempMatches = cell(1, NumTemperatures);

for t_index = 1:NumTemperatures
    TempMatches{t_index} = find((this.Temp_sps == temperatures(t_index)) & UseSet);
end


%% Load relevant parameters into memory
[PlottedParams, PlottedParamSEs,R2s, ylab,OutputString,GlobalPlotYmax,GlobalPlotYmin,LogPlotYmin] = ...
    getPlottingVariables(this, parameter,  TraceType, R2bound);

% Calculate Plot Xlims
ParamMatCopy = PlottedParams;
if ~all(all(all(isnan(PlottedParamSEs))))
    ParamMatCopy((R2s < R2bound) | (PlottedParams./PlottedParamSEs < 1) | isnan(PlottedParams)) = ...
        NaN;
else
    ParamMatCopy((R2s < R2bound) | isnan(PlottedParams)) = NaN;
end
ObservedAPbins =  sum((sum(~isnan(ParamMatCopy),3) > 0), 1) > 0;
PlotXmin = max([0, APResolution*(find(ObservedAPbins > 0, 1)-2)]);
PlotXmax = min([1, APResolution*(find(ObservedAPbins > 0, 1, 'last')+1)]);
if UsePhysicalAPLength
    PlotXmin = PlotXmin*min(this.APLengths);
    PlotXmax = PlotXmax*max(this.APLengths);
end


[BinnedParams, BinnedSEParams, Counts, ParamTemperatures, ParamSETemperatures] = ...
    getBinnedPlottingVariables(this, PlottedParams, PlottedParamSEs,R2s, R2bound);
%%

outdir2 = [outdir,filesep,OutputString];
if ~exist(outdir2, 'dir')
    mkdir(outdir2)
end
for TemperatureIndex =1:NumTemperatures
    CurrentTemperature = temperatures(TemperatureIndex);

    if isempty(PhysicalAPString)
    outdir4 = [outdir2,filesep, TraceType];
    if ~exist(outdir4, 'dir')
        mkdir(outdir4)
    end
    else
    outdir4 = [outdir3,filesep, TraceType,'_', PhysicalAPString];
    if ~exist(outdir4, 'dir')
        mkdir(outdir4)
    end
    end
        
    outdir5 = [outdir4,filesep, datestr(now, 'yyyymmdd')];
    if ~exist(outdir5, 'dir')
        mkdir(outdir5)
    end
    Subset = TempMatches{TemperatureIndex};
    SubNumSets = length(Subset);
    eb = cell(1, SubNumSets);
    prof = cell(1, SubNumSets);
    close all
    FrameProfFig = figure(1);
    set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .9, .7]);
    set(gcf,'color','w');
    FrameProfAx = axes(FrameProfFig);
    for SetIndex =1:SubNumSets
        if UseDifferentColors
            ColIndex = SetIndex;
        else
            ColIndex = TemperatureIndex;
        end
        
        MarkerIndex = SetIndex;
        if isempty(MarkerIndex)
            MarkerIndex = length(MarkerStyles);
        end
        
        eb{SetIndex} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)),...
            'vertical', 'LineStyle', 'none');
        hold on
        if strcmp(lower(PlottingColors), 'gradient')
            set(eb{SetIndex}, 'color', colors(ColIndex,:), 'LineWidth', 1, 'CapSize', 0);
        else
            set(eb{SetIndex}, 'color', colors(ColIndex,:), 'LineWidth', 1, 'CapSize', 0);
        end
        set(get(get(eb{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        if UseLines
            prof{SetIndex} = plot(APbins, ones(1, length(APbins)), [MarkerStyles{MarkerIndex}, '-'],...
                'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:), 'Color', colors(ColIndex,:),...
                'MarkerSize', 10);
            
        else
            prof{SetIndex} = plot(APbins, ones(1, length(APbins)), MarkerStyles{MarkerIndex},...
                'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:),...
                'MarkerSize', 10, 'linestyle', 'none');
            
        end
        
        set(eb{SetIndex},'Visible','off'); %'off' or 'on'
        set(prof{SetIndex},'Visible','off'); %'off' or 'on'
        if ~ismember(SetIndex, this.ProcessedExperiments)
            set(get(get(prof{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        end
        
    end
    
    grid on
    hold off
    if ~UsePhysicalAPLength
        xlabel('Fraction Embryo Length')
    else
        xlabel('Distance from the Anterior Pole (\mum)')
    end
    xlim([PlotXmin, PlotXmax])
    
    ylabel(ylab)
    ylim([GlobalPlotYmin, GlobalPlotYmax])
    
    FrameProfAx.YAxis.FontSize = 14;
    FrameProfAx.XAxis.FontSize = 14;
    title(FrameProfAx, '', 'FontSize', 14);
    
    %%
    
    for NCIndex=1:length(this.IncludedNCs)
        NC = this.IncludedNCs(NCIndex);
        
        
        % Prepare Traces for plotting
        
        NCMaxParams = NaN(1, SubNumSets);
        AllNCParams = NaN(SubNumSets, NumAPbins);
        AllNCParamSEs = NaN(SubNumSets, NumAPbins);
        AllR2s = NaN(SubNumSets, NumAPbins);
        for SetIndex=1:SubNumSets
            SetParams = PlottedParams(Subset(SetIndex),:,NC-8).';
            SetSEParams = PlottedParamSEs(Subset(SetIndex),:,NC-8).';
            SetR2s = R2s(Subset(SetIndex),:,NC-8).';
            if ~all(isnan(SetSEParams))
                IncludedBins = find(~isnan(SetParams) & (SetR2s > R2bound) & (SetParams./SetSEParams >= 1)) ;
            else
                IncludedBins = find(~isnan(SetParams) & (SetR2s > R2bound)) ;
            end
            if ~isempty(IncludedBins)
                AllNCParams(SetIndex,:) = SetParams;
                AllNCParamSEs(SetIndex,:) = SetSEParams;
                TempSEParams = SetSEParams;
                TempSEParams(isnan(TempSEParams)) = 0;
                NCMaxParams(SetIndex) = max(SetParams+TempSEParams);
                AllR2s(SetIndex,:) = SetR2s;
            end
        end
        
        if all(isnan(NCMaxParams))
            continue
        end
        
        
        
        
        
        PlottedSets = zeros(1, NumSets, 'logical');
        for SetIndex=1:SubNumSets
            if ~ismember(Subset(SetIndex), this.ProcessedExperiments)
                continue
            end
            if UsePhysicalAPLength
                APLength = this.APLengths(Subset(SetIndex));
            else
                APLength = 1;
            end
            if ~all(isnan(AllNCParamSEs(SetIndex, :)))
                UseIndex = ~isnan(AllNCParams(SetIndex,:)) &  (AllNCParams(SetIndex, :)./AllNCParamSEs(SetIndex, :) >= 1) & ...
                    (AllR2s(SetIndex,:)>= R2bound);
            else
                UseIndex = ~isnan(AllNCParams(SetIndex,:)) & ...
                    (AllR2s(SetIndex,:)>= R2bound);
            end
            
            if sum(UseIndex) == 0 %| sum(DiffMeanFluoMat(i, use_SetIndex, j) == 0)
                FrameProfAx.Children(end-(2*(SetIndex-1)+1)).XData = APLength*APbins;
                FrameProfAx.Children(end-(2*(SetIndex-1)+1)).YData = zeros(1, length(APbins));
                FrameProfAx.Children(end-(2*(SetIndex-1))).XData = APLength*APbins;
                FrameProfAx.Children(end-(2*(SetIndex-1))).YData = zeros(1, length(APbins));
                FrameProfAx.Children(end-(2*(SetIndex-1))).YPositiveDelta = .1*ones(1, length(APbins));
                FrameProfAx.Children(end-(2*(SetIndex-1))).YNegativeDelta = .1*ones(1, length(APbins));
                set(FrameProfAx.Children(end-(2*(SetIndex-1)+1)),'Visible','off'); %'off' or 'on'
                set(FrameProfAx.Children(end-(2*(SetIndex-1))),'Visible','off'); %'off' or 'on'
                set(get(get(prof{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            else
                PlottedSets(Subset(SetIndex)) = true;
                FrameProfAx.Children(end-(2*(SetIndex-1)+1)).YData = AllNCParams(SetIndex,UseIndex);
                FrameProfAx.Children(end-(2*(SetIndex-1)+1)).XData = APbins(UseIndex);
                FrameProfAx.Children(end-(2*(SetIndex-1))).YData = AllNCParams(SetIndex,UseIndex);
                FrameProfAx.Children(end-(2*(SetIndex-1))).XData = APbins(UseIndex);
                FrameProfAx.Children(end-(2*(SetIndex-1))).YPositiveDelta = AllNCParamSEs(SetIndex,UseIndex);
                FrameProfAx.Children(end-(2*(SetIndex-1))).YNegativeDelta  = AllNCParamSEs(SetIndex,UseIndex);
                
                set(FrameProfAx.Children(end-(2*(SetIndex-1)+1)),'Visible','on'); %'off' or 'on'
                set(FrameProfAx.Children(end-(2*(SetIndex-1))),'Visible','on'); %'off' or 'on'
                set(get(get(prof{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
                
                
            end
        end
        
        if all(~PlottedSets)
            continue
        end
        %try
        if exist('PlotTitle', 'var')
            
            title(FrameProfAx, {PlotTitle,...
                ['Nuclear Cycle ',num2str(NC), ', ', num2str(CurrentTemperature), ' °C']})
            
        else
            
            title(FrameProfAx,  ['Nuclear Cycle ',num2str(NC), ', ', num2str(CurrentTemperature), ' °C'])
            
            
        end
        
        hlegend = legend(legend_labels(PlottedSets), 'Location', 'eastoutside',...
            'FontSize', 14);
        
        saveas(FrameProfFig,[outdir5, filesep,...
                OutputString,'_NC',num2str(NC),'_T', strrep(num2str(CurrentTemperature), '.', '_'), 'C.png']);
        
        
    end
end