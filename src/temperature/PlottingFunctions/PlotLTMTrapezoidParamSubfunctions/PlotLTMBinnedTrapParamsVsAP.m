function PlotLTMBinnedTrapParamsVsAP(this, parameter, outdir, varargin)
%%

% PlotTitle, PlottingColors, UseDifferentColors,
% UseDiffProfiles, UsePhysicalAPLength
UsePhysicalAPLength = false;
UseLines = false;
UseRescaledFluo = false;
UseRescaledParamTiming = false;

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
    elseif strcmpi(varargin{x}, 'userescaledtime') | strcmpi(varargin{x}, 'rescaletime') | ...
            strcmpi(varargin{x}, 'rescaletiming') | strcmpi(varargin{x}, 'userescaledtiming') | ...
            strcmpi(varargin{x}, 'userescaledparamtime') | strcmpi(varargin{x}, 'rescaleparamtime') | ...
            strcmpi(varargin{x}, 'rescaleparamtiming') | strcmpi(varargin{x}, 'userescaledparamtiming')
        UseRescaledParamTiming = true;
    elseif strcmpi(varargin{x}, 'rescalefluo') | strcmpi(varargin{x}, 'userescaledfluo')
        UseRescaledFluo = true;
    end
    x = x+1;
end


if ~exist('PlottingColors', 'var')
    PlottingColors = 'default';
    PlottingColors = 'spectral';
elseif ~strcmpi(PlottingColors, 'gradient') &~strcmp(lower(PlottingColors), 'default')  & ~strcmp(lower(PlottingColors), 'pboc')& ~strcmp(lower(PlottingColors), 'spectral')
    error('Invalid choice of plotting colors. Can use either "default", "pboc", or "gradient".') % change to error
end
if ~exist('TraceType', 'var')
    TraceType = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'Fluo3D') | strcmpi(TraceType, 'Unaligned3D')
    TraceType = 'Unaligned3D';
elseif strcmpi(TraceType, 'Fluo')| strcmpi(TraceType, 'Unaligned')
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
if strcmpi(parameter, 'fractionons') | strcmpi(parameter, 'fractionon') 
    UseRescaledFluo = false;
    UseRescaledParamTiming = false;
elseif strcmpi(parameter, 'meanspotfluos') | strcmpi(parameter, 'meanspotfluo')  | ...
        strcmpi(parameter, 'maxfluos') | strcmpi(parameter, 'maxfluo') | ...
        strcmpi(parameter, 'plateauheights') | strcmpi(parameter, 'plateauheight') | ...
        strcmpi(parameter, 'mrnaproductions') | strcmpi(parameter, 'mrnaproduction') |...
        strcmpi(parameter, 'totalmrnaproductions') | strcmpi(parameter, 'totalmrnaproduction') | ...
        strcmpi(parameter, 'unweightedmrnaproductions') | strcmpi(parameter, 'unweightedmrnaproduction') | ...
        strcmpi(parameter, 'unweightedtotalmrnaproductions') | strcmpi(parameter, 'unweightedtotalmrnaproduction')
    
    UseRescaledParamTiming = false;
elseif strcmpi(parameter, 'timeons') | strcmpi(parameter, 'timeon') | ...
        strcmpi(parameter, 'timeoffs') | strcmpi(parameter, 'timeoff') | ...
        strcmpi(parameter, 'elongationtimes') | strcmpi(parameter, 'elongationtime') | ...
        strcmpi(parameter, 'elongationrates') | strcmpi(parameter, 'elongationrate')
    UseRescaledFluo = false;


elseif ~(strcmpi(parameter, 'loadingrates') | strcmpi(parameter, 'loadingrate') | ...
        strcmpi(parameter, 'initiationrates') | strcmpi(parameter, 'initiationrate') |...
        strcmpi(parameter, 'meaninitiationrates') | strcmpi(parameter, 'meaninitiationrate') | ...
        strcmpi(parameter, 'unloadingrates') | strcmpi(parameter, 'unloadingrate'))
    error('Invalid choice of parameter.')
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

if UseRescaledFluo
    FluoString = 'RescaledFluo';
else
    FluoString = '';
end

if UseRescaledParamTiming
    TimingString = 'DevTime';
else
    TimingString = '';
end
Temp_obs = this.Temp_obs;
Temp_sp = this.Temp_sps;
temperatures = flip(unique(this.Temp_sps(this.ProcessedExperiments)));
NumTemperatures = length(temperatures);

if strcmp(lower(PlottingColors), 'spectral')
    colors = brewermap(10,'Spectral');
    colors = [colors(2:3,:);colors(8:10,:)];

    %colormap(cmap);
    %colors =  brewermap(NumTemperatures,'Spectral');
    
    GradString = '';
elseif strcmp(lower(PlottingColors), 'default') | strcmp(lower(PlottingColors), 'gradient')
    [~, colors] = getColorPalettes();
    GradString = '';
elseif strcmp(lower(PlottingColors), 'pboc')
    [colors, ~] = getColorPalettes();
    GradString = '';
end

SpecificDirString = FluoString;
if ~isempty(TimingString)
    if ~isempty(SpecificDirString)
        SpecificDirString = [SpecificDirString, '_', TimingString];
    else
        SpecificDirString = TimingString;
    end
end

if ~isempty(PhysicalAPString)
    if ~isempty(SpecificDirString)
        SpecificDirString = [SpecificDirString, '_', PhysicalAPString];
    else
        SpecificDirString = PhysicalAPString;
    end
end

if ~isempty(GradString)
    if ~isempty(SpecificDirString)
        SpecificDirString = [SpecificDirString, '_', GradString];
    else
        SpecificDirString = GradString;
    end
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
legend_labels = {};
for j =1:NumTemperatures
    legend_labels = [legend_labels, [num2str(temperatures(j)), ' °C']];
end


MarkerStyles = {'o', 'd', 's', '>', '^','p', 'h', '*', 'x'};

TempMatches = cell(1, NumTemperatures);

for t_index = 1:NumTemperatures
    TempMatches{t_index} = find((this.Temp_sps == temperatures(t_index)) & UseSet);
end


%% Load relevant parameters into memory
[PlottedParams, PlottedParamSEs,R2s, ylab,OutputString,GlobalPlotYmax,GlobalPlotYmin,LogPlotYmin] = ...
    getPlottingVariables(this, parameter,  TraceType, R2bound, UseRescaledFluo, UseRescaledParamTiming);

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

outdir2 = [outdir,filesep,'Binned', OutputString];
if ~exist(outdir2, 'dir')
    mkdir(outdir2)
end

if isempty(SpecificDirString)
    outdir3 = [outdir2, filesep, TraceType];
    if ~exist(outdir3, 'dir')
        mkdir(outdir3)
    end
else
    outdir3 = [outdir2, filesep, TraceType, filesep, SpecificDirString];
    if ~exist(outdir3, 'dir')
        mkdir(outdir3)
    end
end
outdir4 = [outdir3, filesep, datestr(now, 'yyyymmdd')];
if ~exist(outdir4, 'dir')
    mkdir(outdir4)
end


eb = cell(1, NumTemperatures);
prof = cell(1, NumTemperatures);
close all
FrameProfFig = figure(1);
set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .9, .7]);
set(gcf,'color','w');
FrameProfAx = axes(FrameProfFig);
for SetIndex =1:NumTemperatures
    ColIndex = SetIndex;
    MarkerIndex = 1;

  
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
box on 
set(FrameProfAx ,'Fontsize',18)
hold off
if ~UsePhysicalAPLength
    xlabel('Fraction Embryo Length (x/L)')
else
    xlabel('Distance from the Anterior Pole (\mum)')
end
xlim([PlotXmin, PlotXmax])
xlim([.1, .475])
ylabel(ylab)
ylim([GlobalPlotYmin, GlobalPlotYmax])
ylim([2 12])
FrameProfAx.YAxis.FontSize = 14; 
FrameProfAx.XAxis.FontSize = 14; 
%title(FrameProfAx, ['Nuclear Cycle '], 'FontSize', 14);

%%

for NCIndex=4:4%1:length(this.IncludedNCs)
    NC = this.IncludedNCs(NCIndex);
    
    
    % Prepare Traces for plotting
    
    NCMaxParams = NaN(1, NumTemperatures);
    AllNCParams = NaN(NumTemperatures, NumAPbins);
    AllNCParamSEs = NaN(NumTemperatures, NumAPbins);
    AllCounts = NaN(NumTemperatures, NumAPbins);
    for SetIndex=1:NumTemperatures
        SetParams = BinnedParams(SetIndex,:,NC-8).';
        SetSEParams = BinnedSEParams(SetIndex,:,NC-8).';
        SetTemperatures = ParamTemperatures(SetIndex,:,NC-8).';
        SetSETemperatures= ParamSETemperatures(SetIndex,:,NC-8).';
        SetCounts = Counts(SetIndex,:,NC-8).';
        if ~all(isnan(SetSEParams))
            IncludedBins = find(~isnan(SetParams) & (SetCounts >= this.MinimumBinCount) & (SetParams./SetSEParams >= 1)) ;
        else
            IncludedBins = find(~isnan(SetParams) & (SetCounts >= this.MinimumBinCount));
        end
        if ~isempty(IncludedBins)
            AllNCParams(SetIndex,:) = SetParams;
            AllNCParamSEs(SetIndex,:) = SetSEParams;
            TempSEParams = SetSEParams;
            TempSEParams(isnan(TempSEParams)) = 0;
            NCMaxParams(SetIndex) = max(SetParams+TempSEParams);
            AllCounts(SetIndex,:) = SetCounts;
        end
    end
    
    if all(isnan(NCMaxParams))
        continue
    end

    
    
    
    
    PlottedSets = zeros(1, NumTemperatures, 'logical');
    for SetIndex=1:NumTemperatures
        if UsePhysicalAPLength
            APLength = this.APLengths(SetIndex);
        else
            APLength = 1;
        end
        if ~all(isnan(AllNCParamSEs(SetIndex, :)))
            UseIndex = ~isnan(AllNCParams(SetIndex,:)) &  (AllNCParams(SetIndex, :)./AllNCParamSEs(SetIndex, :) >= 1) & ...
                (AllCounts(SetIndex,:)>= this.MinimumBinCount);
        else
            UseIndex = ~isnan(AllNCParams(SetIndex,:)) & ...
                (AllCounts(SetIndex,:)>= this.MinimumBinCount);
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
            PlottedSets(SetIndex) = true;
            FrameProfAx.Children(end-(2*(SetIndex-1)+1)).YData = AllNCParams(SetIndex,UseIndex);
            FrameProfAx.Children(end-(2*(SetIndex-1)+1)).XData = APLength*APbins(UseIndex);
            FrameProfAx.Children(end-(2*(SetIndex-1))).YData = AllNCParams(SetIndex,UseIndex);
            FrameProfAx.Children(end-(2*(SetIndex-1))).XData =  APLength*APbins(UseIndex);
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
%     if exist('PlotTitle', 'var')
%         title(FrameProfAx, {PlotTitle,['Nuclear Cycle ',num2str(NC)]}, 'FontSize', 14)
%     else
%         title(FrameProfAx,  ['Nuclear Cycle ',num2str(NC)], 'FontSize', 14)
%     end
    
        
    hlegend = legend(legend_labels(PlottedSets), 'Location', 'eastoutside',...
        'FontSize', 14);

    outpath = [outdir4, filesep,'Binned', OutputString, '_NC',num2str(NC),'.png'];
    saveas(FrameProfFig,outpath);
    
    
end