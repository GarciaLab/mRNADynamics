function PlotLTMvsHMMLoadingRates(this, CompiledParameters, outdir, varargin)
%%

% PlotTitle, PlottingColors, UseDifferentColors,
% UseDiffProfiles, UsePhysicalAPLength
UsePhysicalAPLength = false;
UseLines = true;
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
elseif ~strcmpi(PlottingColors, 'gradient') &~strcmp(lower(PlottingColors), 'default')  & ~strcmp(lower(PlottingColors), 'pboc') & ~strcmp(lower(PlottingColors), 'spectral')
    error('Invalid choice of plotting colors. Can use either "default", "pboc", "spectral" or "gradient".') % change to error
end
if ~exist('TraceType', 'var')
    TraceType = 'Tbinned';
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
parameter = 'loadingrates';

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

if strcmp(lower(PlottingColors), 'default') | strcmp(lower(PlottingColors), 'gradient')
    [~, colors] = getColorPalettes();
    GradString = '';
elseif strcmp(lower(PlottingColors), 'pboc')
    [colors, ~] = getColorPalettes();
    GradString = '';
elseif strcmp(lower(PlottingColors), 'spectral')
    colors = flipud(brewermap(5,'Spectral'));
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

outdir2 = [outdir,filesep,'Binned', OutputString, '_hmmComps'];
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

if ~all(all(all(isnan(BinnedSEParams))))
    BinnedSEParams((Counts < this.MinimumBinCount) | (isnan(BinnedParams)) | (BinnedParams./BinnedSEParams < 1)) = NaN;
    BinnedParams((Counts < this.MinimumBinCount) | (isnan(BinnedParams)) | (BinnedParams./BinnedSEParams < 1)) = NaN;
else
    BinnedSEParams((Counts < this.MinimumBinCount) | (isnan(BinnedParams))) = NaN;
    BinnedParams((Counts < this.MinimumBinCount) | (isnan(BinnedParams))) = NaN;
end




WhereValid = squeeze(sum(~isnan(BinnedParams), 1));
WhereValidAP = sum(WhereValid, 2).';
ValidAPIndices = find(WhereValidAP);
if size(CompiledParameters.InitiationRates, 2) == 1
    cpAPs = find(sum(~isnan(squeeze(CompiledParameters.InitiationRates)), 1));
else
    cpAPs = find(sum(~isnan(squeeze(CompiledParameters.InitiationRates(:,1,:))), 1));
end
ValidAPIndices = ValidAPIndices(ismember(ValidAPIndices, cpAPs));

WhereValidNC = sum(WhereValid, 1);
ValidNCIndices = find(WhereValidNC);
cpNCs = unique(CompiledParameters.NC);
ValidNCIndices = ValidNCIndices(ismember(ValidNCIndices, cpNCs-8));

SubplotPositionInfo = GetSubplotPositioningParameters(length(ValidAPIndices));




[HmmMeanValues, HmmStdErrors, HmmTimeVector, HmmYmax, HmmYmin, HmmYLabel, HmmOutString, HmmLogYmax, HmmLogYmin] = ...
    GetPlottingMats(CompiledParameters, 'InitiationRates',false,...
    UseRescaledFluo, UseRescaledParamTiming,SubplotPositionInfo.SubplotDims(1));

%%
for nc_idx=1:length(ValidNCIndices)
    
    NC = ValidNCIndices(nc_idx)+8;
    % Prepare Traces for plotting
    TempSetPoints = this.UniqueTemperatures;
    
    NCMaxParams = NaN(1,NumTemperatures);
    NCMinParams = NaN(1,NumTemperatures);
    AllNCParams = NaN(NumTemperatures, length(ValidAPIndices));
    AllNCParamSEs = NaN(NumTemperatures, length(ValidAPIndices));
    AllNCTemperatures = NaN(NumTemperatures, length(ValidAPIndices));
    
    hmmNCMaxParams = NaN(1,NumTemperatures);
    hmmNCMinParams = NaN(1,NumTemperatures);
    hmmAllNCParams = NaN(NumTemperatures, length(ValidAPIndices));
    hmmAllNCParamSEs = NaN(NumTemperatures, length(ValidAPIndices));
    hmmAllNCTemperatures = NaN(NumTemperatures, length(ValidAPIndices));
    for i=1:NumTemperatures
        SetParams = BinnedParams(i,ValidAPIndices,NC-8);
        SetSEParams = BinnedSEParams(i,ValidAPIndices,NC-8);
        SetTemps = ones(size(SetParams))*this.Temp_sps(i);
        
        IncludedBins = find(~isnan(SetParams)) ;
        
        if ~isempty(IncludedBins)
            TempParams = SetParams;
            TempSEParams = SetSEParams;
            TempSEParams(TempParams > GlobalPlotYmax) = 0;
            TempParams(TempParams > GlobalPlotYmax) = NaN;
            AllNCParams(i,:) = SetParams;
            AllNCParamSEs(i,:) = SetSEParams;
            TempSEParams(isnan(TempSEParams )) = 0;
            SumParams = TempParams +TempSEParams;
            SumParams(SumParams > GlobalPlotYmax)= 0;
            NCMaxParams(i) = max(SumParams);
            
            LessParams = TempParams -TempSEParams;
            LessParams(LessParams < GlobalPlotYmin)= min(LessParams(~isnan(LessParams)));
            NCMinParams(i) = min(LessParams);
            
            AllNCTemperatures(i,:) =SetTemps;
        end
        
        hmmSetIndex = find(single(CompiledParameters.UniqueTemperatures) == single(TempSetPoints(i)));
        
        hmmSetParams = NaN(1, length(ValidAPIndices));
        hmmSetSEParams = NaN(1, length(ValidAPIndices));
        for l = 1:length(ValidAPIndices)
            timeIndex = find(~isnan(HmmMeanValues(hmmSetIndex,:,ValidAPIndices(l))), 1);
            if ~isempty(timeIndex)
                hmmSetParams(l) = HmmMeanValues(hmmSetIndex,timeIndex,ValidAPIndices(l));
                hmmSetSEParams(l) = HmmStdErrors(hmmSetIndex,timeIndex,ValidAPIndices(l));
            end
        end
        hmmIncludedBins = find(~isnan(hmmSetParams)) ;
        if ~isempty(hmmIncludedBins)
            hmmTempParams = hmmSetParams;
            hmmTempSEParams = hmmSetSEParams;
            hmmTempSEParams(hmmTempParams > HmmYmax) = 0;
            hmmTempParams(hmmTempParams > HmmYmax) = NaN;
            hmmAllNCParams(i,:) = hmmSetParams;
            hmmAllNCParamSEs(i,:) = hmmSetSEParams;
            hmmTempSEParams(isnan(hmmTempSEParams )) = 0;
            hmmSumParams = hmmTempParams +hmmTempSEParams;
            hmmSumParams(hmmSumParams > HmmYmax)= 0;
            hmmNCMaxParams(i) = max(hmmSumParams);
            hmmLessParams = hmmTempParams -hmmTempSEParams;
            hmmLessParams(hmmLessParams < HmmYmin)= min(hmmLessParams(~isnan(hmmLessParams)));
            hmmNCMinParams(i) = min(hmmLessParams);
            
        end
    end
    
    if all(isnan(NCMaxParams)) | all(isnan(hmmNCMaxParams))
        continue
    end
    
    NCPlotXmax = max(NCMaxParams)*1.05;
    NCPlotYmax = max(hmmNCMaxParams)*1.05;
    NCPlotXmin = min(NCMinParams)*0.95;
    NCPlotYmin = min(hmmNCMinParams)*0.95;
    
    NCPlotMin = floor(min(NCPlotXmin, NCPlotYmin)/100)*100;
    NCPlotMax = ceil(max(NCPlotXmax, NCPlotYmax)/100)*100;
    
    eb = cell(NumTemperatures, length(ValidAPIndices));
    prof = cell(NumTemperatures, length(ValidAPIndices));
    FrameProfAx = cell(1,  length(ValidAPIndices));
    hlegend2 = cell(1,  length(ValidAPIndices));
    PlottedSets = zeros(NumTemperatures,length(ValidAPIndices), 'logical');
    close all
    FrameProfFig = figure(1);
    set(FrameProfFig,'units', 'normalized', 'position',[0.05, 0.05, SubplotPositionInfo.SubFigDims(1), SubplotPositionInfo.SubFigDims(2)]);
    set(gcf,'color','w');
    
    colormap(colors);
    hold on
    
    
    % Initiatialize all subplots before messing with any positioning.
    for SubplotIndex = 1:length(ValidAPIndices)
        if SubplotIndex == 1
            FrameProfAx{SubplotIndex} = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotPositionInfo.SubplotIndexList(SubplotIndex));
        else
            FrameProfAx{SubplotIndex} = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotPositionInfo.SubplotIndexList(SubplotIndex));
        end
    end
    
    LegendAx = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotPositionInfo.LegendSubplots);
    
    for l = 1:length(ValidAPIndices)
        APindex = ValidAPIndices(l);
        APbin = APbins(APindex);
        
        SubplotIndex = SubplotPositionInfo.SubplotIndexList(l );
        set(FrameProfFig, 'CurrentAxes', FrameProfAx{l });
        plot([NCPlotMin NCPlotMax], [NCPlotMin NCPlotMax], 'k-')
        hold(FrameProfAx{l}, 'on')
        for idx=1:NumTemperatures
            ColIndex = length(TempSetPoints)+1- idx;
            
            
            if  isnan(AllNCParams(idx, l)) | isnan(hmmAllNCParams(idx, l))
                
                eb{idx, l} = errorbar(AllNCParams(idx, l), hmmAllNCParams(idx, l),...
                    hmmAllNCParamSEs(idx, l),hmmAllNCParamSEs(idx, l),...
                    AllNCParamSEs(idx, l),AllNCParamSEs(idx, l),...
                    'LineStyle', 'none', 'Color', 'k', 'Capsize', 0);
                
                hold(FrameProfAx{l}, 'on')
                set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
                
                prof{idx,l} = scatter(AllNCParams(idx, l), hmmAllNCParams(idx, l),...
                    'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:));
                
                % 'MarkerSize', 8);
                
                set(eb{idx,l},'Visible','off'); %'off' or 'on'
                set(prof{idx,l},'Visible','off'); %'off' or 'on'
                
                
            else
                PlottedSets(idx, l) = true;
                eb{idx, l} = errorbar(AllNCParams(idx, l), hmmAllNCParams(idx, l),...
                    hmmAllNCParamSEs(idx, l),hmmAllNCParamSEs(idx, l),...
                    AllNCParamSEs(idx, l),AllNCParamSEs(idx, l),...
                    'LineStyle', 'none', 'Color', 'k', 'Capsize', 0);
                
                hold(FrameProfAx{l}, 'on')
                set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
                
                prof{idx,l} = scatter(AllNCParams(idx, l), hmmAllNCParams(idx, l),75,...
                    'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:));
                
                
                
                
                
                
                set(eb{idx,l},'Visible','on'); %'off' or 'on'
                set(prof{idx,l},'Visible','on'); %'off' or 'on'
                
            end
        end
        
        
        grid on
        
        SubplotIndex = SubplotPositionInfo.SubplotIndexList(l);
        pos = get(FrameProfAx{l}, 'position');
        pos(1) = SubplotPositionInfo.SubplotXPositions(l);
        pos(3) = SubplotPositionInfo.SubplotWidth;
        pos(2) = SubplotPositionInfo.SubplotYPositions(l);
        pos(4) = SubplotPositionInfo.SubplotHeight;
        set(FrameProfAx{l}, 'position', pos);
        
        
        
        hold(FrameProfAx{l}, 'off')
        
        xlim([NCPlotMin NCPlotMax])
        ylim([NCPlotMin NCPlotMax])
        if ~UseRescaledFluo & ~UseRescaledParamTiming
            xlabel({'Trapezoid Loading','Rates (AU/min)'})
            ylabel({'HMM Loading','Rates (AU/min)'})
        elseif ~UseRescaledFluo
            xlabel({'Trapezoid Loading Rates','(AU/min at 25 ºC)'})
            ylabel({'HMM Loading Rates','(AU/min at 25 ºC)'})
        elseif ~UseRescaledParamTiming
            xlabel({'Fluo-Rescaled Trapezoid','Loading Rates (AU/min)'})
            ylabel({'Fluo-Rescaled HMM','Loading Rates (AU/min)'})
        else
            xlabel({'Fluo-Rescaled Trapezoid','Loading Rates (AU/min at 25 ºC)'})
            ylabel({'Fluo-Rescaled HMM','Loading Rates (AU/min at 25 ºC)'})
            
        end
        
        
        
        
        FrameProfAx{l}.YAxis.FontSize = 14;
        FrameProfAx{l}.XAxis.FontSize = 14;
        
        title(['AP: ', num2str(round(APbin, 3))])
        FrameProfAx{l}.Title.FontSize = 16;
    end
    
    
    
    set(FrameProfFig, 'CurrentAxes', LegendAx);
    hold on
    axis off
    
    pos = get(LegendAx, 'Position');
    pos(1) = SubplotPositionInfo.LegendXPosition;
    set(LegendAx, 'Position', pos);
    
    set(LegendAx, 'Clim', [1, NumTemperatures+1])
    h = colorbar('west');
    h.Ticks = 1.5:NumTemperatures+0.5;
    h.TickLabels = round(fliplr(TempSetPoints), 1);
    ylabel(h,' Temperature (ºC)')
    
    
    
    hold off
    set(LegendAx,'Fontsize',16)
    
    %
    
    outpath = [outdir4, filesep,'LoadingRateCompSubplots', '_NC',num2str(NC),'.png'];
    saveas(FrameProfFig,outpath);
    
    eb = cell(NumTemperatures, length(ValidAPIndices));
    prof = cell(NumTemperatures, length(ValidAPIndices));
    PlottedSets = zeros(NumTemperatures,length(ValidAPIndices), 'logical');
    close all
    FrameProfFig = figure(1);
    FrameProfAx = axes(FrameProfFig);
    set(FrameProfFig,'units', 'normalized', 'position',[0.05, 0.05, SubplotPositionInfo.SubFigDims(1), SubplotPositionInfo.SubFigDims(2)]);
    set(gcf,'color','w');
    
    colormap(colors);
    
    set(FrameProfFig, 'CurrentAxes', FrameProfAx);
    plot([NCPlotMin NCPlotMax], [NCPlotMin NCPlotMax], 'k-')
    hold(FrameProfAx, 'on')
    
    % Initiatialize all subplots before messing with any positioning.
    
    
    for l = 1:length(ValidAPIndices)
        APindex = ValidAPIndices(l);
        APbin = APbins(APindex);
        
        
        
        for idx=1:NumTemperatures
            ColIndex = length(TempSetPoints)+1- idx;
            
            
            if  isnan(AllNCParams(idx, l)) | isnan(hmmAllNCParams(idx, l))
                
                eb{idx, l} = errorbar(AllNCParams(idx, l), hmmAllNCParams(idx, l),...
                    hmmAllNCParamSEs(idx, l),hmmAllNCParamSEs(idx, l),...
                    AllNCParamSEs(idx, l),AllNCParamSEs(idx, l),...
                    'LineStyle', 'none', 'Color', 'k', 'Capsize', 0);
                
                hold(FrameProfAx, 'on')
                set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
                
                prof{idx,l} = scatter(AllNCParams(idx, l), hmmAllNCParams(idx, l),...
                    'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:));
                
                % 'MarkerSize', 8);
                
                set(eb{idx,l},'Visible','off'); %'off' or 'on'
                set(prof{idx,l},'Visible','off'); %'off' or 'on'
                
                
            else
                PlottedSets(idx, l) = true;
                eb{idx, l} = errorbar(AllNCParams(idx, l), hmmAllNCParams(idx, l),...
                    hmmAllNCParamSEs(idx, l),hmmAllNCParamSEs(idx, l),...
                    AllNCParamSEs(idx, l),AllNCParamSEs(idx, l),...
                    'LineStyle', 'none', 'Color', 'k', 'Capsize', 0);
                
                hold(FrameProfAx, 'on')
                set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
                
                prof{idx,l} = scatter(AllNCParams(idx, l), hmmAllNCParams(idx, l),150,...
                    'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:));
                
                
                
                
                
                
                set(eb{idx,l},'Visible','on'); %'off' or 'on'
                set(prof{idx,l},'Visible','on'); %'off' or 'on'
                
            end
        end
    end
    
    grid on
    
    
    
    
    hold(FrameProfAx, 'off')
    
    xlim([NCPlotMin NCPlotMax])
    ylim([NCPlotMin NCPlotMax])
    if ~UseRescaledFluo & ~UseRescaledParamTiming
        xlabel('Trapezoid Loading Rates (AU/min)')
        ylabel('HMM Loading Rates (AU/min)')
    elseif ~UseRescaledFluo
        xlabel('Trapezoid Loading Rates (AU/min at 25 ºC)')
        ylabel('HMM Loading Rates (AU/min at 25 ºC)')
    elseif ~UseRescaledParamTiming
        xlabel('Fluo-Rescaled Trapezoid Loading Rates (AU/min)')
        ylabel('Fluo-Rescaled HMM Loading Rates (AU/min)')
    else
        xlabel('Fluo-Rescaled Trapezoid Loading Rates (AU/min at 25 ºC)')
        ylabel('Fluo-Rescaled HMM Loading Rates (AU/min at 25 ºC)')
        
    end
    
    
    
    
    FrameProfAx.YAxis.FontSize = 16;
    FrameProfAx.XAxis.FontSize = 16;
    
    
    
    
    
  
    
    set(FrameProfAx, 'Clim', [1, NumTemperatures+1])
    h = colorbar;
    h.Ticks = 1.5:NumTemperatures+0.5;
    h.TickLabels = round(fliplr(TempSetPoints), 1);
    ylabel(h,' Temperature (ºC)')
    
    
    
    hold off
    set(FrameProfAx,'Fontsize',18)
    
    outpath = [outdir4, filesep,'LoadingRateComp', '_NC',num2str(NC),'.png'];
    saveas(FrameProfFig,outpath);
end
close all
%%
[PlottedParams, PlottedParamSEs,R2s, ylab,OutputString,GlobalPlotYmax,GlobalPlotYmin,LogPlotYmin] = ...
    getPlottingVariables(this, 'ElongationTimes',  TraceType, R2bound, UseRescaledFluo, UseRescaledParamTiming);

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

outdir2 = [outdir,filesep,'Binned', OutputString, '_hmmComps'];
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

if ~all(all(all(isnan(BinnedSEParams))))
    BinnedSEParams((Counts < this.MinimumBinCount) | (isnan(BinnedParams)) | (BinnedParams./BinnedSEParams < 1)) = NaN;
    BinnedParams((Counts < this.MinimumBinCount) | (isnan(BinnedParams)) | (BinnedParams./BinnedSEParams < 1)) = NaN;
else
    BinnedSEParams((Counts < this.MinimumBinCount) | (isnan(BinnedParams))) = NaN;
    BinnedParams((Counts < this.MinimumBinCount) | (isnan(BinnedParams))) = NaN;
end




WhereValid = squeeze(sum(~isnan(BinnedParams), 1));
WhereValidAP = sum(WhereValid, 2).';
ValidAPIndices = find(WhereValidAP);


WhereValidNC = sum(WhereValid, 1);
ValidNCIndices = find(WhereValidNC);
cpNCs = unique(CompiledParameters.NC);
ValidNCIndices = ValidNCIndices(ismember(ValidNCIndices, cpNCs-8));

SubplotPositionInfo = GetSubplotPositioningParameters(length(ValidAPIndices));

HmmMeanValues = CompiledParameters.ElongationTimes;
HmmStdErrors = NaN(1, length(HmmMeanValues));
HmmYmax = max(HmmMeanValues);
HmmYmin = min(HmmMeanValues);



%%
for nc_idx=1:length(ValidNCIndices)
    
    NC = ValidNCIndices(nc_idx)+8;
    % Prepare Traces for plotting
    TempSetPoints = this.UniqueTemperatures;
    
    NCMaxParams = NaN(1,NumTemperatures);
    NCMinParams = NaN(1,NumTemperatures);
    AllNCParams = NaN(NumTemperatures, length(ValidAPIndices));
    AllNCParamSEs = NaN(NumTemperatures, length(ValidAPIndices));
    AllNCTemperatures = NaN(NumTemperatures, length(ValidAPIndices));
    
    hmmNCMaxParams = NaN(1,NumTemperatures);
    hmmNCMinParams = NaN(1,NumTemperatures);
    hmmAllNCParams = NaN(NumTemperatures, length(ValidAPIndices));
    hmmAllNCParamSEs = NaN(NumTemperatures, length(ValidAPIndices));
    hmmAllNCTemperatures = NaN(NumTemperatures, length(ValidAPIndices));
    for i=1:NumTemperatures
        SetParams = BinnedParams(i,ValidAPIndices,NC-8);
        SetSEParams = BinnedSEParams(i,ValidAPIndices,NC-8);
        SetTemps = ones(size(SetParams))*this.Temp_sps(i);
        
        IncludedBins = find(~isnan(SetParams)) ;
        
        if ~isempty(IncludedBins)
            TempParams = SetParams;
            TempSEParams = SetSEParams;
            TempSEParams(TempParams > GlobalPlotYmax) = 0;
            TempParams(TempParams > GlobalPlotYmax) = NaN;
            AllNCParams(i,:) = SetParams;
            AllNCParamSEs(i,:) = SetSEParams;
            TempSEParams(isnan(TempSEParams )) = 0;
            SumParams = TempParams +TempSEParams;
            SumParams(SumParams > GlobalPlotYmax)= 0;
            NCMaxParams(i) = max(SumParams);
            
            LessParams = TempParams -TempSEParams;
            LessParams(LessParams < GlobalPlotYmin)= min(LessParams(~isnan(LessParams)));
            NCMinParams(i) = min(LessParams);
            
            AllNCTemperatures(i,:) =SetTemps;
        end
        
        hmmSetIndex = find(single(CompiledParameters.UniqueTemperatures) == single(TempSetPoints(i)));
        
        hmmSetParams = ones(1, length(ValidAPIndices))*HmmMeanValues(hmmSetIndex);
        hmmSetSEParams = NaN(1, length(ValidAPIndices));
       
        hmmIncludedBins = find(~isnan(hmmSetParams)) ;
        if ~isempty(hmmIncludedBins)
            hmmTempParams = hmmSetParams;
            hmmTempSEParams = hmmSetSEParams;
            hmmTempSEParams(hmmTempParams > HmmYmax) = 0;
            hmmTempParams(hmmTempParams > HmmYmax) = NaN;
            hmmAllNCParams(i,:) = hmmSetParams;
            hmmAllNCParamSEs(i,:) = hmmSetSEParams;
            hmmTempSEParams(isnan(hmmTempSEParams )) = 0;
            hmmSumParams = hmmTempParams +hmmTempSEParams;
            hmmSumParams(hmmSumParams > HmmYmax)= 0;
            hmmNCMaxParams(i) = max(hmmSumParams);
            hmmLessParams = hmmTempParams -hmmTempSEParams;
            hmmLessParams(hmmLessParams < HmmYmin)= min(hmmLessParams(~isnan(hmmLessParams)));
            hmmNCMinParams(i) = min(hmmLessParams);
            
        end
    end
    
    if all(isnan(NCMaxParams)) | all(isnan(hmmNCMaxParams))
        continue
    end
    
    NCPlotXmax = max(NCMaxParams)*1.05;
    NCPlotYmax = max(hmmNCMaxParams)*1.05;
    NCPlotXmin = min(NCMinParams)*0.95;
    NCPlotYmin = min(hmmNCMinParams)*0.95;
    
    NCPlotMin = floor(min(NCPlotXmin, NCPlotYmin)/2)*2;
    NCPlotMax = ceil(max(NCPlotXmax, NCPlotYmax)/2)*2;
    
    eb = cell(NumTemperatures, length(ValidAPIndices));
    prof = cell(NumTemperatures, length(ValidAPIndices));
    FrameProfAx = cell(1,  length(ValidAPIndices));
    hlegend2 = cell(1,  length(ValidAPIndices));
    PlottedSets = zeros(NumTemperatures,length(ValidAPIndices), 'logical');
    close all
    FrameProfFig = figure(1);
    set(FrameProfFig,'units', 'normalized', 'position',[0.05, 0.05, SubplotPositionInfo.SubFigDims(1), SubplotPositionInfo.SubFigDims(2)]);
    set(gcf,'color','w');
    
    colormap(colors);
    hold on
    
    
    % Initiatialize all subplots before messing with any positioning.
    for SubplotIndex = 1:length(ValidAPIndices)
        if SubplotIndex == 1
            FrameProfAx{SubplotIndex} = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotPositionInfo.SubplotIndexList(SubplotIndex));
        else
            FrameProfAx{SubplotIndex} = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotPositionInfo.SubplotIndexList(SubplotIndex));
        end
    end
    
    LegendAx = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotPositionInfo.LegendSubplots);
    
    for l = 1:length(ValidAPIndices)
        APindex = ValidAPIndices(l);
        APbin = APbins(APindex);
        
        SubplotIndex = SubplotPositionInfo.SubplotIndexList(l );
        set(FrameProfFig, 'CurrentAxes', FrameProfAx{l });
        plot([NCPlotMin NCPlotMax], [NCPlotMin NCPlotMax], 'k-')
        hold(FrameProfAx{l}, 'on')
        for idx=1:NumTemperatures
            ColIndex = length(TempSetPoints)+1- idx;
            
            
            if  isnan(AllNCParams(idx, l)) | isnan(hmmAllNCParams(idx, l))
                
                eb{idx, l} = errorbar(AllNCParams(idx, l), hmmAllNCParams(idx, l),...
                    hmmAllNCParamSEs(idx, l),hmmAllNCParamSEs(idx, l),...
                    AllNCParamSEs(idx, l),AllNCParamSEs(idx, l),...
                    'LineStyle', 'none', 'Color', 'k', 'Capsize', 0);
                
                hold(FrameProfAx{l}, 'on')
                set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
                
                prof{idx,l} = scatter(AllNCParams(idx, l), hmmAllNCParams(idx, l),75,...
                    'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:));
                
                % 'MarkerSize', 8);
                
                set(eb{idx,l},'Visible','off'); %'off' or 'on'
                set(prof{idx,l},'Visible','off'); %'off' or 'on'
                
                
            else
                PlottedSets(idx, l) = true;
                eb{idx, l} = errorbar(AllNCParams(idx, l), hmmAllNCParams(idx, l),...
                    hmmAllNCParamSEs(idx, l),hmmAllNCParamSEs(idx, l),...
                    AllNCParamSEs(idx, l),AllNCParamSEs(idx, l),...
                    'LineStyle', 'none', 'Color', 'k', 'Capsize', 0);
                
                hold(FrameProfAx{l}, 'on')
                set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
                
                prof{idx,l} = scatter(AllNCParams(idx, l), hmmAllNCParams(idx, l),75,...
                    'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:));
                
                
                
                
                
                
                set(eb{idx,l},'Visible','on'); %'off' or 'on'
                set(prof{idx,l},'Visible','on'); %'off' or 'on'
                
            end
        end
        
        
        grid on
        
        SubplotIndex = SubplotPositionInfo.SubplotIndexList(l);
        pos = get(FrameProfAx{l}, 'position');
        pos(1) = SubplotPositionInfo.SubplotXPositions(l);
        pos(3) = SubplotPositionInfo.SubplotWidth;
        pos(2) = SubplotPositionInfo.SubplotYPositions(l);
        pos(4) = SubplotPositionInfo.SubplotHeight;
        set(FrameProfAx{l}, 'position', pos);
        
        
        
        hold(FrameProfAx{l}, 'off')
        
        xlim([NCPlotMin NCPlotMax])
        ylim([NCPlotMin NCPlotMax])
        if ~UseRescaledParamTiming
            xlabel({'Trapezoid Elongation','Times (min)'})
            ylabel({'HMM Elongation','Times (min)'})
        else
            xlabel({'Trapezoid Elongation','Times (min at 25 ºC)'})
            ylabel({'HMM Elongation','Times (min at 25 ºC)'})
            
        end
        
        
        
        
        FrameProfAx{l}.YAxis.FontSize = 14;
        FrameProfAx{l}.XAxis.FontSize = 14;
        
        title(['AP: ', num2str(round(APbin, 3))])
        FrameProfAx{l}.Title.FontSize = 16;
    end
    
    
    
    set(FrameProfFig, 'CurrentAxes', LegendAx);
    hold on
    axis off
    
    pos = get(LegendAx, 'Position');
    pos(1) = SubplotPositionInfo.LegendXPosition;
    set(LegendAx, 'Position', pos);
    
    set(LegendAx, 'Clim', [1, NumTemperatures+1])
    h = colorbar('west');
    h.Ticks = 1.5:NumTemperatures+0.5;
    h.TickLabels = round(fliplr(TempSetPoints), 1);
    ylabel(h,' Temperature (ºC)')
    
    
    
    hold off
    set(LegendAx,'Fontsize',16)
    
    %
    
    outpath = [outdir4, filesep,'ElongationTimesCompSubplots', '_NC',num2str(NC),'.png'];
    saveas(FrameProfFig,outpath);
    
    eb = cell(NumTemperatures, length(ValidAPIndices));
    prof = cell(NumTemperatures, length(ValidAPIndices));
    PlottedSets = zeros(NumTemperatures,length(ValidAPIndices), 'logical');
    close all
    FrameProfFig = figure(1);
    FrameProfAx = axes(FrameProfFig);
    set(FrameProfFig,'units', 'normalized', 'position',[0.05, 0.05, SubplotPositionInfo.SubFigDims(1), SubplotPositionInfo.SubFigDims(2)]);
    set(gcf,'color','w');
    
    colormap(colors);
    
    set(FrameProfFig, 'CurrentAxes', FrameProfAx);
    plot([NCPlotMin NCPlotMax], [NCPlotMin NCPlotMax], 'k-')
    hold(FrameProfAx, 'on')
    
    % Initiatialize all subplots before messing with any positioning.
    
    
    
    for idx=1:NumTemperatures
        ColIndex = length(TempSetPoints)+1- idx;
        AllNCParamsVector = AllNCParams(idx,:);
        AllNCParamSEsVector = AllNCParamSEs(idx,:);
        AllNCParamSEsVector = AllNCParamSEsVector(~isnan(AllNCParamsVector));
        AllNCParamSEsVector(isnan(AllNCParamSEsVector)) = 0;
        AllNCParamsVector = AllNCParamsVector(~isnan(AllNCParamsVector));
        if ~isempty(AllNCParamsVector)
            CurrentNCParam = mean(AllNCParamsVector);
        else
            CurrentNCParam = NaN;
        end
        CurrentNCParamSE = sqrt(sum(AllNCParamSEsVector.^2))/length(AllNCParamSEsVector);
        if  isnan(CurrentNCParam) | isnan(hmmAllNCParams(idx, 1))
            
            eb{idx, 1} = errorbar(CurrentNCParam, hmmAllNCParams(idx, 1),...
                hmmAllNCParamSEs(idx, 1),hmmAllNCParamSEs(idx, 1),...
                CurrentNCParamSE,CurrentNCParamSE,...
                'LineStyle', 'none', 'Color', 'k', 'Capsize', 0);
            
            hold(FrameProfAx, 'on')
            set(get(get(eb{idx,1}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            
            
            prof{idx,1} = scatter(CurrentNCParam, hmmAllNCParams(idx, 1),...
                'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:));
            
            % 'MarkerSize', 8);
            
            set(eb{idx,1},'Visible','off'); %'off' or 'on'
            set(prof{idx,1},'Visible','off'); %'off' or 'on'
            
            
        else
            eb{idx, 1} = errorbar(CurrentNCParam, hmmAllNCParams(idx, 1),...
                hmmAllNCParamSEs(idx, 1),hmmAllNCParamSEs(idx, 1),...
                CurrentNCParamSE,CurrentNCParamSE,...
                'LineStyle', 'none', 'Color', 'k', 'Capsize', 0);
            
            hold(FrameProfAx, 'on')
            set(get(get(eb{idx,1}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            
            
            prof{idx,1} = scatter(CurrentNCParam, hmmAllNCParams(idx, 1),150,...
                'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:));
            
            
            
            
            
            
            set(eb{idx,1},'Visible','on'); %'off' or 'on'
            set(prof{idx,1},'Visible','on'); %'off' or 'on'
            
        end
    end
    
    grid on
    
    
    
    
    hold(FrameProfAx, 'off')
    
    xlim([NCPlotMin NCPlotMax])
    ylim([NCPlotMin NCPlotMax])
    if ~UseRescaledFluo 
        xlabel('Trapezoid Elongation Times (min)')
        ylabel('HMM Elongation Times (min)')
    else
        xlabel('Trapezoid Elongation Times (min at 25 ºC)')
        ylabel('HMM Elongation Times (min at 25 ºC)')

    end
    
    
    
    
    FrameProfAx.YAxis.FontSize = 16;
    FrameProfAx.XAxis.FontSize = 16;
    
    
    
    
    
  
    
    set(FrameProfAx, 'Clim', [1, NumTemperatures+1])
    h = colorbar;
    h.Ticks = 1.5:NumTemperatures+0.5;
    h.TickLabels = round(fliplr(TempSetPoints), 1);
    ylabel(h,' Temperature (ºC)')
    
    
    
    hold off
    set(FrameProfAx,'Fontsize',18)
    
    outpath = [outdir4, filesep,'ElongationTimesComp', '_NC',num2str(NC),'.png'];
    saveas(FrameProfFig,outpath);
end




%%
close all
[PlottedParams, PlottedParamSEs,R2s, ylab,OutputString,GlobalPlotYmax,GlobalPlotYmin,LogPlotYmin] = ...
    getPlottingVariables(this, 'InitiationRates',  TraceType, R2bound, UseRescaledFluo, UseRescaledParamTiming);

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

outdir2 = [outdir,filesep,'Binned', OutputString, '_hmmComps'];
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

if ~all(all(all(isnan(BinnedSEParams))))
    BinnedSEParams((Counts < this.MinimumBinCount) | (isnan(BinnedParams)) | (BinnedParams./BinnedSEParams < 1)) = NaN;
    BinnedParams((Counts < this.MinimumBinCount) | (isnan(BinnedParams)) | (BinnedParams./BinnedSEParams < 1)) = NaN;
else
    BinnedSEParams((Counts < this.MinimumBinCount) | (isnan(BinnedParams))) = NaN;
    BinnedParams((Counts < this.MinimumBinCount) | (isnan(BinnedParams))) = NaN;
end




WhereValid = squeeze(sum(~isnan(BinnedParams), 1));
WhereValidAP = sum(WhereValid, 2).';
ValidAPIndices = find(WhereValidAP);
if size(CompiledParameters.InitiationRates, 2) == 1
    cpAPs = find(sum(~isnan(squeeze(CompiledParameters.InitiationRates)), 1));
else
    cpAPs = find(sum(~isnan(squeeze(CompiledParameters.InitiationRates(:,1,:))), 1));
end
ValidAPIndices = ValidAPIndices(ismember(ValidAPIndices, cpAPs));

WhereValidNC = sum(WhereValid, 1);
ValidNCIndices = find(WhereValidNC);
cpNCs = unique(CompiledParameters.NC);
ValidNCIndices = ValidNCIndices(ismember(ValidNCIndices, cpNCs-8));

SubplotPositionInfo = GetSubplotPositioningParameters(length(ValidAPIndices));


[HmmMeanValues, HmmStdErrors, HmmTimeVector, HmmYmax, HmmYmin, HmmYLabel, HmmOutString, HmmLogYmax, HmmLogYmin] = ...
    GetPlottingMats(CompiledParameters, 'MeanInitiationRates',false,...
    UseRescaledFluo, UseRescaledParamTiming,SubplotPositionInfo.SubplotDims(1));

%%
for nc_idx=1:length(ValidNCIndices)
    
    NC = ValidNCIndices(nc_idx)+8;
    % Prepare Traces for plotting
    TempSetPoints = this.UniqueTemperatures;
    
    NCMaxParams = NaN(1,NumTemperatures);
    NCMinParams = NaN(1,NumTemperatures);
    AllNCParams = NaN(NumTemperatures, length(ValidAPIndices));
    AllNCParamSEs = NaN(NumTemperatures, length(ValidAPIndices));
    AllNCTemperatures = NaN(NumTemperatures, length(ValidAPIndices));
    
    hmmNCMaxParams = NaN(1,NumTemperatures);
    hmmNCMinParams = NaN(1,NumTemperatures);
    hmmAllNCParams = NaN(NumTemperatures, length(ValidAPIndices));
    hmmAllNCParamSEs = NaN(NumTemperatures, length(ValidAPIndices));
    hmmAllNCTemperatures = NaN(NumTemperatures, length(ValidAPIndices));
    for i=1:NumTemperatures
        SetParams = BinnedParams(i,ValidAPIndices,NC-8);
        SetSEParams = BinnedSEParams(i,ValidAPIndices,NC-8);
        SetTemps = ones(size(SetParams))*this.Temp_sps(i);
        
        IncludedBins = find(~isnan(SetParams)) ;
        
        if ~isempty(IncludedBins)
            TempParams = SetParams;
            TempSEParams = SetSEParams;
            TempSEParams(TempParams > GlobalPlotYmax) = 0;
            TempParams(TempParams > GlobalPlotYmax) = NaN;
            AllNCParams(i,:) = SetParams;
            AllNCParamSEs(i,:) = SetSEParams;
            TempSEParams(isnan(TempSEParams )) = 0;
            SumParams = TempParams +TempSEParams;
            SumParams(SumParams > GlobalPlotYmax)= 0;
            NCMaxParams(i) = max(SumParams);
            
            LessParams = TempParams -TempSEParams;
            LessParams(LessParams < GlobalPlotYmin)= min(LessParams(~isnan(LessParams)));
            NCMinParams(i) = min(LessParams);
            
            AllNCTemperatures(i,:) =SetTemps;
        end
        
        hmmSetIndex = find(single(CompiledParameters.UniqueTemperatures) == single(TempSetPoints(i)));
        
        hmmSetParams = NaN(1, length(ValidAPIndices));
        hmmSetSEParams = NaN(1, length(ValidAPIndices));
        for l = 1:length(ValidAPIndices)
            timeIndex = find(~isnan(HmmMeanValues(hmmSetIndex,:,ValidAPIndices(l))), 1);
            if ~isempty(timeIndex)
                hmmSetParams(l) = HmmMeanValues(hmmSetIndex,timeIndex,ValidAPIndices(l));
                hmmSetSEParams(l) = HmmStdErrors(hmmSetIndex,timeIndex,ValidAPIndices(l));
            end
        end
        hmmIncludedBins = find(~isnan(hmmSetParams)) ;
        if ~isempty(hmmIncludedBins)
            hmmTempParams = hmmSetParams;
            hmmTempSEParams = hmmSetSEParams;
            hmmTempSEParams(hmmTempParams > HmmYmax) = 0;
            hmmTempParams(hmmTempParams > HmmYmax) = NaN;
            hmmAllNCParams(i,:) = hmmSetParams;
            hmmAllNCParamSEs(i,:) = hmmSetSEParams;
            hmmTempSEParams(isnan(hmmTempSEParams )) = 0;
            hmmSumParams = hmmTempParams +hmmTempSEParams;
            hmmSumParams(hmmSumParams > HmmYmax)= 0;
            hmmNCMaxParams(i) = max(hmmSumParams);
            hmmLessParams = hmmTempParams -hmmTempSEParams;
            hmmLessParams(hmmLessParams < HmmYmin)= min(hmmLessParams(~isnan(hmmLessParams)));
            hmmNCMinParams(i) = min(hmmLessParams);
            
        end
    end
    
    if all(isnan(NCMaxParams)) | all(isnan(hmmNCMaxParams))
        continue
    end
    
    NCPlotXmax = max(NCMaxParams)*1.05;
    NCPlotYmax = max(hmmNCMaxParams)*1.05;
    NCPlotXmin = min(NCMinParams)*0.95;
    NCPlotYmin = min(hmmNCMinParams)*0.95;
    
    NCPlotMin = floor(min(NCPlotXmin, NCPlotYmin)/100)*100;
    NCPlotMax = ceil(max(NCPlotXmax, NCPlotYmax)/100)*100;
    
    eb = cell(NumTemperatures, length(ValidAPIndices));
    prof = cell(NumTemperatures, length(ValidAPIndices));
    FrameProfAx = cell(1,  length(ValidAPIndices));
    hlegend2 = cell(1,  length(ValidAPIndices));
    PlottedSets = zeros(NumTemperatures,length(ValidAPIndices), 'logical');
    close all
    FrameProfFig = figure(1);
    set(FrameProfFig,'units', 'normalized', 'position',[0.05, 0.05, SubplotPositionInfo.SubFigDims(1), SubplotPositionInfo.SubFigDims(2)]);
    set(gcf,'color','w');
    
    colormap(colors);
    hold on
    
    
    % Initiatialize all subplots before messing with any positioning.
    for SubplotIndex = 1:length(ValidAPIndices)
        if SubplotIndex == 1
            FrameProfAx{SubplotIndex} = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotPositionInfo.SubplotIndexList(SubplotIndex));
        else
            FrameProfAx{SubplotIndex} = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotPositionInfo.SubplotIndexList(SubplotIndex));
        end
    end
    
    LegendAx = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotPositionInfo.LegendSubplots);
    
    for l = 1:length(ValidAPIndices)
        APindex = ValidAPIndices(l);
        APbin = APbins(APindex);
        
        SubplotIndex = SubplotPositionInfo.SubplotIndexList(l );
        set(FrameProfFig, 'CurrentAxes', FrameProfAx{l });
        plot([NCPlotMin NCPlotMax], [NCPlotMin NCPlotMax], 'k-')
        hold(FrameProfAx{l}, 'on')
        for idx=1:NumTemperatures
            ColIndex = length(TempSetPoints)+1- idx;
            
            
            if  isnan(AllNCParams(idx, l)) | isnan(hmmAllNCParams(idx, l))
                
                eb{idx, l} = errorbar(AllNCParams(idx, l), hmmAllNCParams(idx, l),...
                    hmmAllNCParamSEs(idx, l),hmmAllNCParamSEs(idx, l),...
                    AllNCParamSEs(idx, l),AllNCParamSEs(idx, l),...
                    'LineStyle', 'none', 'Color', 'k', 'Capsize', 0);
                
                hold(FrameProfAx{l}, 'on')
                set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
                
                prof{idx,l} = scatter(AllNCParams(idx, l), hmmAllNCParams(idx, l),...
                    'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:));
                
                % 'MarkerSize', 8);
                
                set(eb{idx,l},'Visible','off'); %'off' or 'on'
                set(prof{idx,l},'Visible','off'); %'off' or 'on'
                
                
            else
                PlottedSets(idx, l) = true;
                eb{idx, l} = errorbar(AllNCParams(idx, l), hmmAllNCParams(idx, l),...
                    hmmAllNCParamSEs(idx, l),hmmAllNCParamSEs(idx, l),...
                    AllNCParamSEs(idx, l),AllNCParamSEs(idx, l),...
                    'LineStyle', 'none', 'Color', 'k', 'Capsize', 0);
                
                hold(FrameProfAx{l}, 'on')
                set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
                
                prof{idx,l} = scatter(AllNCParams(idx, l), hmmAllNCParams(idx, l),75,...
                    'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:));
                
                
                
                
                
                
                set(eb{idx,l},'Visible','on'); %'off' or 'on'
                set(prof{idx,l},'Visible','on'); %'off' or 'on'
                
            end
        end
        
        
        grid on
        
        SubplotIndex = SubplotPositionInfo.SubplotIndexList(l);
        pos = get(FrameProfAx{l}, 'position');
        pos(1) = SubplotPositionInfo.SubplotXPositions(l);
        pos(3) = SubplotPositionInfo.SubplotWidth;
        pos(2) = SubplotPositionInfo.SubplotYPositions(l);
        pos(4) = SubplotPositionInfo.SubplotHeight;
        set(FrameProfAx{l}, 'position', pos);
        
        
        
        hold(FrameProfAx{l}, 'off')
        
        xlim([NCPlotMin NCPlotMax])
        ylim([NCPlotMin NCPlotMax])
        if ~UseRescaledFluo & ~UseRescaledParamTiming
            xlabel({'Trapezoid Loading','Rates (AU/min)'})
            ylabel({'HMM Mean Loading','Rates (AU/min)'})
        elseif ~UseRescaledFluo
            xlabel({'Trapezoid Loading Rates','(AU/min at 25 ºC)'})
            ylabel({'HMM Mean Loading Rates','(AU/min at 25 ºC)'})
        elseif ~UseRescaledParamTiming
            xlabel({'Fluo-Rescaled Trapezoid','Loading Rates (AU/min)'})
            ylabel({'Fluo-Rescaled Mean HMM','Loading Rates (AU/min)'})
        else
            xlabel({'Fluo-Rescaled Trapezoid','Loading Rates (AU/min at 25 ºC)'})
            ylabel({'Fluo-Rescaled Mean HMM','Loading Rates (AU/min at 25 ºC)'})
            
        end
        
        
        
        
        FrameProfAx{l}.YAxis.FontSize = 14;
        FrameProfAx{l}.XAxis.FontSize = 14;
        
        title(['AP: ', num2str(round(APbin, 3))])
        FrameProfAx{l}.Title.FontSize = 16;
    end
    
    
    
    set(FrameProfFig, 'CurrentAxes', LegendAx);
    hold on
    axis off
    
    pos = get(LegendAx, 'Position');
    pos(1) = SubplotPositionInfo.LegendXPosition;
    set(LegendAx, 'Position', pos);
    
    set(LegendAx, 'Clim', [1, NumTemperatures+1])
    h = colorbar('west');
    h.Ticks = 1.5:NumTemperatures+0.5;
    h.TickLabels = round(fliplr(TempSetPoints), 1);
    ylabel(h,' Temperature (ºC)')
    
    
    
    hold off
    set(LegendAx,'Fontsize',16)
    
    %
    
    outpath = [outdir4, filesep,'MeanLoadingRateCompSubplots', '_NC',num2str(NC),'.png'];
    saveas(FrameProfFig,outpath);
    
    eb = cell(NumTemperatures, length(ValidAPIndices));
    prof = cell(NumTemperatures, length(ValidAPIndices));
    PlottedSets = zeros(NumTemperatures,length(ValidAPIndices), 'logical');
    close all
    FrameProfFig = figure(1);
    FrameProfAx = axes(FrameProfFig);
    set(FrameProfFig,'units', 'normalized', 'position',[0.05, 0.05, SubplotPositionInfo.SubFigDims(1), SubplotPositionInfo.SubFigDims(2)]);
    set(gcf,'color','w');
    
    colormap(colors);
    
    set(FrameProfFig, 'CurrentAxes', FrameProfAx);
    plot([NCPlotMin NCPlotMax], [NCPlotMin NCPlotMax], 'k-')
    hold(FrameProfAx, 'on')
    
    % Initiatialize all subplots before messing with any positioning.
    
    
    for l = 1:length(ValidAPIndices)
        APindex = ValidAPIndices(l);
        APbin = APbins(APindex);
        
        
        
        for idx=1:NumTemperatures
            ColIndex = length(TempSetPoints)+1- idx;
            
            
            if  isnan(AllNCParams(idx, l)) | isnan(hmmAllNCParams(idx, l))
                
                eb{idx, l} = errorbar(AllNCParams(idx, l), hmmAllNCParams(idx, l),...
                    hmmAllNCParamSEs(idx, l),hmmAllNCParamSEs(idx, l),...
                    AllNCParamSEs(idx, l),AllNCParamSEs(idx, l),...
                    'LineStyle', 'none', 'Color', 'k', 'Capsize', 0);
                
                hold(FrameProfAx, 'on')
                set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
                
                prof{idx,l} = scatter(AllNCParams(idx, l), hmmAllNCParams(idx, l),...
                    'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:));
                
                % 'MarkerSize', 8);
                
                set(eb{idx,l},'Visible','off'); %'off' or 'on'
                set(prof{idx,l},'Visible','off'); %'off' or 'on'
                
                
            else
                PlottedSets(idx, l) = true;
                eb{idx, l} = errorbar(AllNCParams(idx, l), hmmAllNCParams(idx, l),...
                    hmmAllNCParamSEs(idx, l),hmmAllNCParamSEs(idx, l),...
                    AllNCParamSEs(idx, l),AllNCParamSEs(idx, l),...
                    'LineStyle', 'none', 'Color', 'k', 'Capsize', 0);
                
                hold(FrameProfAx, 'on')
                set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
                
                prof{idx,l} = scatter(AllNCParams(idx, l), hmmAllNCParams(idx, l),150,...
                    'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:));
                
                
                
                
                
                
                set(eb{idx,l},'Visible','on'); %'off' or 'on'
                set(prof{idx,l},'Visible','on'); %'off' or 'on'
                
            end
        end
    end
    
    grid on
    
    
    
    
    hold(FrameProfAx, 'off')
    
    xlim([NCPlotMin NCPlotMax])
    ylim([NCPlotMin NCPlotMax])
    if ~UseRescaledFluo & ~UseRescaledParamTiming
        xlabel('Trapezoid Loading Rates (AU/min)')
        ylabel('HMM Mean Loading Rates (AU/min)')
    elseif ~UseRescaledFluo
        xlabel('Trapezoid Loading Rates (AU/min at 25 ºC)')
        ylabel('HMM Mean Loading Rates (AU/min at 25 ºC)')
    elseif ~UseRescaledParamTiming
        xlabel('Fluo-Rescaled Trapezoid Loading Rates (AU/min)')
        ylabel('Fluo-Rescaled Mean HMM Loading Rates (AU/min)')
    else
        xlabel('Fluo-Rescaled Trapezoid Loading Rates (AU/min at 25 ºC)')
        ylabel('Fluo-Rescaled Mean HMM Loading Rates (AU/min at 25 ºC)')
        
    end
    
    
    
    
    FrameProfAx.YAxis.FontSize = 16;
    FrameProfAx.XAxis.FontSize = 16;
    
    
    
    
    
  
    
    set(FrameProfAx, 'Clim', [1, NumTemperatures+1])
    h = colorbar;
    h.Ticks = 1.5:NumTemperatures+0.5;
    h.TickLabels = round(fliplr(TempSetPoints), 1);
    ylabel(h,' Temperature (ºC)')
    
    
    
    hold off
    set(FrameProfAx,'Fontsize',18)
    
    outpath = [outdir4, filesep,'MeanLoadingRateComp', '_NC',num2str(NC),'.png'];
    saveas(FrameProfFig,outpath);
end
close all
