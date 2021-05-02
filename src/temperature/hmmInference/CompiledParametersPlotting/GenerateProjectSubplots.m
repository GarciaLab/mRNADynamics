function GenerateProjectSubplots(ResultsPaths, outdir, varargin)
close all

useRescaledTiming = false;
useRescaledFluo = false;
useRescaledParamTiming = false;
includeRescaling = true;


x = 1;
while x <= length(varargin)
    if strcmpi(varargin{x}, 'userescaledtime') | strcmpi(varargin{x}, 'rescaletime') | ...
            strcmpi(varargin{x}, 'rescaletiming') | strcmpi(varargin{x}, 'userescaledtiming')
        useRescaledTiming = true;

    elseif strcmpi(varargin{x}, 'rescalefluo') | strcmpi(varargin{x}, 'userescaledfluo')
        useRescaledFluo = true;

    elseif strcmpi(varargin{x}, 'rescaleparamtiming') | strcmpi(varargin{x}, 'rescaleparametertiming') |...
        strcmpi(varargin{x}, 'userescaledparamtiming') | strcmpi(varargin{x}, 'userescaledparametertiming')   
        useRescaledParamTiming = true;
    elseif strcmpi(varargin{x}, 'excluderescaling') 
        includeRescaling = false;
    end
    x = x+1;
end



if ~exist(outdir, 'dir')
    mkdir(outdir)
end

R = 8.314*10^(-3); % kJ mol^-1 K^-1



APResolution = 0.025;
APbins = (0:APResolution:1)*100;

CompiledParameters = CompileParameters(ResultsPaths, includeRescaling);
NumSets = size(CompiledParameters.InitiationRates, 1);
HasSingleTimeBin = all(CompiledParameters.nTimebins == 1);

if ~HasSingleTimeBin
    if ~useRescaledTiming & ~useRescaledParamTiming
        NumTimeBins  =  length(CompiledParameters.TimeVector);
        IncludedAPbins =  find((sum(squeeze(sum(~isnan(CompiledParameters.InitiationRates), 1)), 1)) > 0);
        IncludedTimeBins =  find((sum(squeeze(sum(~isnan(CompiledParameters.InitiationRates), 1)), 2)) > 0).';
        RescaledString = '';
    else
        NumTimeBins  =  length(CompiledParameters.ScaledTimeVector);
        IncludedAPbins =  find((sum(squeeze(sum(~isnan(CompiledParameters.ScaledInitiationRates), 1)), 1)) > 0);
        IncludedTimeBins =  find((sum(squeeze(sum(~isnan(CompiledParameters.ScaledInitiationRates), 1)), 2)) > 0).';
        RescaledString = 'DevTime';
    end
else
    NumTimeBins = 1;
    IncludedTimeBins = 1;
    IncludedAPbins =  find((sum(squeeze(sum(~isnan(CompiledParameters.InitiationRates), 1)), 2)) > 0).';
    RescaledString = '';
    useRescaledTiming = false;
end

if useRescaledFluo
    FluoString = 'RescaledFluo';
else
    FluoString = '';
end

if useRescaledParamTiming
    ParamString = 'RescaledParams';
    RescaledString = '';
else
    ParamString = '';
end
MinAPbin = min(IncludedAPbins);
MaxAPbin = max(IncludedAPbins);





if max(IncludedTimeBins) > 1
    cmap = brewermap(max(IncludedTimeBins),'Spectral');
else
    cmap = brewermap(2,'Spectral');
end

SubplotPositionInfo = GetSubplotPositioningParameters(NumSets);
hmmVarStrings = { 'InitiationRates','Durations', 'Frequencies','BurstCycleTimes',...
    'InitiationRates','Durations', 'Frequencies','BurstCycleTimes'};
useLogVector = [false, false, false, false, true, true, true, true];


%%
for i = 1:length(hmmVarStrings)
    if useLogVector(i)
        LogString = 'Log';
    else
        LogString = '';
    end
    close all
    [MeanValues, StdErrors, TimeVector, Ymax, Ymin, YLabel, OutString, LogYmax, LogYmin] = ...
        GetPlottingMats(CompiledParameters, hmmVarStrings{i},useRescaledTiming,...
        useRescaledFluo, useRescaledParamTiming, SubplotPositionInfo.SubplotDims(1));
    
    FigAx = cell(1, NumSets);
    
    FigHandle = figure(1);
    
    set(FigHandle,'units', 'normalized', 'position',[0.05, 0.05, SubplotPositionInfo.SubFigDims(1), SubplotPositionInfo.SubFigDims(2)]);
    set(gcf,'color','w');
    
    colormap(cmap);
    hold on
    
    
    % Initiatialize all subplots before messing with any positioning.
    for SetIndex = 1:NumSets
        if SetIndex == 1
            FigAx{SetIndex} = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotPositionInfo.SubplotIndexList(SetIndex));
        else
            FigAx{SetIndex} = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotPositionInfo.SubplotIndexList(SetIndex));
        end
    end
    
    LegendAx = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotPositionInfo.LegendSubplots);
    
    
    
    
    
    
    for SetIndex = 1:NumSets
        SubplotIndex = SubplotPositionInfo.SubplotIndexList(SetIndex);
        set(FigHandle, 'CurrentAxes', FigAx{SetIndex});
        
        PlottedTimes = zeros(1, length(TimeVector), 'logical');
        
        y_vals = squeeze(MeanValues(SetIndex,:,:));
        y_val_errs = squeeze(StdErrors(SetIndex,:,:));
        if ~all(size(y_vals) > 1)
            y_vals = y_vals.';
            y_val_errs = y_val_errs.';
        end
        for t = 1:length(TimeVector)
            
            if all(size(y_vals) > 1)
                if all(isnan(y_vals(t,:)))
                    continue
                end
                set_ids = ~isnan(y_vals(t,:));
                
                x = APbins(set_ids);
                y = y_vals(t, set_ids);
                y_err = y_val_errs(t,set_ids);
            else
                if all(isnan(y_vals))
                    continue
                end
                set_ids = ~isnan(y_vals);
                
                x = APbins(set_ids);
                y = y_vals(set_ids);
                y_err = y_val_errs(set_ids);
            end
                
            
            errorbar(FigAx{SetIndex}, x,y,y_err,'Color','k','Capsize',0)
            hold on
            scatter(FigAx{SetIndex}, x,y,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
            
        end
        
        
        
        
        
        pos = get(FigAx{SetIndex}, 'position');
        pos(1) = SubplotPositionInfo.SubplotXPositions(SetIndex);
        pos(3) = SubplotPositionInfo.SubplotWidth;
        pos(2) = SubplotPositionInfo.SubplotYPositions(SetIndex);
        pos(4) = SubplotPositionInfo.SubplotHeight;
        set(FigAx{SetIndex}, 'position', pos);
        
        xlabel('AP position')
        ylabel(YLabel)
        set(FigAx{SetIndex} ,'Fontsize',14)
        xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
        %xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
        title({[CompiledParameters.ReporterLabels{SetIndex}, ' ',...
            num2str(CompiledParameters.SetTemperatures(SetIndex)), '°C w', num2str(CompiledParameters.nSteps(SetIndex))],...
            ['dt=',...
            num2str(CompiledParameters.dt(SetIndex)), 's, t_{elong}= ',...
            num2str(round(CompiledParameters.ElongationTimes(SetIndex), 2)), ' m']});
        
        
        if useLogVector(i)
            
            newylim = get(FigAx{SetIndex}, 'ylim');
            newylim(1) = LogYmin;
            newylim(2) = LogYmax;
            set(FigAx{SetIndex}, 'ylim', newylim)
            set(FigAx{SetIndex}, 'YScale', 'log');
        else
            newylim = get(FigAx{SetIndex}, 'ylim');
            newylim(1) = Ymin;
            newylim(2) = Ymax;
            set(FigAx{SetIndex}, 'ylim', newylim)
        end
        %     ylim([DurYmin, DurYmax])
        
        grid on
        
        
    end
    
    set(FigHandle, 'CurrentAxes', LegendAx);
    hold on
    axis off
    
    pos = get(LegendAx, 'Position');
    pos(1) = SubplotPositionInfo.LegendXPosition;
    set(LegendAx, 'Position', pos);
    if max(IncludedTimeBins) > 1
        set(LegendAx, 'Clim', [1, max(IncludedTimeBins)+1])
        h = colorbar('west');
        h.Ticks = 1.5:max(IncludedTimeBins)+0.5;
        h.TickLabels = TimeVector(1:max(IncludedTimeBins));
        if useRescaledTiming | useRescaledParamTiming
            ylabel(h,'dev time cohort (normalized minutes into nc14)')
        else
            ylabel(h,'time cohort (minutes into nc14)')
        end
    end
    hold off
    set(LegendAx,'Fontsize',16)
    
    
    if strcmpi(hmmVarStrings{i}, 'initiationrates')
        OutFileName = [LogString, FluoString, RescaledString,ParamString, 'ProjectSubplots_', OutString,'.png'];
    else
        OutFileName = [LogString, RescaledString, ParamString,'ProjectSubplots_', OutString,'.png'];
    end
    
    
    saveas(FigHandle,[outdir, OutFileName])
    
    
    
end
%%

close all

%%










