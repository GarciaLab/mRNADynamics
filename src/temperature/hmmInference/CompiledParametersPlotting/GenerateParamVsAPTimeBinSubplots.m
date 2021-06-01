function GenerateParamVsAPTimeBinSubplots(ResultsPaths, outdir, varargin)
close all



if ~exist(outdir, 'dir')
    mkdir(outdir)
end

useRescaledTiming = false;
useRescaledFluo = false;
useRescaledParamTiming = false;
PrescaledBins = false;


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
    elseif strcmpi(varargin{x}, 'prescaledbins')
        PrescaledBins = true;
    end
    x = x+1;
end




R = 8.314*10^(-3); % kJ mol^-1 K^-1



APResolution = 0.025;
APbins = (0:APResolution:1)*100;

CompiledParameters = CompileParameters(ResultsPaths, true, PrescaledBins);
ReferenceTemperature = CompiledParameters.RefTemperature;
if single(round(ReferenceTemperature, 0)) == single(ReferenceTemperature)
    TemperatureString = sprintf('%.0f ', ReferenceTemperature);
else
    TemperatureString = sprintf('%.1f ', ReferenceTemperature);
end
NumSets = size(CompiledParameters.InitiationRates, 1);
NumTemperatures = length(CompiledParameters.UniqueTemperatures);

if ~useRescaledTiming & ~useRescaledParamTiming
    NumTimeBins  =  length(CompiledParameters.TimeVector);
    IncludedAPbins =  find((sum(squeeze(sum(~isnan(CompiledParameters.InitiationRates), 1)), 1)) > 0);
    IncludedTimeBins =  find((sum(squeeze(sum(~isnan(CompiledParameters.InitiationRates), 1)), 2)) > 0).';
    NumTimeBins  =  length(IncludedTimeBins);
    RescaledString = '';
else
    NumTimeBins  =  length(CompiledParameters.ScaledTimeVector);
    IncludedAPbins =  find((sum(squeeze(sum(~isnan(CompiledParameters.ScaledInitiationRates), 1)), 1)) > 0);
    IncludedTimeBins =  find((sum(squeeze(sum(~isnan(CompiledParameters.ScaledInitiationRates), 1)), 2)) > 0).';
    NumTimeBins  =  length(IncludedTimeBins);
    RescaledString = 'DevTime';
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
NumAPbins = MaxAPbin-MinAPbin + 1;

if max(IncludedTimeBins) > 1
    cmap = flipud(brewermap(NumTemperatures,'Spectral'));
else
    cmap = flipud(brewermap(2,'Spectral'));
end

cmap = brewermap(10,'Spectral');
cmap = [cmap(2:3,:);cmap(8:10,:)];
SubplotPositionInfo = GetSubplotPositioningParameters(NumTimeBins);
hmmVarStrings = { 'InitiationRates','Durations', 'Frequencies','BurstCycleTimes','MeanInitiationRates'};
%%
for i = 1:length(hmmVarStrings)
    close all
    
    [MeanValues, StdErrors, TimeVector, Ymax, Ymin, YLabel, OutString, LogYmax, LogYmin] = ...
        GetPlottingMats(CompiledParameters, hmmVarStrings{i}, useRescaledTiming, useRescaledFluo,...
        useRescaledParamTiming, SubplotPositionInfo.SubplotDims(1),ReferenceTemperature);
    if strcmpi(hmmVarStrings{i}, 'initiationrates') | strcmpi(hmmVarStrings{i}, 'initiationrates')
        paramYinterval =50;
    elseif strcmpi(hmmVarStrings{i}, 'durations')
        paramYinterval =2;
    elseif strcmpi(hmmVarStrings{i}, 'Frequencies')
        paramYinterval =0.2;
    elseif strcmpi(hmmVarStrings{i}, 'BurstCycleTimes')
        paramYinterval =2;
    end
    dt = TimeVector(2)-TimeVector(1);
    
    ValidTimePoints = TimeVector(find(sum((squeeze(sum(~isnan(MeanValues), 1))), 2).' > 0));
    FigAx = cell(1, NumTimeBins);
    DataPlots = cell(NumTimeBins, NumSets);
    ErrorPlots = cell(NumTimeBins, NumSets);
    
    FigHandle = figure(1);
    
    set(FigHandle,'units', 'normalized', 'position',[0.05, 0.05, SubplotPositionInfo.SubFigDims(1), SubplotPositionInfo.SubFigDims(2)]);
    set(gcf,'color','w');
    
    colormap(cmap);
    hold on
    
    
    % Initiatialize all subplots before messing with any positioning.
    for TimeIndex = 1:NumTimeBins
        TimeBin = IncludedTimeBins(TimeIndex);
        if TimeIndex == 1
            FigAx{TimeIndex} = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotPositionInfo.SubplotIndexList(TimeIndex));
        else
            FigAx{TimeIndex} = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotPositionInfo.SubplotIndexList(TimeIndex));
        end
    end
    
    LegendAx = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotPositionInfo.LegendSubplots);
    
    for TimeIndex = 1:NumTimeBins
        pos = get(FigAx{TimeIndex}, 'position');
        pos(1) = SubplotPositionInfo.SubplotXPositions(TimeIndex);
        pos(3) = SubplotPositionInfo.SubplotWidth;
        pos(2) = SubplotPositionInfo.SubplotYPositions(TimeIndex);
        pos(4) = SubplotPositionInfo.SubplotHeight;
        set(FigAx{TimeIndex}, 'position', pos);
    end
    
    
    
    
    for TimeIndex = 1:NumTimeBins
        TimeBin = IncludedTimeBins(TimeIndex);
        SubplotIndex = SubplotPositionInfo.SubplotIndexList(TimeIndex);
        set(FigHandle, 'CurrentAxes', FigAx{TimeIndex});
        
        TimeBin_y_vals = squeeze(MeanValues(:,TimeBin,:));
        TimeBin_y_val_errs = squeeze(StdErrors(:,TimeBin,:));
        MaxTimeBinValue = max(TimeBin_y_vals(~isnan(TimeBin_y_vals)));
        MinTimeBinValue = min(TimeBin_y_vals(~isnan(TimeBin_y_vals)));
        
        MaxTimeBinValue = max([MaxTimeBinValue*1.3, MaxTimeBinValue*0.7]);
        MinTimeBinValue = min([MinTimeBinValue*1.1, MinTimeBinValue*0.9]);
        
        TimeBinYmax = paramYinterval*ceil(MaxTimeBinValue/paramYinterval);
        TimeBinYmin = paramYinterval*floor(MinTimeBinValue/paramYinterval);
       
        if strcmpi(hmmVarStrings{i}, 'initiationrates')
            LogTimeBinYmax = 10^(ceil(log10(TimeBinYmax)));
            LogTimeBinYmin =TimeBinYmin;
        else
            LogTimeBinYmin =TimeBinYmin;
            LogTimeBinYmax =TimeBinYmax;
        end
        for SetIndex = 1:NumSets
            t = find(round(CompiledParameters.UniqueTemperatures, 1) == round(CompiledParameters.SetTemperatures(SetIndex), 1));
            if all(isnan(TimeBin_y_vals(SetIndex,:)))
                continue
            end
            
            
            set_ids = ~isnan(TimeBin_y_vals(SetIndex,:));
            
            x = APbins(set_ids)/100;

            
            y = TimeBin_y_vals(SetIndex, set_ids);
            y_err = TimeBin_y_val_errs(SetIndex, set_ids);
            

            ErrorPlots{TimeIndex, SetIndex} =  errorbar(FigAx{TimeIndex}, x,y,y_err,'Color','k','Capsize',0);
            hold on
            
            DataPlots{TimeIndex, SetIndex} = scatter(FigAx{TimeIndex}, x,y,'MarkerFaceColor',cmap(end-t+1,:),'MarkerEdgeColor','k');
        end
        
        
        
        
        
        
        
        if ~useRescaledTiming & ~useRescaledParamTiming
            xlabel('Time (min)')
        else
            xlabel({'Rescaled Time', [' (min. at ',TemperatureString,' ºC)']})
        end
        
        ylabel(YLabel)
        set(FigAx{TimeIndex} ,'Fontsize',18)
        xlim([APbins(MinAPbin-1)/100 APbins(MaxAPbin+1)/100])
        xlabel('Position (x/L)')
        
        
        if NumTimeBins > 1
        TimeInterval = TimeVector(1)*2;
        HalfTimeInterval = TimeVector(1);
        if ~useRescaledTiming
            subtitle = strcat('$t =  ',sprintf('%.1f ', TimeVector(TimeBin)-HalfTimeInterval),'\,\textup{--}\,',sprintf('%.1f ', TimeVector(TimeBin)+HalfTimeInterval),'\, \textrm{m}$');
            title(subtitle, 'Interpreter', 'latex');
        else
            subtitle = strcat('$\bar{t} =  ',sprintf('%.1f ', TimeVector(TimeBin)-HalfTimeInterval),'\,\textup{--}\,',sprintf('%.1f ', TimeVector(TimeBin)+HalfTimeInterval),'\, \textrm{m}$');
            title(subtitle, 'Interpreter', 'latex');
        end
        end
        
        
        
        newylim = get(FigAx{TimeIndex}, 'ylim');
        newylim(1) = Ymin;
        newylim(2) = Ymax;
        set(FigAx{TimeIndex}, 'ylim', newylim)
        
        %     ylim([DurYmin, DurYmax])
        
        grid on
        box on 
        
        
    end
    
    set(FigHandle, 'CurrentAxes', LegendAx);
    hold on
    axis off
    
    pos = get(LegendAx, 'Position');
    pos(1) = SubplotPositionInfo.LegendXPosition;
    set(LegendAx, 'Position', pos);
    if max(IncludedTimeBins) > 1
        set(LegendAx, 'Clim', [1, NumTemperatures+1])
        h = colorbar('west');
        h.Ticks = 1.5:NumTemperatures+0.5;
        h.TickLabels = fliplr(CompiledParameters.UniqueTemperatures);
        
        ylabel(h,'Temperature (ºC)')
        
    end
    hold off
    set(LegendAx,'Fontsize',18)
    
    if strcmpi(hmmVarStrings{i}, 'initiationrates') | strcmpi(hmmVarStrings{i}, 'meaninitiationrates')
        OutFileName = [FluoString, RescaledString,ParamString, 'TimeSubplots_', OutString,'VsAP.png'];
    else
        OutFileName = [RescaledString,ParamString, 'TimeSubplots_', OutString,'VsAP.png'];
    end
    
    
    
    saveas(FigHandle,[outdir, OutFileName])
    
    
    
end
%%

close all
end
%%










