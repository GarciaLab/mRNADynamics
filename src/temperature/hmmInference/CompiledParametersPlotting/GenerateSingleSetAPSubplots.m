function GenerateSingleSetAPSubplots(ResultsPaths, varargin)
close all

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
    cmap = brewermap(max(IncludedTimeBins),'Spectral');
else
    cmap = brewermap(2,'Spectral');
end

SubplotPositionInfo = GetSubplotPositioningParameters(NumAPbins);
hmmVarStrings = { 'InitiationRates','Durations', 'Frequencies','BurstCycleTimes','MeanInitiationRates'};
%%

for SetIndex = 1:NumSets
    fig_path = CompiledParameters.FigurePaths{SetIndex};
    
    %%
    for i = 1:length(hmmVarStrings)
        close all
        [MeanValues, StdErrors, TimeVector, Ymax, Ymin, YLabel, OutString, LogYmax, LogYmin] = ...
            GetPlottingMats(CompiledParameters, hmmVarStrings{i}, useRescaledTiming,...
            useRescaledFluo, useRescaledParamTiming, SubplotPositionInfo.SubplotDims(1),ReferenceTemperature);
        if length(TimeVector) > 1
            dt = TimeVector(2)-TimeVector(1);
        end
        
        ValidTimePoints = TimeVector(find(sum((squeeze(sum(~isnan(MeanValues), 1))), 2).' > 0));
        FigAx = cell(1, NumAPbins);
        
        FigHandle = figure(1);
        
        set(FigHandle,'units', 'normalized', 'position',[0.05, 0.05, SubplotPositionInfo.SubFigDims(1), SubplotPositionInfo.SubFigDims(2)]);
        set(gcf,'color','w');
        
        colormap(cmap);
        hold on
        
        
        % Initiatialize all subplots before messing with any positioning.
        for APindex = 1:NumAPbins
            if APindex == 1
                FigAx{APindex} = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotPositionInfo.SubplotIndexList(APindex));
            else
                FigAx{APindex} = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotPositionInfo.SubplotIndexList(APindex));
            end
        end
        
        LegendAx = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotPositionInfo.LegendSubplots);
        
        
        
        
        
        
        for APindex = 1:NumAPbins
            APbin = APindex + MinAPbin -1;
            SubplotIndex = SubplotPositionInfo.SubplotIndexList(APindex);
            set(FigHandle, 'CurrentAxes', FigAx{APindex});
            
            
            
            y_vals = squeeze(MeanValues(:,:,APbin));
            y_val_errs = squeeze(StdErrors(:,:,APbin));
            for t = 1:length(TimeVector)
                if isnan(y_vals(SetIndex,t))
                    continue
                end
                
                
                x = TimeVector(t);
                
                y = y_vals(SetIndex, t);
                y_err = y_val_errs(SetIndex, t);
                
                errorbar(FigAx{APindex}, x,y,y_err,'Color','k','Capsize',0, 'LineStyle', 'None')
                hold on
                scatter(FigAx{APindex}, x,y,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
                
            end
            
            
            
            
            
            pos = get(FigAx{APindex}, 'position');
            pos(1) = SubplotPositionInfo.SubplotXPositions(APindex);
            pos(3) = SubplotPositionInfo.SubplotWidth;
            pos(2) = SubplotPositionInfo.SubplotYPositions(APindex);
            pos(4) = SubplotPositionInfo.SubplotHeight;
            set(FigAx{APindex}, 'position', pos);
            
            if ~useRescaledTiming & ~useRescaledParamTiming
                xlabel('Time (min)')
            else
                xlabel({'Rescaled Time', ['(min at ',TemperatureString,'ºC)']})
            end
            
            ylabel(YLabel)
            set(FigAx{APindex} ,'Fontsize',14)
            
            xlim([min(ValidTimePoints)-dt, max(ValidTimePoints)+dt])
            
            %xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
            title({['AP: ', num2str((APbin-1)*APResolution*100), '%']});
            
            
            
            newylim = get(FigAx{APindex}, 'ylim');
            newylim(1) = Ymin;
            newylim(2) = Ymax;
            set(FigAx{APindex}, 'ylim', newylim)
            
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
            if useRescaledTiming
                ylabel(h,['dev time cohort (',TemperatureString,'ºC minutes into nc14)'])
            else
                ylabel(h,'time cohort (minutes into nc14)')
            end
        end
        hold off
        set(LegendAx,'Fontsize',16)
        
        if strcmpi(hmmVarStrings{i}, 'initiationrates') | strcmpi(hmmVarStrings{i}, 'meaninitiationrates')
            OutFileName = [FluoString, RescaledString, ParamString, 'APSubplots_', OutString,'.png'];
        else
            OutFileName = [RescaledString, ParamString, 'APSubplots_', OutString,'.png'];
        end

        sgtitle([CompiledParameters.ReporterLabels{SetIndex}, ' ',...
            num2str(CompiledParameters.SetTemperatures(SetIndex)), '°C w',...
            num2str(CompiledParameters.nSteps(SetIndex)), ', dt = ',...
            num2str(CompiledParameters.dt(SetIndex)), 'seconds, Elongation time = ',...
            num2str(round(CompiledParameters.ElongationTimes(SetIndex), 2)), ' minutes'],...
            'FontWeight', 'bold', 'FontSize', 18);
        
        saveas(FigHandle,[fig_path, OutFileName])
        
        
        
        
        
    end
    %%
    
    close all
end
%%










