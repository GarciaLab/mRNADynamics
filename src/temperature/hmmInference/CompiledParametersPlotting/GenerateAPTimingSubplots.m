function GenerateAPTimingSubplots(ResultsPaths, outdir, varargin)
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
    cmap = flipud(brewermap(max(NumTemperatures),'Spectral'));
else
    cmap = flipud(brewermap(2,'Spectral'));
end
cmap = brewermap(10,'Spectral');
cmap = [cmap(2:3,:);cmap(8:10,:)];
SubplotPositionInfo = GetSubplotPositioningParameters(NumAPbins);
hmmVarStrings = { 'InitiationRates','Durations', 'Frequencies','BurstCycleTimes','MeanInitiationRates'};
%%
for i = 1:length(hmmVarStrings)
    close all

        
        
        %%
        
        [MeanValues, StdErrors, TimeVector, Ymax, Ymin, YLabel, OutString, LogYmax, LogYmin] = ...
            GetPlottingMats(CompiledParameters, hmmVarStrings{i}, useRescaledTiming, useRescaledFluo,...
            useRescaledParamTiming, SubplotPositionInfo.SubplotDims(1),ReferenceTemperature);
        dt = TimeVector(2)-TimeVector(1);
        
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
            for SetIndex = 1:NumSets
                t = find(round(CompiledParameters.UniqueTemperatures, 1) == round(CompiledParameters.SetTemperatures(SetIndex), 1));
                if all(isnan(y_vals(SetIndex,:)))
                    continue
                end
                
                
                
                
                y = y_vals(SetIndex, :);
                y_err = y_val_errs(SetIndex, :);
                
                x = TimeVector(~isnan(y));
                y_err = y_err(~isnan(y));
                y = y(~isnan(y));
                
                errorbar(FigAx{APindex}, x,y,y_err,'Color','k','Capsize',0)
                hold on
                
                scatter(FigAx{APindex}, x,y,'MarkerFaceColor',cmap(end-t+1,:),'MarkerEdgeColor','k')
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
                xlabel({'Rescaled Time', [' (min. at ',TemperatureString,' ºC)']})
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
            set(LegendAx, 'Clim', [1, NumTemperatures+1])
            h = colorbar('west');
            h.Ticks = 1.5:NumTemperatures+0.5;
            h.TickLabels = fliplr(CompiledParameters.UniqueTemperatures);
        
            ylabel(h,'Temperature (ºC)')
           
        end
        hold off
        set(LegendAx,'Fontsize',16)
        
        if strcmpi(hmmVarStrings{i}, 'initiationrates') | strcmpi(hmmVarStrings{i}, 'meaninitiationrates') 
            OutFileName = [FluoString, RescaledString,ParamString, 'APTimingSubplots_', OutString,'VsTime.png'];
        else
            OutFileName = [RescaledString,ParamString, 'APTimingSubplots_', OutString,'VsTime.png'];
        end
        

        
        saveas(FigHandle,[outdir, OutFileName])
        
        
        
    end
    %%
    
    close all
end
%%










