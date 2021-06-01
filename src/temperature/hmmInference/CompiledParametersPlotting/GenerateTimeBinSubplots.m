function GenerateTimeBinSubplots(ResultsPaths, outdir, varargin)
close all


IncludeFits = true;
useRescaledTiming = false;
useRescaledFluo = false;
useRescaledParamTiming = false;
includeRescaling = true;
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
    elseif strcmpi(varargin{x}, 'excluderescaling')
        includeRescaling = false;
    elseif strcmpi(varargin{x}, 'ExcludeFits')
        IncludeFits = false;
    elseif strcmpi(varargin{x}, 'prescaledbins')
        PrescaledBins = true;
    end
    x = x+1;
end

if useRescaledParamTiming
    IncludeFits = false;
end




if ~exist(outdir, 'dir')
    mkdir(outdir)
end

R = 8.314*10^(-3); % kJ mol^-1 K^-1



APResolution = 0.025;
APbins = (0:APResolution:1)*100;

CompiledParameters = CompileParameters(ResultsPaths, includeRescaling, PrescaledBins);
ReferenceTemperature = CompiledParameters.RefTemperature;
if single(round(ReferenceTemperature, 0)) == single(ReferenceTemperature)
    TemperatureString = sprintf('%.0f ', ReferenceTemperature);
else
    TemperatureString = sprintf('%.1f ', ReferenceTemperature);
end
NumSets = size(CompiledParameters.InitiationRates, 1);
HasSingleTimeBin = all(CompiledParameters.nTimebins == 1);

if ~HasSingleTimeBin
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
if IncludeFits
    FitString = 'IncludedFits';
else
    FitString = '';
end
MinAPbin = min(IncludedAPbins);
MaxAPbin = max(IncludedAPbins);
NumAPbins = MaxAPbin-MinAPbin + 1;

cmap = brewermap(NumAPbins,'Spectral');


SubplotPositionInfo = GetSubplotPositioningParameters(NumTimeBins, true, false, false);
hmmVarStrings = { 'InitiationRates','Durations', 'Frequencies','BurstCycleTimes','MeanInitiationRates',...
    'InitiationRates','Durations', 'Frequencies','BurstCycleTimes','MeanInitiationRates'};
useLogVector = [false, false, false, false,false, true, true, true, true,true];


%%
for i = 1:length(hmmVarStrings)
    if useLogVector(i)
        LogString = 'Log';
    else
        LogString = '';
    end
    close all
    [MeanValues, StdErrors, TimeVector, Ymax, Ymin, YLabel, OutString, LogYmax, LogYmin] = ...
        GetPlottingMats(CompiledParameters, hmmVarStrings{i}, useRescaledTiming,...
        useRescaledFluo,useRescaledParamTiming, SubplotPositionInfo.SubplotDims(1),ReferenceTemperature);
    
    if IncludeFits & ~(strcmpi(hmmVarStrings{i}, 'BurstCycleTimes') | strcmpi(hmmVarStrings{i}, 'MeanInitiationRates'))
        [ActivationEnergies, SEActivationEnergies, LogAs, SELogAs, R2s, PlotTitle] =...
            getHMMActivationEnergyMatrices(CompiledParameters, hmmVarStrings{i}, useRescaledFluo, useRescaledTiming);
    end
    FigAx = cell(1, NumTimeBins);
    
    FigHandle = figure(1);
    
    set(FigHandle,'units', 'normalized', 'position',[0.05, 0.05, SubplotPositionInfo.SubFigDims(1), SubplotPositionInfo.SubFigDims(2)]);
    set(gcf,'color','w');
    
    colormap(cmap);
    hold on
    
    
    % Initiatialize all subplots before messing with any positioning.
    for TimeIndex = 1:NumTimeBins
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
        hold(FigAx{TimeIndex}, 'on');
        
        if IncludeFits & ~(strcmpi(hmmVarStrings{i}, 'BurstCycleTimes') | strcmpi(hmmVarStrings{i}, 'MeanInitiationRates'))
            Ea_vals = squeeze(ActivationEnergies(TimeBin,:));
            Ea_val_errs = squeeze(SEActivationEnergies(TimeBin,:));
            LogA_vals = squeeze(LogAs(TimeBin,:));
            LogA_val_errs = squeeze(SELogAs(TimeBin,:));
            
            
            for APindex = 1:length(IncludedAPbins)
                APbin = IncludedAPbins(APindex);
                if isnan(Ea_vals(APbin))
                    continue
                end
                if useLogVector(i)
                    FitXLims = [0.3975 .415];
                    x = linspace(FitXLims(1), FitXLims(2), 50);
                    log_fity = LogA_vals(APbin)-(Ea_vals(APbin))*x;
                    fity = exp(log_fity);
                else
                    x = linspace(16, 29, 50);
                    log_fity = LogA_vals(APbin)-(Ea_vals(APbin)/R)*(1./(x+273.15));
                    fity = exp(log_fity);

                end
                
                plot( FigAx{TimeIndex},x, fity,'-', 'Color', cmap(APindex,:),'LineWidth', 1);
                
            end
        end
        
        
        y_vals = squeeze(MeanValues(:,TimeBin,:));
        y_val_errs = squeeze(StdErrors(:,TimeBin,:));
        
        for APindex = 1:length(IncludedAPbins)
            APbin = IncludedAPbins(APindex);

            if all(isnan(y_vals(:,APbin)))
                continue
            end
            set_ids = ~isnan(y_vals(:,APbin));
            
            x = CompiledParameters.SetTemperatures(set_ids);
            if useLogVector(i)
                x = 1./(R*(x+273.15));
            end
            y = y_vals(set_ids, APbin);
            y_err = y_val_errs(set_ids, APbin);
            

            errorbar(FigAx{TimeIndex}, x,y,y_err,'Color','k','Capsize',0, 'LineStyle', 'None')
            hold on
            scatter(FigAx{TimeIndex}, x,y,'MarkerFaceColor',cmap(APindex,:),'MarkerEdgeColor','k')
            
        end
        
        
        
        
        
        
        
        ylabel(YLabel)
        set(FigAx{TimeIndex} ,'Fontsize',14)
        if ~useLogVector(i)
            xlim([16, 29])
            xlabel('Temperature (ºC)')
            xticks([17.5 20 22.5 25 27.5])
            if SubplotPositionInfo.SubplotDims(2) > 5
                xtickangle(45)
            end
        else
            xlim([.3975 .415])
            xlabel('1/RT (mol \cdot kJ^{-1})')
        end
        %xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
        TimeInterval = TimeVector(1)*2;
        HalfTimeInterval = TimeVector(1);
        if NumTimeBins > 1
            if ~useRescaledTiming
                subtitle = strcat('$t =  ',sprintf('%.1f ', TimeVector(TimeBin)-HalfTimeInterval),'\,\textup{--}\,',sprintf('%.1f ', TimeVector(TimeBin)+HalfTimeInterval),'\, \textrm{m}$');
                title(subtitle, 'Interpreter', 'latex', 'FontSize', 14);
            else
                subtitle = strcat('$\bar{t} =  ',sprintf('%.1f ', TimeVector(TimeBin)-HalfTimeInterval),'\,\textup{--}\,',sprintf('%.1f ', TimeVector(TimeBin)+HalfTimeInterval),'\, \textrm{m}$');
                title(subtitle, 'Interpreter', 'latex', 'FontSize', 14);
            end
        end
        
        
        if useLogVector(i)
            
            newylim = get(FigAx{TimeIndex}, 'ylim');
            newylim(1) = LogYmin;
            newylim(2) = LogYmax;
            set(FigAx{TimeIndex}, 'ylim', newylim)
            set(FigAx{TimeIndex}, 'YScale', 'log');
        else
            newylim = get(FigAx{TimeIndex}, 'ylim');
            newylim(1) = Ymin;
            newylim(2) = Ymax;
            set(FigAx{TimeIndex}, 'ylim', newylim)
        end
        %     ylim([DurYmin, DurYmax])
        
        grid on
        hold off
        box on 
    end
    
    set(FigHandle, 'CurrentAxes', LegendAx);
    hold on
    axis off
    
    pos = get(LegendAx, 'Position');
    pos(1) = SubplotPositionInfo.LegendXPosition;
    set(LegendAx, 'Position', pos);
    set(LegendAx, 'Clim', [1,NumAPbins+1])
    h = colorbar('west');
    h.Ticks = 1.5:length(IncludedAPbins)+0.5;
    h.TickLabels = round(APbins(IncludedAPbins)/100, 3);
 
    ylabel(h,'Fraction of Embryo Length (x/L)')
    

    hold off
    set(LegendAx,'Fontsize',16)
    
    if strcmpi(hmmVarStrings{i}, 'initiationrates') | strcmpi(hmmVarStrings{i}, 'meaninitiationrates')
        OutFileName = [LogString, FluoString, RescaledString,FitString,ParamString, 'TimeBinSubplots_', OutString,'VsTemp.png'];
    else
        OutFileName = [LogString, RescaledString,FitString,ParamString, 'TimeBinSubplots_', OutString,'VsTemp.png'];
    end
    
    saveas(FigHandle,[outdir, OutFileName])
    
    
    
end
%%

close all

%%










