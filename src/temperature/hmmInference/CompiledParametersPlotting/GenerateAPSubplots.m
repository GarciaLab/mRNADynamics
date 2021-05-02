function GenerateAPSubplots(ResultsPaths, outdir, varargin)
close all


IncludeFits = true;
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
    elseif strcmpi(varargin{x}, 'ExcludeFits')
        IncludeFits = false;
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
if IncludeFits
    FitString = 'IncludedFits';
else
    FitString = '';
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
        GetPlottingMats(CompiledParameters, hmmVarStrings{i}, useRescaledTiming,...
        useRescaledFluo,useRescaledParamTiming, SubplotPositionInfo.SubplotDims(1));
    
    if IncludeFits & ~strcmpi(hmmVarStrings{i}, 'BurstCycleTimes')
        [ActivationEnergies, SEActivationEnergies, LogAs, SELogAs, R2s, PlotTitle] =...
            getHMMActivationEnergyMatrices(CompiledParameters, hmmVarStrings{i}, useRescaledFluo, useRescaledTiming);
    end
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
        hold on
        PlottedTimes = zeros(1, length(TimeVector), 'logical');
        
        if IncludeFits & ~strcmpi(hmmVarStrings{i}, 'BurstCycleTimes')
            Ea_vals = squeeze(ActivationEnergies(:,APbin)).';
            Ea_val_errs = squeeze(SEActivationEnergies(:,APbin)).';
            LogA_vals = squeeze(LogAs(:,APbin)).';
            LogA_val_errs = squeeze(SELogAs(:,APbin)).';
            
            
            for t = 1:length(IncludedTimeBins)
                TimeIndex = IncludedTimeBins(t);
                if isnan(Ea_vals(TimeIndex))
                    continue
                end
                if useLogVector(i)
                    FitXLims = [0.3975 .415];
                    x = linspace(FitXLims(1), FitXLims(2), 50);
                    log_fity = LogA_vals(TimeIndex)-(Ea_vals(TimeIndex))*x;
                    fity = exp(log_fity);
                else
                    x = linspace(16, 29, 50);
                    log_fity = LogA_vals(TimeIndex)-(Ea_vals(TimeIndex)/R)*(1./(x+273.15));
                    fity = exp(log_fity);

                end
                
                plot( FigAx{APindex},x, fity,'-', 'Color', cmap(TimeIndex,:),'LineWidth', 1);
                
            end
        end
        
        
        y_vals = squeeze(MeanValues(:,:,APbin));
        y_val_errs = squeeze(StdErrors(:,:,APbin));
        if ~all(size(y_vals) > 1)
            y_vals = y_vals.';
            y_val_errs = y_val_errs.';
        end
        for t = 1:length(IncludedTimeBins)
            TimeIndex = IncludedTimeBins(t);
            if all(size(y_vals) > 1)
                if all(isnan(y_vals(:,TimeIndex)))
                    continue
                end
                set_ids = ~isnan(y_vals(:,TimeIndex));
                
                x = CompiledParameters.SetTemperatures(set_ids);
                if useLogVector(i)
                    x = 1./(R*(x+273.15));
                end
                y = y_vals(set_ids, TimeIndex);
                y_err = y_val_errs(set_ids, TimeIndex);
            else
                if all(isnan(y_vals))
                    continue
                end
                set_ids = ~isnan(y_vals);
                
                x = CompiledParameters.SetTemperatures(set_ids);
                if useLogVector(i)
                    x = 1./(R*(x+273.15));
                end
                y = y_vals(set_ids);
                y_err = y_val_errs(set_ids);
                
            end
            errorbar(FigAx{APindex}, x,y,y_err,'Color','k','Capsize',0, 'LineStyle', 'None')
            hold on
            scatter(FigAx{APindex}, x,y,'MarkerFaceColor',cmap(TimeIndex,:),'MarkerEdgeColor','k')
            
        end
        
        
        
        
        
        pos = get(FigAx{APindex}, 'position');
        pos(1) = SubplotPositionInfo.SubplotXPositions(APindex);
        pos(3) = SubplotPositionInfo.SubplotWidth;
        pos(2) = SubplotPositionInfo.SubplotYPositions(APindex);
        pos(4) = SubplotPositionInfo.SubplotHeight;
        set(FigAx{APindex}, 'position', pos);
        
        
        ylabel(YLabel)
        set(FigAx{APindex} ,'Fontsize',14)
        if ~useLogVector(i)
            xlim([16, 29])
            xlabel('Temperature (ºC)')
        else
            xlim([.3975 .415])
            xlabel('1/RT (mol \cdot kJ^{-1})')
        end
        %xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
        title({['AP: ', num2str((APbin-1)*APResolution*100), '%']});
        
        
        if useLogVector(i)
            
            newylim = get(FigAx{APindex}, 'ylim');
            newylim(1) = LogYmin;
            newylim(2) = LogYmax;
            set(FigAx{APindex}, 'ylim', newylim)
            set(FigAx{APindex}, 'YScale', 'log');
        else
            newylim = get(FigAx{APindex}, 'ylim');
            newylim(1) = Ymin;
            newylim(2) = Ymax;
            set(FigAx{APindex}, 'ylim', newylim)
        end
        %     ylim([DurYmin, DurYmax])
        
        grid on
        hold off
        
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
        OutFileName = [LogString, FluoString, RescaledString,FitString,ParamString, 'APSubplots_', OutString,'.png'];
    else
        OutFileName = [LogString, RescaledString,FitString,ParamString, 'APSubplots_', OutString,'.png'];
    end
    
    saveas(FigHandle,[outdir, OutFileName])
    
    
    
end
%%

close all

%%










