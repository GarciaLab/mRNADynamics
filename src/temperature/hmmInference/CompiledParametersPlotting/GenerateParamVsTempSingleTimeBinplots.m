function GenerateParamVsTempSingleTimeBinplots(ResultsPaths, outdir, varargin)
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
        useRescaledFluo,useRescaledParamTiming, 1);
    
    if IncludeFits & ~strcmpi(hmmVarStrings{i}, 'BurstCycleTimes')
        [ActivationEnergies, SEActivationEnergies, LogAs, SELogAs, R2s, PlotTitle] =...
            getHMMActivationEnergyMatrices(CompiledParameters, hmmVarStrings{i}, useRescaledFluo, useRescaledTiming);
    end
    
    
    
    
    
    % Initiatialize all subplots before messing with any positioning.

    
    if strcmpi(hmmVarStrings{i}, 'initiationrates')
        paramYinterval =50;
    elseif strcmpi(hmmVarStrings{i}, 'durations')
        paramYinterval =2;
    elseif strcmpi(hmmVarStrings{i}, 'Frequencies')
        paramYinterval =0.2;
    elseif strcmpi(hmmVarStrings{i}, 'BurstCycleTimes')
        paramYinterval =2;
    end
    
    
    
    for TimeIndex = 1:NumTimeBins
        TimeBin = IncludedTimeBins(TimeIndex);
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
        close all
        FigHandle = figure(1);
        FitPlots = cell(1, NumAPbins);
        DataPlots = cell(1, NumAPbins);
        ErrorPlots = cell(1, NumAPbins);
        set(FigHandle,'units', 'normalized', 'position',[0.05, 0.05, 0.9, 0.8]);
        set(gcf,'color','w');
        FigAx = axes(FigHandle);
        colormap(cmap);
        hold on
       
        legend_labels = {};
        if IncludeFits & ~strcmpi(hmmVarStrings{i}, 'BurstCycleTimes')
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
                    x = linspace(FitXLims(1), FitXLims(2), 1000);
                    log_fity = LogA_vals(APbin)-(Ea_vals(APbin))*x;
                    fity = exp(log_fity);
                else
                    x = linspace(16, 29, 1000);
                    log_fity = LogA_vals(APbin)-(Ea_vals(APbin)/R)*(1./(x+273.15));
                    fity = exp(log_fity);

                end
                
                if ~strcmpi(hmmVarStrings{i}, 'durations')
                    legend_labels{end+1} = strcat('$AP = ',sprintf('%.1f ', APbins(APbin)),'\%,\,\,E_A = ', sprintf('%.1f ', Ea_vals(APbin)), ' \pm ', sprintf('%.1f ', Ea_val_errs(APbin)), '\,\,\textrm{kJ}\cdot\textrm{mol}^{-1}$');
                else
                     legend_labels{end+1} = strcat('$AP = ',sprintf('%.1f ', APbins(APbin)),'\%,\,\,E_A = ', sprintf('%.1f ', -1*Ea_vals(APbin)), ' \pm ', sprintf('%.1f ', Ea_val_errs(APbin)), '\,\,\textrm{kJ}\cdot\textrm{mol}^{-1}$');
                end
                
                FitPlots{APindex} = plot( FigAx,x, fity,'-', 'Color', cmap(APindex,:),'LineWidth', 2,...
                    'DisplayName', legend_labels{end});   
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
            
       
            ErrorPlots{APindex} = errorbar(FigAx, x,y,y_err,'Color','k','Capsize',0, 'LineStyle', 'None');
            hold on
            DataPlots{APindex} = scatter(FigAx, x,y,80, 'MarkerFaceColor',cmap(APindex,:),'MarkerEdgeColor','k');
            set(get(get(ErrorPlots{APindex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            set(get(get(DataPlots{APindex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            
        end
        
        
        

        
        ylabel(YLabel)
        set(FigAx ,'Fontsize',18)
        if ~useLogVector(i)
            xlim([16, 29])
            xlabel('Temperature (ºC)')
            xticks([17.5 20 22.5 25 27.5])
        else
            xlim([.3975 .415])
            xlabel('1/RT (mol \cdot kJ^{-1})')
        end
        %xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
        TimeInterval = TimeVector(1)*2;
        HalfTimeInterval = TimeVector(1);
        if ~useRescaledTiming
            subtitle = strcat('$t =  ',sprintf('%.1f ', TimeVector(TimeBin)-HalfTimeInterval),'\,\textup{--}\,',sprintf('%.1f ', TimeVector(TimeBin)+HalfTimeInterval),'\, \textrm{m}$');
            title(subtitle, 'Interpreter', 'latex');
        else
            subtitle = strcat('$\bar{t} =  ',sprintf('%.1f ', TimeVector(TimeBin)-HalfTimeInterval),'\,\textup{--}\,',sprintf('%.1f ', TimeVector(TimeBin)+HalfTimeInterval),'\, \textrm{m}$');
            title(subtitle, 'Interpreter', 'latex');
        end
        
        
        if useLogVector(i)
            
            newylim = get(FigAx, 'ylim');
            newylim(1) = LogTimeBinYmin;
            newylim(2) = LogTimeBinYmax;
            set(FigAx, 'ylim', newylim)
            set(FigAx, 'YScale', 'log');
        else
            newylim = get(FigAx, 'ylim');
            newylim(1) = TimeBinYmin;
            newylim(2) = TimeBinYmax;
            set(FigAx, 'ylim', newylim)
        end
        %     ylim([DurYmin, DurYmax])
        
        grid on
        box on
        hold off
        
        if IncludeFits & ~strcmpi(hmmVarStrings{i}, 'BurstCycleTimes') & ~isempty(legend_labels)
            if (~useLogVector(i) & ~strcmpi(hmmVarStrings{i}, 'durations')) | (useLogVector(i) & strcmpi(hmmVarStrings{i}, 'durations'))
              
                LegendAxes{TimeIndex} = legend('show', 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'northwest');
            else
                LegendAxes{TimeIndex} = legend('show', 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'northeast');
            end
        end

        
        
        set(FigAx, 'Clim', [1,NumAPbins+1]);
        h = colorbar();
        h.Ticks = 1.5:length(IncludedAPbins)+0.5;
        h.TickLabels = round(APbins(IncludedAPbins)/100, 3);
        
        ylabel(h,'Fraction of Embryo Length (x/L)')
        
        hold off
        if strcmpi(hmmVarStrings{i}, 'initiationrates')
            OutFileName = [LogString, FluoString, RescaledString,FitString,ParamString, 'TimeBin',num2str(TimeBin),'_', OutString,'.png'];
        else
            OutFileName = [LogString, RescaledString,FitString,ParamString, 'TimeBin',num2str(TimeBin),'_', OutString,'.png'];
        end
        
        saveas(FigHandle,[outdir, OutFileName])
        
        
    end


    
    
    
end
%%

close all

%%










