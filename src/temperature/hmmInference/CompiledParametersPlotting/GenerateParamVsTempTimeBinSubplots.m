function GenerateParamVsTempTimeBinSubplots(ResultsPaths, outdir, varargin)
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

CompiledParameters = CompileParameters(ResultsPaths, true, PrescaledBins);
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

if max(IncludedAPbins) > 1
    cmap = brewermap(length(IncludedAPbins),'Spectral');
else
    cmap = brewermap(2,'Spectral');
end


hmmVarStrings = { 'InitiationRates','Durations', 'Frequencies','BurstCycleTimes','MeanInitiationRates',...
    'InitiationRates','Durations', 'Frequencies','BurstCycleTimes','MeanInitiationRates'};
useLogVector = [false, false, false, false,false, true, true, true, true,true];



%%
for i = 1:length(hmmVarStrings)
    SubplotPositionInfo = GetSubplotPositioningParameters(NumTimeBins, false);
    if useLogVector(i)
        LogString = 'Log';
    else
        LogString = '';
    end
    close all
    [MeanValues, StdErrors, TimeVector, Ymax, Ymin, YLabel, OutString, LogYmax, LogYmin] = ...
        GetPlottingMats(CompiledParameters, hmmVarStrings{i}, useRescaledTiming,...
        useRescaledFluo,useRescaledParamTiming, SubplotPositionInfo.SubplotDims(1), ReferenceTemperature);
    
    

    
    if IncludeFits & ~(strcmpi(hmmVarStrings{i}, 'BurstCycleTimes') | strcmpi(hmmVarStrings{i}, 'MeanInitiationRates'))
        [ActivationEnergies, SEActivationEnergies, LogAs, SELogAs, R2s, PlotTitle] =...
            getHMMActivationEnergyMatrices(CompiledParameters, hmmVarStrings{i}, useRescaledFluo, useRescaledTiming);
    end
    
    if strcmpi(hmmVarStrings{i}, 'initiationrates')
        paramYinterval =50;
    elseif strcmpi(hmmVarStrings{i}, 'durations')
        paramYinterval =2;
    elseif strcmpi(hmmVarStrings{i}, 'Frequencies')
        paramYinterval =0.2;
    elseif strcmpi(hmmVarStrings{i}, 'BurstCycleTimes')
        paramYinterval =2;
    end
    
    % Initiatialize all subplots before messing with any positioning.
    for APindex = 1:NumAPbins
        APbin = APindex + MinAPbin -1;
        
        if ~HasSingleTimeBin
            if ~useRescaledTiming & ~useRescaledParamTiming
                
                InitiationRateMat = squeeze(CompiledParameters.InitiationRates(:,:,APbin));
                IncludedSubplotTimeBins = find(sum(~isnan(InitiationRateMat), 1) > 0);
                NumSubplotTimeBins  =  length(IncludedSubplotTimeBins);
                TimeLabels = CompiledParameters.TimeVector;
            else
                InitiationRateMat = squeeze(CompiledParameters.ScaledInitiationRates(:,:,APbin));
                IncludedSubplotTimeBins = find(sum(~isnan(InitiationRateMat), 1) > 0);
                NumSubplotTimeBins  =  length(IncludedSubplotTimeBins);
                TimeLabels = CompiledParameters.ScaledTimeVector;
            end
        end
        close all
        APbin_y_vals = squeeze(MeanValues(:,:,APbin));
        APbin_y_val_errs = squeeze(StdErrors(:,:,APbin));
        MaxAPbinValue = max(APbin_y_vals(~isnan(APbin_y_vals)));
        MinAPbinValue = min(APbin_y_vals(~isnan(APbin_y_vals)));
        
        
        MaxAPbinValue = max([MaxAPbinValue*1.1, MaxAPbinValue*0.9]);
        MinAPbinValue = min([MinAPbinValue*1.1, MinAPbinValue*0.9]);
        
        APbinYmax = paramYinterval*ceil(MaxAPbinValue/paramYinterval);
        APbinYmin = paramYinterval*floor(MinAPbinValue/paramYinterval);
       
        if strcmpi(hmmVarStrings{i}, 'initiationrates')
            LogAPbinYmax = 10^(ceil(log10(APbinYmax)));
            LogAPbinYmin =APbinYmin;
        else
            LogAPbinYmin =APbinYmin;
            LogAPbinYmax =APbinYmax;
        end
        
        SubplotPositionInfo = GetSubplotPositioningParameters(NumSubplotTimeBins, false);
        FigAx = cell(1, NumSubplotTimeBins);
        FitPlots = cell(1, NumSubplotTimeBins);
        DataPlots = cell(1, NumSubplotTimeBins);
        ErrorPlots = cell(1, NumSubplotTimeBins);
        LegendAxes = cell(1, NumSubplotTimeBins);
        FigHandle = figure(1);
        
        set(FigHandle,'units', 'normalized', 'position',[0.05, 0.05, SubplotPositionInfo.SubFigDims(1), SubplotPositionInfo.SubFigDims(2)]);
        set(gcf,'color','w');
        
        colormap(cmap);
        hold on
        for TimeIndex=1:NumSubplotTimeBins
            if TimeIndex == 1
                FigAx{TimeIndex} = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotPositionInfo.SubplotIndexList(TimeIndex));
            else
                FigAx{TimeIndex} = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotPositionInfo.SubplotIndexList(TimeIndex));
            end
        end
        
        for TimeIndex=1:NumSubplotTimeBins
            pos = get(FigAx{TimeIndex}, 'position');
            pos(1) = SubplotPositionInfo.SubplotXPositions(TimeIndex);
            pos(3) = SubplotPositionInfo.SubplotWidth;
            pos(2) = SubplotPositionInfo.SubplotYPositions(TimeIndex);
            pos(4) = SubplotPositionInfo.SubplotHeight;
            set(FigAx{TimeIndex}, 'position', pos);
        end
       
        for TimeIndex=1:NumSubplotTimeBins
            TimeBin = IncludedSubplotTimeBins(TimeIndex);
            
            SubplotIndex = SubplotPositionInfo.SubplotIndexList(TimeIndex);
            set(FigHandle, 'CurrentAxes', FigAx{TimeIndex});
            hold(FigAx{TimeIndex}, 'on')
            hold on
            
            if IncludeFits & ~(strcmpi(hmmVarStrings{i}, 'BurstCycleTimes') | strcmpi(hmmVarStrings{i}, 'MeanInitiationRates'))
                Ea_val = ActivationEnergies(TimeBin,APbin);
                Ea_val_err = SEActivationEnergies(TimeBin,APbin);
                LogA_val = LogAs(TimeBin,APbin);
                LogA_val_errs = SELogAs(TimeBin,APbin);
                
                if ~isnan(Ea_val)
                    if useLogVector(i)
                        FitXLims = [0.3975 .415];
                        x = linspace(FitXLims(1), FitXLims(2), 50);
                        log_fity = LogA_val-(Ea_val)*x;
                        fity = exp(log_fity);
                    else
                        x = linspace(16, 29, 50);
                        log_fity = LogA_val-(Ea_val/R)*(1./(x+273.15));
                        fity = exp(log_fity);  
                    end
                    
                    FitPlots{1, TimeIndex} = plot( FigAx{TimeIndex},x, fity,'-', 'Color', cmap(APindex,:),'LineWidth', 1);
                    %a = ['$\it{\bf{(f_{ohm}P)^{',rats(m), '}}}$'];
                    legend_label = strcat('$E_A = ', sprintf('%.0f ', Ea_val), ' \pm ', sprintf('%.0f ', Ea_val_err), '\,\,\frac{\textrm{kJ}}{\textrm{mol}}$');
                    
                else
                    legend_label = '';
                end
                
            end
            
            
            y_vals = squeeze(MeanValues(:,TimeBin,APbin)).';
            y_val_errs = squeeze(StdErrors(:,TimeBin,APbin)).';
            
            
            if ~all(isnan(y_vals))
                
                set_ids = ~isnan(y_vals);
                
                x = CompiledParameters.SetTemperatures(set_ids);
                if useLogVector(i)
                    x = 1./(R*(x+273.15));
                end
                y = y_vals(set_ids);
                y_err = y_val_errs(set_ids);
                
                ErrorPlots{TimeIndex} = errorbar(FigAx{TimeIndex}, x,y,y_err,'Color','k','Capsize',0, 'LineStyle', 'None');
                hold on
                DataPlots{TimeIndex} = scatter(FigAx{TimeIndex}, x,y,'MarkerFaceColor',cmap(APindex,:),'MarkerEdgeColor','k');
                set(get(get(ErrorPlots{TimeIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                set(get(get(DataPlots{TimeIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    
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
            title({[num2str(round(TimeLabels(TimeBin), 1)), ' min.']});
            
            if useLogVector(i)
                newylim = get(FigAx{TimeIndex}, 'ylim');
                newylim(1) = LogAPbinYmin;
                newylim(2) = LogAPbinYmax;
                set(FigAx{TimeIndex}, 'ylim', newylim)
                set(FigAx{TimeIndex}, 'YScale', 'log');
            else
                newylim = get(FigAx{TimeIndex}, 'ylim');
                newylim(1) = APbinYmin;
                newylim(2) = APbinYmax;
                set(FigAx{TimeIndex}, 'ylim', newylim)
            end
            grid on
            box on
            hold off
            
            if IncludeFits & ~(strcmpi(hmmVarStrings{i}, 'BurstCycleTimes') | strcmpi(hmmVarStrings{i}, 'MeanInitiationRates'))
                if (~useLogVector(i) & ~strcmpi(hmmVarStrings{i}, 'durations')) | (useLogVector(i) & strcmpi(hmmVarStrings{i}, 'durations'))
                    LegendAxes{TimeIndex} = legend(legend_label, 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'northwest');
                else
                    LegendAxes{TimeIndex} = legend(legend_label, 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'northeast');
                end
            end
        end
        if strcmpi(hmmVarStrings{i}, 'initiationrates') | strcmpi(hmmVarStrings{i}, 'meaninitiationrates') 
            OutFileName = [LogString, FluoString, RescaledString,FitString,ParamString, 'TimeSubplotsAP',num2str(APbin),'_', OutString,'VsTemp.png'];
        else
            OutFileName = [LogString, RescaledString,FitString,ParamString, 'TimeSubplotsAP',num2str(APbin),'_', OutString,'VsTemp.png'];
        end
        
        sgtitle(['AP: ', num2str(round(APbins(APbin)/100, 2))], 'FontSize', 18, 'FontWeight', 'bold')
        
        
        
        saveas(FigHandle,[outdir, OutFileName])
        
        
    end
    
    
    
end



%%

close all

%%










