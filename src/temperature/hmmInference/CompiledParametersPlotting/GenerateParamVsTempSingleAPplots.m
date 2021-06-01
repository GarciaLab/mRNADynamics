function GenerateParamVsTempSingleAPplots(ResultsPaths, outdir, varargin)
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


cmap = brewermap(max(IncludedTimeBins),'Spectral');

hmmVarStrings = { 'InitiationRates','Durations', 'Frequencies','BurstCycleTimes','MeanInitiationRates',...
    'InitiationRates','Durations', 'Frequencies','BurstCycleTimes','MeanInitiationRates'};
useLogVector = [false, false, false, false, false,true, true, true, true,true];


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
        useRescaledFluo,useRescaledParamTiming, 1,ReferenceTemperature);
    
    if IncludeFits & ~(strcmpi(hmmVarStrings{i}, 'BurstCycleTimes') | strcmpi(hmmVarStrings{i}, 'MeanInitiationRates'))
        [ActivationEnergies, SEActivationEnergies, LogAs, SELogAs, R2s, PlotTitle] =...
            getHMMActivationEnergyMatrices(CompiledParameters, hmmVarStrings{i}, useRescaledFluo, useRescaledTiming);
    end
    
    
    
    
    
    % Initiatialize all subplots before messing with any positioning.

    
    if strcmpi(hmmVarStrings{i}, 'initiationrates') | strcmpi(hmmVarStrings{i}, 'meaninitiationrates')
        paramYinterval =50;
    elseif strcmpi(hmmVarStrings{i}, 'durations')
        paramYinterval =2;
    elseif strcmpi(hmmVarStrings{i}, 'Frequencies')
        paramYinterval =0.2;
    elseif strcmpi(hmmVarStrings{i}, 'BurstCycleTimes')
        paramYinterval =2;
    end
    
    
    
    for APindex = 1:NumAPbins
        APbin = APindex + MinAPbin -1;
        APbin_y_vals = squeeze(MeanValues(:,:,APbin));
        APbin_y_val_errs = squeeze(StdErrors(:,:,APbin));
        MaxAPbinValue = max(APbin_y_vals(~isnan(APbin_y_vals)));
        MinAPbinValue = min(APbin_y_vals(~isnan(APbin_y_vals)));
        
        if NumTimeBins > 1
            MaxAPbinValue = max([MaxAPbinValue*1.3, MaxAPbinValue*0.7]);
        else
            MaxAPbinValue = max([MaxAPbinValue*1.1, MaxAPbinValue*0.9]);
        end
        MinAPbinValue = min([MinAPbinValue*1.1, MinAPbinValue*0.9]);
        
        APbinYmax = paramYinterval*ceil(MaxAPbinValue/paramYinterval);
        APbinYmin = paramYinterval*floor(MinAPbinValue/paramYinterval);
        
        Ymax = paramYinterval*ceil(Ymax/paramYinterval);
        Ymin = paramYinterval*floor(Ymin/paramYinterval);
        LogYmin = Ymin;
        LogYmax = Ymax;
       
        if strcmpi(hmmVarStrings{i}, 'initiationrates') | strcmpi(hmmVarStrings{i}, 'meaninitiationrates') 
            LogAPbinYmax = 600;
            LogYmax = 700;
            LogYmin = 200;
       
            LogAPbinYmin =APbinYmin;
        else
            LogAPbinYmin =APbinYmin;
            LogAPbinYmax =APbinYmax;
        end
        close all
        FigHandle = figure(1);
        FitPlots = cell(1, NumTimeBins);
        DataPlots = cell(NumSets, NumTimeBins);
        ErrorPlots = cell(NumSets, NumTimeBins);
        set(FigHandle,'units', 'normalized', 'position',[0.05, 0.05, 0.35, 0.4]);
        set(gcf,'color','w');
        FigAx = axes(FigHandle);
        colormap(cmap);
        hold on
       
        legend_labels = {};
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
                    x = linspace(FitXLims(1), FitXLims(2), 1000);
                    log_fity = LogA_vals(TimeIndex)-(Ea_vals(TimeIndex))*x;
                    fity = exp(log_fity);
                else
                    x = linspace(16, 29, 1000);
                    log_fity = LogA_vals(TimeIndex)-(Ea_vals(TimeIndex)/R)*(1./(x+273.15));
                    fity = exp(log_fity);

                end
                
                if ~strcmpi(hmmVarStrings{i}, 'durations')
                    if NumTimeBins > 1
                        legend_labels{end+1} = strcat('$t = ',sprintf('%.1f ', TimeVector(TimeIndex)),'\,\, \textrm{min},\,\,E_A = ', sprintf('%.1f ', Ea_vals(TimeIndex)), ' \pm ', sprintf('%.1f ', Ea_val_errs(TimeIndex)), '\,\,\textrm{kJ}\cdot\textrm{mol}^{-1}$');
                
                    else
                        legend_labels{end+1} = strcat('$E_A = ', sprintf('%.1f ', Ea_vals(TimeIndex)), ' \pm ', sprintf('%.1f ', Ea_val_errs(TimeIndex)), '\,\,\textrm{kJ}\cdot\textrm{mol}^{-1}$');
                    end
                else
                    if NumTimeBins > 1
                        legend_labels{end+1} = strcat('$t = ',sprintf('%.1f ', TimeVector(TimeIndex)),'\,\, \textrm{min},\,\,E_A = ', sprintf('%.1f ', -1*Ea_vals(TimeIndex)), ' \pm ', sprintf('%.1f ', Ea_val_errs(TimeIndex)), '\,\,\textrm{kJ}\cdot\textrm{mol}^{-1}$');
                    else
                        legend_labels{end+1} = strcat('$E_A = ', sprintf('%.1f ', -1*Ea_vals(TimeIndex)), ' \pm ', sprintf('%.1f ', Ea_val_errs(TimeIndex)), '\,\,\textrm{kJ}\cdot\textrm{mol}^{-1}$');
                    end
                end
                
                FitPlots{t} = plot( FigAx,x, fity,'-', 'Color', 'k',...%cmap(TimeIndex,:),
                    'LineWidth', 2,...
                    'DisplayName', legend_labels{end});   
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
                SetIndices = find(~isnan(y_vals(:,TimeIndex)));
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
                SetIndices = find(~isnan(y_vals));
                x = CompiledParameters.SetTemperatures(set_ids);
                if useLogVector(i)
                    x = 1./(R*(x+273.15));
                end
                y = y_vals(set_ids);
                y_err = y_val_errs(set_ids);
                
            end
            for xIndex = 1:length(x)
                ErrorPlots{SetIndices(xIndex),t} = errorbar(FigAx, x(xIndex),y(xIndex),y_err(xIndex),'Color','k','Capsize',0, 'LineStyle', 'None');
                hold on
                DataPlots{SetIndices(xIndex),t} = scatter(FigAx, x(xIndex),y(xIndex),80, 'MarkerFaceColor',cmap(TimeIndex,:),'MarkerEdgeColor','k');
                set(get(get(ErrorPlots{SetIndices(xIndex),t}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                set(get(get(DataPlots{SetIndices(xIndex),t}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            end
        end
        
        
        

        
        %ylabel(YLabel)
        if strcmpi(hmmVarStrings{i}, 'initiationrates')
            ylabel('$r(x) \textrm{ (AU/min)}$', 'Interpreter', 'latex');
        elseif strcmpi(hmmVarStrings{i}, 'meaninitiationrates')
            ylabel('$R(x) \textrm{ (AU/min)}$', 'Interpreter', 'latex');
        elseif strcmpi(hmmVarStrings{i}, 'Durations')
            ylabel('$k^{-1}_{\textrm{off}}(x)\textrm{ (min)}$', 'Interpreter', 'latex');
        elseif strcmpi(hmmVarStrings{i}, 'frequencies')
            ylabel('$k_{\textrm{on}}(x)\textrm{ (1/min)}$', 'Interpreter', 'latex');
        else
            ylabel('$k^{-1}_{\textrm{off}}(x)+k^{-1}_{\textrm{on}}(x)\textrm{ (min)}$', 'Interpreter', 'latex');
        end
        set(FigAx ,'Fontsize',16)
        if ~useLogVector(i)
            xlim([16, 29])
            xlabel('Temperature (ºC)')
            xticks([17.5 20 22.5 25 27.5])
        else
            xlim([.3975 .415])
            xlabel('1/RT (mol \cdot kJ^{-1})')
        end
        %xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
        title({['AP: ', num2str((APbin-1)*APResolution*100), '%']});
        
        
        if useLogVector(i)
            
            newylim = get(FigAx, 'ylim');
            newylim(1) = LogYmin;%LogAPbinYmin;
            newylim(2) = LogYmax;%LogAPbinYmax;
            set(FigAx, 'ylim', newylim)
            set(FigAx, 'YScale', 'log');
        else
            newylim = get(FigAx, 'ylim');
            newylim(1) = Ymin;% APbinYmin;
            newylim(2) = Ymax;%APbinYmax;
            set(FigAx, 'ylim', newylim)
        end
        %     ylim([DurYmin, DurYmax])
        
        grid on
        box on
        hold off
        if IncludeFits & ~(strcmpi(hmmVarStrings{i}, 'BurstCycleTimes') | strcmpi(hmmVarStrings{i}, 'MeanInitiationRates')) & ~isempty(legend_labels)
            if (~useLogVector(i) & strcmpi(hmmVarStrings{i}, 'initiationrates')) | (useLogVector(i) & ~strcmpi(hmmVarStrings{i}, 'initiationrates'))
              
                LegendAxes{TimeIndex} = legend('show', 'Interpreter', 'latex', 'FontSize', 16,...
                    'Location', 'northwest');
            else
                LegendAxes{TimeIndex} = legend('show', 'Interpreter', 'latex', 'FontSize', 16,...
                    'Location', 'northeast');
            end
        end

        
%         if max(IncludedTimeBins) > 1
%             set(FigAx, 'Clim', [1, max(IncludedTimeBins)+1])
%             h = colorbar();
%             h.Ticks = 1.5:max(IncludedTimeBins)+0.5;
%             h.TickLabels = TimeVector(1:max(IncludedTimeBins));
%             if useRescaledTiming | useRescaledParamTiming
%                 ylabel(h,['dev time cohort (',TemperatureString,'ºC minutes into nc14)'], 'FontSize', 18)
%             else
%                 ylabel(h,'time cohort (minutes into nc14)', 'FontSize', 18)
%             end
%         end
        hold off
        if strcmpi(hmmVarStrings{i}, 'initiationrates')  | strcmpi(hmmVarStrings{i}, 'meaninitiationrates') 
            OutFileName = [LogString, FluoString, RescaledString,FitString,ParamString, 'APbin',num2str(APbin),'_', OutString,'VsTemp.png'];
        else
            OutFileName = [LogString, RescaledString,FitString,ParamString, 'APbin',num2str(APbin),'_', OutString,'VsTemp.png'];
        end
        
        saveas(FigHandle,[outdir, OutFileName])
        
        
    end


    
    
    
end
%%

close all

%%



close all

%%

LegHandles = figure(2);
set(LegHandles,'units', 'normalized', 'position',[0.05, 0.05, 0.075, 0.4]);
set(LegHandles,'color','w');
ax2 = axes(LegHandles);
colormap(cmap)
set(ax2, 'Clim',  [1, max(IncludedTimeBins)+1]);

h = colorbar(ax2, 'west');% 'YDir', 'reverse');
h.Ticks = 1.5:max(IncludedTimeBins)+0.5;
 h.TickLabels = TimeVector(1:max(IncludedTimeBins));


                ylabel(h,'dev. time cohort (minutes into nc14)', 'FontSize', 18)


axis off
x1=get(ax2,'position');
x=get(h,'Position');
x(1) = .1;
x(2) = .025;
x(3)=0.2;
x(4)=0.95;
set(h,'Position',x)
set(ax2,'position',x1)
%h.Ticks = 1.5:5.5;
%h.TickLabels = T_vector(T_matches)-273.15;
set(ax2,'Fontsize',16)
saveas(LegHandles, [outdir, filesep, 'TemperatureLegend.png'])


















