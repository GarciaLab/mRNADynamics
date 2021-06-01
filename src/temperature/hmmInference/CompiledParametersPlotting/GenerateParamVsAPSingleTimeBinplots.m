function GenerateParamVsAPSingleTimeBinplots(ResultsPaths, outdir, varargin)
close all



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

MinAPbin = min(IncludedAPbins);
MaxAPbin = max(IncludedAPbins);
NumAPbins = MaxAPbin-MinAPbin + 1;

cmap = brewermap(NumSets,'Spectral');

cmap = brewermap(10,'Spectral');
cmap = [cmap(2:3,:);cmap(8:10,:)];
hmmVarStrings = { 'InitiationRates','Durations', 'Frequencies','BurstCycleTimes','MeanInitiationRates',};



%%
for i = 1:length(hmmVarStrings)

    close all
    [MeanValues, StdErrors, TimeVector, Ymax, Ymin, YLabel, OutString, LogYmax, LogYmin] = ...
        GetPlottingMats(CompiledParameters, hmmVarStrings{i}, useRescaledTiming,...
        useRescaledFluo,useRescaledParamTiming, 1,ReferenceTemperature);

    
    
    
    
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
    
    
    
    for TimeIndex = 1:NumTimeBins
        TimeBin = IncludedTimeBins(TimeIndex);
        TimeBin_y_vals = squeeze(MeanValues(:,TimeBin,:));
        TimeBin_y_val_errs = squeeze(StdErrors(:,TimeBin,:));
        MaxTimeBinValue = max(TimeBin_y_vals(~isnan(TimeBin_y_vals)));
        MinTimeBinValue = min(TimeBin_y_vals(~isnan(TimeBin_y_vals)));
        
        MaxTimeBinValue = max([MaxTimeBinValue*1.1, MaxTimeBinValue*0.9]);
        MinTimeBinValue = min([MinTimeBinValue*1.1, MinTimeBinValue*0.9]);
        
        TimeBinYmax = paramYinterval*ceil(MaxTimeBinValue/paramYinterval);
        TimeBinYmin = paramYinterval*floor(MinTimeBinValue/paramYinterval);
       
        if strcmpi(hmmVarStrings{i}, 'initiationrates') | strcmpi(hmmVarStrings{i}, 'meaninitiationrates')
            LogTimeBinYmax = 10^(ceil(log10(TimeBinYmax)));
            LogTimeBinYmax= TimeBinYmax;
            LogTimeBinYmin =TimeBinYmin;
        else
            LogTimeBinYmin =TimeBinYmin;
            LogTimeBinYmax =TimeBinYmax;
        end
        close all
        FigHandle = figure(1);
        FitPlots = cell(1, NumSets);
        DataPlots = cell(1, NumSets);
        ErrorPlots = cell(1, NumSets);
        set(FigHandle,'units', 'normalized', 'position',[0.05, 0.05, 0.35, 0.4]);
        set(gcf,'color','w');
        FigAx = axes(FigHandle);
        colormap(cmap);
        hold on
       
       
        
        y_vals = squeeze(MeanValues(:,TimeBin,:));
        y_val_errs = squeeze(StdErrors(:,TimeBin,:));

        for SetIndex = 1:NumSets
            
            if all(isnan(y_vals(SetIndex,:)))
                continue
            end
            set_ids = ~isnan(y_vals(SetIndex,:));
            
            x = APbins(set_ids)/100;

            y = y_vals(SetIndex, set_ids);
            y_err = y_val_errs(SetIndex, set_ids);
            
       
            ErrorPlots{SetIndex} = errorbar(FigAx, x,y,y_err,'Color','k','Capsize',0, 'LineStyle', 'None');
            hold on
            DataPlots{SetIndex} = scatter(FigAx, x,y,80, 'MarkerFaceColor',cmap(SetIndex,:),'MarkerEdgeColor','k');
            set(get(get(ErrorPlots{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            set(get(get(DataPlots{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            
        end
        
        
        

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
        xlim([APbins(MinAPbin-1)/100 APbins(MaxAPbin+1)/100])
        xlabel('Position (x/L)')
        
        %xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
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
        
        
        newylim = get(FigAx, 'ylim');
        newylim(1) = TimeBinYmin;
        newylim(2) = TimeBinYmax;
        set(FigAx, 'ylim', newylim)
        
        %     ylim([DurYmin, DurYmax])
        
        grid on
        box on
        hold off
        


        
        
        set(FigAx, 'Clim', [1,NumSets+1]);

        
        hold off
        if strcmpi(hmmVarStrings{i}, 'initiationrates') | strcmpi(hmmVarStrings{i}, 'meaninitiationrates')
            OutFileName = [FluoString, RescaledString,ParamString, 'TimeBin',num2str(TimeBin),'_', OutString,'VsAP.png'];
        else
            OutFileName = [RescaledString,ParamString, 'TimeBin',num2str(TimeBin),'_', OutString,'VsAP.png'];
        end
        
        saveas(FigHandle,[outdir, OutFileName])
        
        
    end


    
    
    
end
%%

close all

%%

LegHandles = figure(2);
set(LegHandles,'units', 'normalized', 'position',[0.05, 0.05, 0.075, 0.4]);
set(LegHandles,'color','w');
ax2 = axes(LegHandles);
colormap(cmap)
set(ax2, 'Clim', [1,NumSets+1]);

h = colorbar(ax2, 'west', 'YDir', 'reverse');
h.Ticks = 1.5:NumSets+0.5;
h.TickLabels = CompiledParameters.UniqueTemperatures;

ylabel(h,'Temperature (ºC)','Fontsize',18)


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










