function GenerateStandardPlots(ResultsPaths, outdir)


if ~exist(outdir, 'dir')
    mkdir(outdir)
end
%%
close all

UseSharedYAxis = false;
UseSharedYAxisFrequency = false;
UseSharedYAxisInit = false;

APResolution = 0.025;
APbins = (0:APResolution:1)*100;

CompiledParameters = CompileParameters(ResultsPaths);
NumSets = size(CompiledParameters.InitiationRates, 1);



NumTimeBins  =  max(CompiledParameters.nTimebins);
if NumTimeBins > 1
    IncludedAPbins =  find((sum(squeeze(sum(~isnan(CompiledParameters.InitiationRates), 1)), 1)) > 0);
    IncludedTimeBins =  find((sum(squeeze(sum(~isnan(CompiledParameters.InitiationRates), 1)), 2)) > 0).';
else
    IncludedTimeBins = 1;
    IncludedAPbins =  find((sum(squeeze(sum(~isnan(CompiledParameters.InitiationRates), 1)), 2)) > 0).';
end
MinAPbin = min(IncludedAPbins);
MaxAPbin = max(IncludedAPbins);

InitYmin = nanmin(nanmin(nanmin(CompiledParameters.InitiationRates+CompiledParameters.InitiationRatesStdErr)));
InitYmax = nanmax(nanmax(nanmax(CompiledParameters.InitiationRates+CompiledParameters.InitiationRatesStdErr)));

InitYmin = max(floor(InitYmin/50)*50, 0);
InitYmax = ceil(InitYmax/50)*50;

DurYmin = nanmin(nanmin(nanmin(CompiledParameters.Durations+CompiledParameters.DurationsStdErr)));
DurYmax = nanmax(nanmax(nanmax(CompiledParameters.Durations+CompiledParameters.DurationsStdErr)));

DurYmin = max(floor(DurYmin/2)*2, 0);
DurYmax = min(ceil(DurYmax/2)*2, 20);


FreqYmin = nanmin(nanmin(nanmin(CompiledParameters.Frequencies+CompiledParameters.FrequenciesStdErr)));
FreqYmax = nanmax(nanmax(nanmax(CompiledParameters.Frequencies+CompiledParameters.FrequenciesStdErr)));

FreqYmin = max(floor(FreqYmin/0.1)*0.1, 0);
FreqYmax = min(ceil(FreqYmax/0.1)*0.1, 4);


for SetIndex = 1:NumSets
    fig_path = CompiledParameters.FigurePaths{SetIndex};
    
    
    dur_fig = figure;
    if max(IncludedTimeBins) > 1
        cmap = brewermap(max(IncludedTimeBins),'Spectral');
    else
        cmap = brewermap(2,'Spectral');
    end

    colormap(cmap);
    hold on
    
    for t = 1:max(CompiledParameters.nTimebins)
        if all(isnan(CompiledParameters.Durations(SetIndex, t, :)))
            continue
        end
        
        apbin_ids = ~isnan(CompiledParameters.Durations(SetIndex, t, :));
        x = APbins(apbin_ids);
        y = squeeze(CompiledParameters.Durations(SetIndex, t, apbin_ids)).';
        y_err = squeeze(CompiledParameters.DurationsStdErr(SetIndex, t, apbin_ids)).';
        
        
        errorbar(x,y,y_err,'Color','k','Capsize',0)
        hold on
        scatter(x,y,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
    end
    ax = gca;
    if max(IncludedTimeBins) > 1
        set(ax, 'Clim', [1, max(IncludedTimeBins)])
        h = colorbar;
        h.Ticks = 1.5:max(IncludedTimeBins)+0.5;
        h.TickLabels = CompiledParameters.TimeVector(1:max(IncludedTimeBins));
        ylabel(h,'time cohort (minutes into nc14)')
    end
    
    grid on
    xlabel('AP position')
    ylabel('burst duration (minutes)')
    set(gca,'Fontsize',14)
    xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
    title([CompiledParameters.ReporterLabels{SetIndex}, ' ', num2str(CompiledParameters.SetTemperatures(SetIndex)), '°C, Cycle ',num2str(CompiledParameters.NC(SetIndex)) ]);
    if UseSharedYAxis
        ylim([DurYmin, DurYmax])
        
    else
        ax = gca;
        newylim = get(ax, 'ylim');
        newylim(1) = DurYmin;
        set(ax, 'ylim', newylim);
    end
    
    if UseSharedYAxis
        saveas(dur_fig,[fig_path, 'burst_duration_sharedy.png'])
    else
        saveas(dur_fig,[fig_path, 'burst_duration.png'])
    end
    
    freq_fig = figure;
    if max(IncludedTimeBins) > 1
        cmap = brewermap(max(IncludedTimeBins),'Spectral');
    else
        cmap = brewermap(2,'Spectral');
    end
    colormap(cmap);
    hold on
    
    for t = 1:NumTimeBins
        if all(isnan(CompiledParameters.Frequencies(SetIndex, t, :)))
            continue
        end
        
        apbin_ids = ~isnan(CompiledParameters.Frequencies(SetIndex, t, :));
        x = APbins(apbin_ids);
        y = squeeze(CompiledParameters.Frequencies(SetIndex, t, apbin_ids)).';
        y_err = squeeze(CompiledParameters.FrequenciesStdErr(SetIndex, t, apbin_ids)).';
        
        
        errorbar(x,y,y_err,'Color','k','Capsize',0)
        scatter(x,y,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
    end
    ax = gca;
    if max(IncludedTimeBins) > 1
        set(ax, 'Clim', [1, max(IncludedTimeBins)])
        h = colorbar;
        h.Ticks = 1.5:max(IncludedTimeBins)+0.5;
        h.TickLabels = CompiledParameters.TimeVector(1:max(IncludedTimeBins));
        
        ylabel(h,'time cohort (minutes into nc14)')
    end
    grid on
    xlabel('AP position')
    ylabel('burst frequency (1/min)')
    set(gca,'Fontsize',14)
    xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
    title([CompiledParameters.ReporterLabels{SetIndex}, ' ', num2str(CompiledParameters.SetTemperatures(SetIndex)), '°C, Cycle ',num2str(CompiledParameters.NC(SetIndex)) ]);
    if UseSharedYAxisFrequency
        ylim([FreqYmin, FreqYmax])
        
    else
        ax = gca;
        newylim = get(ax, 'ylim');
        newylim(1) = FreqYmin;
        set(ax, 'ylim', newylim);
    end
    
    if UseSharedYAxisFrequency
        saveas(freq_fig,[fig_path, 'burst_frequency_sharedy.png'])
    else
        saveas(freq_fig,[fig_path, 'burst_frequency.png'])
    end
    
    init_fig = figure;
    if max(IncludedTimeBins) > 1
        cmap = brewermap(max(IncludedTimeBins),'Spectral');
    else
        cmap = brewermap(2,'Spectral');
    end
    colormap(cmap);
    hold on
    
    for t = 1:NumTimeBins
        if all(isnan(CompiledParameters.InitiationRates(SetIndex, t, :)))
            continue
        end
        
        apbin_ids = ~isnan(CompiledParameters.InitiationRates(SetIndex, t, :));
        x = APbins(apbin_ids);
        y = squeeze(CompiledParameters.InitiationRates(SetIndex, t, apbin_ids)).';
        y_err = squeeze(CompiledParameters.InitiationRatesStdErr(SetIndex, t, apbin_ids)).';
        
        
        errorbar(x,y,y_err,'Color','k','Capsize',0)
        scatter(x,y,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
    end
    
    ax = gca;
    if max(IncludedTimeBins) > 1
        set(ax, 'Clim', [1, max(IncludedTimeBins)])
        h = colorbar;
        h.Ticks = 1.5:max(IncludedTimeBins)+0.5;
        h.TickLabels = CompiledParameters.TimeVector(1:max(IncludedTimeBins));
        ylabel(h,'time cohort (minutes into nc14)')
    end
    grid on
    xlabel('AP position')
    ylabel('initiation rate (au/min)')
    set(gca,'Fontsize',14)
    xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
    title([CompiledParameters.ReporterLabels{SetIndex}, ' ', num2str(CompiledParameters.SetTemperatures(SetIndex)), '°C, Cycle ',num2str(CompiledParameters.NC(SetIndex)) ]);
    if UseSharedYAxisInit
        ylim([InitYmin, InitYmax])
        
    else
        ax = gca;
        newylim = get(ax, 'ylim');
        newylim(1) = InitYmin;
        set(ax, 'ylim', newylim);
    end
    
    if UseSharedYAxisInit
        saveas(init_fig,[fig_path, 'burst_initiation_sharedy.png'])
    else
        saveas(init_fig,[fig_path, 'burst_initiation.png'])
    end
    
end

%%


close all

if all(CompiledParameters.nTimebins == 1) & length(CompiledParameters.nTimebins)  > 1
    t = 1;
    temperatures = sort(unique(CompiledParameters.SetTemperatures));
    fig_path = outdir;
    dur_fig = figure;
    cmap = flipud(brewermap(length(temperatures),'Spectral'));
    colormap(cmap);
    hold on
    for SetIndex = 1:NumSets
        t_index = find(round(temperatures, 2) == round(CompiledParameters.SetTemperatures(SetIndex), 2));
        
        
        if all(isnan(CompiledParameters.Durations(SetIndex, t, :)))
            continue
        end
        
        apbin_ids = ~isnan(CompiledParameters.Durations(SetIndex, t, :));
        x = APbins(apbin_ids);
        y = squeeze(CompiledParameters.Durations(SetIndex, t, apbin_ids)).';
        y_err = squeeze(CompiledParameters.DurationsStdErr(SetIndex, t, apbin_ids)).';
        
        
        errorbar(x,y,y_err,'Color','k','Capsize',0)
        hold on
        scatter(x,y,'MarkerFaceColor',cmap(t_index,:),'MarkerEdgeColor','k')
    end
    ax = gca;

        set(ax, 'Clim', [1, length(temperatures)])
        h = colorbar;
        h.Ticks = 1:length(temperatures);
        h.TickLabels = temperatures;
        ylabel(h,' Temperature (°C)')
 
    
    grid on
    xlabel('AP position')
    ylabel('burst duration (minutes)')
    set(gca,'Fontsize',14)
    xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
    
    if UseSharedYAxis
        ylim([DurYmin, DurYmax])
        
    else
        ax = gca;
        newylim = get(ax, 'ylim');
        newylim(1) = DurYmin;
        set(ax, 'ylim', newylim);
    end
    hold off
    
    if UseSharedYAxis
        saveas(dur_fig,[outdir, 'burst_duration_sharedy.png'])
    else
        saveas(dur_fig,[outdir, 'burst_duration.png'])
    end
    
    freq_fig = figure(2);
    cmap = flipud(brewermap(length(temperatures),'Spectral'));
    colormap(cmap);
    hold on
    for SetIndex = 1:NumSets
        t_index = find(round(temperatures, 2) == round(CompiledParameters.SetTemperatures(SetIndex), 2));
        
        
        if all(isnan(CompiledParameters.Durations(SetIndex, t, :)))
            continue
        end
        
        apbin_ids = ~isnan(CompiledParameters.Frequencies(SetIndex, t, :));
        x = APbins(apbin_ids);
        y = squeeze(CompiledParameters.Frequencies(SetIndex, t, apbin_ids)).';
        y_err = squeeze(CompiledParameters.FrequenciesStdErr(SetIndex, t, apbin_ids)).';
        
        
        errorbar(x,y,y_err,'Color','k','Capsize',0)
        hold on
        scatter(x,y,'MarkerFaceColor',cmap(t_index,:),'MarkerEdgeColor','k')
    end
    ax = gca;

        set(ax, 'Clim', [1, length(temperatures)])
        h = colorbar;
        h.Ticks = 1:length(temperatures);
        h.TickLabels = temperatures;
        ylabel(h,' Temperature (°C)')
 
    
    grid on
    xlabel('AP position')
    ylabel('burst frequency (1/min)')
    set(gca,'Fontsize',14)
    xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
    
    if UseSharedYAxis
        ylim([FreqYmin, FreqYmax])
        
    else
        ax = gca;
        newylim = get(ax, 'ylim');
        newylim(1) = FreqYmin;
        set(ax, 'ylim', newylim);
    end
    
    hold off
    saveas(freq_fig,[outdir, 'all_temps_burst_frequency.png'])
    
    init_fig = figure(3);
    cmap = flipud(brewermap(length(temperatures),'Spectral'));
    colormap(cmap);
    hold on
    for SetIndex = 1:NumSets
        t_index = find(round(temperatures, 2) == round(CompiledParameters.SetTemperatures(SetIndex), 2));
        
        
        if all(isnan(CompiledParameters.Durations(SetIndex, t, :)))
            continue
        end
        
        apbin_ids = ~isnan(CompiledParameters.InitiationRates(SetIndex, t, :));
        x = APbins(apbin_ids);
        y = squeeze(CompiledParameters.InitiationRates(SetIndex, t, apbin_ids)).';
        y_err = squeeze(CompiledParameters.InitiationRatesStdErr(SetIndex, t, apbin_ids)).';
        
        
        errorbar(x,y,y_err,'Color','k','Capsize',0)
        hold on
        scatter(x,y,'MarkerFaceColor',cmap(t_index,:),'MarkerEdgeColor','k')
    end
    ax = gca;

        set(ax, 'Clim', [1, length(temperatures)])
        h = colorbar;
        h.Ticks = 1:length(temperatures);
        h.TickLabels = temperatures;
        ylabel(h,' Temperature (°C)')
 
    
    grid on
    xlabel('AP position')
    ylabel('initiation rate (au/min)')
    set(gca,'Fontsize',14)
    xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
    
    if UseSharedYAxis
        ylim([InitYmin, InitYmax])
        
    else
        ax = gca;
        newylim = get(ax, 'ylim');
        newylim(1) = InitYmin;
        set(ax, 'ylim', newylim);
    end
    

    saveas(init_fig,[outdir, 'all_temps_burst_initiation.png'])
    
end





