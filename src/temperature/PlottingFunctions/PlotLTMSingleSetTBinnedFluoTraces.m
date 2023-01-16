function PlotLTMSingleSetTBinnedFluoTraces(this, outdir, varargin)
%%

% PlotTitle, PlottingColors, UseDifferentColors,
% UseDiffProfiles, UsePhysicalAPLength

NormedCycleTimes = false;
useRescaledTiming = false;
useRescaledFluo = false;
SuppressErrorbars = false;
SuppressMarkers = false;
DownsampleTraces = false;
DownsamplingRate = 1;
UsePerNucleusTraces = false;
UseFractionOns = false;
AltFractionOns = false;
IncludeFits = true;

x = 1;
while x <= length(varargin)
    if strcmp(lower(varargin{x}), 'plottitle')
        PlotTitle = varargin{x+1};
        x = x+1;
    elseif strcmp(lower(varargin{x}), 'plottingcolors')
        PlottingColors = varargin{x+1};
        x = x+1;
    elseif strcmp(lower(varargin{x}), 'tracetype')
        TraceType = lower(varargin{x+1});
        x = x+1;
    elseif strcmpi(varargin{x}, 'NormalizeCycleTimes')
        NormedCycleTimes = true;
    elseif strcmpi(varargin{x}, 'SuppressErrorbars')
        SuppressErrorbars = true;
    elseif strcmpi(varargin{x}, 'SuppressMarkers')
        SuppressMarkers = true;
    elseif strcmpi(varargin{x}, 'UsePerNucleusTraces')
        UsePerNucleusTraces = true;
    elseif strcmpi(varargin{x}, 'UseFractionOns')
        UseFractionOns = true;
    elseif strcmpi(varargin{x}, 'UseAltFractionOns')
        UseFractionOns = true;
        AltFractionOns = true;
    elseif strcmpi(varargin{x}, 'DownsampleTraces')
        DownsampleTraces = true;
        DownsamplingRate = 2;
    elseif strcmpi(varargin{x}, 'DownsamplingRate')
        DownsampleTraces = true;
        DownsamplingRate = varargin{x+1};
        x = x+1;
    elseif strcmpi(varargin{x}, 'userescaledtime') | strcmpi(varargin{x}, 'rescaletime') | ...
            strcmpi(varargin{x}, 'rescaletiming') | strcmpi(varargin{x}, 'userescaledtiming')
        useRescaledTiming = true;
    elseif strcmpi(varargin{x}, 'rescalefluo') | strcmpi(varargin{x}, 'userescaledfluo')
        useRescaledFluo = true;
    elseif strcmp(lower(varargin{x}), 'excludefits')
        IncludeFits = false;
    else
        error(['Invalid option ',varargin{x},' passed.'])
    end
    x = x+1;
end

if ~exist('PlottingColors', 'var')
    PlottingColors = 'thesis';
elseif ~strcmp(lower(PlottingColors), 'default')  & ~strcmp(lower(PlottingColors), 'pboc') & ~strcmp(lower(PlottingColors), 'thesis')
    error('Invalid choice of plotting colors. Can use either "default", "pboc", or "thesis".') % change to error
end
if ~exist('TraceType', 'var')
    TraceType = 'anaphasealigned';
elseif ~strcmp(lower(TraceType), 'tbinned3d') & ~strcmp(lower(TraceType), 'tbinned')  & ~strcmp(lower(TraceType), 'anaphasealigned')& ~strcmp(lower(TraceType), 'anaphasealigned3d')
    error('Invalid choice of trace type. Can use either "tbinned", "tbinned3d", "anaphasealigned", or "anaphasealigned3d".') % change to error
end

if ~exist(outdir, 'dir')
    mkdir(outdir)
end

if NormedCycleTimes & useRescaledTiming
    error('Cannot use NormedCycleTimes & useRescaledTiming options simulatenously.')
end
%%
if strcmpi(TraceType, 'anaphasealigned')
    traceName = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'anaphasealigned3d')
    traceName = 'AnaphaseAligned3D';
elseif strcmpi(TraceType, 'tbinned')
    traceName = 'Tbinned';
elseif strcmpi(TraceType, 'tbinned3d')
    traceName = 'Tbinned3D';
end

%%


if ~exist(outdir, 'dir')
    mkdir(outdir)
end



if NormedCycleTimes
    NormString = 'Normalized';
else
    NormString = '';
end


if useRescaledTiming
    TimingString = 'DevTime';
else
    TimingString = '';
end

if useRescaledFluo
    FluoString = 'RescaledFluo';
else
    FluoString = '';
end

if SuppressErrorbars
    ErrorbarString = 'NoErrorbars';
else
    ErrorbarString = '';
end

if UsePerNucleusTraces
    PerNucleusString = 'PerNucleus';
else
    PerNucleusString = '';
end


if SuppressMarkers
    MarkerString = 'NoMarkers';
else
    MarkerString = '';
end

if DownsampleTraces
    DownsamplingString = ['Downsampled', num2str(DownsamplingRate),'x'];
else
    DownsamplingString = '';
end

SpecificDirString = PerNucleusString;
if ~isempty(FluoString)
    if ~isempty(SpecificDirString)
        SpecificDirString = [SpecificDirString, '_', FluoString];
    else
        SpecificDirString = FluoString;
    end
end
if ~isempty(TimingString)
    if ~isempty(SpecificDirString)
        SpecificDirString = [SpecificDirString, '_', TimingString];
    else
        SpecificDirString = TimingString;
    end
elseif ~isempty(NormString)
    if ~isempty(SpecificDirString)
        SpecificDirString = [SpecificDirString, '_', NormString];
    else
        SpecificDirString = NormString;
    end
end

if ~isempty(MarkerString)
    if ~isempty(SpecificDirString)
        SpecificDirString = [SpecificDirString, '_', MarkerString];
    else
        SpecificDirString = MarkerString;
    end
elseif ~isempty(ErrorbarString)
    if ~isempty(SpecificDirString)
        SpecificDirString = [SpecificDirString, '_', ErrorbarString];
    else
        SpecificDirString = ErrorbarString;
    end
end

if ~isempty(DownsamplingString)
    if ~isempty(SpecificDirString)
        SpecificDirString = [SpecificDirString, '_', DownsamplingString];
    else
        SpecificDirString = DownsamplingString;
    end
end


%%

NumSets = length(this.ExperimentPrefixes);
temperatures = flip(unique(this.Temp_sps));
NumTemperatures = length(temperatures);

APResolution = this.Experiments{1}.APResolution;
APbins = 0:APResolution:1;

if strcmp(lower(PlottingColors), 'default')
    [~, colors, ~] = getColorPalettes();
elseif strcmp(lower(PlottingColors), 'pboc')
    [colors, ~, ~] = getColorPalettes();
elseif strcmp(lower(PlottingColors), 'thesis')
    [~, ~, colors] = getColorPalettes();
end

legend_labels = {};
for i = 1:NumTemperatures
    legend_labels{end+1} = [num2str(temperatures(i)), 'ºC'];
end
MinimumEmbryoCount = this.MinimumEmbryos;


MarkerStyles = {'o', 'd', 's', '>', '^','p', 'h', '*', 'x'};

Nsigfigs = 3;
%%

if ~UseFractionOns
    outdir2 = [outdir,filesep,'TemperatureBinnedFluoTraces'];
elseif ~AltFractionOns
    outdir2 = [outdir,filesep,'TemperatureBinnedFractionOnTraces'];
else
    outdir2 = [outdir,filesep,'TemperatureBinnedAltFractionOnTraces'];
end
if ~exist(outdir2, 'dir')
    mkdir(outdir2)
end

if isempty(SpecificDirString)
    outdir3 = [outdir2, filesep, TraceType];
    if ~exist(outdir3, 'dir')
        mkdir(outdir3)
    end
else
    outdir3 = [outdir2, filesep, TraceType, filesep, SpecificDirString];
    if ~exist(outdir3, 'dir')
        mkdir(outdir3)
    end
end



outdir4 =  [outdir3, filesep, datestr(now, 'yyyymmdd')];
if ~exist(outdir4, 'dir')
    mkdir(outdir4)
end
%%
for nc_idx=4:4%1:length(this.IncludedNCs)
    
    NC = this.IncludedNCs(nc_idx);
    if NormedCycleTimes & NC == 14
        continue
    end
    if NC < 13
        continue
    end
    
    
    
    
    MeanFluoMats = cell(1, NumTemperatures);
    StdFluoMats = cell(1, NumTemperatures);
    NumEmbryoMats = cell(1, NumTemperatures);
    NCTimes = cell(1, NumTemperatures);
    MaximumNCTimes = NaN(1, NumTemperatures);
    MaxFluos = NaN(1, NumTemperatures);
    MinFluos = NaN(1, NumTemperatures);
    NumFrames = NaN(1, NumTemperatures);
    MinAPs = NaN(1, NumTemperatures);
    MaxAPs = NaN(1, NumTemperatures);
    
    for idx=1:NumTemperatures
        ExpNCTimes = this.BinnedMeanProfiles.([traceName, 'CycleTimes'])/60;
        IncludedRows = 1:DownsamplingRate:length(ExpNCTimes);
        if ~ismember(length(ExpNCTimes), IncludedRows)
            IncludedRows(end+1) = length(ExpNCTimes);
        end
        
        ExpNCTimes = ExpNCTimes(IncludedRows);
        if ~UseFractionOns & ~AltFractionOns
            if UsePerNucleusTraces
                ExpFluoMat = this.BinnedMeanProfiles.([traceName, 'CycleNuclearMeanTraces'])(IncludedRows,:,NC-8,idx);
                ExpStdMat =  this.BinnedMeanProfiles.([traceName, 'CycleNuclearStdErrors'])(IncludedRows,:,NC-8,idx);
            else
                ExpFluoMat = this.BinnedMeanProfiles.([traceName, 'CycleMeanTraces'])(IncludedRows,:,NC-8,idx);
                ExpStdMat =  this.BinnedMeanProfiles.([traceName, 'CycleStdErrors'])(IncludedRows,:,NC-8,idx);
            end
        elseif ~AltFractionOns
            ExpFluoMat = this.BinnedMeanProfiles.([traceName, 'CycleTotalOnNuclei'])(IncludedRows,:,NC-8,idx)./...
                this.BinnedMeanProfiles.([traceName, 'CycleTotalNuclei'])(IncludedRows,:,NC-8,idx);
            ExpStdMat = NaN(size(ExpFluoMat));
        else
            ExpFluoMat = this.BinnedMeanProfiles.([traceName, 'CycleFractionOn'])(IncludedRows,:,NC-8,idx);
            ExpStdMat = NaN(size(ExpFluoMat));
        end
        ExpNumEmbryoMat = this.BinnedMeanProfiles.([traceName, 'CycleNumEmbryos'])(IncludedRows,:,NC-8,idx);
        NCLength = this.TimeScalingInfo.MeanNCDivisionInfo(idx,NC-8);
        
        
        IncludedColumns = find(sum(~isnan(ExpFluoMat),1).' > 0);
        if ~isempty(IncludedColumns)
            MinAPs(idx) = min(IncludedColumns);
            MaxAPs(idx) = max(IncludedColumns);
        end
        
        IncludedRows = find(sum(~isnan(ExpFluoMat),2).' > 0);
        if isempty(IncludedRows)
            MeanFluoMats{idx} = [];
            StdFluoMats{idx} = [];
            NumEmbryoMats{idx} = [];
            NCTimes{idx} = [];
        else
            
            MeanFluoMats{idx} = ExpFluoMat(IncludedRows,:);
            StdFluoMats{idx} = ExpStdMat(IncludedRows,:);
            NumEmbryoMats{idx} = ExpNumEmbryoMat(IncludedRows,:);
            
            if useRescaledFluo
                MeanFluoMats{idx} = MeanFluoMats{idx}*this.FluoCoeffs(idx);
                StdFluoMats{idx} = StdFluoMats{idx}*this.FluoCoeffs(idx);
            end
            
            
            if NormedCycleTimes
                NCTimes{idx} = ExpNCTimes(IncludedRows)/NCLength;
            elseif useRescaledTiming
                if NC == 14
                    if ~isnan(this.TimeScalingInfo.PropNCDivisionInfo(idx, NC-9))
                        NCTimes{idx} = ExpNCTimes(IncludedRows)/this.TimeScalingInfo.PropNCDivisionInfo(idx, NC-9);
                    else
                        NCTimes{idx} = this.TimingCoeffs(idx)*ExpNCTimes(IncludedRows);
                    end
                elseif ~isnan(this.TimeScalingInfo.PropNCDivisionInfo(idx, NC-8))
                    NCTimes{idx} = ExpNCTimes(IncludedRows)/this.TimeScalingInfo.PropNCDivisionInfo(idx, NC-8);
                else
                    NCTimes{idx} = this.TimingCoeffs(idx)*ExpNCTimes(IncludedRows);
                end
            else
                NCTimes{idx} = ExpNCTimes(IncludedRows);
            end
            
            MaximumNCTimes(idx) = max(NCTimes{idx});
            if ~SuppressErrorbars & ~SuppressMarkers
                MaxFluos(idx) = max(max(MeanFluoMats{idx}+StdFluoMats{idx}));
                MinFluos(idx) = min(min(MeanFluoMats{idx}-StdFluoMats{idx}));
            else
                MaxFluos(idx) = max(max(MeanFluoMats{idx}));
                MinFluos(idx) = min(min(MeanFluoMats{idx}));
            end
            NumFrames(idx) = length(NCTimes{idx});
        end
    end
    
    if all(isnan(NumFrames))
        continue
    end
    
    MinAPbin = nanmin(MinAPs);
    MaxAPbin = nanmax(MaxAPs);
    
    
    close all
    
    FrameProfFig = figure(1);
    set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .5, .4]);
    set(gcf,'color','w');
    FrameProfAx = axes(FrameProfFig);
    
    eb = cell(1, NumTemperatures);
    prof = cell(1, NumTemperatures);
    
    if IncludeFits
        ci_plotline = fill([0, 1, 1, 0], [0, 0, max([MaxFluos*1.2,1]), max([MaxFluos*1.2,1])], 'k');%colors(temp_idx,:));
        ci_plotline.FaceAlpha = 0.2;
        set(ci_plotline, 'EdgeColor', 'none');
        hold on
        t_on_plotline = xline(0,'-', 'color', [0,0,0]);
        t_off_plotline = xline(0, '--',  'color', [0,0,0]);
        t_peak_plotline = xline(0, '-.',  'color', [0,0,0]);
        %t_peak_plotline = xline(0, '--',  'color', [0,0,0]+0.5);
        
        
        fitted_prof = plot(APbins, ones(1, length(APbins)), '-', 'Color','k', 'Linewidth', 1.5); %colors(temp_idx,:));
        set(fitted_prof,'Visible','off'); %'off' or 'on'
        set(t_on_plotline,'Visible','off'); %'off' or 'on'
        set(t_off_plotline,'Visible','off'); %'off' or 'on'
        set(t_peak_plotline,'Visible','off'); %'off' or 'on'
        set(ci_plotline,'Visible','off'); %'off' or 'on'
    end
    
    
    for idx =1:NumTemperatures
        temp_idx = idx;
        marker_idx = 1;
        
        eb{idx} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none',...
            'CapSize', 0);
        hold on
        set(eb{idx}, 'color', colors(temp_idx,:), 'LineWidth', 1.5);
        set(get(get(eb{idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        set(eb{idx},'Visible','off'); %'off' or 'on'
    end
    
    for idx =1:NumTemperatures
        temp_idx = idx;
        marker_idx = 1;
        if ~SuppressErrorbars & ~SuppressMarkers
            prof{idx} = plot(APbins, ones(1, length(APbins)), MarkerStyles{marker_idx},...
                'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colors(temp_idx,:),...
                'MarkerSize', 5, 'LineStyle', 'none', 'Color', colors(temp_idx,:),...
                'LineWidth', 1.5);
        elseif ~SuppressMarkers
            prof{idx} = plot(APbins, ones(1, length(APbins)), MarkerStyles{marker_idx},...
                'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colors(temp_idx,:),...
                'MarkerSize', 5, 'LineStyle', '-', 'Color', colors(temp_idx,:),...
                'LineWidth', 1.5);
        else
            prof{idx} = plot(APbins, ones(1, length(APbins)), '-',...
                'Color', colors(temp_idx,:),'LineWidth', 3);
        end
        
        
        
        set(prof{idx},'Visible','off'); %'off' or 'on'
        set(get(get(prof{idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    end
    
    
    hold off
    if NormedCycleTimes
        xlabel('Fraction of Nuclear Cycle')
    elseif useRescaledTiming & strcmpi(TraceType, 'anaphasealigned')
        xlabel('Re-scaled time since anaphase (25 °C min)')
    elseif useRescaledTiming
        xlabel('Re-scaled time since NC start (25 °C min)')
    elseif strcmpi(TraceType, 'anaphasealigned')
        xlabel('Time since anaphase (min)')
    else
        xlabel('Time since NC start (min)')
    end
    
    
    if ~UseFractionOns
        if ~useRescaledFluo
            if strcmp(lower(TraceType), 'fluo') | strcmp(lower(TraceType), 'anaphasealigned') | strcmp(lower(TraceType), 'tbinned')
                ylabel('Fluo (AU)')
            else
                ylabel('3D Fluo (AU)')
            end
        else
            if strcmp(lower(TraceType), 'fluo') | strcmp(lower(TraceType), 'anaphasealigned') | strcmp(lower(TraceType), 'tbinned')
                ylabel('Re-scaled fluo (AU)')
            else
                ylabel('Re-scaled 3D fluo (AU)')
            end
            
        end
        ylim([floor(min([min(MinFluos),0])/100)*100, ceil(max(max(MaxFluos)*1.1,1)/200)*200])
    else
        ylim([ -0.05 1.05])
        ylabel('Fraction Active Nuclei')
    end
    %
    
    FrameProfAx.XAxis.FontSize = 10;
    FrameProfAx.YAxis.FontSize = 10;
    title(FrameProfAx, {'',...
        ['Nuclear Cycle ',num2str(NC),', Fraction Embryo Length: ',num2str(-1) ]}, 'FontSize', 10)
    for idx=1:NumTemperatures
        ci_plotline.FaceColor = colors(idx,:);
        %fitted_prof.Color = colors(idx,:);
        xlim([0, ceil(MaximumNCTimes(idx)/5)*5])
        current_prof = prof{idx};
        for i = MinAPbin:MaxAPbin
            PlottedSets = false;
            APBinHasData = false;
            
            
            set(get(get(FrameProfAx.Children(end-NumTemperatures-idx-4), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            set(FrameProfAx.Children(end-NumTemperatures-idx-4),'Visible','off'); %'off' or 'on'
            set(FrameProfAx.Children(end-idx-4),'Visible','off'); %'off' or 'on'
            if isempty(NCTimes{idx})
                continue
            end
            use_idx = NumEmbryoMats{idx}(:,i) >= MinimumEmbryoCount;
            
            if sum(use_idx) == 0 %| sum(DiffMeanFluoMat(i, use_idx, j) == 0)
                FrameProfAx.Children(end-NumTemperatures-idx-4).XData = NCTimes{idx}(use_idx);
                FrameProfAx.Children(end-NumTemperatures-idx-4).YData = zeros(1, length(NCTimes{idx}(use_idx)));
                FrameProfAx.Children(end-idx-4).XData = NCTimes{idx}(use_idx);
                FrameProfAx.Children(end-idx-4).YData = zeros(1, length(NCTimes{idx}(use_idx)));
                FrameProfAx.Children(end-idx-4).YPositiveDelta = zeros(1, length(NCTimes{idx}(use_idx)));
                FrameProfAx.Children(end-idx-4).YNegativeDelta = zeros(1, length(NCTimes{idx}(use_idx)));
                set(FrameProfAx.Children(end-NumTemperatures-idx-4),'Visible','off'); %'off' or 'on'
                set(FrameProfAx.Children(end-idx-4),'Visible','off'); %'off' or 'on'
                set(get(get(FrameProfAx.Children(end-NumTemperatures-idx-4), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
            else
                APBinHasData = true;
                PlottedSets= true;
                FrameProfAx.Children(end-NumTemperatures-idx-4).YData = MeanFluoMats{idx}(use_idx, i);
                FrameProfAx.Children(end-NumTemperatures-idx-4).XData =NCTimes{idx}(use_idx);
                FrameProfAx.Children(end-idx-4).YData = MeanFluoMats{idx}(use_idx, i);
                FrameProfAx.Children(end-idx-4).XData = NCTimes{idx}(use_idx);
                FrameProfAx.Children(end-idx-4).YPositiveDelta = StdFluoMats{idx}(use_idx, i);
                FrameProfAx.Children(end-idx-4).YNegativeDelta  = StdFluoMats{idx}(use_idx, i);
                
                set(FrameProfAx.Children(end-NumTemperatures-idx-4),'Visible','on'); %'off' or 'on'
                
                if ~SuppressErrorbars & ~SuppressMarkers
                    %set(FrameProfAx.Children(end-idx-4),'Visible','on'); %'off' or 'on'
                end
                set(get(get(FrameProfAx.Children(end-NumTemperatures-idx-4), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
                if IncludeFits
                    plot_legend_labels = cell(1, 5);
                    
                    plot_legend_labels{1} = legend_labels{idx};
                    
                    [t_vector, fit_solution, ci, pos_slope, se_pos_slope, neg_slope,...
                        se_neg_slope, time_on, se_time_on, time_off,...
                        se_time_off, time_peak, se_time_peak, R2] =...
                        getFittedTrace(this, idx, i, NC, TraceType, max(MaximumNCTimes), true);
                    if ~isempty(t_vector)
                        
                        FrameProfAx.Children(11).YData = fit_solution;
                        FrameProfAx.Children(11).XData =t_vector;
                        set(FrameProfAx.Children(11),'Visible','on'); %'off' or 'on'
                        if isnan(neg_slope) | (NC == 14)
                            lab2 = MeanSE_num2str(pos_slope, se_pos_slope, Nsigfigs);
                            plot_legend_labels{2} = ['Loading Rate: ', lab2.m,...
                                ' \pm ', lab2.se, ' AU/min, R^2: ', num2str(R2, 3)];
                            
                            
                        else
                            lab2a = MeanSE_num2str(pos_slope, se_pos_slope, Nsigfigs);
                            lab2b = MeanSE_num2str(neg_slope, se_neg_slope, Nsigfigs);
                            plot_legend_labels{2} = ['Loading Rate: ', lab2a.m,...
                                ' \pm ', lab2a.se, ' AU/min, R^2: ', num2str(R2, 3),'\newlineUnloading Rate: ',...
                                lab2b.m,' \pm ', lab2b.se, ' AU/min'];
                        end
                        if ~isempty(ci)
                            curve1 = ci(:,1).';
                            curve2 = ci(:,2).';
                            curve1(curve1 > max([MaxFluos*1.05,1])) = max([MaxFluos*1.05,1]);
                            curve1(curve1 < 0) = 0;
                            curve2(curve2 > max([MaxFluos*1.05,1])) = max([MaxFluos*1.05,1]);
                            curve2(curve2 < 0) = 0;
                            t_vector2 = [t_vector, fliplr(t_vector)];
                            inBetween = [curve1, fliplr(curve2)];
                            FrameProfAx.Children(15).Faces= 1:length(t_vector2);
                            FrameProfAx.Children(15).Vertices = [t_vector2; inBetween].';
                            set(FrameProfAx.Children(15),'Visible','on');
                            plot_legend_labels{3} = 'Prediction 95% Confidence Interval';
                        else
                            set(FrameProfAx.Children(15),'Visible','off');
                            plot_legend_labels{3} = '';
                        end
                        
                        
                        
                        
                    else
                        FrameProfAx.Children(11).XData = NCTimes{idx};
                        FrameProfAx.Children(11).YData = zeros(1, length(NCTimes{idx}));
                        set(FrameProfAx.Children(11),'Visible','off'); %'off' or oon'
                        plot_legend_labels{2} = 'No Fitting Info';
                        set(FrameProfAx.Children(15),'Visible','off');
                        plot_legend_labels{3} = '';
                    end
                    if ~isnan(time_peak) & ~isnan(time_on)
                        %                         FrameProfAx.Children(4).Vertices(:, 1) = [time_on; time_peak; time_peak; time_on];
                        %                         set(FrameProfAx.Children(4),'Visible','on');
                        %                         se_time_elongation = sqrt(se_time_peak^2+se_time_on^2);
                        %                         lab3 = MeanSE_num2str(time_peak-time_on, se_time_elongation, Nsigfigs);
                        %                         plot_legend_labels{5} = ['t_{elongation}: ', lab3.m, ' \pm ', lab3.se,...
                        %                             ' min'];%, t_{elongation}: ', num2str(time_peak-time_on, 3), ' min'];
                        FrameProfAx.Children(12).Value = time_peak;
                        set(FrameProfAx.Children(12),'Visible','on');
                        lab3 = MeanSE_num2str(time_peak, se_time_peak, Nsigfigs);
                        time_elongation = time_peak -time_on;
                        se_time_elongation = sqrt(se_time_peak^2 + se_time_on^2);
                        lab8 = MeanSE_num2str(time_elongation, se_time_elongation, Nsigfigs);
                        plot_legend_labels{6} = ['t_{peak}: ', lab3.m,' \pm ', lab3.se,...
                            ' min\newlinet_{elongation}: ', lab8.m,' \pm ', lab8.se, ' min'];
                        
                    else
                        FrameProfAx.Children(12).Value = 0;
                        set(FrameProfAx.Children(12),'Visible','off');
                        plot_legend_labels{6} = 'No t_{peak} info';
                    end
                    
                    if ~isnan(time_off)
                        FrameProfAx.Children(13).Value = time_off;
                        set(FrameProfAx.Children(13),'Visible','on');
                        lab4 = MeanSE_num2str(time_off, se_time_off, Nsigfigs);
                        plot_legend_labels{5} = ['t_{off}: ', lab4.m,' \pm ', lab4.se,  ' min'];
                    else
                        FrameProfAx.Children(13).Value = 0;
                        set(FrameProfAx.Children(13),'Visible','off');
                        plot_legend_labels{5} = 'No t_{off} info';
                    end
                    
                    if ~isnan(time_on)
                        FrameProfAx.Children(14).Value = time_on;
                        set(FrameProfAx.Children(14),'Visible','on');
                        lab5 = MeanSE_num2str(time_on, se_time_on, Nsigfigs);
                        plot_legend_labels{4} =  ['t_{on}: ', lab5.m,...
                            ' \pm ', lab5.se, ' min'];
                    else
                        FrameProfAx.Children(14).Value = 0;
                        set(FrameProfAx.Children(14),'Visible','off');
                        plot_legend_labels{4} = 'No t_{on} info';
                    end
                    
                    
                    
                    
                end
                
                
            end

            
            %try
            %legend(FrameProfAx, legend_labels(PlottedSets), 'Location', 'northeast', 'FontSize', 18)
            if exist('PlotTitle', 'var')
                title(FrameProfAx, {PlotTitle,...
                    ['Nuclear Cycle ',num2str(NC),', Fraction Embryo Length: ',num2str(APbins(i)) ]})
                
                
            else
                title(FrameProfAx, ['Nuclear Cycle ',num2str(NC),', Fraction Embryo Length: ',num2str(APbins(i)) ])
                
            end
            
            
            %end
            if PlottedSets
                if (NC == 14) & IncludeFits
                    plot_handles ={current_prof, fitted_prof, ci_plotline, t_on_plotline, t_off_plotline,  t_peak_plotline};
                    hlegend = legend([plot_handles{:}], plot_legend_labels, 'Location', 'eastoutside',...
                        'FontSize', 8);
                    
                elseif (NC < 14) & IncludeFits
                    plot_handles ={current_prof, fitted_prof, ci_plotline, t_on_plotline, t_off_plotline,  t_peak_plotline};
                    hlegend = legend([plot_handles{:}], plot_legend_labels, 'Location', 'eastoutside',...
                        'FontSize', 8);
                end
                
                if ~UseFractionOns
                    outpath=[outdir4, filesep, 'FluoTrace_NC',num2str(NC),'T',strrep(num2str(temperatures(idx)), '.', '_'),...
                        'C_AP', num2str(i),'.pdf'];
                    FrameProfFig.Position = [0.01, 0.05, .5, .4];
                    exportgraphics(FrameProfFig,outpath, 'ContentType','vector');
                end
                
                
                
            end
        end
        set(FrameProfAx.Children(end-NumTemperatures-idx-4),'Visible','off'); %'off' or 'on'
        set(FrameProfAx.Children(end-idx-4),'Visible','off'); %'off' or 'on'
    end
    
end
close all
end