function PlotLTMFluoTraces(this, outdir, varargin)
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
    elseif strcmpi(varargin{x}, 'DownsampleTraces')
        DownsampleTraces = true;
        DownsamplingRate = 2;
    elseif strcmpi(varargin{x}, 'UsePerNucleusTraces')
        UsePerNucleusTraces = true;
        SuppressErrorbars = true;
    elseif strcmpi(varargin{x}, 'DownsamplingRate')
        DownsampleTraces = true;
        DownsamplingRate = varargin{x+1};
        x = x+1;
    elseif strcmpi(varargin{x}, 'userescaledtime') | strcmpi(varargin{x}, 'rescaletime') | ...
            strcmpi(varargin{x}, 'rescaletiming') | strcmpi(varargin{x}, 'userescaledtiming')
        useRescaledTiming = true;
    elseif strcmpi(varargin{x}, 'rescalefluo') | strcmpi(varargin{x}, 'userescaledfluo')
        useRescaledFluo = true;
    end
    x = x+1;
end

if ~exist('PlottingColors', 'var')
    PlottingColors = 'default';
elseif ~strcmp(lower(PlottingColors), 'gradient') & ~strcmp(lower(PlottingColors), 'default')  & ~strcmp(lower(PlottingColors), 'pboc')
    error('Invalid choice of plotting colors. Can use either "default", "pboc", or "gradient".') % change to error
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
elseif strcmpi(TraceType, 'fluo') | strcmpi(TraceType, 'unaligned') 
    traceName = 'Unaligned';
elseif strcmpi(TraceType, 'fluo3d') | strcmpi(TraceType, 'unaligned3d') 
    traceName = 'Unaligned3D';
elseif strcmpi(TraceType, 'tbinned')
    traceName = 'Tbinned';
elseif strcmpi(TraceType, 'tbinned3d')
    traceName = 'Tbinned3D';
end

%%
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

Temp_obs = this.Temp_obs;
Temp_sp = this.Temp_sps;
if strcmp(lower(PlottingColors), 'default')
    [~, colors] = getColorPalettes();
    GradString = '';
elseif strcmp(lower(PlottingColors), 'pboc')
    [colors, ~] = getColorPalettes();
    GradString = '';
else
    
    Temp_range = 15:0.1:max(Temp_obs);
    colors = jet(length(Temp_range));
    FractionalTempRange = (Temp_range-min(Temp_range))/(max(Temp_range)-min(Temp_range));
    GradString = 'Gradient';
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

if ~isempty(GradString)
    if ~isempty(SpecificDirString)
        SpecificDirString = [SpecificDirString, '_', GradString];
    else
        SpecificDirString = GradString;
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

legend_labels = this.LegendLabels;
MinimumTraceCount = this.MinimumTraceCount;


MarkerStyles = {'o', 'd', 's', '>', '^','p', 'h', '*', 'x'};

UseSet = ismember(1:NumSets, this.ProcessedExperiments);
TempMatches = cell(1, NumTemperatures);
LegendOrder = zeros(1, NumSets);

%%


counter = 0;
for t_index = 1:NumTemperatures
    TempMatches{t_index} = find((this.Temp_sps == temperatures(t_index)) & UseSet);
    TempMatches2 =  find(this.Temp_sps == temperatures(t_index));
    LegendOrder(counter + 1:counter+length(TempMatches2)) = TempMatches2;
    counter = counter + length(TempMatches);
end

outdir2 = [outdir,filesep,'FluoTraces'];
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
for nc_idx=1:length(this.IncludedNCs)
    
    NC = this.IncludedNCs(nc_idx);
    if NormedCycleTimes & NC == 14
        continue
    end
    
 
    
    
    MeanFluoMats = cell(1, NumSets);
    StdFluoMats = cell(1, NumSets);
    NumNucMats = cell(1, NumSets);
    NCTimes = cell(1, NumSets);
    MaximumNCTimes = NaN(1, NumSets);
    MaxFluos = NaN(1, NumSets);
    MinFluos = NaN(1, NumSets);
    NumFrames = NaN(1, NumSets);
    MinAPs = NaN(1, NumSets);
    MaxAPs = NaN(1, NumSets);
    
    for idx=1:NumSets
        if ~ismember(idx, this.ProcessedExperiments)
            continue
        end
        ExpNCTimes = this.MeanProfiles{idx}.([traceName, 'CycleFrameTimes']){NC-8}/60;
        if isempty(ExpNCTimes)
             MeanFluoMats{idx} = [];
            StdFluoMats{idx} = [];
            NumNucMats{idx} = [];
            NCTimes{idx} = [];
            continue
        end
        IncludedRows = 1:DownsamplingRate:length(ExpNCTimes);
        if ~ismember(length(ExpNCTimes), IncludedRows)
            IncludedRows(end+1) = length(ExpNCTimes);
        end
        
        ExpNCTimes = ExpNCTimes(IncludedRows);
        if UsePerNucleusTraces
            ExpFluoMat = squeeze(this.MeanProfiles{idx}.([traceName, 'CycleMeanPerNucleusTraces'])(IncludedRows,:,NC-8));
        else
            ExpFluoMat = squeeze(this.MeanProfiles{idx}.([traceName, 'CycleMeanTraces'])(IncludedRows,:,NC-8));
        end
        ExpStdMat = squeeze(this.MeanProfiles{idx}.([traceName, 'CycleTraceStdErrors'])(IncludedRows,:,NC-8));
        ExpNumNucMat = squeeze(this.MeanProfiles{idx}.([traceName, 'CycleNumOnNuclei'])(IncludedRows,:,NC-8));
        NCLength = this.EmbryoStats.NCDivisionInfo(idx,NC-8);
        
        
        IncludedColumns = find(sum(~isnan(ExpFluoMat),1).' > 0);
        if ~isempty(IncludedColumns)
            MinAPs(idx) = min(IncludedColumns);
            MaxAPs(idx) = max(IncludedColumns);
        end
        
        IncludedRows = find(sum(~isnan(ExpFluoMat),2).' > 0);
        if isempty(IncludedRows)
            MeanFluoMats{idx} = [];
            StdFluoMats{idx} = [];
            NumNucMats{idx} = [];
            NCTimes{idx} = [];
        else
            
            MeanFluoMats{idx} = ExpFluoMat(IncludedRows,:);
            StdFluoMats{idx} = ExpStdMat(IncludedRows,:);
            NumNucMats{idx} = ExpNumNucMat(IncludedRows,:);
            
            if useRescaledFluo
                temp_match = find(round(temperatures, 2) ==  round(Temp_sp(idx), 2));
                MeanFluoMats{idx} = MeanFluoMats{idx}*this.FluoCoeffs(temp_match);
                StdFluoMats{idx} = StdFluoMats{idx}*this.FluoCoeffs(temp_match);
            end
            
            
            if NormedCycleTimes
                NCTimes{idx} = ExpNCTimes(IncludedRows)/NCLength;
            elseif useRescaledTiming
                temp_match = find(round(temperatures, 2) ==  round(Temp_sp(idx), 2));
                if NC == 14
                    if ~isnan(this.TimeScalingInfo.PropNCDivisionInfo(temp_match, NC-9))
                        NCTimes{idx} = ExpNCTimes(IncludedRows)/this.TimeScalingInfo.PropNCDivisionInfo(temp_match, NC-9);
                    else
                        NCTimes{idx} = this.TimingCoeffs(temp_match)*ExpNCTimes(IncludedRows);
                    end
                elseif ~isnan(this.TimeScalingInfo.PropNCDivisionInfo(temp_match, NC-8))
                    NCTimes{idx} = ExpNCTimes(IncludedRows)/this.TimeScalingInfo.PropNCDivisionInfo(temp_match, NC-8);
                else
                    NCTimes{idx} = this.TimingCoeffs(temp_match)*ExpNCTimes(IncludedRows);
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
    if NumSets > 25
        set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .9, .6]);
    else
        set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    end
    set(gcf,'color','w');
    FrameProfAx = axes(FrameProfFig);
    eb = cell(1, NumSets);
    prof = cell(1, NumSets);
    
    
    
    for idx =1:NumSets
        if strcmpi(PlottingColors, 'gradient')
            temp_idx = find(abs(Temp_range - Temp_obs(idx)) == min(abs(Temp_range - Temp_obs(idx))), 1);
            temp_idx2 =  find(temperatures == Temp_sp(idx));
            if ~isempty(TempMatches{temp_idx2})
                marker_idx = find(TempMatches{temp_idx2} == idx);
            else
                marker_idx = [];
            end
           
        else
            temp_idx = find(temperatures == Temp_sp(idx));
            if ~isempty(TempMatches{temp_idx})
                marker_idx = find(TempMatches{temp_idx} == idx);
            else
                marker_idx = [];
            end
        end
        
        if isempty(marker_idx)
            marker_idx = 1;
        end
     
        eb{idx} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical',...
            'LineStyle', 'none', 'CapSize', 0);
        hold on
   
        set(eb{idx}, 'color', colors(temp_idx,:), 'LineWidth', 1);
     
        set(get(get(eb{idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        if ~SuppressErrorbars & ~SuppressMarkers
            if marker_idx <= length(MarkerStyles)
                prof{idx} = plot(APbins, ones(1, length(APbins)), MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                    'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colors(temp_idx,:),...
                    'MarkerSize', 6, 'LineStyle', '-', 'Color', colors(temp_idx,:), 'LineWidth', 1.5);
            else
                prof{idx} = plot(APbins, ones(1, length(APbins)), MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                    'MarkerEdgeColor', colors(temp_idx,:), 'MarkerFaceColor', colors(temp_idx,:),...
                    'MarkerSize', 6, 'LineStyle', '-', 'Color', colors(temp_idx,:), 'LineWidth', 1.5);
            end
        elseif ~SuppressMarkers
            if marker_idx <= length(MarkerStyles)
                prof{idx} = plot(APbins, ones(1, length(APbins)), MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                    'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colors(temp_idx,:),...
                    'MarkerSize', 6, 'LineStyle', '-', 'Color', colors(temp_idx,:), 'LineWidth', 1.5);
            else
                prof{idx} = plot(APbins, ones(1, length(APbins)), MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                    'MarkerEdgeColor', colors(temp_idx,:), 'MarkerFaceColor', colors(temp_idx,:),...
                    'MarkerSize', 6, 'LineStyle', '-', 'Color', colors(temp_idx,:), 'LineWidth', 1.5);
            end
        else
            prof{idx} = plot(APbins, ones(1, length(APbins)), '-',...
                'Color', colors(temp_idx,:), 'LineWidth', 3);
        end
        
  
 
        set(eb{idx},'Visible','off'); %'off' or 'on'
        set(prof{idx},'Visible','off'); %'off' or 'on'
        set(get(get(prof{idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    end
        grid on 
    if strcmp(lower(PlottingColors), 'gradient')
        map = colormap(colors);
        h = colorbar;
        % %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
        hold off
        colorTitleHandle = get(h,'Title');
        titleString = 'Temperature (°C)';
        set(colorTitleHandle ,'String',titleString);
        h.Ticks =  FractionalTempRange(1:25:126); %Create 8 ticks from zero to 1
        h.TickLabels = {'15','17.5','20','22.5','25','27.5'} ;
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
    
     xlim([0, ceil(max(MaximumNCTimes)*1.1/5)*5])
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
    ylim([floor(min([0 min(MinFluos)])/100)*100, ceil(max([1 max(MaxFluos)])/200)*200])
    %
    FrameProfAx.XAxis.FontSize = 18;
    FrameProfAx.YAxis.FontSize = 18;
    title(FrameProfAx, {'',...
        ['Nuclear Cycle ',num2str(NC),', Fraction Embryo Length: ',num2str(-1) ]}, 'FontSize', 20)

    for i = MinAPbin:MaxAPbin
        APBinHasData = false;
        PlottedSets = zeros(1, NumSets, 'logical');
        for idx=1:NumSets
            set(get(get(FrameProfAx.Children(end-(2*(idx-1)+1)), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','off'); %'off' or 'on'
            set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','off'); %'off' or 'on'    
            if ~ismember(idx, this.ProcessedExperiments)
                continue
            end
            if isempty(NCTimes{idx})
                continue
            end
            use_idx = NumNucMats{idx}(:,i) >= MinimumTraceCount;
            
            if sum(use_idx) == 0 %| sum(DiffMeanFluoMat(i, use_idx, j) == 0)
                FrameProfAx.Children(end-(2*(idx-1)+1)).XData = NCTimes{idx}(use_idx);
                FrameProfAx.Children(end-(2*(idx-1)+1)).YData = zeros(1, length(NCTimes{idx}(use_idx)));
                FrameProfAx.Children(end-(2*(idx-1))).XData = NCTimes{idx}(use_idx);
                FrameProfAx.Children(end-(2*(idx-1))).YData = zeros(1, length(NCTimes{idx}(use_idx)));
                FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = zeros(1, length(NCTimes{idx}(use_idx)));
                FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta = zeros(1, length(NCTimes{idx}(use_idx)));
                set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','off'); %'off' or 'on'
                set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','off'); %'off' or 'on'
                set(get(get(FrameProfAx.Children(end-(2*(idx-1)+1)), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
            else
                APBinHasData = true;
                PlottedSets(idx) = true;
                FrameProfAx.Children(end-(2*(idx-1)+1)).YData = MeanFluoMats{idx}(use_idx, i);
                FrameProfAx.Children(end-(2*(idx-1)+1)).XData =NCTimes{idx}(use_idx);
                FrameProfAx.Children(end-(2*(idx-1))).YData = MeanFluoMats{idx}(use_idx, i);
                FrameProfAx.Children(end-(2*(idx-1))).XData = NCTimes{idx}(use_idx);
                FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = StdFluoMats{idx}(use_idx, i);
                FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta  = StdFluoMats{idx}(use_idx, i);
                
                set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','on'); %'off' or 'on'
                if ~SuppressErrorbars
                set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','on'); %'off' or 'on'
                end
                set(get(get(FrameProfAx.Children(end-(2*(idx-1)+1)), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
                
            end
        end
        %try
        if NumSets > 25
            legend(FrameProfAx, legend_labels(PlottedSets),'NumColumns', 2, 'Location', 'eastoutside')
        else
            legend(FrameProfAx, legend_labels(PlottedSets), 'Location', 'eastoutside')
        end
        if exist('PlotTitle', 'var')
            title(FrameProfAx, {PlotTitle,...
                ['Nuclear Cycle ',num2str(NC),', Fraction Embryo Length: ',num2str(APbins(i)) ]})
            
            
        else
            title(FrameProfAx, ['Nuclear Cycle ',num2str(NC),', Fraction Embryo Length: ',num2str(APbins(i)) ])
            
            
        end
        %end
    if sum(PlottedSets) > 0 
        outpath=[outdir4, filesep,'FluoTrace_NC',num2str(NC), '_AP', num2str(i),'.png'];
        saveas(FrameProfFig,outpath);
    end
        
    end
    
end
close all
end