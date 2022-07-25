function PlotLTMTBinnedFluoProfiles(this, outdir, varargin)
%%

% PlotTitle, PlottingColors, UseDifferentColors,
% UseDiffProfiles, UsePhysicalAPLength

useRescaledTiming = false;
useRescaledFluo = false;
SuppressErrorbars = false; 
SuppressMarkers = false; 
DownsampleTraces = false;
DownsamplingRate = 1;
UsePhysicalAPLength = false;
UsePerNucleusTraces = false;
UseFractionOns = false;
AltFractionOns = false;
UseThesisStyle = false;

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
    elseif strcmpi(varargin{x}, 'SuppressErrorbars')
        SuppressErrorbars = true;
    elseif strcmpi(varargin{x}, 'SuppressMarkers')
        SuppressMarkers = true;
    elseif strcmpi(varargin{x}, 'UsePerNucleusTraces')
        UsePerNucleusTraces = true;
    elseif strcmpi(varargin{x}, 'DownsampleTraces')
        DownsampleTraces = true;
        DownsamplingRate = 2;
    elseif strcmpi(varargin{x}, 'DownsamplingRate')
        DownsampleTraces = true;
        DownsamplingRate = varargin{x+1};
        x = x+1;
    elseif strcmpi(varargin{x}, 'UseFractionOns')
        UseFractionOns = true;
    elseif strcmpi(varargin{x}, 'UseThesisStyle')
        UseThesisStyle = true;
        PlottingColors = 'thesis';
    elseif strcmpi(varargin{x}, 'UseAltFractionOns')
        UseFractionOns = true;
        AltFractionOns = true;
    elseif strcmp(lower(varargin{x}), 'usephysicalaplength')
        UsePhysicalAPLength = true;
    elseif strcmpi(varargin{x}, 'userescaledtime') | strcmpi(varargin{x}, 'rescaletime') | ...
            strcmpi(varargin{x}, 'rescaletiming') | strcmpi(varargin{x}, 'userescaledtiming')
        useRescaledTiming = true;
    elseif strcmpi(varargin{x}, 'rescalefluo') | strcmpi(varargin{x}, 'userescaledfluo')
        useRescaledFluo = true;
    else
        error('Invalid option passed to function.')
    end
    x = x+1;
end

if ~exist('PlottingColors', 'var')
    PlottingColors = 'default';
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

if UsePhysicalAPLength
    PhysicalAPString = 'PhysAP';
else
    PhysicalAPString = '';
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
end

if ~isempty(PhysicalAPString)
    if ~isempty(SpecificDirString)
        SpecificDirString = [SpecificDirString, '_', PhysicalAPString];
    else
        SpecificDirString = PhysicalAPString;
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
NumAPbins = length(APbins);

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


%%

if ~UseFractionOns
    outdir2 = [outdir,filesep,'TemperatureBinnedFluoProfiles'];
elseif ~AltFractionOns
    outdir2 = [outdir,filesep,'TemperatureBinnedFractionOnProfiles'];
else
    outdir2 = [outdir,filesep,'TemperatureBinnedAltFractionOnProfiles'];
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
for nc_idx=1:length(this.IncludedNCs)
    NC = this.IncludedNCs(nc_idx);
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
        
        if ~useRescaledTiming
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
                
                NCTimes{idx} = ExpNCTimes(IncludedRows);
                
                
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
        
        else
        [MeanFluoMats, StdFluoMats, NumEmbryoMats, NCTimes, MaximumNCTimes,...
            MaxFluos,MinFluos, NumFrames, MinAPs, MaxAPs] = getLTMRescaledTBinnedTimingMats(this, traceName, NC,...
            useRescaledFluo, DownsamplingRate, UseFractionOns, AltFractionOns);
        end
           
        

    end
    SupersetTimes = NCTimes{1};
    for idx = 2:NumTemperatures
        SupersetTimes = unique([SupersetTimes NCTimes{idx}]);
    end
    
    if all(isnan(NumFrames))
        continue
    end
    
    MinAPbin = 1;
    MaxAPbin = NumAPbins;

    if ~UsePhysicalAPLength
        MinAPLimit = max((MinAPbin-3), 0)*APResolution;
        MaxAPLimit = min((MaxAPbin+1), NumAPbins)*APResolution;
    else
        MinAPLimit = max((MinAPbin-3), 0)*APResolution*min(this.APLengths);
        MaxAPLimit = min((MaxAPbin+1), NumAPbins)*APResolution*max(this.APLengths);
    end
    
    close all
    
    FrameProfFig = figure(1);
    set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    FrameProfAx = axes(FrameProfFig);
    
    eb = cell(1, NumTemperatures);
    prof = cell(1, NumTemperatures);
    
    
    
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
                'MarkerSize', 10, 'LineStyle', '-', 'Color', colors(temp_idx,:),...
                'LineWidth', 1.5);
        elseif ~SuppressMarkers
            prof{idx} = plot(APbins, ones(1, length(APbins)), MarkerStyles{marker_idx},...
                'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colors(temp_idx,:),...
                'MarkerSize', 8, 'LineStyle', '-', 'Color', colors(temp_idx,:),...
                'LineWidth', 1.5);
        else
            prof{idx} = plot(APbins, ones(1, length(APbins)), '-',...
                'Color', colors(temp_idx,:),'LineWidth', 3);
        end

 
      
        set(prof{idx},'Visible','off'); %'off' or 'on'
        set(get(get(prof{idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    end
    
    grid on 
    hold off
    
    if ~UsePhysicalAPLength
        xlabel('Fraction Embryo Length', 'FontSize', 16)
        xlim([MinAPLimit, MaxAPLimit])
    else
        xlabel('Distance from the Anterior Pole (\mum)', 'FontSize', 16)
        xlim([MinAPLimit, MaxAPLimit])
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
    
    FrameProfAx.XAxis.FontSize = 18;
    FrameProfAx.YAxis.FontSize = 18;
    if strcmp(lower(TraceType), 'anaphasealigned') | strcmp(lower(TraceType), 'anaphasealigned3d')
        title(FrameProfAx, {'',...
            ['Nuclear Cycle ',num2str(NC),', ',num2str(-1), ' min since anaphase' ]}, 'FontSize', 16)
    else
        title(FrameProfAx, {'',...
            ['Nuclear Cycle ',num2str(NC),', ',num2str(-1), ' min' ]}, 'FontSize', 16)
    end
    
    for i = 1:length(SupersetTimes)
        time_value = SupersetTimes(i);

        PlottedSets = zeros(1, NumTemperatures, 'logical');
        for idx=1:NumTemperatures

            set(get(get(FrameProfAx.Children(end-NumTemperatures-idx+1), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            set(FrameProfAx.Children(end-NumTemperatures-idx+1),'Visible','off'); %'off' or 'on'
            set(FrameProfAx.Children(end-idx+1),'Visible','off'); %'off' or 'on'    
            if ismember(SupersetTimes(i), NCTimes{idx})
                match_idx = find(NCTimes{idx} == SupersetTimes(i), 1);
            else
                continue
            end

            if match_idx > size(NumEmbryoMats{idx}, 1)
                continue
            end
            if UsePhysicalAPLength
                APLength = this.APLengths(idx);
            else
                APLength = 1;
            end
            use_idx = NumEmbryoMats{idx}(match_idx,:) >= MinimumEmbryoCount;
            
            if sum(use_idx) == 0 %| sum(DiffMeanFluoMat(i, use_idx, j) == 0)
                FrameProfAx.Children(end-NumTemperatures-idx+1).XData = APLength*APbins;
                FrameProfAx.Children(end-NumTemperatures-idx+1).YData = zeros(1, length(APbins));
                FrameProfAx.Children(end-idx+1).XData = APLength*APbins;
                FrameProfAx.Children(end-idx+1).YData = zeros(1, length(APbins));
                FrameProfAx.Children(end-idx+1).YPositiveDelta = .1*ones(1, length(APbins));
                FrameProfAx.Children(end-idx+1).YNegativeDelta = .1*ones(1, length(APbins));
                set(FrameProfAx.Children(end-NumTemperatures-idx+1),'Visible','off'); %'off' or 'on'
                set(FrameProfAx.Children(end-idx+1),'Visible','off'); %'off' or 'on'
                set(get(get(FrameProfAx.Children(end-NumTemperatures-idx+1), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
            else
                APBinHasData = true;
                PlottedSets(idx) = true;
                FrameProfAx.Children(end-NumTemperatures-idx+1).YData = MeanFluoMats{idx}(match_idx, use_idx);
                FrameProfAx.Children(end-NumTemperatures-idx+1).XData =APLength*APbins(use_idx);
                FrameProfAx.Children(end-idx+1).YData = MeanFluoMats{idx}(match_idx, use_idx);
                FrameProfAx.Children(end-idx+1).XData = APLength*APbins(use_idx);
                FrameProfAx.Children(end-idx+1).YPositiveDelta = StdFluoMats{idx}(match_idx, use_idx);
                FrameProfAx.Children(end-idx+1).YNegativeDelta  = StdFluoMats{idx}(match_idx, use_idx);
                
                set(FrameProfAx.Children(end-NumTemperatures-idx+1),'Visible','on'); %'off' or 'on'
                if ~SuppressErrorbars & ~SuppressMarkers
                    set(FrameProfAx.Children(end-idx+1),'Visible','on'); %'off' or 'on'
                end
                set(get(get(FrameProfAx.Children(end-NumTemperatures-idx+1), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
                
            end
        end
        %try
        legend(FrameProfAx, legend_labels(PlottedSets), 'Location', 'northeast', 'FontSize', 18)
         if ~useRescaledTiming
            if exist('PlotTitle', 'var')
                
                if strcmp(lower(TraceType), 'anaphasealigned') | strcmp(lower(TraceType), 'anaphasealigned3d')
                    title(FrameProfAx, {PlotTitle,...
                        ['Nuclear Cycle ',num2str(NC),', ',num2str(SupersetTimes(i)), ' min. since anaphase' ]}, 'FontSize', 16)
                else
                    title(FrameProfAx, {PlotTitle,...
                        ['Nuclear Cycle ',num2str(NC),', ',num2str(SupersetTimes(i)), ' min.' ]}, 'FontSize', 16)
                end
            else
                if strcmp(lower(TraceType), 'anaphasealigned') | strcmp(lower(TraceType), 'anaphasealigned3d')
                    title(FrameProfAx,  ['Nuclear Cycle ',num2str(NC),', ',num2str(SupersetTimes(i)), ' min. since anaphase' ], 'FontSize', 16)
                else
                    title(FrameProfAx,  ['Nuclear Cycle ',num2str(NC),', ',num2str(SupersetTimes(i)), ' min.' ], 'FontSize', 16)
                    
                end
            end
        else
            if exist('PlotTitle', 'var')
                
                if strcmp(lower(TraceType), 'anaphasealigned') | strcmp(lower(TraceType), 'anaphasealigned3d')
                    title(FrameProfAx, {PlotTitle,...
                        ['Nuclear Cycle ',num2str(NC),', ',num2str(SupersetTimes(i)), ' temperature-scaled min. since anaphase' ]}, 'FontSize', 16)
                else
                    title(FrameProfAx, {PlotTitle,...
                        ['Nuclear Cycle ',num2str(NC),', ',num2str(SupersetTimes(i)), ' temperature-scaled min.' ]}, 'FontSize', 16)
                end
            else
                if strcmp(lower(TraceType), 'anaphasealigned') | strcmp(lower(TraceType), 'anaphasealigned3d')
                    title(FrameProfAx,  ['Nuclear Cycle ',num2str(NC),', ',num2str(SupersetTimes(i)), ' temperature-scaled min. since anaphase' ], 'FontSize', 16)
                else
                    title(FrameProfAx,  ['Nuclear Cycle ',num2str(NC),', ',num2str(SupersetTimes(i)), ' temperature-scaled min.' ], 'FontSize', 16)
                    
                end
            end
            
            
         end
         %end
         
         if sum(PlottedSets) > 0
             if ~UseFractionOns
                 outpath=[outdir4, filesep, 'FluoProfile_NC',num2str(NC), '_AP', num2str(i),'.png'];
             else
                 outpath=[outdir4, filesep, 'FractionOnProfile_NC',num2str(NC), '_AP', num2str(i),'.png'];
             end
             saveas(FrameProfFig,outpath);
         end
    end
    
end
close all
end