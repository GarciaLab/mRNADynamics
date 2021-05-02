function Plot2LTMSingleTempFluoTraces(project1, project2, outdir, varargin)
%%

% PlotTitle, PlottingColors, UseDifferentColors,
% UseDiffProfiles, UsePhysicalAPLength

NormedCycleTimes = false;



x = 1;
while x <= length(varargin)
    if strcmp(lower(varargin{x}), 'plottitle')
        PlotTitle = varargin{x+1};
        x = x+1;
    elseif strcmp(lower(varargin{x}), 'plottingcolors')
        PlottingColors = varargin{x+1};
        x = x+1;
    elseif strcmpi(varargin{x}, 'NormalizeCycleTimes')
        NormedCycleTimes = true;
    elseif strcmp(lower(varargin{x}), 'tracetype')
        TraceType = lower(varargin{x+1});
        x = x+1;
    end
    x = x+1;
end

if ~exist('PlottingColors', 'var')
    PlottingColors = 'default';
elseif  ~strcmp(lower(PlottingColors), 'default')  & ~strcmp(lower(PlottingColors), 'pboc')
    error('Invalid choice of plotting colors. Can use either "default", "pboc", or "gradient".') % change to error
end
if ~exist('TraceType', 'var')
    TraceType = 'AnaphaseAligned';
elseif  ~strcmp(lower(TraceType), 'tbinned')  & ~strcmp(lower(TraceType), 'anaphasealigned') & ~strcmp(lower(TraceType), 'anaphasealigned3d') &  ~strcmp(lower(TraceType), 'tbinned3d')
    error('Invalid choice of trace type. Can use either "tbinned" or "anaphasealigned".') % change to error
end

if strcmpi(TraceType, 'anaphasealigned')
    TraceType = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'anaphasealigned3d')
    TraceType = 'AnaphaseAligned3D';
elseif strcmpi(TraceType, 'fluo')
    TraceType = 'Unaligned';
elseif strcmpi(TraceType, 'fluo3d')
    TraceType = 'Unaligned3D';
elseif strcmpi(TraceType, 'tbinned')
    TraceType = 'Tbinned';
elseif strcmpi(TraceType, 'tbinned3d')
    TraceType = 'Tbinned3D';
end

%%

if ~exist(outdir, 'dir')
    mkdir(outdir)
end


Temp_obs{1} = project1.Temp_obs;
Temp_sp{1} = project1.Temp_sps;

Temp_obs{2} = project2.Temp_obs;
Temp_sp{2} = project2.Temp_sps;



if strcmp(lower(PlottingColors), 'default')
    [~, colors] = getColorPalettes();
elseif strcmp(lower(PlottingColors), 'pboc')
    [colors, ~] = getColorPalettes();
end

if NormedCycleTimes
    NormString = 'Normalized';
else
    NormString = '';
end



%%

NumSets(1) = length(project1.ExperimentPrefixes);
temperatures{1} = flip(unique(project1.Temp_sps(project1.ProcessedExperiments)));
NumTemperatures(1) = length(temperatures{1});


NumSets(2) = length(project2.ExperimentPrefixes);
temperatures{2} = flip(unique(project2.Temp_sps(project2.ProcessedExperiments)));
NumTemperatures(2)= length(temperatures{2});


APResolution = project1.Experiments{1}.APResolution;
APbins = 0:APResolution:1;
NumAPbins = length(APbins);


legend_labels{1} = project1.LegendLabels;
legend_labels{2} = project2.LegendLabels;
MinimumTraceCount = project1.MinimumTraceCount;


MarkerStyles = {'o', 'd', 's', '>', '^','<', 'h', '*', 'x'};


TempMatches = cell(2, NumTemperatures(1));

UseSet{1} = ismember(1:NumSets(1), project1.ProcessedExperiments);

for t_index = 1:NumTemperatures(1)
    TempMatches{1, t_index} = find((project1.Temp_sps == temperatures{1}(t_index)) & UseSet{1});
end

UseSet{2} = ismember(1:NumSets(2), project2.ProcessedExperiments);

for t_index = 1:NumTemperatures(1)
    TempMatches{2, t_index} = find((project2.Temp_sps == temperatures{1}(t_index)) & UseSet{2});
end
%%

NC = 14;
TemperatureIndex = 2;

CurrentTemperature = temperatures{1}(TemperatureIndex);
TemperatureIndex2 = find( temperatures{2} == CurrentTemperature);

outdir2 = [outdir, filesep, 'T', strrep(num2str(CurrentTemperature), '.', '_'), 'C_NC', num2str(NC)];

if ~exist(outdir2, 'dir')
    mkdir(outdir2)
end

if ~NormedCycleTimes
    outdir3 = [outdir2, filesep, TraceType];
else
    outdir3 = [outdir2, filesep, TraceType, '_NormalizedNCTimes'];
end

if ~exist(outdir3, 'dir')
    mkdir(outdir3)
end

outdir4 = [outdir3, filesep, datestr(now, 'yyyymmdd')];

if ~exist(outdir4, 'dir')
    mkdir(outdir4)
end


% Prepare Traces for plotting
% Prepare Traces for plotting
MaximumNCTimes = NaN(2,max( length(TempMatches{1, TemperatureIndex}), length(TempMatches{2, TemperatureIndex})));
MaxFluos =NaN(2,max( length(TempMatches{1, TemperatureIndex}), length(TempMatches{2, TemperatureIndex})));
MinFluos = NaN(2,max( length(TempMatches{1, TemperatureIndex}), length(TempMatches{2, TemperatureIndex})));
NumTimePoints = NaN(2,max( length(TempMatches{1, TemperatureIndex}), length(TempMatches{2, TemperatureIndex})));
MeanFluoMats = cell(2,max( length(TempMatches{1, TemperatureIndex}), length(TempMatches{2, TemperatureIndex})));
StdFluoMats = cell(2,max( length(TempMatches{1, TemperatureIndex}), length(TempMatches{2, TemperatureIndex})));
NumNucMats = cell(2,max( length(TempMatches{1, TemperatureIndex}), length(TempMatches{2, TemperatureIndex})));
NCTimes =cell(2,max( length(TempMatches{1, TemperatureIndex}), length(TempMatches{2, TemperatureIndex})));
MinAPs = NaN(2,max( length(TempMatches{1, TemperatureIndex}), length(TempMatches{2, TemperatureIndex})));
MaxAPs = NaN(2,max( length(TempMatches{1, TemperatureIndex}), length(TempMatches{2, TemperatureIndex})));

for idx2=1:length(TempMatches{1, TemperatureIndex})
    idx = TempMatches{1, TemperatureIndex}(idx2);
    if ~ismember(idx, project1.ProcessedExperiments)
        continue
    end
    ExpNCTimes = project1.MeanProfiles{idx}.([TraceType, 'CycleFrameTimes']){NC-8};
    IncludedRows = 1:length(ExpNCTimes);
    ExpFluoMat = squeeze(project1.MeanProfiles{idx}.([TraceType, 'CycleMeanTraces'])(IncludedRows,:,NC-8));
    ExpStdMat = squeeze(project1.MeanProfiles{idx}.([TraceType, 'CycleTraceStdErrors'])(IncludedRows,:,NC-8));
    ExpNumNucMat = squeeze(project1.MeanProfiles{idx}.([TraceType, 'CycleNumOnNuclei'])(IncludedRows,:,NC-8));
    NCLength = project1.EmbryoStats.NCDivisionInfo(idx,NC-8);
    
    IncludedColumns = find(sum(~isnan(ExpFluoMat),1).' > 0);
    if ~isempty(IncludedColumns)
        MinAPs(1, idx2) = min(IncludedColumns);
        MaxAPs(1, idx2) = max(IncludedColumns);
    end
    
    IncludedRows = find(sum(~isnan(ExpFluoMat),2).' > 0);
    if isempty(IncludedRows)
        MeanFluoMats{1, idx2} = [];
        StdFluoMats{1, idx2} = [];
        NumNucMats{1, idx2} = [];
        NCTimes{1, idx2} = [];
    else
        MeanFluoMats{1, idx2} = ExpFluoMat(IncludedRows,:);
        StdFluoMats{1, idx2} = ExpStdMat(IncludedRows,:);
        NumNucMats{1, idx2} = ExpNumNucMat(IncludedRows,:);
        
        if ~NormedCycleTimes
            NCTimes{1, idx2} = ExpNCTimes(IncludedRows)/60;
        else
            NCTimes{1, idx2} = (ExpNCTimes(IncludedRows)/60)/NCLength;
        end
        
        MaximumNCTimes(1, idx2) = max(NCTimes{1, idx2});
        MaxFluos(1, idx2) = max(max(MeanFluoMats{1, idx2}+StdFluoMats{1, idx2}));
        MinFluos(1, idx2) = min(min(MeanFluoMats{1, idx2}+StdFluoMats{1, idx2}));
        NumTimePoints(1, idx2) = uint16(max(NCTimes{1, idx2})/(project1.time_delta/60)+1);
    end
    
    
end

for idx2=1:length(TempMatches{2, TemperatureIndex})
    idx = TempMatches{2, TemperatureIndex}(idx2);
    if ~ismember(idx, project2.ProcessedExperiments)
        continue
    end
    ExpNCTimes = project2.MeanProfiles{idx}.([TraceType, 'CycleFrameTimes']){NC-8};
    IncludedRows = 1:length(ExpNCTimes);
    ExpFluoMat = squeeze(project2.MeanProfiles{idx}.([TraceType, 'CycleMeanTraces'])(IncludedRows,:,NC-8));
    ExpStdMat = squeeze(project2.MeanProfiles{idx}.([TraceType, 'CycleTraceStdErrors'])(IncludedRows,:,NC-8));
    ExpNumNucMat = squeeze(project2.MeanProfiles{idx}.([TraceType, 'CycleNumOnNuclei'])(IncludedRows,:,NC-8));
    NCLength = project2.EmbryoStats.NCDivisionInfo(idx,NC-8);
    
    IncludedColumns = find(sum(~isnan(ExpFluoMat),1).' > 0);
    if ~isempty(IncludedColumns)
        MinAPs(2, idx2) = min(IncludedColumns);
        MaxAPs(2, idx2) = max(IncludedColumns);
    end
    
    IncludedRows = find(sum(~isnan(ExpFluoMat),2).' > 0);
    if isempty(IncludedRows)
        MeanFluoMats{2, idx2} = [];
        StdFluoMats{2, idx2} = [];
        NumNucMats{2, idx2} = [];
        NCTimes{2, idx2} = [];
    else
        MeanFluoMats{2, idx2} = ExpFluoMat(IncludedRows,:);
        StdFluoMats{2, idx2} = ExpStdMat(IncludedRows,:);
        NumNucMats{2, idx2} = ExpNumNucMat(IncludedRows,:);
        
        if ~NormedCycleTimes
            NCTimes{2, idx2} = ExpNCTimes(IncludedRows)/60;
        else
            NCTimes{2, idx2} = (ExpNCTimes(IncludedRows)/60)/NCLength;
        end
        
        MaximumNCTimes(2, idx2) = max(NCTimes{2, idx2});
        MaxFluos(2, idx2) = max(max(MeanFluoMats{2, idx2}+StdFluoMats{2, idx2}));
        MinFluos(2, idx2) = min(min(MeanFluoMats{2, idx2}+StdFluoMats{2, idx2}));
        NumTimePoints(2, idx2) = uint16(max(NCTimes{2, idx2})/(project2.time_delta/60)+1);
    end
    
    
end


MinAPbin = nanmin(nanmin(MinAPs));
MaxAPbin = nanmax(nanmax(MaxAPs));



close all


FrameProfFig = figure(1);
set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .8, .6]);
set(gcf,'color','w');
FrameProfAx = axes(FrameProfFig);
eb = cell(2,max( length(TempMatches{1, TemperatureIndex}), length(TempMatches{2, TemperatureIndex})));
prof = cell(2,max( length(TempMatches{1, TemperatureIndex}), length(TempMatches{2, TemperatureIndex})));

yyaxis(FrameProfAx, 'left');
for idx =1:length(TempMatches{1, TemperatureIndex})
    
    
    temp_idx = 1;
    

    
    marker_idx = idx;
    
    
    
    
    
    eb{1, idx} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
    hold on
    
    set(eb{1, idx}, 'color', colors(temp_idx,:), 'LineWidth', 1);
    
    set(get(get(eb{1, idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    
    prof{1, idx} = plot(APbins, ones(1, length(APbins)), MarkerStyles{marker_idx},...
        'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colors(temp_idx,:),...
        'MarkerSize', 4, 'LineStyle', '-', 'Color', colors(temp_idx,:));
    
    
    set(eb{1, idx},'Visible','off'); %'off' or 'on'
    set(prof{1, idx},'Visible','off'); %'off' or 'on'
    set(get(get(prof{1, idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
end



if NormedCycleTimes
    xlabel('Fraction of Nuclear Cycle')
elseif strcmpi(TraceType, 'anaphasealigned') | strcmp(lower(TraceType), 'anaphasealigned3d')
    xlabel('Time since anaphase (min)')
else
    xlabel('Time since NC start (min)')
end


if ~NormedCycleTimes
    xlim([0-0.01* max(max(MaximumNCTimes)), max(max(MaximumNCTimes))*1.1]);
else
    xlim([-0.01, 1.01]);
end

if strcmp(lower(TraceType), 'tbinned') | strcmp(lower(TraceType), 'anaphasealigned') | strcmpi(TraceType, 'unaligned')
    ylabel('Fluo (AU)')
else
    ylabel('3D Fluo (AU)')
end

ylim([0, max(max(MaxFluos(1,:))*1.1,1)])

%%


yyaxis(FrameProfAx, 'right');
for idx =1:length(TempMatches{2, TemperatureIndex})
    
    
    temp_idx = 2;
    

    
    marker_idx = idx;
    
    
    
    
    
    eb{2, idx} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
    hold on
    
    set(eb{2, idx}, 'color', colors(temp_idx,:), 'LineWidth', 1);
    
    set(get(get(eb{2, idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    
    prof{2, idx} = plot(APbins, ones(1, length(APbins)), MarkerStyles{marker_idx},...
        'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colors(temp_idx,:),...
        'MarkerSize', 4, 'LineStyle', '-', 'Color', colors(temp_idx,:));
    
    
    set(eb{2, idx},'Visible','off'); %'off' or 'on'
    set(prof{2, idx},'Visible','off'); %'off' or 'on'
    set(get(get(prof{2, idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
end





if strcmp(lower(TraceType), 'tbinned') | strcmp(lower(TraceType), 'anaphasealigned') | strcmpi(TraceType, 'unaligned')
    ylabel('Fluo (AU)')
else
    ylabel('3D Fluo (AU)')
end


ylim([0, max(max(MaxFluos(2,:))*1.1,1)])
%%

%

title(FrameProfAx, {'',...
    ['T = ', num2str(temperatures{1}(TemperatureIndex)), '°C, Nuclear Cycle ',num2str(NC),', Fraction Embryo Length: ',num2str(-1) ]})



for i = MinAPbin:MaxAPbin
    PlottedSets = zeros(2, max(NumSets(1), NumSets(2)), 'logical');
    yyaxis(FrameProfAx, 'left');
    for idx=1:length(TempMatches{1, TemperatureIndex})
        set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','off'); %'off' or 'on'
        set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','off'); %'off' or 'on'
        set(get(get(FrameProfAx.Children(end-(2*(idx-1)+1)), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        
        if ~ismember(TempMatches{1, TemperatureIndex}(idx), project1.ProcessedExperiments)
            continue
        end
        if isempty(NCTimes{1, idx})
            continue
        end
        
        use_idx = NumNucMats{1, idx}(:,i) >= MinimumTraceCount;
        
        if sum(use_idx) == 0 %| sum(DiffMeanFluoMat(i, use_idx, j) == 0)
            FrameProfAx.Children(end-(2*(idx-1)+1)).XData = NCTimes{1, idx}(use_idx);
            FrameProfAx.Children(end-(2*(idx-1)+1)).YData = zeros(1, length(NCTimes{1, idx}(use_idx)));
            FrameProfAx.Children(end-(2*(idx-1))).XData = NCTimes{1, idx}(use_idx);
            FrameProfAx.Children(end-(2*(idx-1))).YData = zeros(1, length(NCTimes{1, idx}(use_idx)));
            FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = zeros(1, length(NCTimes{1, idx}(use_idx)));
            FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta = zeros(1, length(NCTimes{1, idx}(use_idx)));
            set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','off'); %'off' or 'on'
            set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','off'); %'off' or 'on'
            set(get(get(FrameProfAx.Children(end-(2*(idx-1)+1)), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        else
            PlottedSets(1, TempMatches{1, TemperatureIndex}(idx)) = true;
            FrameProfAx.Children(end-(2*(idx-1)+1)).YData = MeanFluoMats{1, idx}(use_idx, i);
            FrameProfAx.Children(end-(2*(idx-1)+1)).XData =NCTimes{1, idx}(use_idx);
            FrameProfAx.Children(end-(2*(idx-1))).YData = MeanFluoMats{1, idx}(use_idx, i);
            FrameProfAx.Children(end-(2*(idx-1))).XData = NCTimes{1, idx}(use_idx);
            FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = StdFluoMats{1, idx}(use_idx, i);
            FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta  = StdFluoMats{1, idx}(use_idx, i);
            
            set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','on'); %'off' or 'on'
            set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','on'); %'off' or 'on'
            set(get(get(FrameProfAx.Children(end-(2*(idx-1)+1)), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
            
        end
    end
    
    yyaxis(FrameProfAx, 'right');
    for idx=1:length(TempMatches{2, TemperatureIndex})
        set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','off'); %'off' or 'on'
        set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','off'); %'off' or 'on'
        set(get(get(FrameProfAx.Children(end-(2*(idx-1)+1)), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        
        if ~ismember(TempMatches{2, TemperatureIndex}(idx), project2.ProcessedExperiments)
            continue
        end
        if isempty(NCTimes{2, idx})
            continue
        end
        
        use_idx = NumNucMats{2, idx}(:,i) >= MinimumTraceCount;
        
        if sum(use_idx) == 0 %| sum(DiffMeanFluoMat(i, use_idx, j) == 0)
            FrameProfAx.Children(end-(2*(idx-1)+1)).XData = NCTimes{2, idx}(use_idx);
            FrameProfAx.Children(end-(2*(idx-1)+1)).YData = zeros(1, length(NCTimes{2, idx}(use_idx)));
            FrameProfAx.Children(end-(2*(idx-1))).XData = NCTimes{2, idx}(use_idx);
            FrameProfAx.Children(end-(2*(idx-1))).YData = zeros(1, length(NCTimes{2, idx}(use_idx)));
            FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = zeros(1, length(NCTimes{2, idx}(use_idx)));
            FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta = zeros(1, length(NCTimes{2, idx}(use_idx)));
            set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','off'); %'off' or 'on'
            set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','off'); %'off' or 'on'
            set(get(get(FrameProfAx.Children(end-(2*(idx-1)+1)), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        else
            PlottedSets(2, TempMatches{2, TemperatureIndex}(idx)) = true;
            FrameProfAx.Children(end-(2*(idx-1)+1)).YData = MeanFluoMats{2, idx}(use_idx, i);
            FrameProfAx.Children(end-(2*(idx-1)+1)).XData =NCTimes{2, idx}(use_idx);
            FrameProfAx.Children(end-(2*(idx-1))).YData = MeanFluoMats{2, idx}(use_idx, i);
            FrameProfAx.Children(end-(2*(idx-1))).XData = NCTimes{2, idx}(use_idx);
            FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = StdFluoMats{2, idx}(use_idx, i);
            FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta  = StdFluoMats{2, idx}(use_idx, i);
            
            set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','on'); %'off' or 'on'
            set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','on'); %'off' or 'on'
            set(get(get(FrameProfAx.Children(end-(2*(idx-1)+1)), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
            
        end
    end
    %try
    if ~sum(sum(PlottedSets))
        continue
    end
    legend(FrameProfAx, [legend_labels{1}(PlottedSets(1,1:length(legend_labels{1}))) legend_labels{2}(PlottedSets(2,1:length(legend_labels{2})))], 'Location', 'eastoutside')
    if exist('PlotTitle', 'var')
        title(FrameProfAx, {PlotTitle,...
            ['T = ', num2str(temperatures(TemperatureIndex)), '°C, Nuclear Cycle ',num2str(NC),', Fraction Embryo Length: ',num2str(APbins(i)) ]})
    else
        title(FrameProfAx, ['T = ', num2str(temperatures{1}(TemperatureIndex)), '°C, Nuclear Cycle ',num2str(NC),', Fraction Embryo Length: ',num2str(APbins(i)) ])
        
    end
    
    saveas(FrameProfFig,[outdir4, filesep, NormString,...
        'FluoTrace_NC',num2str(NC),'_T',strrep(num2str(CurrentTemperature), '.', '_'),'_Bin', num2str(i),'.png']);
    
end

close all
end