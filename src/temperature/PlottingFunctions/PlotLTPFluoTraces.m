function PlotLTPFluoTraces(this, outdir, varargin)
%%

% PlotTitle, PlottingColors, UseDifferentColors,
% UseDiffProfiles, UsePhysicalAPLength



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
elseif ~strcmp(lower(TraceType), 'unaligned') &  ~strcmp(lower(TraceType), 'anaphasealigned')& ~strcmp(lower(TraceType), 'tbinned')
    error('Invalid choice of trace type. Can use either "unaligned", "anaphasealigned", or "tbinned".') % change to error
end

if ~exist(outdir, 'dir')
    mkdir(outdir)
end
%%
if strcmpi(TraceType, 'anaphasealigned')
    traceName = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'unaligned')
    traceName = 'Unaligned';
elseif strcmpi(TraceType, 'tbinned')
    traceName = 'Tbinned';
end

%%


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



%%

NumSets = length(this.ExperimentPrefixes);
temperatures = flip(unique(this.Temp_sps(this.ProcessedExperiments)));
NumTemperatures = length(temperatures);


APResolution = this.Experiments{1}.APResolution;

APbins = 0:APResolution:1;

legend_labels = this.LegendLabels;
MinimumTraceCount = this.MinimumNuclearCount;


MarkerStyles = {'o', 'd', 's', '>', '^','p', 'h', '*', 'x'};

UseSet = ismember(1:NumSets, this.ProcessedExperiments);
TempMatches = cell(1, NumTemperatures);
LegendOrder = zeros(1, NumSets);

counter = 0;
for t_index = 1:NumTemperatures
    TempMatches{t_index} = find((this.Temp_sps == temperatures(t_index)) & UseSet);
    TempMatches2 =  find(this.Temp_sps == temperatures(t_index));
    LegendOrder(counter + 1:counter+length(TempMatches2)) = TempMatches2;
    counter = counter + length(TempMatches);
end



for NC=10:14

    
    
    
    
    outdir3 = [outdir, filesep, datestr(now, 'yyyymmdd')];
    if ~exist(outdir3, 'dir')
        mkdir(outdir3)
    end
    
    % Prepare Traces for plotting
    MaximumNCTimes = NaN(1, NumSets);
    MaxFluos = NaN(1, NumSets);
    NumFrames = NaN(1, NumSets);
    MeanFluoMats = cell(1, NumSets);
    StdFluoMats = cell(1, NumSets);
    NumNucMats = cell(1, NumSets);
    NCTimes = cell(1, NumSets);
    for idx=1:NumSets
        if ~ismember(idx, this.ProcessedExperiments)
            continue
        end
        ExpNCTimes = this.TFProfiles{idx}.([traceName, 'CycleFrameTimes']){NC-8};
        IncludedRows = 1:length(ExpNCTimes);
        ExpFluoMat = squeeze(this.TFProfiles{idx}.([traceName, 'CycleMeanTraces'])(IncludedRows,:,NC-8));
        ExpStdMat = squeeze(this.TFProfiles{idx}.([traceName, 'CycleTraceStdErrors'])(IncludedRows,:,NC-8));
        ExpNumNucMat = squeeze(this.TFProfiles{idx}.([traceName, 'CycleNumOnNuclei'])(IncludedRows,:,NC-8));
       
        
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
            
            NCTimes{idx} = ExpNCTimes(IncludedRows)/60;
           
            MaximumNCTimes(idx) = max(NCTimes{idx});
            MaxFluos(idx) = max(max(MeanFluoMats{idx}+StdFluoMats{idx}));
            NumFrames(idx) = length(NCTimes{idx});
        end
    end
    
    if all(isnan(NumFrames))
        continue
    end
    close all
    
    
    FrameProfFig = figure(1);
    set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    FrameProfAx = axes(FrameProfFig);
    for idx =1:NumSets
        
        if strcmp(lower(PlottingColors), 'gradient')
            temp_idx = find(abs(Temp_range - Temp_obs(idx)) == min(abs(Temp_range - Temp_obs(idx))), 1);
            temp_idx2 =  find(temperatures == Temp_sp(idx));
            marker_idx = find(TempMatches{temp_idx2} == idx);
           
        else
            temp_idx = find(temperatures == Temp_sp(idx));
            marker_idx = find(TempMatches{temp_idx} == idx);
        end
        
        if isempty(marker_idx)
            marker_idx = 1;
        end
        
        

        eb{idx} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
        hold on
   
        set(eb{idx}, 'color', colors(temp_idx,:), 'LineWidth', 1);
     
        set(get(get(eb{idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    
        prof{idx} = plot(APbins, ones(1, length(APbins)), MarkerStyles{marker_idx},...
        'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colors(temp_idx,:),...
            'MarkerSize', 8, 'LineStyle', '-', 'Color', colors(temp_idx,:));
  
 
        set(eb{idx},'Visible','off'); %'off' or 'on'
        set(prof{idx},'Visible','off'); %'off' or 'on'
        set(get(get(prof{idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    end
    
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
    if strcmp(lower(TraceType), 'anaphasealigned')
        xlabel('Time since anaphase (min)')
    else
        xlabel('Time since NC start (min)')
    end
    xlim([0, nanmax(MaximumNCTimes)])
    

    ylabel('Fluo (AU)')
    
    ylim([0, nanmax(nanmax(MaxFluos)*1.1,1)])
    %
    
    title(FrameProfAx, {'',...
        ['Nuclear Cycle ',num2str(NC),', Fraction Embryo Length: ',num2str(-1) ]})
    
    
    
    

    for i = 1:length(APbins)
        APBinHasData = false;
        PlottedSets = zeros(1, NumSets, 'logical');
        for idx=1:NumSets
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
                set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','on'); %'off' or 'on'
                 set(get(get(FrameProfAx.Children(end-(2*(idx-1)+1)), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
                
            end
        end
        %try
        legend(FrameProfAx, legend_labels(PlottedSets), 'Location', 'eastoutside')
        if exist('PlotTitle', 'var')
            title(FrameProfAx, {PlotTitle,...
                ['Nuclear Cycle ',num2str(NC),', Fraction Embryo Length: ',num2str(APbins(i)) ]})
            
            
        else
            title(FrameProfAx, ['Nuclear Cycle ',num2str(NC),', Fraction Embryo Length: ',num2str(APbins(i)) ])
            
            
        end
        %end
        if ~APBinHasData
            continue
        end
        
        saveas(FrameProfFig,[outdir3, filesep,...
            'NC',num2str(NC), GradString, 'FluoTrace',TraceType,'_', num2str(i),'.png']);
        
    end
    
end
close all
end