function PlotLTPSingleSetFluoTraces(this, varargin)
%%

% PlotTitle, PlottingColors, UseDifferentColors,
% UseDiffProfiles, UsePhysicalAPLength

IncludedSets = this.ProcessedExperiments;
UseOffsets = false;

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
    elseif strcmpi(varargin{x}, 'subtractoffsets') | strcmpi(varargin{x}, 'useoffsets')
        UseOffsets = true;
    elseif strcmp(lower(varargin{x}), 'includedsets')
        IncludedSets = lower(varargin{x+1});
        x = x+1;
    end
    x = x+1;
end

if ~exist('PlottingColors', 'var')
    PlottingColors = 'default';
elseif ~strcmp(lower(PlottingColors), 'default')  & ~strcmp(lower(PlottingColors), 'pboc')
    error('Invalid choice of plotting colors. Can use either "default", "pboc", or "gradient".') % change to error
end
if ~exist('TraceType', 'var')
    TraceType = 'AnaphaseAligned';
elseif  ~strcmp(lower(TraceType), 'tbinned')  & ~strcmp(lower(TraceType), 'anaphasealigned')
    error('Invalid choice of trace type. Can use either "unaligned", "tbinned" or "anaphasealigned".') % change to error
end

if strcmpi(TraceType, 'anaphasealigned')
    TraceType = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'tbinned')
    TraceType = 'Tbinned';
end

%%


Temp_obs = this.Temp_obs;
Temp_sp = this.Temp_sps;


if strcmp(lower(PlottingColors), 'default')
    [~, colors] = getColorPalettes();
elseif strcmp(lower(PlottingColors), 'pboc')
    [colors, ~] = getColorPalettes();
end


if UseOffsets
    OffsetString = 'Offset';
else
    OffsetString = '';
end
%%
timeSubDir = 'MeanTraces';
Nsigfigs = 3;
%%
NumSets = length(this.ExperimentPrefixes);
temperatures = flip(unique(this.Temp_sps(this.ProcessedExperiments)));
NumTemperatures = length(temperatures);


APResolution = this.Experiments{1}.APResolution;
APbins = 0:APResolution:1;
NumAPbins = length(APbins);


legend_labels = this.LegendLabels;
MinimumTraceCount = this.MinimumNuclearCount;

for SetIndex=1:NumSets%length(temperatures)
    if ~ismember(SetIndex, this.ProcessedExperiments)
        continue
    end
    resultsFolder = this.Experiments{SetIndex}.resultsFolder;
    
    outdir2 = [resultsFolder, filesep, timeSubDir];
    
    if ~exist(outdir2, 'dir')
        mkdir(outdir2)
    end
    
    
    
    
    CurrentTemperature = this.Temp_sps(SetIndex);
    
    for NC = this.IncludedNCs
        
        % Prepare Traces for plotting
        
        
        
        NCTimes = this.TFProfiles{SetIndex}.([TraceType, 'CycleFrameTimes']){NC-8};
        IncludedRows = 1:length(NCTimes);
        MeanFluoMat = squeeze(this.TFProfiles{SetIndex}.([TraceType, 'CycleMeanTraces'])(IncludedRows,:,NC-8));
        StdFluoMat = squeeze(this.TFProfiles{SetIndex}.([TraceType, 'CycleTraceStdErrors'])(IncludedRows,:,NC-8));
        NumNucMat = squeeze(this.TFProfiles{SetIndex}.([TraceType, 'CycleNumOnNuclei'])(IncludedRows,:,NC-8));
        if UseOffsets
            IncludedRows2 = uint16(NCTimes/this.time_delta+1);
            
            ExpOffsets = squeeze(this.SetFluoOffsets.(TraceType).mean(SetIndex, IncludedRows2,NC-8)).';
            ExpOffsetStd = squeeze(this.SetFluoOffsets.(TraceType).std(SetIndex, IncludedRows2,NC-8)).';
            ExpOffsetCount = squeeze(this.SetFluoOffsets.(TraceType).count(SetIndex, IncludedRows2,NC-8)).';
        end
        
        IncludedRows = find(sum(~isnan(MeanFluoMat),2).' > 0);
        if isempty(IncludedRows)
            continue
        end
        
        IncludedColumns = find(sum(~isnan(MeanFluoMat),1).' > 0);
        if ~isempty(IncludedColumns)
            MinAPbin = min(IncludedColumns);
            MaxAPbin = max(IncludedColumns);
        end
        
        if ~UseOffsets
            MeanFluoMat = MeanFluoMat(IncludedRows,:);
            StdFluoMat = StdFluoMat(IncludedRows,:);
        else
            StdErrorOffsetVector = ExpOffsetStd(IncludedRows)./sqrt(ExpOffsetCount(IncludedRows));
            StdErrorOffsetMat = repmat(StdErrorOffsetVector, 1, NumAPbins);
            FluoOffsetVector = ExpOffsets(IncludedRows);
            FluoOffsetMat = repmat(FluoOffsetVector, 1, NumAPbins);
            MeanMat = MeanFluoMat(IncludedRows,:)-FluoOffsetMat;
            StdMat = sqrt(StdFluoMat(IncludedRows,:).^2 + (StdErrorOffsetMat).^2);
            MeanFluoMat = MeanMat;
            StdFluoMat =  StdMat;
        end
        NumNucMat = NumNucMat(IncludedRows,:);
        NCTimes = NCTimes(IncludedRows)/60;
        
        
        
        
        
        MaximumNCTime = max(NCTimes);
        MaxFluo = max(max(MeanFluoMat+StdFluoMat));
        MinFluo = min(min(MeanFluoMat+StdFluoMat));
        NumFrames = length(NCTimes);
        
        
        
        
        
        
        
        outdir3 = [outdir2, filesep, 'NC', num2str(NC)];
        
        if ~exist(outdir3, 'dir')
            mkdir(outdir3)
        end
        
        
        outdir4 = [outdir3, filesep, TraceType];
        
        
        if ~exist(outdir4, 'dir')
            mkdir(outdir4)
        end
        
        
        outdir5 = [outdir4, filesep, datestr(now, 'yyyymmdd')];
        
        if ~exist(outdir5, 'dir')
            mkdir(outdir5)
        end
        close all
        
        
        FrameProfFig = figure(1);
        set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
        set(gcf,'color','w');
        FrameProfAx = axes(FrameProfFig);
        
        
        
        temp_idx = find(temperatures == CurrentTemperature);
        
        
        eb = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
        hold on
        set(eb, 'color', colors(temp_idx,:), 'LineWidth', 1);
        set(get(get(eb, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        %set(get(get(eb, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        
        
        prof = plot(APbins, ones(1, length(APbins)),'o',...
            'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colors(temp_idx,:),...
            'MarkerSize', 6, 'LineStyle', '-', 'Color', colors(temp_idx,:));
        
        set(eb,'Visible','off'); %'off' or 'on'
        set(prof,'Visible','off'); %'off' or 'on'
        
        
        hold off
        if strcmp(lower(TraceType), 'anaphasealigned')
            xlabel('Time since anaphase (min)')
        else
            xlabel('Time since NC start (min)')
        end
        
        xlim([0-0.01* max(MaximumNCTime), max(MaximumNCTime)*1.1]);
        
        
        if ~UseOffsets
            ylim([max(0, min(MinFluo*0.95)), max(MaxFluo*1.05,1)])
        else
            ylim([min(min(min(min(MinFluo))*1.1,min(min(MinFluo))*0.95), 0) , max(max(MaxFluo)*1.1,1)]);
        end
        
        
        
        if ~UseOffsets
            ylabel('Fluo (AU)')
        else
            ylabel('Offset Fluo (AU)')
        end
        
        title(FrameProfAx, {'',...
            ['Nuclear Cycle ',num2str(NC),', Fraction Embryo Length: ',num2str(-1) ]})
        
        
        
        
        legend(FrameProfAx, legend_labels(SetIndex), 'Location', 'northwest')
   
        for i = MinAPbin:MaxAPbin
            APBinHasData = false;
            set(FrameProfAx.Children(2),'Visible','off'); %'off' or 'on'
            set(FrameProfAx.Children(1),'Visible','off'); %'off' or 'on'
            set(get(get(FrameProfAx.Children(1), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            
            
            use_idx = NumNucMat(:,i) >= this.MinimumNuclearCount;
            
            
            if sum(use_idx) == 0 %| sum(DiffMeanFluoMat(i, use_idx, j) == 0)
                FrameProfAx.Children(1).XData = NCTimes(use_idx);
                FrameProfAx.Children(1).YData = zeros(1, length(NCTimes(use_idx)));
                FrameProfAx.Children(2).XData = NCTimes(use_idx);
                FrameProfAx.Children(2).YData = zeros(1, length(NCTimes(use_idx)));
                FrameProfAx.Children(2).YPositiveDelta = zeros(1, length(NCTimes(use_idx)));
                FrameProfAx.Children(2).YNegativeDelta = zeros(1, length(NCTimes(use_idx)));
                set(FrameProfAx.Children(1),'Visible','off'); %'off' or 'on'
                set(FrameProfAx.Children(2),'Visible','off'); %'off' or 'on'
                set(get(get(FrameProfAx.Children(1), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            else
                APBinHasData = true;
                set(get(get(FrameProfAx.Children(1), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
                set(FrameProfAx.Children(1),'Visible','on'); %'off' or 'on'
                set(FrameProfAx.Children(2),'Visible','on'); %'off' or 'on'
               
                FrameProfAx.Children(1).YData = MeanFluoMat(use_idx, i);
                FrameProfAx.Children(1).XData = NCTimes(use_idx);
                FrameProfAx.Children(2).YData = MeanFluoMat(use_idx, i);
                FrameProfAx.Children(2).XData = NCTimes(use_idx);
                FrameProfAx.Children(2).YPositiveDelta = StdFluoMat(use_idx, i);
                FrameProfAx.Children(2).YNegativeDelta  = StdFluoMat(use_idx, i);
                
                set(FrameProfAx.Children(1),'Visible','on'); %'off' or 'on'
                set(FrameProfAx.Children(2),'Visible','on'); %'off' or 'on'
                
            end
            
            %try
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
            
            
            saveas(FrameProfFig,[outdir5, filesep,...
                OffsetString,'FluoTrace_NC',num2str(NC),'_Bin', num2str(i),'.png']);
            
        end
    end
end
close all
end