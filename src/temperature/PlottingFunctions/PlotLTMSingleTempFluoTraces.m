function PlotLTMSingleTempFluoTraces(this, outdir, varargin)
%%

% PlotTitle, PlottingColors, UseDifferentColors,
% UseDiffProfiles, UsePhysicalAPLength

NormedCycleTimes = false;

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
elseif ~strcmp(lower(PlottingColors), 'gradient') & ~strcmp(lower(PlottingColors), 'default')  & ~strcmp(lower(PlottingColors), 'pboc')
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


Temp_obs = this.Temp_obs;
Temp_sp = this.Temp_sps;



if strcmp(lower(PlottingColors), 'default')
    [~, colors] = getColorPalettes();
    GradString = '';
elseif strcmp(lower(PlottingColors), 'pboc')
    [colors, ~] = getColorPalettes();
    GradString = '';
else
    Temp_obs = this.Temp_obs(this.ProcessedExperiments);
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

if ~isempty(NormString)
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
temperatures = flip(unique(this.Temp_sps(this.ProcessedExperiments)));
NumTemperatures = length(temperatures);


APResolution = this.Experiments{1}.APResolution;
APbins = 0:APResolution:1;
NumAPbins = length(APbins);


legend_labels = this.LegendLabels;
MinimumTraceCount = this.MinimumTraceCount;


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


outdir2 = [outdir,filesep,'SingleTempFluoTraces'];
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




%%
for NC=this.IncludedNCs
    if NormedCycleTimes & NC == 14
        continue
    end
    
   for TemperatureIndex = 1:length(temperatures)

        CurrentTemperature = temperatures(TemperatureIndex);
        if isempty(TempMatches{TemperatureIndex})
            continue
        end
        outdir4 = [outdir3, filesep, 'T', strrep(num2str(CurrentTemperature), '.', '_'), 'C_NC', num2str(NC)];
        
        if ~exist(outdir4, 'dir')
            mkdir(outdir4)
        end
      
        
        outdir5 = [outdir4, filesep, datestr(now, 'yyyymmdd')];
        
        if ~exist(outdir5, 'dir')
            mkdir(outdir5)
        end
      
        
        % Prepare Traces for plotting
        % Prepare Traces for plotting
        MaximumNCTimes = NaN(1, length(TempMatches{TemperatureIndex}));
        MaxFluos = NaN(1, length(TempMatches{TemperatureIndex}));
        MinFluos = NaN(1, length(TempMatches{TemperatureIndex}));
        NumTimePoints = NaN(1,length(TempMatches{TemperatureIndex}));
        MeanFluoMats = cell(1, length(TempMatches{TemperatureIndex}));
        StdFluoMats = cell(1, length(TempMatches{TemperatureIndex}));
        NumNucMats = cell(1, length(TempMatches{TemperatureIndex}));
        NCTimes = cell(1, length(TempMatches{TemperatureIndex}));
        MinAPs = NaN(1, length(TempMatches{TemperatureIndex}));
        MaxAPs = NaN(1,length(TempMatches{TemperatureIndex}));
        
        for idx2=1:length(TempMatches{TemperatureIndex})
            idx = TempMatches{TemperatureIndex}(idx2);
            if ~ismember(idx, this.ProcessedExperiments)
                continue
            end
            ExpNCTimes = this.MeanProfiles{idx}.([TraceType, 'CycleFrameTimes']){NC-8};
            IncludedRows = 1:length(ExpNCTimes);
            if UsePerNucleusTraces
                ExpFluoMat = squeeze(this.MeanProfiles{idx}.([TraceType, 'CycleMeanPerNucleusTraces'])(IncludedRows,:,NC-8));
            else
                ExpFluoMat = squeeze(this.MeanProfiles{idx}.([TraceType, 'CycleMeanTraces'])(IncludedRows,:,NC-8));
            end
            ExpStdMat = squeeze(this.MeanProfiles{idx}.([TraceType, 'CycleTraceStdErrors'])(IncludedRows,:,NC-8));
            ExpNumNucMat = squeeze(this.MeanProfiles{idx}.([TraceType, 'CycleNumOnNuclei'])(IncludedRows,:,NC-8));
            NCLength = this.EmbryoStats.NCDivisionInfo(idx,NC-8);
            
            IncludedColumns = find(sum(~isnan(ExpFluoMat),1).' > 0);
            if ~isempty(IncludedColumns)
                MinAPs(idx2) = min(IncludedColumns);
                MaxAPs(idx2) = max(IncludedColumns);
            end
            
            IncludedRows = find(sum(~isnan(ExpFluoMat),2).' > 0);
            if isempty(IncludedRows)
                MeanFluoMats{idx2} = [];
                StdFluoMats{idx2} = [];
                NumNucMats{idx2} = [];
                NCTimes{idx2} = [];
            else
                MeanFluoMats{idx2} = ExpFluoMat(IncludedRows,:);
                StdFluoMats{idx2} = ExpStdMat(IncludedRows,:);
                NumNucMats{idx2} = ExpNumNucMat(IncludedRows,:);
                
                if ~NormedCycleTimes
                    NCTimes{idx2} = ExpNCTimes(IncludedRows)/60;
                else
                    NCTimes{idx2} = (ExpNCTimes(IncludedRows)/60)/NCLength;
                end
                
                MaximumNCTimes(idx2) = max(NCTimes{idx2});
                MaxFluos(idx2) = max(max(MeanFluoMats{idx2}+StdFluoMats{idx2}));
                MinFluos(idx2) = min(min(MeanFluoMats{idx2}+StdFluoMats{idx2}));
                NumTimePoints(idx2) = uint16(max(NCTimes{idx2})/(this.time_delta/60)+1);
            end
            
    
        end
        
        if all(isnan(NumTimePoints))
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
        eb = cell(1, length(TempMatches{TemperatureIndex}));
        prof = cell(1, length(TempMatches{TemperatureIndex}));
    
        for idx =1:length(TempMatches{TemperatureIndex})
            
            if strcmp(lower(PlottingColors), 'gradient')
                temp_idx = find(abs(Temp_range - Temp_obs(idx)) == min(abs(Temp_range - Temp_obs(TempMatches{TemperatureIndex}(idx)))), 1);

            else
                temp_idx = TemperatureIndex;
                
            end
            
            marker_idx = idx;
            
            
            
            
            
            eb{idx} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
            hold on
            
            set(eb{idx}, 'color', colors(temp_idx,:), 'LineWidth', 1);
            
            set(get(get(eb{idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            
            if marker_idx <= length(MarkerStyles) 
            prof{idx} = plot(APbins, ones(1, length(APbins)), MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colors(temp_idx,:),...
                'MarkerSize', 6, 'LineStyle', '-', 'Color', colors(temp_idx,:));
            else
                prof{idx} = plot(APbins, ones(1, length(APbins)), MarkerStyles{mod(marker_idx, length(MarkerStyles))+1},...
                'MarkerEdgeColor', colors(temp_idx,:), 'MarkerFaceColor', colors(temp_idx,:),...
                'MarkerSize', 6, 'LineStyle', '-', 'Color', colors(temp_idx,:));
            end
            
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
        grid on 
        hold off
        
        if NormedCycleTimes
            xlabel('Fraction of Nuclear Cycle')
        elseif strcmpi(TraceType, 'anaphasealigned') | strcmp(lower(TraceType), 'anaphasealigned3d')
            xlabel('Time since anaphase (min)')
        else
            xlabel('Time since NC start (min)')
        end
        FrameProfAx.XAxis.FontSize = 18;
        FrameProfAx.YAxis.FontSize = 18;
        if ~NormedCycleTimes
            xlim([0-0.01* max(MaximumNCTimes), max(MaximumNCTimes)*1.1]);
        else
            xlim([-0.01, 1.01]);
        end
        
        if strcmp(lower(TraceType), 'tbinned') | strcmp(lower(TraceType), 'anaphasealigned') | strcmpi(TraceType, 'unaligned')
            ylabel('Fluo (AU)')
        else
            ylabel('3D Fluo (AU)')
        end
   
        ylim([0, max(max(MaxFluos)*1.1,1)])
      
        %
        
        title(FrameProfAx, {'',...
            ['T = ', num2str(temperatures(TemperatureIndex)), '°C, Nuclear Cycle ',num2str(NC),', Fraction Embryo Length: ',num2str(-1) ]})
        
        

        for i = MinAPbin:MaxAPbin
            PlottedSets = zeros(1, NumSets, 'logical');
            for idx=1:length(TempMatches{TemperatureIndex})
                set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','off'); %'off' or 'on'
                set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','off'); %'off' or 'on'
                set(get(get(FrameProfAx.Children(end-(2*(idx-1)+1)), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');

                if ~ismember(TempMatches{TemperatureIndex}(idx), this.ProcessedExperiments)
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
                    PlottedSets(TempMatches{TemperatureIndex}(idx)) = true;
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
            if ~sum(PlottedSets)
                continue
            end
            if NumSets > 25
                lgd = legend(FrameProfAx, legend_labels(PlottedSets),'NumColumns', 2, 'Location', 'eastoutside');
            else
                lgd = legend(FrameProfAx, legend_labels(PlottedSets), 'Location', 'eastoutside');
            end
            if exist('PlotTitle', 'var')
                title(FrameProfAx, {PlotTitle,...
                    ['T = ', num2str(temperatures(TemperatureIndex)), '°C, Nuclear Cycle ',num2str(NC),', Fraction Embryo Length: ',num2str(APbins(i)) ]}, 'FontSize', 20)
            else
                title(FrameProfAx, ['T = ', num2str(temperatures(TemperatureIndex)), '°C, Nuclear Cycle ',num2str(NC),', Fraction Embryo Length: ',num2str(APbins(i)) ], 'FontSize', 20)
  
            end
       
             saveas(FrameProfFig,[outdir4, filesep, GradString, NormString,...
                'FluoTrace_NC',num2str(NC),'_T',strrep(num2str(CurrentTemperature), '.', '_'),'_Bin', num2str(i),'.png']);
            
        end
    end
end
close all
end