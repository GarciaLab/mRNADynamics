function PlotLTMFluoTraces(this, outdir, varargin)
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
    TraceType = 'fluo3d';
elseif ~strcmp(lower(TraceType), 'fluo3d') & ~strcmp(lower(TraceType), 'fluo')  & ~strcmp(lower(TraceType), 'anaphasealigned')& ~strcmp(lower(TraceType), 'anaphasealigned3d')
    error('Invalid choice of trace type. Can use either "fluo", "fluo3d", "anaphasealigned", or "anaphasealigned3d".') % change to error
end

if ~exist(outdir, 'dir')
    mkdir(outdir)
end
%%
if strcmpi(TraceType, 'anaphasealigned')
    traceName = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'anaphasealigned3d')
    traceName = 'AnaphaseAligned3D';
elseif strcmpi(TraceType, 'fluo')
    traceName = 'Unaligned';
elseif strcmpi(TraceType, 'fluo3d')
    traceName = 'Unaligned3D';
elseif strcmpi(TraceType, 'tbinned')
    traceName = 'Tbinned';
elseif strcmpi(TraceType, 'tbinned3d')
    traceName = 'Tbinned3D';
end

%%


Temp_obs = this.Temp_obs(this.ProcessedExperiments);
Temp_sp = this.Temp_sps(this.ProcessedExperiments);
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

temperatures = flip(unique(this.Temp_sps(this.ProcessedExperiments)));
APResolution = this.Experiments{1}.APResolution;

APbins = 0:APResolution:1;

legend_labels = this.LegendLabels(this.ProcessedExperiments);
MinimumTraceCount = this.MinimumTraceCount;
for nc_idx=1:length(this.IncludedNCs)
    NC = this.IncludedNCs(nc_idx);
    
    
    
    
    
    outdir3 = [outdir, filesep, datestr(now, 'yyyymmdd')];
    if ~exist(outdir3, 'dir')
        mkdir(outdir3)
    end
    
    % Prepare Traces for plotting
    MaximumNCTimes = NaN(1, length(Temp_sp));
    MaxFluos = NaN(1, length(Temp_sp));
    NumFrames = NaN(1, length(Temp_sp));
    MeanFluoMats = {};
    StdFluoMats = {};
    NumNucMats = {};
    NCTimes = {};
    for i=1:length(Temp_sp)
        idx = this.ProcessedExperiments(i);
        ExpNCTimes = this.MeanProfiles{idx}.([traceName, 'CycleFrameTimes']){NC-8};
        IncludedRows = 1:length(ExpNCTimes);
        ExpFluoMat = squeeze(this.MeanProfiles{idx}.([traceName, 'CycleMeanTraces'])(IncludedRows,:,NC-8));
        ExpStdMat = squeeze(this.MeanProfiles{idx}.([traceName, 'CycleTraceStdErrors'])(IncludedRows,:,NC-8));
        ExpNumNucMat = squeeze(this.MeanProfiles{idx}.([traceName, 'CycleNumOnNuclei'])(IncludedRows,:,NC-8));
       
        
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
    set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .4]);
    set(gcf,'color','w');
    FrameProfAx = axes(FrameProfFig);
    for idx =1:length(Temp_sp)
        if strcmp(lower(PlottingColors), 'gradient')
            temp_idx = find(abs(Temp_range - Temp_obs(idx)) == min(abs(Temp_range - Temp_obs(idx))), 1);
        else
            temp_idx = find(temperatures == Temp_sp(idx));
        end
        %temp_idx = mod(i, 7)+1;
        eb{idx} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
        hold on
        if strcmp(lower(PlottingColors), 'gradient')
            set(eb{idx}, 'color', colors(temp_idx,:), 'LineWidth', 1);
        else
            set(eb{idx}, 'color', colors(temp_idx,:), 'LineWidth', 1);
        end
        set(get(get(eb{idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        if strcmp(lower(PlottingColors), 'gradient')
            prof{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(temp_idx,:));
        else
            prof{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(temp_idx,:));
        end
        set(eb{idx},'Visible','off'); %'off' or 'on'
        set(prof{idx},'Visible','off'); %'off' or 'on'
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
    if strcmp(lower(TraceType), 'anaphasealigned') | strcmp(lower(TraceType), 'anaphasealigned3d')
        xlabel('Time since anaphase (min)')
    else
        xlabel('Time since NC start (min)')
    end
    xlim([0, max(MaximumNCTimes)])
    
    if strcmp(lower(TraceType), 'fluo') | strcmp(lower(TraceType), 'anaphasealigned')
        ylabel('Fluo (AU)')
    else
        ylabel('3D Fluo (AU)')
    end
    ylim([0, max(max(MaxFluos)*1.1,1)])
    %
    
    title(FrameProfAx, {'',...
        ['Nuclear Cycle ',num2str(NC),', Fraction Embryo Length: ',num2str(-1) ]})
    
    
    
    
    legend(FrameProfAx, legend_labels, 'Location', 'eastoutside')
    for i = 1:length(APbins)
        APBinHasData = false;
        for idx=1:length(Temp_sp)
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
            else
                APBinHasData = true;
                FrameProfAx.Children(end-(2*(idx-1)+1)).YData = MeanFluoMats{idx}(use_idx, i);
                FrameProfAx.Children(end-(2*(idx-1)+1)).XData =NCTimes{idx}(use_idx);
                FrameProfAx.Children(end-(2*(idx-1))).YData = MeanFluoMats{idx}(use_idx, i);
                FrameProfAx.Children(end-(2*(idx-1))).XData = NCTimes{idx}(use_idx);
                FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = StdFluoMats{idx}(use_idx, i);
                FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta  = StdFluoMats{idx}(use_idx, i);
                
                set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','on'); %'off' or 'on'
                set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','on'); %'off' or 'on'
                
            end
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
        
        saveas(FrameProfFig,[outdir3, filesep,...
            'NC',num2str(NC), GradString, 'FluoTrace',TraceType,'_', num2str(i),'.png']);
        
    end
    
end
close all
end