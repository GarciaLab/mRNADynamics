function PlotLTPFluoProfiles(this, outdir, varargin)
%%

% PlotTitle, PlottingColors, UseDifferentColors,
% UseDiffProfiles, UsePhysicalAPLength
UsePhysicalAPLength = false;
UseOffsets = false;
useRescaledFluo = false;
useRescaledTiming = false;

x = 1;
while x <= length(varargin)
    if strcmp(lower(varargin{x}), 'plottitle')
        PlotTitle = varargin{x+1};
        x = x+1;
    elseif strcmp(lower(varargin{x}), 'plottingcolors')
        PlottingColors = varargin{x+1};
        x = x+1;
    elseif strcmpi(varargin{x}, 'subtractoffsets') | strcmpi(varargin{x}, 'useoffsets') 
        UseOffsets = true;
    elseif strcmp(lower(varargin{x}), 'usephysicalaplength')
        UsePhysicalAPLength = true;
    elseif strcmp(lower(varargin{x}), 'tracetype')
        TraceType = lower(varargin{x+1});
        x = x+1;
    elseif strcmpi(varargin{x}, 'rescalefluo') | strcmpi(varargin{x}, 'userescaledfluo')
        useRescaledFluo = true;
    elseif strcmpi(varargin{x}, 'userescaledtime') | strcmpi(varargin{x}, 'rescaletime') | ...
            strcmpi(varargin{x}, 'rescaletiming') | strcmpi(varargin{x}, 'userescaledtiming')
        useRescaledTiming = true;
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
elseif  ~strcmp(lower(TraceType), 'tbinned')  & ~strcmp(lower(TraceType), 'anaphasealigned')
    error('Invalid choice of trace type. Can use either "tbinned" or "anaphasealigned".') % change to error
end

if strcmpi(TraceType, 'anaphasealigned')
    TraceType = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'tbinned')
    TraceType = 'Tbinned';
end


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
if UsePhysicalAPLength
    PhysicalAPString = 'PhysAP';
else
    PhysicalAPString = '';
end

if UseOffsets
    OffsetString = 'Offset';
else
    OffsetString = '';
end

if useRescaledFluo
    FluoString = 'RescaledFluo';
else
    FluoString = '';
end

if useRescaledTiming
    TimingString = 'DevTime';
else
    TimingString = '';
end


SpecificDirString = FluoString;
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

if ~isempty(OffsetString)
    if ~isempty(SpecificDirString)
        SpecificDirString = [SpecificDirString, '_', OffsetString];
    else
        SpecificDirString = OffsetString;
    end
end

if ~isempty(GradString)
    if ~isempty(SpecificDirString)
        SpecificDirString = [SpecificDirString, '_', GradString];
    else
        SpecificDirString = GradString;
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

%%

outdir2 = [outdir,filesep,'FluoProfiles'];
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

for NC=this.IncludedNCs
    
    if ~useRescaledTiming
    % Prepare Traces for plotting
    MaximumNCTimes = NaN(1, NumSets);
    MaxFluos = NaN(1, NumSets);
    MinFluos = NaN(1, NumSets);
    NumTimePoints = NaN(1, NumSets);
    MeanFluoMats = cell(1, NumSets);
    StdFluoMats =cell(1, NumSets);
    NumNucMats = cell(1, NumSets);
    NCTimes = cell(1, NumSets);
    MinAPs = NaN(1, NumSets);
    MaxAPs = NaN(1, NumSets);
    
    
    for idx=1:NumSets
        if ~ismember(idx, this.ProcessedExperiments)
            continue
        end
        ExpNCTimes = this.TFProfiles{idx}.([TraceType, 'CycleFrameTimes']){NC-8};
        IncludedRows = 1:length(ExpNCTimes);
        ExpFluoMat = squeeze(this.TFProfiles{idx}.([TraceType, 'CycleMeanTraces'])(IncludedRows,:,NC-8));
        ExpStdMat = squeeze(this.TFProfiles{idx}.([TraceType, 'CycleTraceStdErrors'])(IncludedRows,:,NC-8));
        ExpNumNucMat = squeeze(this.TFProfiles{idx}.([TraceType, 'CycleNumOnNuclei'])(IncludedRows,:,NC-8));
        if UseOffsets 
            IncludedRows2 = uint16(ExpNCTimes/this.time_delta+1);
            
            ExpOffsets = squeeze(this.SetFluoOffsets.(TraceType).mean(idx, IncludedRows2,NC-8)).';
            ExpOffsetStd = squeeze(this.SetFluoOffsets.(TraceType).std(idx, IncludedRows2,NC-8)).';
            ExpOffsetCount = squeeze(this.SetFluoOffsets.(TraceType).count(idx, IncludedRows2,NC-8)).';
        end
        
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
            if ~UseOffsets 
                MeanFluoMats{idx} = ExpFluoMat(IncludedRows,:);
                StdFluoMats{idx} = ExpStdMat(IncludedRows,:);
                
            else
                StdErrorOffsetVector = ExpOffsetStd(IncludedRows)./sqrt(ExpOffsetCount(IncludedRows));
                StdErrorOffsetMat = repmat(StdErrorOffsetVector, 1, NumAPbins);
                FluoOffsetVector = ExpOffsets(IncludedRows);
                FluoOffsetMat = repmat(FluoOffsetVector, 1, NumAPbins);
                MeanMat = ExpFluoMat(IncludedRows,:)-FluoOffsetMat;
                StdMat = sqrt(ExpStdMat(IncludedRows,:).^2 + (StdErrorOffsetMat).^2);
                
                
                MeanFluoMats{idx} = MeanMat;
                StdFluoMats{idx} = StdMat;
                
                
                
            end
            
            if useRescaledFluo
                temp_match = find(round(temperatures, 2) ==  round(Temp_sp(idx), 2));
                MeanFluoMats{idx} = MeanFluoMats{idx}*this.FluoCoeffs(temp_match);
                StdFluoMats{idx} = StdFluoMats{idx}*this.FluoCoeffs(temp_match);
            end
            NumNucMats{idx} = ExpNumNucMat(IncludedRows,:);
                
            NCTimes{idx} = ExpNCTimes(IncludedRows)/60;
          
            MaximumNCTimes(idx) = max(NCTimes{idx});
            MaxFluos(idx) = max(max(MeanFluoMats{idx}+StdFluoMats{idx}));
            MinFluos(idx) = min(min(MeanFluoMats{idx}-StdFluoMats{idx}));
            NumTimePoints(idx) = uint16(max(NCTimes{idx})/(this.time_delta/60)+1);
        end
    end
    
    else
        [MeanFluoMats, StdFluoMats, NumNucMats, NCTimes, MaximumNCTimes,...
            MaxFluos, MinFluos, NumTimePoints, MinAPs, MaxAPs] = getLTPRescaledTimingMats(this, TraceType, NC, useRescaledFluo, UseOffsets);
    end
    
    if all(isnan(NumTimePoints))
        continue
    end
    
    MinAPbin = nanmin(MinAPs);
    MaxAPbin = nanmax(MaxAPs);
    
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
    eb = cell(1, NumSets);
    prof = cell(1, NumSets);
    
    
    for idx =1:NumSets
        if strcmp(lower(PlottingColors), 'gradient')
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
        
        %temp_idx = mod(i, 7)+1;
        eb{idx} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
        hold on
        set(eb{idx}, 'color', colors(temp_idx,:), 'LineWidth', 1);
        set(get(get(eb{idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        
        prof{idx} = plot(APbins, ones(1, length(APbins)), MarkerStyles{marker_idx},...
        'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colors(temp_idx,:),...
            'MarkerSize', 6, 'LineStyle', '-', 'Color', colors(temp_idx,:));
        
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
    
    if ~UsePhysicalAPLength
        xlabel('Fraction Embryo Length')
        xlim([MinAPLimit, MaxAPLimit])
    else
        xlabel('Distance from the Anterior Pole (\mum)')
        xlim([MinAPLimit, MaxAPLimit])
    end
    

    
    if ~useRescaledFluo
        if ~UseOffsets
        if strcmp(lower(TraceType), 'fluo') | strcmp(lower(TraceType), 'anaphasealigned') | strcmp(lower(TraceType), 'tbinned')
            ylabel('Fluo (AU)')
        else
            ylabel('3D Fluo (AU)')
        end
        else
            if strcmp(lower(TraceType), 'fluo') | strcmp(lower(TraceType), 'anaphasealigned') | strcmp(lower(TraceType), 'tbinned')
            ylabel('Offset Fluo (AU)')
        else
            ylabel('Offset 3D Fluo (AU)')
        end 
        end
    else
        if ~UseOffsets
        if strcmp(lower(TraceType), 'fluo') | strcmp(lower(TraceType), 'anaphasealigned') | strcmp(lower(TraceType), 'tbinned')
            ylabel('Re-scaled fluo (AU)')
        else
            ylabel('Re-scaled 3D fluo (AU)')
        end
        else
        if strcmp(lower(TraceType), 'fluo') | strcmp(lower(TraceType), 'anaphasealigned') | strcmp(lower(TraceType), 'tbinned')
            ylabel('Re-scaled offset fluo (AU)')
        else
            ylabel('Re-scaled offset 3D fluo (AU)')
        end
        end
    end
    
    if ~UseOffsets
        ylim([0, max(max(MaxFluos)*1.1,1)])
    else
        ylim([min(min(min(min(MinFluos))*1.1,min(min(MinFluos))*0.95), 0) , max(max(MaxFluos)*1.1,1)]);
    end
    %
    if strcmp(lower(TraceType), 'anaphasealigned')
        title(FrameProfAx, {'',...
            ['Nuclear Cycle ',num2str(NC),', ',num2str(-1), ' min since anaphase' ]})
    else
        title(FrameProfAx, {'',...
            ['Nuclear Cycle ',num2str(NC),', ',num2str(-1), ' min' ]})
    end
    
    
    for i = 1:nanmax(NumTimePoints)
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
            if i > size(NumNucMats{idx}, 1)
                continue
            end
            
            if UsePhysicalAPLength
                APLength = this.APLengths(idx);
            else
                APLength = 1;
            end
            
            time_index = find(NCTimes{idx}/(this.time_delta/60)+1 == i, 1);
            if isempty(time_index)
                continue 
            end
            use_idx = NumNucMats{idx}(time_index,:) >= MinimumTraceCount;
            
            if sum(use_idx) == 0 %| sum(DiffMeanFluoMat(i, use_idx, j) == 0)
                FrameProfAx.Children(end-(2*(idx-1)+1)).XData = APLength*APbins;
                FrameProfAx.Children(end-(2*(idx-1)+1)).YData = zeros(1, length(APbins));
                FrameProfAx.Children(end-(2*(idx-1))).XData = APLength*APbins;
                FrameProfAx.Children(end-(2*(idx-1))).YData = zeros(1, length(APbins));
                FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = .1*ones(1, length(APbins));
                FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta = .1*ones(1, length(APbins));
                set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','off'); %'off' or 'on'
                set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','off'); %'off' or 'on'
                set(get(get(FrameProfAx.Children(end-(2*(idx-1)+1)), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            else
                PlottedSets(idx) = true;
                FrameProfAx.Children(end-(2*(idx-1)+1)).YData = MeanFluoMats{idx}(time_index, use_idx);
                FrameProfAx.Children(end-(2*(idx-1)+1)).XData =APLength*APbins(use_idx);
                FrameProfAx.Children(end-(2*(idx-1))).YData = MeanFluoMats{idx}(time_index, use_idx);
                FrameProfAx.Children(end-(2*(idx-1))).XData = APLength*APbins(use_idx);
                FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = StdFluoMats{idx}(time_index, use_idx);
                FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta  = StdFluoMats{idx}(time_index, use_idx);
                
                set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','on'); %'off' or 'on'
                set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','on'); %'off' or 'on'
                set(get(get(FrameProfAx.Children(end-(2*(idx-1)+1)), 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
                
            end
        end
        %try
        if sum(PlottedSets) == 0
            continue
        end
        
        legend(FrameProfAx, legend_labels(PlottedSets), 'Location', 'eastoutside')
        if ~useRescaledTiming
        if exist('PlotTitle', 'var')
            
            if strcmp(lower(TraceType), 'anaphasealigned') 
                title(FrameProfAx, {PlotTitle,...
                    ['Nuclear Cycle ',num2str(NC),', ',num2str((i-1)*this.time_delta/60), ' min. since anaphase' ]})
            else
                title(FrameProfAx, {PlotTitle,...
                    ['Nuclear Cycle ',num2str(NC),', ',num2str((i-1)*this.time_delta/60), ' min.' ]})
            end
        else
            if strcmp(lower(TraceType), 'anaphasealigned') 
                title(FrameProfAx,  ['Nuclear Cycle ',num2str(NC),', ',num2str((i-1)*this.time_delta/60), ' min. since anaphase' ])
            else
                title(FrameProfAx,  ['Nuclear Cycle ',num2str(NC),', ',num2str((i-1)*this.time_delta/60), ' min.' ])
                
            end
        end
        else
           if exist('PlotTitle', 'var')
            
            if strcmp(lower(TraceType), 'anaphasealigned') 
                title(FrameProfAx, {PlotTitle,...
                    ['Nuclear Cycle ',num2str(NC),', ',num2str((i-1)*this.time_delta/60), ' temperature-scaled min. since anaphase' ]})
            else
                title(FrameProfAx, {PlotTitle,...
                    ['Nuclear Cycle ',num2str(NC),', ',num2str((i-1)*this.time_delta/60), ' temperature-scaled min.' ]})
            end
        else
            if strcmp(lower(TraceType), 'anaphasealigned') 
                title(FrameProfAx,  ['Nuclear Cycle ',num2str(NC),', ',num2str((i-1)*this.time_delta/60), ' temperature-scaled min. since anaphase' ])
            else
                title(FrameProfAx,  ['Nuclear Cycle ',num2str(NC),', ',num2str((i-1)*this.time_delta/60), ' temperature-scaled min.' ])
                
            end
        end 
            
        end
        
    
        %end
        
        outpath=[outdir3, filesep,...
            PhysicalAPString, GradString, OffsetString, FluoString, TimingString,'FluoProfile_NC',num2str(NC),'_Frame', num2str(i),'.png'];
        saveas(FrameProfFig,outpath);
        
    end
    
end
close all
end