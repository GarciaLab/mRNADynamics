classdef LTP
    %liveProject Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        ProjectName = '';
        ProjectType = '';
        ExperimentTypes = {};
        ExperimentPrefixes = {};
        Experiments = {};
        ExperimentStatuses = {};
        IncludedExperiments = [];
        ProcessedExperiments = [];
        
        IncludedNCs = [];
        Regions = {};
        Temp_sps = [];
        Temp_obs = [];
        MinAPs = [];
        MaxAPs = [];
        
        
        MinOffsetAP = [];
        MaxOffsetAP = [];
        time_delta = [];
        MinimumNuclearCount = [];
        MinimumTimePoints = [];
        
        TFProfiles = {};
        TempFluoOffsets = {};
        SetFluoOffsets = {};

        
        LegendLabels = {};
        
        APLengths = [];
        
        EmbryoStats = {};
        
        UniqueTemperatures = [];
        TimingCoeffs = [];
        FluoCoeffs = [];
        
        FluoScalingInfo = {};
        TimeScalingInfo = {};
        
        
    end
    
    
    
    methods
        %% Constructors
        
        
        function this = LTP(ProjectName, ProjectType, IncludedNCs, time_delta, MinimumNuclearCount)
            %liveProject Construct an instance of this class
            %   Detailed explanation goes here
            this.ProjectName = ProjectName;
            
            if exist('ProjectType', 'var')
                this.ProjectType = ProjectType;
            end
            
            if exist('IncludedNCs', 'var')
                this.IncludedNCs = IncludedNCs;
            else
                this.IncludedNCs =[12, 13, 14];
            end
            
            if exist('time_delta', 'var')
                this.time_delta = time_delta;
            else
                this.time_delta = 60; % unit: seconds
            end
            
            if exist('MinimumNuclearCount', 'var')
                this.MinimumNuclearCount = MinimumNuclearCount;
            else
                this.MinimumNuclearCount = 5; % unit: minutes
            end
            
            AllPrefixes = getProjectPrefixes(ProjectName);
            TempVector = zeros(1, length(AllPrefixes));
            for i = 1:length(AllPrefixes)
                TempVector(i) = pullTempSpFromDataStatus(ProjectName, AllPrefixes{i});
            end
             [~, indexorder] = sort(TempVector, 'descend');
            this.ExperimentPrefixes = AllPrefixes(indexorder);
            
            
            
            for i = 1:length(this.ExperimentPrefixes)
                Prefix = this.ExperimentPrefixes{i};
                this.Experiments{i} = LiveExperiment(Prefix);
                this.ExperimentTypes{i} = this.Experiments{i}.experimentType;
                this.ExperimentStatuses{i} = TemperatureExpStatus(ProjectName, Prefix);
                include_set(i) = this.ExperimentStatuses{i}.include_set;
                if exist('ProjectType', 'var')
                    if lower(this.ExperimentTypes{i}) ~= lower(ProjectType)
                        include_set(i) = 0;
                    end
                end
                if this.ExperimentStatuses{1}.hasCompiledNuclearProtein
                    finished_set(i) = 1;
                else
                    finished_set(i) = 0;
                end
                this.Regions{i} = this.ExperimentStatuses{i}.Region;
                this.Temp_sps(i) = this.ExperimentStatuses{i}.Temp_sp;
                this.Temp_obs(i) = this.ExperimentStatuses{i}.Temp_obs;
                if ~isempty(this.ExperimentStatuses{i}.MinAP)
                    this.MinAPs(i) = this.ExperimentStatuses{i}.MinAP;
                else
                    this.MinAPs(i) = NaN;
                end
                if ~isempty(this.ExperimentStatuses{i}.MaxAP)
                    this.MaxAPs(i) = this.ExperimentStatuses{i}.MaxAP;
                else
                    this.MaxAPs(i) = NaN;
                end
            end
            
            this.IncludedExperiments = find(include_set > 0);
            this.ProcessedExperiments = find((include_set > 0) & (finished_set > 0));
            
            if isempty(this.MinOffsetAP)
                if contains(lower(this.ProjectName), 'hb')
                    this.MinOffsetAP = 0.575;
                end
            end
            if isempty(this.MaxOffsetAP)
                if contains(lower(this.ProjectName), 'hb')
                    this.MaxOffsetAP = 0.675;
                end
            end
            

    
            this = AddTFProfiles(this);
            this = AddNBFluoOffsets(this);
            this.APLengths = AddAPLengths(this);
            this = AddHealthInfo(this);
        end
        
        
        
        %% Methods
        
        
        
        
      
        
        function [FluoOffset, FluoOffsetStdError] = getTFoffset(this, NC)
            
            temperatures = flip(unique(this.Temp_sps(this.ProcessedExperiments)));
            nc_idx = find(this.IncludedNCs == NC, 1);
            MeanFluoMat = this.MeanFluoMats{nc_idx};
            StdFluoMat = this.StdFluoMats{nc_idx};
            NumNucMat = this.NumNucMats{nc_idx};
            
            MinimumNucleiCount = this.MinimumNuclearCount;
            APResolution = 1/(size(MeanFluoMat, 2)-1);
            numAPbins = size(MeanFluoMat, 2);
            numTimePoints = size(MeanFluoMat, 1);
            numTemperatures = length(temperatures);
            numSets = size(MeanFluoMat, 3);
            
            FluoOffset = NaN(numTimePoints, numSets);
            FluoOffsetStdError = NaN(numTimePoints,numSets);
            MinOffsetAPbin =  uint16(ceil(this.MinOffsetAP/APResolution))+1;
            MaxOffsetAPbin =  uint16(floor(this.MaxOffsetAP/APResolution))+1;
            OffsetAPbins = MinOffsetAPbin:MaxOffsetAPbin;
            for t_idx=1:numTemperatures
                temp_matches =  find(this.Temp_sps(this.ProcessedExperiments) == temperatures(t_idx));
                for i=1:numTimePoints
                    if sum(sum(NumNucMat(i, OffsetAPbins, temp_matches) > 5)) == 0
                        continue
                    end
                    OffsetNumNucMatrix = squeeze(NumNucMat(i, OffsetAPbins, temp_matches));
                    OffsetMatrix = squeeze(MeanFluoMat(i, OffsetAPbins, temp_matches));
                    OffsetMatrix(OffsetNumNucMatrix <= MinimumNucleiCount) = NaN;
                    OffsetStdMatrix = squeeze(StdFluoMat(i, OffsetAPbins, temp_matches));
                    OffsetStdMatrix(OffsetNumNucMatrix < MinimumNucleiCount) = NaN;
                    Offset = nanmean(nanmean(OffsetMatrix, 1));
                    TimePointVariances = (OffsetStdMatrix.^2);
                    EmbryoTimePointStdErrors = sqrt(nansum(TimePointVariances./sum(~isnan(OffsetStdMatrix), 1), 1));
                    OffsetStdError = sqrt(sum(EmbryoTimePointStdErrors.^2))/sqrt(sum(~isnan(EmbryoTimePointStdErrors)));
                    for prefix_idx=temp_matches
                        if sum(NumNucMat(i, OffsetAPbins, prefix_idx)) == 0
                            FluoOffset(i, prefix_idx) = Offset;
                            FluoOffsetStdError(i, prefix_idx) = OffsetStdError;
                        else
                            SingleEmbryoOffsetNumNucMatrix = squeeze(NumNucMat(i, OffsetAPbins, prefix_idx));
                            SingleEmbryoOffsetMatrix = squeeze(MeanFluoMat(i, OffsetAPbins, prefix_idx));
                            SingleEmbryoOffsetStdMatrix = squeeze(StdFluoMat(i, OffsetAPbins, prefix_idx));
                            SingleEmbryoOffset = nanmean(SingleEmbryoOffsetMatrix(SingleEmbryoOffsetNumNucMatrix >= MinimumNucleiCount));
                            SingleEmbryoVariances = SingleEmbryoOffsetStdMatrix.^2;
                            SingleEmbryoOffsetStdError = sqrt(sum(SingleEmbryoVariances/sum(~isnan(SingleEmbryoVariances))));
                            FluoOffset(i, prefix_idx) = SingleEmbryoOffset;
                            FluoOffsetStdError(i, prefix_idx) = SingleEmbryoOffsetStdError;
                        end
                    end
                end
                
            end
        end
        
        function PlotSingleTempFluoProfiles(this, outdir, varargin)
            % PlotTitle, PlottingColors, UseDifferentColors,
            % UseDiffProfiles, UsePhysicalAPLength
            UsePhysicalAPLength = false;
            UseDiffProfiles = false;
            UseDifferentColors = true;
            
            x = 1;
            while x <= length(varargin)
                if strcmp(lower(varargin{x}), 'plottitle')
                    PlotTitle = varargin{x+1};
                    x = x+1;
                elseif strcmp(lower(varargin{x}), 'plottingcolors')
                    PlottingColors = varargin{x+1};
                    x = x+1;
                elseif strcmp(lower(varargin{x}), 'useuniformcolor')
                    UseDifferentColors = false;
                elseif strcmp(lower(varargin{x}), 'usediffprofiles')
                    UseDiffProfiles = true;
                elseif strcmp(lower(varargin{x}), 'usephysicalaplength')
                    UsePhysicalAPLength = true;
                end
                x = x+1;
            end
            if ~exist('PlottingColors', 'var')
                PlottingColors = 'default';
            elseif ~strcmp(lower(PlottingColors), 'gradient') & ~strcmp(lower(PlottingColors), 'default')  & ~strcmp(lower(PlottingColors), 'pboc')
                error('Invalid choice of plotting colors. Can use either "default", "pboc", or "gradient".') % change to error
            end
            if ~exist(outdir, 'dir')
                mkdir(outdir)
            end
            if UsePhysicalAPLength
                PhysicalAPString = 'PhysAP';
            else
                PhysicalAPString = '';
            end
            if UseDiffProfiles
                DiffString = 'Diff';
            else
                DiffString = '';
            end
            
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
            temperatures = flip(unique(this.Temp_sps(this.ProcessedExperiments)));
            APResolution = 1/(size(this.MeanFluoMats{length(this.IncludedNCs)}, 2)-1);
            
            APbins = 0:APResolution:1;
            
            legend_labels = this.LegendLabels(this.ProcessedExperiments);
            MinimumNuclearCount = this.MinimumNuclearCount;
            for nc_idx=1:length(this.IncludedNCs)
                NC = this.IncludedNCs(nc_idx);
                MeanAPMat = this.MeanAPMats{nc_idx};
                if ~UseDiffProfiles
                    MeanFluoMat = this.MeanFluoMats{nc_idx};
                    StdFluoMat = this.StdFluoMats{nc_idx};
                else
                    MeanFluoMat = this.DiffMeanFluoMats{nc_idx};
                    StdFluoMat = this.DiffStdErrorMats{nc_idx};
                end
                NumNucMat = this.NumNucMats{nc_idx};
                for t_idx=1:length(temperatures)
                    current_temp = temperatures(t_idx);
                    temp_matches =  find(this.Temp_sps(this.ProcessedExperiments) == current_temp);
                    if isempty(temp_matches)
                        continue
                    end
                    outdir2 = [outdir, filesep, 'T', strrep(num2str(current_temp), '.', '_'), 'C_NC', num2str(NC)];
                    if ~exist(outdir2, 'dir')
                        mkdir(outdir2)
                    end
                    outdir3 = [outdir2, filesep, datestr(now, 'yyyymmdd')];
                    if ~exist(outdir3, 'dir')
                        mkdir(outdir3)
                    end
                    close all
                    
                    
                    FrameProfFig = figure(1);
                    set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .4]);
                    set(gcf,'color','w');
                    FrameProfAx = axes(FrameProfFig);
                    for idx =1:length(temp_matches)
                        i = temp_matches(idx);
                        if strcmp(lower(PlottingColors), 'gradient')
                            temp_idx = find(Temp_range == Temp_obs(i));
                        end
                        %temp_idx = mod(i, 7)+1;
                        eb{idx} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
                        hold on
                        if strcmp(lower(PlottingColors), 'gradient')
                            set(eb{idx}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                        elseif UseDifferentColors
                            set(eb{idx}, 'color', colors(idx,:), 'LineWidth', 1);
                        else
                            set(eb{idx}, 'color', colors(t_idx,:), 'LineWidth', 1);
                        end
                        set(get(get(eb{idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                        if strcmp(lower(PlottingColors), 'gradient')
                            prof{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(temp_idx,:));
                        elseif UseDifferentColors
                            prof{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(idx,:));
                        else
                            prof{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(t_idx,:));
                        end
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
                        xlim([0, 1])
                    else
                        xlabel('Distance from the Anterior Pole (\mum)')
                        xlim([0, max(this.APLengths)])
                    end
                    ylabel('Fluo (AU)')
                    ylim([0, max(max(max(max(MeanFluoMat+StdFluoMat./NumNucMat)))*1.1,1)])
                    %
                    title(FrameProfAx, {'Hunchback Anterior Profiles',...
                        ['T = ', num2str(temperatures(t_idx)), '°C, Nuclear Cycle ',num2str(NC),' (',num2str(-1), ' min)' ]})
                    legend(FrameProfAx, legend_labels(temp_matches), 'Location', 'eastoutside')
                    for i = 1:size(MeanFluoMat, 1)
                        for idx=1:length(temp_matches)
                            j = temp_matches(idx);
                            if UsePhysicalAPLength
                                EmbryoIndex =  this.ProcessedExperiments(j);
                                APLength = this.APLengths(j);
                            else
                                APLength = 1;
                            end
                            use_idx = NumNucMat(i,:,j) >= MinimumNuclearCount;
                            
                            if sum(use_idx) == 0 %| sum(DiffMeanFluoMat(i, use_idx, j) == 0)
                                FrameProfAx.Children(end-(2*(idx-1)+1)).XData = APLength*APbins;
                                FrameProfAx.Children(end-(2*(idx-1)+1)).YData = zeros(1, length(APbins));
                                FrameProfAx.Children(end-(2*(idx-1))).XData = APLength*APbins;
                                FrameProfAx.Children(end-(2*(idx-1))).YData = zeros(1, length(APbins));
                                FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = .1*ones(1, length(APbins));
                                FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta = .1*ones(1, length(APbins));
                                set(eb{idx},'Visible','off'); %'off' or 'on'
                                set(prof{idx},'Visible','off'); %'off' or 'on'
                            else
                                FrameProfAx.Children(end-(2*(idx-1)+1)).YData = MeanFluoMat(i, use_idx, j);
                                FrameProfAx.Children(end-(2*(idx-1)+1)).XData =APLength*MeanAPMat(i, use_idx, j);
                                FrameProfAx.Children(end-(2*(idx-1))).YData = MeanFluoMat(i, use_idx, j);
                                FrameProfAx.Children(end-(2*(idx-1))).XData = APLength*MeanAPMat(i, use_idx, j);
                                FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = StdFluoMat(i, use_idx,j)./NumNucMat(i, use_idx,j);
                                FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta  = StdFluoMat(i, use_idx,j)./NumNucMat(i, use_idx,j);
                                
                                set(eb{idx},'Visible','on'); %'off' or 'on'
                                set(prof{idx},'Visible','on'); %'off' or 'on'
                                
                            end
                        end
                        %try
                        if exist('PlotTitle', 'var')
                            title(FrameProfAx, {PlotTitle,...
                                ['Nuclear Cycle ',num2str(NC),' (',num2str((i-1)*this.time_delta), ' min)' ]})
                        else
                            title(FrameProfAx, ['Nuclear Cycle ',num2str(NC),' (',num2str((i-1)*this.time_delta), ' min)' ])
                        end
                        %end
                        if sum(use_idx) == 0
                            continue
                        end
                        
                        saveas(FrameProfFig,[outdir3, filesep, 'T',strrep(num2str(current_temp), '.', '_'),...
                            'C_NC',num2str(NC),PhysicalAPString, GradString,DiffString, 'FluoProfile_Plot', num2str(i),'.png']);
                        
                    end
                end
            end
            
        end
        
        function PlotFluoProfiles(this, outdir, PlotTitle,varargin)
            UsePhysicalAPLength = false;
            UseDiffProfiles = false;
            UseDifferentColors = true;
            
            x = 1;
            while x <= length(varargin)
                if strcmp(lower(varargin{x}), 'plottitle')
                    PlotTitle = varargin{x+1};
                    x = x+1;
                elseif strcmp(lower(varargin{x}), 'plottingcolors')
                    PlottingColors = varargin{x+1};
                    x = x+1;
                elseif strcmp(lower(varargin{x}), 'useuniformcolor')
                    UseDifferentColors = false;
                elseif strcmp(lower(varargin{x}), 'usediffprofiles')
                    UseDiffProfiles = true;
                elseif strcmp(lower(varargin{x}), 'usephysicalaplength')
                    UsePhysicalAPLength = true;
                end
                x = x+1;
            end
            if ~exist('PlottingColors', 'var')
                PlottingColors = 'default';
            elseif ~strcmp(lower(PlottingColors), 'gradient') & ~strcmp(lower(PlottingColors), 'default')  & ~strcmp(lower(PlottingColors), 'pboc')
                error('Invalid choice of plotting colors. Can use either "default", "pboc", or "gradient".') % change to error
            end
            if ~exist(outdir, 'dir')
                mkdir(outdir)
            end
            if UsePhysicalAPLength
                PhysicalAPString = 'PhysAP';
            else
                PhysicalAPString = '';
            end
            if UseDiffProfiles
                DiffString = 'Diff';
            else
                DiffString = '';
            end
            
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
            
            temp_sps = this.Temp_sps(this.ProcessedExperiments);
            temperatures = flip(unique(this.Temp_sps(this.ProcessedExperiments)));
            APResolution = 1/(size(this.MeanFluoMats{length(this.IncludedNCs)}, 2)-1);
            APbins = 0:APResolution:1;
            
            legend_labels = this.LegendLabels(this.ProcessedExperiments);
            MinimumNuclearCount = this.MinimumNuclearCount;
            numSets = length(this.ProcessedExperiments);
            for nc_idx=1:length(this.IncludedNCs)
                NC = this.IncludedNCs(nc_idx);
                MeanAPMat = this.MeanAPMats{nc_idx};
                if ~UseDiffProfiles
                    MeanFluoMat = this.MeanFluoMats{nc_idx};
                    StdFluoMat = this.StdFluoMats{nc_idx};
                else
                    MeanFluoMat = this.DiffMeanFluoMats{nc_idx};
                    StdFluoMat = this.DiffStdErrorMats{nc_idx};
                end
                NumNucMat = this.NumNucMats{nc_idx};
                
                outdir2 = [outdir, filesep, 'NC', num2str(NC)];
                if ~exist(outdir2, 'dir')
                    mkdir(outdir2)
                end
                outdir3 = [outdir2, filesep, datestr(now, 'yyyymmdd')];
                if ~exist(outdir3, 'dir')
                    mkdir(outdir3)
                end
                close all
                
                
                FrameProfFig = figure(1);
                set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .7, .6]);
                set(gcf,'color','w');
                FrameProfAx = axes(FrameProfFig);
                for i =1:numSets
                    
                    if strcmp(lower(PlottingColors), 'gradient')
                        temp_idx = find(Temp_range == Temp_obs(i));
                    else
                        temp_idx = find(temperatures == temp_sps(i));
                    end
                    %temp_idx = mod(i, 7)+1;
                    eb{i} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
                    hold on
                    if strcmp(lower(PlottingColors), 'gradient')
                        set(eb{i}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                    elseif UseDifferentColors
                        set(eb{i}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                    else
                        set(eb{i}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                    end
                    set(get(get(eb{i}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                    if strcmp(lower(PlottingColors), 'gradient')
                        prof{i} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(temp_idx,:));
                    elseif UseDifferentColors
                        prof{i} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(temp_idx,:));
                    else
                        prof{i} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(temp_idx,:));
                    end
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
                ylabel('Fluo (AU)')
                ylim([0, max(max(max(max(MeanFluoMat+StdFluoMat./NumNucMat)))*1.1,1)])
                if ~UsePhysicalAPLength
                    xlabel('Fraction Embryo Length')
                    xlim([0, 1])
                else
                    xlabel('Distance from the Anterior Pole (\mum)')
                    xlim([0, max(this.APLengths)])
                end
                
                legend(FrameProfAx, legend_labels, 'Location', 'eastoutside')
                for i = 1:size(MeanFluoMat, 1)
                    for j=1:numSets
                        if UsePhysicalAPLength
                            EmbryoIndex =  this.ProcessedExperiments(j);
                            APLength = this.APLengths(j);
                        else
                            APLength = 1;
                        end
                        use_idx = NumNucMat(i,:,j) >= MinimumNuclearCount;
                        if sum(use_idx) == 0 %| sum(DiffMeanFluoMat(i, use_idx, j) == 0)
                            FrameProfAx.Children(end-(2*(j-1)+1)).XData = APLength*APbins;
                            FrameProfAx.Children(end-(2*(j-1)+1)).YData = zeros(1, length(APbins));
                            FrameProfAx.Children(end-(2*(j-1))).XData = APLength*APbins;
                            FrameProfAx.Children(end-(2*(j-1))).YData = zeros(1, length(APbins));
                            FrameProfAx.Children(end-(2*(j-1))).YPositiveDelta = .1*ones(1, length(APbins));
                            FrameProfAx.Children(end-(2*(j-1))).YNegativeDelta = .1*ones(1, length(APbins));
                            set(eb{j},'Visible','off'); %'off' or 'on'
                            set(prof{j},'Visible','off'); %'off' or 'on'
                        else
                            FrameProfAx.Children(end-(2*(j-1)+1)).YData = MeanFluoMat(i, use_idx, j);
                            FrameProfAx.Children(end-(2*(j-1)+1)).XData = APLength*MeanAPMat(i, use_idx, j);
                            FrameProfAx.Children(end-(2*(j-1))).YData = MeanFluoMat(i, use_idx, j);
                            FrameProfAx.Children(end-(2*(j-1))).XData = APLength*MeanAPMat(i, use_idx, j);
                            FrameProfAx.Children(end-(2*(j-1))).YPositiveDelta = StdFluoMat(i, use_idx,j)./NumNucMat(i, use_idx,j);
                            FrameProfAx.Children(end-(2*(j-1))).YNegativeDelta  = StdFluoMat(i, use_idx,j)./NumNucMat(i, use_idx,j);
                            
                            set(eb{j},'Visible','on'); %'off' or 'on'
                            set(prof{j},'Visible','on'); %'off' or 'on'
                            
                        end
                    end
                    %try
                    if exist('PlotTitle', 'var')
                        title(FrameProfAx, {PlotTitle,...
                            ['Nuclear Cycle ',num2str(NC),' (',num2str((i-1)*this.time_delta), ' min)' ]})
                    else
                        title(FrameProfAx, ['Nuclear Cycle ',num2str(NC),' (',num2str((i-1)*this.time_delta), ' min)' ])
                    end
          
                    %end
                    if sum(use_idx) == 0
                        continue
                    end
                    
                    saveas(FrameProfFig,[outdir3, filesep, 'NC',num2str(NC),...
                        PhysicalAPString, GradString,DiffString, 'FluoProfile_Plot', num2str(i),'.png']);
                    
                end
            end
            
        end
        
        function PlotFluoresenceOffsets(this, outdir, varargin)
            
            PlotTitle = 'Fluorescence Offsets halfway through NC';
            PlottingColors = 'default';
            UseFullTrace = false;
            %for
            x = 1;
            while x <= length(varargin)
                if strcmp(lower(varargin{x}), 'plottitle')
                    PlotTitle = varargin{x+1};
                    x = x+1;
                elseif strcmp(lower(varargin{x}), 'plottingcolors')
                    PlottingColors = varargin{x+1};
                    x =x+1;
                elseif strcmp(lower(varargin{x}), 'usefulltrace')
                    UseFullTrace = true;
                end
                x = x+1;
            end
            
            
            % PlotFluoOffsets
            
            if ~strcmp(lower(PlottingColors), 'gradient') & ~strcmp(lower(PlottingColors), 'default')  & ~strcmp(lower(PlottingColors), 'pboc')
                disp('Invalid choice of plotting colors. Can use either "default", "pboc", or "gradient".') % change to error
            end
            if ~exist(outdir, 'dir')
                mkdir(outdir)
            end
            
            
            
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
            if ~UseFullTrace
                temperatures = this.Temp_obs(this.ProcessedExperiments);
                temp_sps = this.Temp_sps(this.ProcessedExperiments);
                for nc_idx = 1:length(this.IncludedNCs)
                    NC = this.IncludedNCs(nc_idx);
                    outdir2 = [outdir, filesep, 'NC', num2str(NC)];
                    if ~exist(outdir2, 'dir')
                        mkdir(outdir2)
                    end
                    outdir3 = [outdir2, filesep, datestr(now, 'yyyymmdd')];
                    if ~exist(outdir3, 'dir')
                        mkdir(outdir3)
                    end
                    FluoOffsets = this.FluoOffsets{nc_idx};
                    FluoErrors= this.FluoOffsetStdErrors{nc_idx};
                    OffsetTimePoints = zeros(1,size(FluoOffsets, 2));
                    FluoOffsetsToPlot = zeros(1,size(FluoOffsets, 2));
                    FluoErrorsToPlot = zeros(1,size(FluoOffsets, 2));
                    close all
                    FrameOffsetFig = figure(1);
                    set(FrameOffsetFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .5]);
                    set(gcf,'color','w');
                    FrameOffsetAx = axes(FrameOffsetFig);
                    temperatures_labeled = zeros(1, length(temperatures));
                    temp_sp_set = flip(unique(temp_sps));
                    for s=1:size(FluoOffsets, 2)
                        if strcmp(PlottingColors, 'gradient')
                            temp_idx = find(Temp_range == temperatures(s));
                        else
                            temp_idx = find(temp_sp_set == temp_sps(s));
                        end
                        SingleSetFluoOffsets = FluoOffsets(:,s);
                        OffsetTimePoints(s) = floor(max(find(~isnan(SingleSetFluoOffsets)))/2);
                        FluoOffsetsToPlot(s) = FluoOffsets(OffsetTimePoints(s),s);
                        FluoErrorsToPlot(s) = FluoErrors(OffsetTimePoints(s),s);
                        
                        
                        
                        
                        
                        eb{s} = errorbar(temperatures(s),FluoOffsetsToPlot(s), FluoErrorsToPlot(s), 'vertical', 'LineStyle', 'none');
                        hold on
                        set(eb{s}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                        set(get(get(eb{s}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                        prof{s} = scatter(temperatures(s), FluoOffsetsToPlot(s), 50, 'MarkerFaceColor', colors(temp_idx,:),...
                            'MarkerEdgeColor', colors(temp_idx,:));
                    end
                    
                    %end
                    
                    hold off
                    xlabel('Temperature (°C)', 'FontSize', 16)
                    ylabel('Nuclear GFP concentration (AU)', 'FontSize', 16)
                    ylim([0, max(max(FluoOffsetsToPlot+FluoErrorsToPlot))*1.1])
                    xlim([15, 29])
                    FrameOffsetAx.XAxis.FontSize = 16;
                    FrameOffsetAx.YAxis.FontSize = 16;
                    
                    
                    title(FrameOffsetAx, {PlotTitle,...
                        ['Nuclear Cycle ',num2str(NC) ]})
                    
                    saveas(FrameOffsetFig,[outdir3, filesep, 'NC',num2str(NC),GradString, 'Offsets.png']);
                    
                end
            else
                temperatures = this.Temp_obs(this.ProcessedExperiments);
                temp_sps = this.Temp_sps(this.ProcessedExperiments);
                temp_sp_set = flip(unique(temp_sps));
                for nc_idx = 1:length(this.IncludedNCs)
                    profs2 = cell(1, length(temp_sp_set));
                    NC = this.IncludedNCs(nc_idx);
                    outdir2 = [outdir, filesep, 'NC', num2str(NC)];
                    if ~exist(outdir2, 'dir')
                        mkdir(outdir2)
                    end
                    outdir3 = [outdir2, filesep, datestr(now, 'yyyymmdd')];
                    if ~exist(outdir3, 'dir')
                        mkdir(outdir3)
                    end
                    FluoOffsets = this.FluoOffsets{nc_idx};
                    FluoErrors= this.FluoOffsetStdErrors{nc_idx};
                    OffsetTimePoints = zeros(1,size(FluoOffsets, 2));
                    FluoOffsetsToPlot = zeros(1,size(FluoOffsets, 2));
                    FluoErrorsToPlot = zeros(1,size(FluoOffsets, 2));
                    close all
                    FrameOffsetFig = figure(1);
                    set(FrameOffsetFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .5]);
                    set(gcf,'color','w');
                    FrameOffsetAx = axes(FrameOffsetFig);
                    temperatures_labeled = zeros(1, length(temperatures));
                    if ~strcmp(PlottingColors, 'gradient')
                        plotted_temps = zeros(1, length(temp_sp_set));
                    end
                    time_vector = 0:(size(FluoOffsets, 1)-1);
                    time_vector = time_vector*this.time_delta;
                    for s=1:size(FluoOffsets, 2)
                        if strcmp(PlottingColors, 'gradient')
                            temp_idx = find(Temp_range == temperatures(s));
                        else
                            temp_idx = find(temp_sp_set == temp_sps(s));
                        end
                        
                        
                        
                        
                        
                        eb{s} = errorbar(time_vector,FluoOffsets(:,s), FluoErrors(:,s), 'vertical', 'LineStyle', 'none');
                        hold on
                        set(eb{s}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                        set(get(get(eb{s}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                        prof{s} = scatter(time_vector,FluoOffsets(:,s), 50, 'MarkerFaceColor', colors(temp_idx,:),...
                            'MarkerEdgeColor', colors(temp_idx,:));
                        if ~strcmp(PlottingColors, 'gradient')
                            if plotted_temps(temp_idx) == 0
                                plotted_temps(temp_idx) = s;
                            end
                        end
                    end
                    
                    %end
                    
                    hold off
                    xlabel('Temperature (°C)', 'FontSize', 16)
                    ylabel('Nuclear GFP concentration (AU)', 'FontSize', 16)
                    ylim([0, max(max(FluoOffsets+FluoErrors))*1.1])
                    xlim([0, max(time_vector)+1])
                    FrameOffsetAx.XAxis.FontSize = 16;
                    FrameOffsetAx.YAxis.FontSize = 16;
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
                    else
                        legend_labels = {};
                        for i= 1:length(plotted_temps)
                            if plotted_temps(i) > 0
                                current_temp = temp_sps(plotted_temps(i));
                                legend_labels{1, length(legend_labels)+1} =[ num2str(current_temp),'°C'];
                                profs2{1, i} = prof{plotted_temps(i)};
                            end
                        end
                        hlegend = legend([profs2{:}], legend_labels, 'Location', 'northeast',...
                            'FontSize', 12);
                    end
                    
                    title(FrameOffsetAx, {PlotTitle,...
                        ['Nuclear Cycle ',num2str(NC) ]})
                    
                    saveas(FrameOffsetFig,[outdir3, filesep, 'NC',num2str(NC),GradString, 'OffsetTraces.png']);
                    
                end
            end
        end
        
        function [AllNCMaxPositions, HalfMaxPositions, AllNCMaxSetTimes, HalfMaxVectors,SetsContainMax, AllMaxFluos] = ...
                CalculateHbProfileMaximumPositions(this)
            for nc_idx=1:length(this.IncludedNCs)
                NC = this.IncludedNCs(nc_idx);
                numTimePoints = size(this.MeanFluoMats{nc_idx}, 1);
                numAPBins = size(this.MeanFluoMats{nc_idx}, 2);
                numSets = size(this.MeanFluoMats{nc_idx}, 3);
                ProfileRegions = this.Regions(1, this.ProcessedExperiments);
                MaxPositionBins = NaN(numSets, numTimePoints);
                MaxPositionSets = NaN(1, numSets);
                MeanFluoMat = this.MeanFluoMats{nc_idx};
                DiffMeanFluoMat = this.DiffMeanFluoMats{nc_idx};
                SetIsAnterior = zeros(1, numSets);
                APResolution = this.Experiments{this.ProcessedExperiments(1)}.APResolution;
                for i=1:numSets
                    for j = 1:numTimePoints
                        % FIND MINIMUM OF DERIVATIVE
                        if strcmp(lower(ProfileRegions{i}), 'anterior')
                            SetIsAnterior(i) = 1;
                        end
                        if sum(~isnan(DiffMeanFluoMat(j,:,i))) > 1
                            mid_bin = find(diff(DiffMeanFluoMat(j,:,i)) == nanmin(diff(DiffMeanFluoMat(j,:,i))), 1);
                            if ~isempty(mid_bin)
                                MaxPositionBins(i,j) = mid_bin -2;
                            end
                        end
                        
                    end
                    MaxPositionSets(i) =  floor(nanmedian(MaxPositionBins(i,:)));
                    
                end
                
                %% Find Hunchback boundary
                HalfMaxPos = NaN(numSets, numTimePoints);
                for i=1:numSets
                    MaxPos = MaxPositionSets(i);
                    if ~SetIsAnterior(i) | isnan(MaxPos)
                        continue
                    end
                    for j = 1:numTimePoints
                        AllowedBins = find(~isnan(DiffMeanFluoMat(j,MaxPos:end,i)));
                        if length(AllowedBins) >= 4
                            interp_points = MaxPos:0.01:(max(AllowedBins)+MaxPos-1);
                            vq = interp1(AllowedBins+MaxPos-1, DiffMeanFluoMat(j,AllowedBins+MaxPos-1,i), interp_points);
                            halfMax_point = find(vq < DiffMeanFluoMat(j,MaxPos,i)/2, 1);
                            if ~isempty(halfMax_point)
                                HalfMaxPos(i,j) = (interp_points(halfMax_point)-1)*APResolution;
                            end
                        end
                    end
                end
                
                
                
                %% FIGURE OUT TIME TO MAXIMUM FLUO
                HalfMaxVector = NaN(1, numSets);
                MaxPrefixFluos = NaN(1, numSets);
                MaxSetTimes = NaN(1, numSets);
                SetContainsMax = ones(1, numSets);
                ValidBins = uint16((0.325:APResolution:0.425)/APResolution + 1);
                for i =1:numSets
                    if ~SetIsAnterior(i) | isnan(MaxPositionSets(i))
                        continue
                    end
                    if sum(find(~isnan(MeanFluoMat(:,MaxPositionSets(i),i)))) > 0
                        [MaxPrefixFluos(i), MaxSetTimes(i)] =...
                            max(smoothdata(MeanFluoMat(:,MaxPositionSets(i),i), 'gaussian', 5));
                    end
                    ValidTimePointsInSet = find(~isnan(MeanFluoMat(:,MaxPositionSets(i),i)));
                    if max(ValidTimePointsInSet) == MaxSetTimes(i)
                        SetContainsMax(i) = 0;
                    end
                    if ~isnan(MaxSetTimes(i))
                        HalfMaxVector(i) = HalfMaxPos(i, MaxSetTimes(i));
                    end
                end
                
                MaxSetTimes = MaxSetTimes-1;
                AllNCMaxPositions{nc_idx} = MaxPositionSets;
                HalfMaxPositions{nc_idx} = HalfMaxPos;
                AllNCMaxSetTimes{nc_idx} = MaxSetTimes;
                HalfMaxVectors{nc_idx} = HalfMaxVector;
                SetsContainMax{nc_idx} = SetContainsMax;
                AllMaxFluos{nc_idx} = MaxPrefixFluos;
                
            end
        end
        function PlotDivisionCycleTimes(this, NCDivisionInfo,DivisionStdErrorInfo,...
                TemperatureInfo, FitParams, ActivationEnergies, outdir, varargin)
            
            x = 1;
            while x < length(varargin)
                if strcmp(lower(varargin{x}), 'plottingcolors')
                    PlottingColors = varargin{x+1};
                    x = x+1;
                end
                x = x+ 1;
            end
            if ~exist('PlottingColors', 'var')
                PlottingColors = 'default';
            elseif ~strcmp(lower(PlottingColors), 'default')  & ~strcmp(lower(PlottingColors), 'pboc')
                disp('Invalid choice of plotting colors. Can use either "default" or "pboc".') % change to error
            end
            if strcmp(lower(PlottingColors), 'default')
                [~, colors] = getColorPalettes();
            elseif strcmp(lower(PlottingColors), 'pboc')
                [colors, ~] = getColorPalettes();
            end
            outdir2 = [outdir, filesep, 'DivisionCycleTimes'];
            if ~exist(outdir2, 'dir')
                mkdir(outdir2)
            end
            outdir3 = [outdir2, filesep, datestr(now, 'yyyymmdd')];
            if ~exist(outdir3, 'dir')
                mkdir(outdir3)
            end
            R = 8.314*10^(-3); % kJ * K^(-1)*mol^(-1)
            
            %% CalculatePlotLimits using info from all times and temperatures
            AllTemperatures = [];
            AllNCTimesPlusErrors = [];
            AllAlphas = [];
            AllFitIntercepts = [];
            AllFitYValues = [];
            for NC=10:13
                nc_idx = NC-9;
                AllNCTimesPlusErrors = [AllNCTimesPlusErrors, NCDivisionInfo{nc_idx}+DivisionStdErrorInfo{nc_idx}];
                AllTemperatures = [AllTemperatures, TemperatureInfo{nc_idx}];
                for j = 1:length( NCDivisionInfo{nc_idx})
                    AllAlphas = [AllAlphas, FitParams{nc_idx}(1)];
                    AllFitIntercepts = [AllFitIntercepts, FitParams{nc_idx}(2)];
                    AllFitYValues = [AllFitYValues,  FitParams{nc_idx}(2)- FitParams{nc_idx}(1)./(TemperatureInfo{nc_idx}(j)+273)];
                end
            end
            TemperatureVector = -1./(AllTemperatures + 273);
            
            %% Calculate x and y limits for both subplots
            
            Subplot1DataYmin = min(AllNCTimesPlusErrors);
            Subplot1DataYmax = max(AllNCTimesPlusErrors);
            if Subplot1DataYmin/5 == floor(Subplot1DataYmin/5)
                Plot1Ymin = Subplot1DataYmin-5;
            else
                Plot1Ymin =  floor(Subplot1DataYmin/5)*5;
            end
            if Subplot1DataYmax/5 == floor(Subplot1DataYmax/5)
                Plot1Ymax = Subplot1DataYmax+5;
            else
                Plot1Ymax =  ceil(Subplot1DataYmax/5)*5;
            end
            
            Subplot2DataXmin = min(TemperatureVector);
            Subplot2DataXmax = max(TemperatureVector);
            Subplot2Xspan = Subplot2DataXmax-Subplot2DataXmin;
            Plot2Xmin = Subplot2DataXmin - Subplot2Xspan*.05;
            Plot2Xmax = Subplot2DataXmax + Subplot2Xspan*.05;
            
            Subplot2DataYmin = min(min(log(AllNCTimesPlusErrors)), min(AllFitYValues));
            Subplot2DataYmax = max(max(log(AllNCTimesPlusErrors)), max(AllFitYValues));
            Subplot2Yspan = Subplot2DataYmax-Subplot2DataYmin;
            Plot2Ymin = Subplot2DataYmin - Subplot2Yspan*.05;
            Plot2Ymax = Subplot2DataYmax + Subplot2Yspan*.2;
            
            %% Set up both subplots
            close all
            
            TimeMaximaFig = figure(1);
            set(TimeMaximaFig,'units', 'normalized', 'position',[0.01, 0.05, .8, .6]);
            set(gcf,'color','w');
            TimeMaximaAx1 = subplot(1, 2, 1);
            xlabel('Temperature (°C)', 'FontSize', 16)
            ylabel('time (min)')
            xlim([15, 29])
            ylim([Plot1Ymin, Plot1Ymax])
            TimeMaximaAx1.XAxis.FontSize = 16;
            TimeMaximaAx1.YAxis.FontSize = 16;
            TimeMaximaAx2 = subplot(1, 2, 2);
            xlabel('log(A) + E_a/(RT)', 'FontSize', 16)
            ylabel('log(time)', 'FontSize', 16)
            TimeMaximaAx2.XAxis.FontSize = 16;
            TimeMaximaAx2.YAxis.FontSize = 16;
            ylim([Plot2Ymin, Plot2Ymax])
            xlim([Plot2Xmin, Plot2Xmax])
            % %sgtitle( 'Measured time maximum Hunchback concentration', 'FontSize', 20)
            %% Add cycle information to plot
            data_legend_labels = {};
            fit_legend_labels = {};
            p1s = {};
            eb1s = {};
            p2s = {};
            eb2s = {};
            p1_idx = 0;
            eb1_idx = 0;
            p2_idx = 0;
            eb2_idx = 0;
            data_legend_handles = [];
            fit_legend_handles = [];
            for nc_idx = 1:4
                NC = nc_idx + 9;
                if isempty(NCDivisionInfo{nc_idx})
                    continue
                else
                    p1_idx = p1_idx + 1;
                    eb1_idx = eb1_idx + 1;
                    p2_idx = p2_idx + 1;
                    eb2_idx = eb2_idx + 1;
                end
                
                % if isempty(nc_idx)
                %     error('NC 14 is not included in these data sets.')
                % end
                
                TemperatureVector = -1./(TemperatureInfo{nc_idx}+273);
                
                
                
                hold(TimeMaximaAx1, 'on')
                eb1s{eb1_idx} = errorbar(TimeMaximaAx1,  TemperatureInfo{nc_idx}, NCDivisionInfo{nc_idx},...
                    DivisionStdErrorInfo{nc_idx}, 'vertical', 'LineStyle', 'none');
                set(eb1s{eb1_idx}, 'color', colors(nc_idx,:), 'LineWidth', 1);
                set(get(get(eb1s{eb1_idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                hold on
                p1s{p1_idx} = scatter(TimeMaximaAx1, TemperatureInfo{nc_idx}, NCDivisionInfo{nc_idx},50,...
                    'MarkerFaceColor', colors(nc_idx,:), 'MarkerEdgeColor', colors(nc_idx,:));
                
                
                %
                
                % p5 =  scatter(TimeMaximaAx2, TemperatureVector, log(MaxSetTimes(SetContainsMax == 1)),50,...
                %     'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', colors(1,:));
                
                % eb{1} = errorbar(TempNC10s, NC10s, StdNC10s, 'vertical', 'LineStyle', 'none');
                % set(eb{1}, 'color', colors(3,:), 'LineWidth', 1);
                % set(get(get(eb{1}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                hold(TimeMaximaAx2, 'on')
                eb2s{eb2_idx} = errorbar(TimeMaximaAx2, -1./(TemperatureInfo{nc_idx}+273),...
                    log(NCDivisionInfo{nc_idx}), DivisionStdErrorInfo{nc_idx}./NCDivisionInfo{nc_idx}, 'vertical', 'LineStyle', 'none');
                set(eb2s{eb2_idx}, 'color', colors(nc_idx,:), 'LineWidth', 1);
                set(get(get(eb2s{eb2_idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
                
                hold on
                p2s{p2_idx} = scatter(TimeMaximaAx2, -1./(TemperatureInfo{nc_idx}+273), log(NCDivisionInfo{nc_idx}),50,...
                    'MarkerFaceColor', colors(nc_idx,:), 'MarkerEdgeColor', colors(nc_idx,:));
                data_legend_handles= [data_legend_handles p2s{p2_idx}];
                data_legend_labels{length(data_legend_labels) + 1} = ['Cycle ', num2str(NC)];
                if ~isnan(FitParams{nc_idx}(1)) & ~isnan(FitParams{nc_idx}(2))
                    p2_idx = p2_idx + 1;
                    p2s{p2_idx}= plot(TimeMaximaAx2, [Subplot2DataXmin, Subplot2DataXmax], FitParams{nc_idx}(2) + ActivationEnergies(nc_idx)/R*[Subplot2DataXmin, Subplot2DataXmax],...
                        '-', 'Color', colors(nc_idx,:));
                    fit_legend_handles = [fit_legend_handles p2s{p2_idx}];
                    fit_legend_labels{length(fit_legend_labels) + 1} = ['Arrehnius fit, E_a = ', num2str(round(ActivationEnergies(nc_idx))),' kJ/mol'];
                    
                end
                
                
                
            end
            %%
            
            legend_handles = [data_legend_handles fit_legend_handles];
            legend_labels = [data_legend_labels, fit_legend_labels];
            
            hlegend = legend(legend_handles, legend_labels, 'Location', 'northeast',...
                'FontSize', 12);
            % hlegend = legend([p6 p7 p8 p9 p11 p12 p13], {'Cycle 10', 'Cycle 11', 'Cycle 12','Cycle 13',...
            %     'Arrehnius fit, E_a = -64 kJ/mol', 'Arrehnius fit, E_a = -58 kJ/mol',...
            %     'Arrehnius fit, E_a = -66 kJ/mol'}, 'Location', 'northeast',...
            %     'FontSize', 12);
            
            hlegend.NumColumns = 2;
            saveas(TimeMaximaFig,[outdir3, filesep, 'AllNuclearCycleTimingPlots.png']);
        end
        
    end
end

