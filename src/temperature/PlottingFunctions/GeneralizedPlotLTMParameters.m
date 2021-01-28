function GeneralizedPlotLTMParameters(this, parameter, outdir, varargin)
%%

% PlotTitle, PlottingColors, UseDifferentColors,
% UseDiffProfiles, UsePhysicalAPLength
UsePhysicalAPLength = false;
UseDifferentColors = true;
UseLines = true;
R2bound = .9;
SkipParamsVsAP = false;
SkipSingleTempParamsVsAP = false;
SkipParamsVsTemp = false;
SkipBinnedParamsVsAP = false;
SkipBinnedParamsVsTemp = false;
SkipAPSubplots = false;
SkipAPBinnedSubplots = false;
IncludeFits = true;

x = 1;
while x <= length(varargin)
    if strcmp(lower(varargin{x}), 'plottitle')
        PlotTitle = varargin{x+1};
        x = x+1;
    elseif strcmp(lower(varargin{x}), 'plottingcolors')
        PlottingColors = varargin{x+1};
        x = x+1;
    elseif strcmp(lower(varargin{x}), 'usephysicalaplength')
        UsePhysicalAPLength = true;
    elseif strcmp(lower(varargin{x}), 'noline')
        UseLines = false;
    elseif strcmpi(varargin{x}, 'SkipParamsVsAP')
        SkipParamsVsAP = true;
    elseif strcmpi(varargin{x}, 'SkipSingleTempParamsVsAP')
        SkipSingleTempParamsVsAP = true;
    elseif strcmpi(varargin{x}, 'SkipParamsVsTemp')
        SkipParamsVsTemp = true;
    elseif strcmpi(varargin{x}, 'SkipBinnedParamsVsAP')
        SkipBinnedParamsVsAP = true;
    elseif strcmpi(varargin{x}, 'SkipBinnedParamsVsTemp')
        SkipBinnedParamsVsTemp = true;
    elseif strcmpi(varargin{x}, 'SkipAPSubplots')
        SkipAPSubplots = true;
    elseif strcmpi(varargin{x}, 'SkipAPBinnedSubplots')
        SkipAPBinnedSubplots = true;
    elseif strcmpi(varargin{x}, 'ExcludeFits')
        IncludeFits = false;
    elseif strcmp(lower(varargin{x}), 'tracetype')
        TraceType = lower(varargin{x+1});
        x = x+1;
    end
    x = x+1;
end


if ~exist('PlottingColors', 'var')
    PlottingColors = 'default';
elseif ~strcmpi(PlottingColors, 'gradient') &~strcmp(lower(PlottingColors), 'default')  & ~strcmp(lower(PlottingColors), 'pboc')
    error('Invalid choice of plotting colors. Can use either "default", "pboc", or "gradient".') % change to error
end
if ~exist('TraceType', 'var')
    TraceType = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'Fluo3D')
    TraceType = 'Fluo3D';
elseif strcmpi(TraceType, 'Fluo')
    TraceType = 'Fluo';
elseif strcmpi(TraceType, 'AnaphaseAligned')
    TraceType = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'AnaphaseAligned3D')
    TraceType = 'AnaphaseAligned3D';
elseif strcmpi(TraceType, 'Tbinned')
    TraceType = 'Tbinned';
elseif strcmpi(TraceType, 'Tbinned3D')
    TraceType = 'Tbinned3D';
else
    error('Invalid choice of trace type. Can use either "fluo", "fluo3d", "anaphasealigned", or "anaphasealigned3d".') % change to error
end

%%
if ~strcmpi(parameter, 'TimeOns') & ~strcmpi(parameter, 'TranscriptionWindows') & ...
        ~strcmpi(parameter, 'ElongationTimes') & ~strcmpi(parameter, 'ElongationRates') & ...
        ~strcmpi(parameter, 'LoadingRates')
    IncludeFits = false;
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

R = 8.314*10^(-3); % kJ * K^(-1)*mol^(-1)
temperatures = flip(unique(this.Temp_sps(this.ProcessedExperiments)));
APResolution = this.Experiments{1}.APResolution;
APbins = 0:APResolution:1;
NumSets = length(this.ProcessedExperiments);
Nsigfigs = 3;
legend_labels = this.LegendLabels(this.ProcessedExperiments);
% MarkerStyles = {'o', 'd', 's', '>', '^', '*', 'x', 'p'};
% NumTemperatures = length(temperatures);
% TempMatches = cell(1, NumTemperatures);
% UseSet = ismember(1:NumSets, this.ProcessedExperiments);
% for t_index = 1:NumTemperatures
%     TempMatches{t_index} = find((this.Temp_sps == temperatures(t_index)) & UseSet);
% end


%% Load relevant parameters into memory
[PlottedParams, PlottedParamSEs,R2s, ylab,OutputString,GlobalPlotYmax,GlobalPlotYmin,LogPlotYmin] = ...
    getPlottingVariables(this, parameter,  TraceType, R2bound);
[BinnedParams, BinnedSEParams, Counts, ParamTemperatures, ParamSETemperatures] = ...
    getBinnedPlottingVariables(this, PlottedParams, PlottedParamSEs,R2s, R2bound);
%%

outdir2 = [outdir,filesep,OutputString];
if ~exist(outdir2, 'dir')
    mkdir(outdir2)
end
outdir3 = [outdir2, filesep, datestr(now, 'yyyymmdd')];
if ~exist(outdir3, 'dir')
    mkdir(outdir3)
end

if ~SkipParamsVsAP
    for nc_idx=1:length(this.IncludedNCs)
        NC = this.IncludedNCs(nc_idx);
        clear prof eb
        plotted_temps = zeros(1, length(temperatures));
        
        
        
        % Prepare Traces for plotting
        
        NCMaxParams = NaN(1, length(Temp_sp));
        AllNCParams = NaN(NumSets, length(APbins));
        AllNCParamSEs = NaN(NumSets, length(APbins));
        AllR2s = NaN(NumSets, length(APbins));
        for i=1:NumSets
            idx = this.ProcessedExperiments(i);
            SetParams = PlottedParams(idx,:,NC-8).';
            SetSEParams = PlottedParamSEs(idx,:,NC-8).';
            SetR2s = R2s(idx,:,NC-8).';
            if ~all(isnan(SetSEParams))
                IncludedBins = find(~isnan(SetParams) & (SetR2s > R2bound) & (SetParams./SetSEParams >= 1)) ;
            else
                IncludedBins = find(~isnan(SetParams) & (SetR2s > R2bound)) ;
            end
            if ~isempty(IncludedBins)
                AllNCParams(idx,:) = SetParams;
                AllNCParamSEs(idx,:) = SetSEParams;
                TempSEParams = SetSEParams;
                TempSEParams(isnan(TempSEParams)) = 0;
                NCMaxParams(idx) = max(SetParams+TempSEParams);
                AllR2s(idx,:) = SetR2s;
            end
        end
        
        if all(isnan(NCMaxParams))
            continue
        end
        close all
        
        
        FrameProfFig = figure(1);
        set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .5, .4]);
        set(gcf,'color','w');
        FrameProfAx = axes(FrameProfFig);
        for idx =1:length(Temp_sp)
            if strcmp(lower(PlottingColors), 'gradient')
                temp_idx = find(abs(Temp_range-Temp_obs(idx)) == min(abs(Temp_range-Temp_obs(idx))));
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
            if UseLines
                if strcmp(lower(PlottingColors), 'gradient')
                    prof{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(temp_idx,:));
                elseif UseDifferentColors
                    prof{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(temp_idx,:));
                else
                    prof{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(temp_idx,:));
                end
            else
                if strcmp(lower(PlottingColors), 'gradient')
                    prof{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(temp_idx,:), 'linestyle', 'none');
                elseif UseDifferentColors
                    prof{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(temp_idx,:), 'linestyle', 'none');
                else
                    prof{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(temp_idx,:), 'linestyle', 'none');
                end
            end
            
            set(eb{idx},'Visible','off'); %'off' or 'on'
            set(prof{idx},'Visible','off'); %'off' or 'on'
        end
        
        
        
        hold off
        if ~UsePhysicalAPLength
            xlabel('Fraction Embryo Length')
            xlim([0, 1])
        else
            xlabel('Distance from the Anterior Pole (\mum)')
            xlim([0, max(this.APLengths)])
        end
        
        ylabel(ylab)
        ylim([GlobalPlotYmin, GlobalPlotYmax])
        
        %ylim([min(0, min(min(AllTimeOns+AllStdErrorTimeOns))*1.2), max(max(max(AllTimeOns+AllStdErrorTimeOns))*1.2, 1)])
        %
        
        title(FrameProfAx, {'',...
            ['Nuclear Cycle ',num2str(NC)]})
        
        
        
        %legend(FrameProfAx, legend_labels, 'Location', 'eastoutside')
        HasDataPlotted = false;
        for idx=1:length(Temp_sp)
            
            if UsePhysicalAPLength
                EmbryoIndex =  this.ProcessedExperiments(idx);
                APLength = this.APLengths(idx);
            else
                APLength = 1;
            end
            if ~all(isnan(AllNCParamSEs(idx, :)))
                use_idx = ~isnan(AllNCParams(idx,:)) &  (AllNCParams(idx, :)./AllNCParamSEs(idx, :) >= 1) & ...
                    (AllR2s(idx,:)>= R2bound);
            else
                use_idx = ~isnan(AllNCParams(idx,:)) & ...
                    (AllR2s(idx,:)>= R2bound);
            end
            
            if sum(use_idx) == 0 %| sum(DiffMeanFluoMat(i, use_idx, j) == 0)
                FrameProfAx.Children(end-(2*(idx-1)+1)).XData = APLength*APbins;
                FrameProfAx.Children(end-(2*(idx-1)+1)).YData = zeros(1, length(APbins));
                FrameProfAx.Children(end-(2*(idx-1))).XData = APLength*APbins;
                FrameProfAx.Children(end-(2*(idx-1))).YData = zeros(1, length(APbins));
                FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = .1*ones(1, length(APbins));
                FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta = .1*ones(1, length(APbins));
                set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','off'); %'off' or 'on'
                set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','off'); %'off' or 'on'
            else
                HasDataPlotted = true;
                FrameProfAx.Children(end-(2*(idx-1)+1)).YData = AllNCParams(idx,use_idx);
                FrameProfAx.Children(end-(2*(idx-1)+1)).XData = APbins(use_idx);
                FrameProfAx.Children(end-(2*(idx-1))).YData = AllNCParams(idx,use_idx);
                FrameProfAx.Children(end-(2*(idx-1))).XData = APbins(use_idx);
                FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = AllNCParamSEs(idx,use_idx);
                FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta  = AllNCParamSEs(idx,use_idx);
                
                set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','on'); %'off' or 'on'
                set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','on'); %'off' or 'on'
                if ~strcmp(lower(PlottingColors), 'gradient')
                    temp_idx = find(temperatures == Temp_sp(idx));
                    if plotted_temps(temp_idx) == 0
                        plotted_temps(temp_idx) = idx;
                    end
                end
                
            end
        end
        %try
        if exist('PlotTitle', 'var')
            
            
            title(FrameProfAx, {PlotTitle,...
                ['Nuclear Cycle ',num2str(NC)]})
            
        else
            
            title(FrameProfAx,  ['Nuclear Cycle ',num2str(NC)])
            
            
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
        else
            legend_labels = {};
            profs2 = {};
            for i= 1:length(plotted_temps)
                if plotted_temps(i) > 0
                    current_temp = Temp_sp(plotted_temps(i));
                    legend_labels{1, length(legend_labels)+1} =[ num2str(current_temp),'°C'];
                    profs2{1, length(profs2)+1} = prof{plotted_temps(i)};
                end
            end
            hlegend = legend([profs2{:}], legend_labels, 'Location', 'northeast',...
                'FontSize', 10);
        end
        
        saveas(FrameProfFig,[outdir3, filesep,...
            TraceType, PhysicalAPString, GradString,'_', OutputString, '_NC',num2str(NC),'.png']);
        
        
    end
end
%% Single Temperature Plots
UseSet = zeros(1, length(this.ExperimentPrefixes));
for i = 1:length(UseSet)
    if ismember(i, this.ProcessedExperiments)
        UseSet(i) = 1;
    end
end
legend_labels2 = this.LegendLabels;
if ~SkipSingleTempParamsVsAP
    for nc_idx=1:length(this.IncludedNCs)
        
        NC = this.IncludedNCs(nc_idx);
        for CurrentTemperature = temperatures
            clear prof eb
            TempMatches = find((this.Temp_sps == CurrentTemperature) & UseSet == 1);
            plotted_temps = zeros(1, length(temperatures));
            outdir2 = [outdir,filesep,OutputString];
            if ~exist(outdir2, 'dir')
                mkdir(outdir2)
            end
            outdir3 = [outdir2,filesep,'T', strrep(num2str(CurrentTemperature), '.', '_'),'C'];
            if ~exist(outdir3, 'dir')
                mkdir(outdir3)
            end
            outdir4 = [outdir3,filesep, datestr(now, 'yyyymmdd')];
            if ~exist(outdir4, 'dir')
                mkdir(outdir4)
            end
            % Prepare Traces for plotting
            
            % Prepare Traces for plotting
            
            NCMaxParams = NaN(1, length(TempMatches));
            AllNCParams = NaN(length(TempMatches), length(APbins));
            AllNCParamSEs = NaN(length(TempMatches), length(APbins));
            AllR2s = NaN(length(TempMatches), length(APbins));
            for i=1:length(TempMatches)
                idx = TempMatches(i);
                SetParams = PlottedParams(idx,:,NC-8).';
                SetSEParams = PlottedParamSEs(idx,:,NC-8).';
                SetR2s = R2s(idx,:,NC-8).';
                if ~all(isnan(SetSEParams))
                    IncludedBins = find(~isnan(SetParams) & (SetR2s > R2bound) & (SetParams./SetSEParams >= 1));
                else
                    IncludedBins = find(~isnan(SetParams) & (SetR2s > R2bound));
                end
                if ~isempty(IncludedBins)
                    AllNCParams(i,:) = SetParams;
                    AllNCParamSEs(i,:) = SetSEParams;
                    TempSEParams = SetSEParams;
                    TempSEParams(isnan(TempSEParams)) = 0;
                    NCMaxParams(i) = max(SetParams +TempSEParams);
                    AllR2s(i,:) = SetR2s;
                end
            end
            
            if all(isnan(NCMaxParams))
                continue
            end
            
            close all
            
            
            FrameProfFig = figure(1);
            set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .4]);
            set(gcf,'color','w');
            FrameProfAx = axes(FrameProfFig);
            for idx =1:length(TempMatches)
                if strcmp(lower(PlottingColors), 'gradient')
                    temp_idx = find(abs(Temp_range-Temp_obs(idx)) == min(abs(Temp_range-Temp_obs(idx))));
                else
                    temp_idx = idx;
                end
                %temp_idx = mod(i, 7)+1;
                eb{idx} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
                hold on
                
                set(eb{idx}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                
                set(get(get(eb{idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                if UseLines
                    
                    prof{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(temp_idx,:));
                    
                else
                    
                    prof{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(temp_idx,:), 'linestyle', 'none');
                    
                end
                
                set(eb{idx},'Visible','off'); %'off' or 'on'
                set(prof{idx},'Visible','off'); %'off' or 'on'
            end
            
            
            
            hold off
            if ~UsePhysicalAPLength
                xlabel('Fraction Embryo Length')
                xlim([0, 1])
            else
                xlabel('Distance from the Anterior Pole (\mum)')
                xlim([0, max(this.APLengths)])
            end
            ylabel(ylab)
            ylim([GlobalPlotYmin, GlobalPlotYmax])
            if strcmp(lower(TraceType), 'anaphasealigned') | strcmp(lower(TraceType), 'anaphasealigned3d')
                title(FrameProfAx, {'',...
                    ['Nuclear Cycle ',num2str(NC) ]})
            else
                title(FrameProfAx, {'',...
                    ['Nuclear Cycle ',num2str(NC),', ',num2str(-1), ' min' ]})
            end
            
            
            %legend(FrameProfAx, legend_labels, 'Location', 'eastoutside')
            HasDataPlotted = false;
            for idx=1:length(TempMatches)
                
                if UsePhysicalAPLength
                    EmbryoIndex =  TempMatches(idx);
                    APLength = this.APLengths(idx);
                else
                    APLength = 1;
                end
                if ~all(isnan(AllNCParamSEs(idx, :)))
                    use_idx = ~isnan(AllNCParams(idx,:)) &  (AllNCParams(idx, :)./AllNCParamSEs(idx, :) >= 1) & ...
                        (AllR2s(idx,:)>= R2bound);
                else
                    use_idx = ~isnan(AllNCParams(idx,:)) & ...
                        (AllR2s(idx,:)>= R2bound);
                end
                
                if (sum(use_idx) == 0)
                    FrameProfAx.Children(end-(2*(idx-1)+1)).XData = APLength*APbins;
                    FrameProfAx.Children(end-(2*(idx-1)+1)).YData = zeros(1, length(APbins));
                    FrameProfAx.Children(end-(2*(idx-1))).XData = APLength*APbins;
                    FrameProfAx.Children(end-(2*(idx-1))).YData = zeros(1, length(APbins));
                    FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = .1*ones(1, length(APbins));
                    FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta = .1*ones(1, length(APbins));
                    set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','off'); %'off' or 'on'
                    set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','off'); %'off' or 'on'
                else
                    HasDataPlotted = true;
                    FrameProfAx.Children(end-(2*(idx-1)+1)).YData = AllNCParams(idx,use_idx);
                    FrameProfAx.Children(end-(2*(idx-1)+1)).XData = APbins(use_idx);
                    FrameProfAx.Children(end-(2*(idx-1))).YData = AllNCParams(idx,use_idx);
                    FrameProfAx.Children(end-(2*(idx-1))).XData = APbins(use_idx);
                    FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = AllNCParamSEs(idx,use_idx);
                    FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta  = AllNCParamSEs(idx,use_idx);
                    
                    set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','on'); %'off' or 'on'
                    set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','on'); %'off' or 'on'
                    
                    
                end
            end
            
            if ~HasDataPlotted
                continue
            end
            %try
            if exist('PlotTitle', 'var')
                
                title(FrameProfAx, {PlotTitle,...
                    ['Nuclear Cycle ',num2str(NC), ', ', num2str(CurrentTemperature), ' °C']})
                
            else
                
                title(FrameProfAx,  ['Nuclear Cycle ',num2str(NC), ', ', num2str(CurrentTemperature), ' °C'])
                
                
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
            else
                %                 legend_labels = {};
                %                 profs2 = {};
                %                 for i= 1:length(plotted_temps)
                %                     if plotted_temps(i) > 0
                %                         current_temp = Temp_sp(plotted_temps(i));
                %                         legend_labels{1, length(legend_labels)+1} =[ num2str(current_temp),'°C'];
                %                         profs2{1, length(profs2)+1} = prof{plotted_temps(i)};
                %                     end
                %                 end
                %                 hlegend = legend([profs2{:}], legend_labels, 'Location', 'northeast',...
                %                     'FontSize', 10);
                hlegend = legend(FrameProfAx, legend_labels2(TempMatches), 'Location', 'eastoutside',...
                    'FontSize', 10);
            end
            
            
            saveas(FrameProfFig,[outdir4, filesep,...
                TraceType, PhysicalAPString,GradString,'_',OutputString,'_NC',num2str(NC),'_T', strrep(num2str(CurrentTemperature), '.', '_'), 'C.png']);
            
        end
        
    end
end
%close all

%% Calculate x  limits for 2nd subplot


TemperatureVector = 1./(R*(temperatures + 273));
Subplot2DataXmin = min(TemperatureVector);
Subplot2DataXmax = max(TemperatureVector);
Subplot2Xspan = Subplot2DataXmax-Subplot2DataXmin;
Plot2Xmin = Subplot2DataXmin - Subplot2Xspan*.05;
Plot2Xmax = Subplot2DataXmax + Subplot2Xspan*.05;

if ~SkipParamsVsTemp
    for nc_idx=1:length(this.IncludedNCs)
        
        NC = this.IncludedNCs(nc_idx);
        % Prepare Traces for plotting
        
        
        NCMaxParams = NaN(1, length(Temp_sp));
        AllNCParams = NaN(NumSets, length(APbins));
        AllNCParamSEs = NaN(NumSets, length(APbins));
        AllR2s = NaN(NumSets, length(APbins));
        for i=1:NumSets
            idx = this.ProcessedExperiments(i);
            SetParams = PlottedParams(idx,:,NC-8).';
            SetSEParams = PlottedParamSEs(idx,:,NC-8).';
            SetR2s = R2s(idx,:,NC-8).';
            if ~all(isnan(SetSEParams))
                IncludedBins = find(~isnan(SetParams) & (SetR2s > R2bound) & (SetParams./SetSEParams >= 1)) ;
            else
                IncludedBins = find(~isnan(SetParams) & (SetR2s > R2bound)) ;
            end
            if ~isempty(IncludedBins)
                AllNCParams(idx,:) = SetParams;
                AllNCParamSEs(idx,:) = SetSEParams;
                TempSEParams = SetSEParams;
                TempSEParams(isnan(TempSEParams)) = 0;
                NCMaxParams(idx) = max(SetParams+TempSEParams);
                AllR2s(idx,:) = SetR2s;
            end
        end
        
        if all(isnan(NCMaxParams))
            continue
        end
        outdir2 = [outdir,filesep,OutputString];
        if ~exist(outdir2, 'dir')
            mkdir(outdir2)
        end
        
        outdir3 = [outdir2,filesep, datestr(now, 'yyyymmdd')];
        if ~exist(outdir3, 'dir')
            mkdir(outdir3)
        end
        for APindex = 1:length(APbins)
            APbin = APbins(APindex);
            clear eb prof eb2 profB
            plotted_temps = zeros(1, length(temperatures));
            %         plotted_temps = zeros(1, length(temperatures));
            
            close all
            
            
            
            
            FrameProfFig = figure(1);
            set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .8, .6]);
            set(gcf,'color','w');
            FrameProfAx = subplot(1, 5, [1,2], gca);
            for idx =1:length(Temp_sp)
                
                temp_idx = find(temperatures == Temp_sp(idx));
                
                eb{idx} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
                hold on
                
                set(eb{idx}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                
                set(get(get(eb{idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
                
                prof{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(temp_idx,:), 'linestyle', 'none');
                
                
                
                set(eb{idx},'Visible','off'); %'off' or 'on'
                set(prof{idx},'Visible','off'); %'off' or 'on'
            end
            
            
            
            hold off
            
            xlabel('Temperature (°C)')
            xlim([15, 30])
            
            ylabel(ylab)
            ylim([GlobalPlotYmin, GlobalPlotYmax])
            
            FrameProfAx2 = subplot(1, 5, [3,4]);
            
            if IncludeFits
                ci_plotline = fill([0, 1, 1, 0], [0, 0, max(GlobalPlotYmax*1.2,1), max(GlobalPlotYmax*1.2,1)], [0, 0, 0] );
                ci_plotline.FaceAlpha = 0.2;
                set(ci_plotline, 'EdgeColor', 'none');
                hold on
                
                fitted_prof = plot(1./(this.R*(Temp_sp+273)), ones(1, length(Temp_sp)), '-', 'Color', [0, 0, 0]);
                set(fitted_prof,'Visible','off'); %'off' or 'on'
                set(ci_plotline,'Visible','off'); %'off' or 'on'
                subplot2_labels = cell(1, 2);
            end
            
            for idx =1:length(Temp_sp)
                
                temp_idx = find(temperatures == Temp_sp(idx));
                
                %temp_idx = mod(i, 7)+1;
                eb2{idx} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
                hold on
                
                set(eb2{idx}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                
                
                set(get(get(eb2{idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
                
                
                profB{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(temp_idx,:), 'linestyle', 'none');
                
                
                set(eb2{idx},'Visible','off'); %'off' or 'on'
                set(profB{idx},'Visible','off'); %'off' or 'on'
            end
            
            
            
            
            hold off
            
            xlabel('1/(RT) (mol/kJ)')
            xlim([Plot2Xmin, Plot2Xmax])
            
            
            
            ylabel(ylab)
            ylim([LogPlotYmin, GlobalPlotYmax])
            
            set(FrameProfAx2, 'YScale', 'log')
            
            %legend(FrameProfAx, legend_labels, 'Location', 'eastoutside')
            PlottedSets = zeros(1, length(Temp_sp), 'logical');
            for idx=1:length(Temp_sp)
                if ~isnan(AllNCParamSEs(idx, APindex))
                    NoDataCondition = isnan(AllNCParams(idx, APindex)) | (AllNCParams(idx, APindex)./AllNCParamSEs(idx, APindex) <  1) | ...
                        (AllR2s(idx,APindex)< R2bound);
                else
                    NoDataCondition = isnan(AllNCParams(idx, APindex)) | (AllR2s(idx,APindex)< R2bound);
                end
                
                if  NoDataCondition
                    FrameProfAx.Children(end-(2*(idx-1)+1)).XData = 25;
                    FrameProfAx.Children(end-(2*(idx-1)+1)).YData = 1;
                    FrameProfAx.Children(end-(2*(idx-1))).XData = 25;
                    FrameProfAx.Children(end-(2*(idx-1))).YData = 1;
                    FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = 1;
                    FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta = 1;
                    set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','off'); %'off' or 'on'
                    set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','off'); %'off' or 'on'
                    
                    if IncludeFits
                        
                        
                        
                        
                        FrameProfAx2.Children(end-(2*(idx-1)+1)-2).XData = 1/(R*(25+273));
                        FrameProfAx2.Children(end-(2*(idx-1)+1)-2).YData = 1;
                        FrameProfAx2.Children(end-(2*(idx-1))-2).XData = 1/(R*(25+273));
                        FrameProfAx2.Children(end-(2*(idx-1))-2).YData = 1;
                        FrameProfAx2.Children(end-(2*(idx-1))-2).YPositiveDelta = 1;
                        FrameProfAx2.Children(end-(2*(idx-1))-2).YNegativeDelta = 1;
                        set(FrameProfAx2.Children(end-(2*(idx-1)+1)-2),'Visible','off'); %'off' or 'on'
                        set(FrameProfAx2.Children(end-(2*(idx-1))-2),'Visible','off'); %'off' or 'on'
                    else
                        FrameProfAx2.Children(end-(2*(idx-1)+1)).XData = 1/(R*(25+273));
                        FrameProfAx2.Children(end-(2*(idx-1)+1)).YData = 1;
                        FrameProfAx2.Children(end-(2*(idx-1))).XData = 1/(R*(25+273));
                        FrameProfAx2.Children(end-(2*(idx-1))).YData = 1;
                        FrameProfAx2.Children(end-(2*(idx-1))).YPositiveDelta = 1;
                        FrameProfAx2.Children(end-(2*(idx-1))).YNegativeDelta = 1;
                        set(FrameProfAx2.Children(end-(2*(idx-1)+1)),'Visible','off'); %'off' or 'on'
                        set(FrameProfAx2.Children(end-(2*(idx-1))),'Visible','off'); %'off' or 'on'
                    end
                else
                    PlottedSets(idx) = 1;
                    FrameProfAx.Children(end-(2*(idx-1)+1)).YData = AllNCParams(idx, APindex);
                    FrameProfAx.Children(end-(2*(idx-1)+1)).XData =Temp_obs(idx);
                    FrameProfAx.Children(end-(2*(idx-1))).YData =  AllNCParams(idx, APindex);
                    FrameProfAx.Children(end-(2*(idx-1))).XData = Temp_obs(idx);
                    FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = AllNCParamSEs(idx,APindex);
                    FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta  = AllNCParamSEs(idx,APindex);
                    
                    set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','on'); %'off' or 'on'
                    set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','on'); %'off' or 'on'
                    
                    if IncludeFits
                        
                        
                        FrameProfAx2.Children(end-(2*(idx-1)+1)-2).XData =1/(R*(Temp_obs(idx)+273));
                        FrameProfAx2.Children(end-(2*(idx-1)+1)-2).YData = AllNCParams(idx, APindex);
                        FrameProfAx2.Children(end-(2*(idx-1))-2).XData = 1/(R*(Temp_obs(idx)+273));
                        FrameProfAx2.Children(end-(2*(idx-1))-2).YData = AllNCParams(idx, APindex);
                        FrameProfAx2.Children(end-(2*(idx-1))-2).YPositiveDelta = AllNCParamSEs(idx,APindex);
                        FrameProfAx2.Children(end-(2*(idx-1))-2).YNegativeDelta = AllNCParamSEs(idx,APindex);
                        set(FrameProfAx2.Children(end-(2*(idx-1)+1)-2),'Visible','on'); %'off' or 'on'
                        set(FrameProfAx2.Children(end-(2*(idx-1))-2),'Visible','on'); %'off' or 'on'
                    else
                        FrameProfAx2.Children(end-(2*(idx-1)+1)).YData = AllNCParams(idx, APindex);
                        FrameProfAx2.Children(end-(2*(idx-1)+1)).XData =1/(R*(Temp_obs(idx)+273));
                        FrameProfAx2.Children(end-(2*(idx-1))).YData =  AllNCParams(idx, APindex);
                        FrameProfAx2.Children(end-(2*(idx-1))).XData = 1/(R*(Temp_obs(idx)+273));
                        FrameProfAx2.Children(end-(2*(idx-1))).YPositiveDelta = AllNCParamSEs(idx,APindex);
                        FrameProfAx2.Children(end-(2*(idx-1))).YNegativeDelta  = AllNCParamSEs(idx,APindex);
                        
                        set(FrameProfAx2.Children(end-(2*(idx-1)+1)),'Visible','on'); %'off' or 'on'
                        set(FrameProfAx2.Children(end-(2*(idx-1))),'Visible','on'); %'off' or 'on'
                    end
                    
                    
                    
                    temp_idx = find(temperatures == Temp_sp(idx));
                    if plotted_temps(temp_idx) == 0
                        plotted_temps(temp_idx) = idx;
                    end
                    
                    
                end
            end
            
            if IncludeFits
                [fitx,fity, ci, Ea, se_Ea, LogA, se_LogA, fitR2] = ...
                    getActivationEnergyFitTraces(this, parameter, NC, APindex, TraceType);
                if ~isempty(fity)
                    curve1 = ci(:,1).';
                    curve2 = ci(:,2).';
                    curve1(curve1 > max(GlobalPlotYmax*1.05,1)) = max(GlobalPlotYmax*1.05,1);
                    curve1(curve1 < LogPlotYmin) = LogPlotYmin;
                    curve2(curve2 > max(GlobalPlotYmax*1.05,1)) = max(GlobalPlotYmax*1.05,1);
                    curve2(curve2 < LogPlotYmin) = LogPlotYmin;
                    t_vector2 = [fitx, fliplr(fitx)];
                    inBetween = [curve1, fliplr(curve2)];
                    FrameProfAx2.Children(end).Faces= 1:length(t_vector2);
                    FrameProfAx2.Children(end).Vertices = [t_vector2; inBetween].';
                    set(FrameProfAx2.Children(end),'Visible','on');
                    
                    
                    FrameProfAx2.Children(end-1).XData = fitx;
                    FrameProfAx2.Children(end-1).YData = fity;
                    set(FrameProfAx2.Children(end-1),'Visible','on');
                    
                    lab1a = MeanSE_num2str(Ea, se_Ea, Nsigfigs);
                    
                    
                    
                    subplot2_labels{2} = '95% Confidence Interval';
                    subplot2_labels{1} = ['E_{A}: ', lab1a.m,...
                        ' \pm ', lab1a.se, ' kJ/mol, R^2: ', num2str(fitR2, 2)];
                    
                    
                    if fity(1) < fity(2)
                        legLoc = 'northeast';
                    else
                        legLoc = 'northwest';
                    end
                    
                    hlegend2 = legend(FrameProfAx2, [fitted_prof, ci_plotline], subplot2_labels, 'Location', legLoc,...
                        'FontSize', 8);
                    
                    
                    
                else
                    set(FrameProfAx2.Children(end),'Visible','off');
                    set(FrameProfAx2.Children(end-1),'Visible','off');
                    subplot2_labels{2} = '';
                    subplot2_labels{1} = 'No Fit Info available';
                end
            else
                set(FrameProfAx2.Children(end),'Visible','off');
                set(FrameProfAx2.Children(end-1),'Visible','off');
                
                subplot2_labels{2} = '';
                subplot2_labels{1} = 'No Fit Info available';
                
                
                
            end
            
            
            if  all(~PlottedSets)
                continue
            end
            %             legend_labels2 = {};
            %             profs2 = {};
            %             for i= 1:length(plotted_temps)
            %                 if plotted_temps(i) > 0
            %                     current_temp = Temp_sp(plotted_temps(i));
            %                     legend_labels2{1, length(legend_labels2)+1} =[ num2str(current_temp),'°C'];
            %                     profs2{1, length(profs2)+1} = profB{plotted_temps(i)};
            %                 end
            %             end
            %try
            if exist('PlotTitle', 'var')
                
                
                sgtitle({PlotTitle,...
                    ['Nuclear Cycle ',num2str(NC), ', Fraction Embryo Length: ', num2str(APbin)]})
                
            else
                
                sgtitle(['Nuclear Cycle ',num2str(NC), ', Fraction Embryo Length: ', num2str(APbin)])
                
            end
            
            LegendAx = subplot(1, 5, 5);
            legend_profs = cell(1, length(Temp_sp));
            legend_labels3 = {};
            for SetIndex = 1:length(Temp_sp)
                temp_idx = find(temperatures == Temp_sp(SetIndex));
                
                legend_profs{SetIndex} = plot([0], [0], '.', 'color', colors(temp_idx,:),...
                    'linestyle', 'none');
                hold on
               
                    legend_labels3{1, length(legend_labels3)+1} = legend_labels{SetIndex};
         
            end
            hold off
            axis off
            
            
            hlegend = legend(LegendAx, legend_labels3,...
                'FontSize', 8);%, 'Location', );
            hlegend.Position(1) = LegendAx.Position(1)+ LegendAx.Position(3)/2 - hlegend.Position(3)/2;
            hlegend.Position(2) = LegendAx.Position(2)+ LegendAx.Position(4)/2 - hlegend.Position(4)/2;
            
            
            
            saveas(FrameProfFig,[outdir3, filesep,...
                TraceType,'_', OutputString,  '_NC',num2str(NC), 'AP', num2str(APindex), '.png']);
        end
        
    end
    close all
end
%%

NumTemperatures = length(temperatures);

if ~SkipBinnedParamsVsAP
    outdir2 = [outdir,filesep,'Binned', OutputString];
    if ~exist(outdir2, 'dir')
        mkdir(outdir2)
    end
    outdir3 = [outdir2, filesep, datestr(now, 'yyyymmdd')];
    if ~exist(outdir3, 'dir')
        mkdir(outdir3)
    end
    
    
    for nc_idx=1:length(this.IncludedNCs)
        NC = this.IncludedNCs(nc_idx);
        clear prof eb
        plotted_temps = zeros(1, length(temperatures));
        
        
        
        % Prepare Traces for plotting
        
        NCMaxParams = NaN(1,NumTemperatures);
        AllNCParams = NaN(NumTemperatures, length(APbins));
        AllNCParamSEs = NaN(NumTemperatures, length(APbins));
        AllNCTemperatures = NaN(NumTemperatures, length(APbins));
        AllNCTemperatureSEs = NaN(NumTemperatures, length(APbins));
        AllCounts = NaN(NumTemperatures, length(APbins));
        
        for i=1:NumTemperatures
            SetParams = BinnedParams(i,:,NC-8).';
            SetSEParams = BinnedSEParams(i,:,NC-8).';
            SetTemps = ParamTemperatures(i,:,NC-8).';
            SetSETemps = ParamSETemperatures(i,:,NC-8).';
            SetCounts = Counts(i,:,NC-8).';
            
            IncludedBins = find(~isnan(SetParams) & (SetCounts >= this.MinimumBinCount) & (SetParams./SetSEParams >= 1)) ;
            if ~isempty(IncludedBins)
                AllNCParams(i,:) = SetParams;
                AllNCParamSEs(i,:) = SetSEParams;
                NCMaxParams(i,:) = max(SetParams(i,:) +SetSEParams(i,:));
                AllCounts(i,:) = SetCounts;
                AllNCTemperatures(i,:) =SetTemps;
                AllNCTemperatureSEs(i,:) =SetSETemps;
            end
        end
        
        if all(isnan(NCMaxParams))
            continue
        end
        close all
        
        
        FrameProfFig = figure(1);
        set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .5, .4]);
        set(gcf,'color','w');
        FrameProfAx = axes(FrameProfFig);
        for idx =1:NumTemperatures
            if strcmp(lower(PlottingColors), 'gradient')
                temp_idx = find(abs(Temp_range-temperatures(idx)) == min(abs(Temp_range-temperatures(idx))));
            else
                temp_idx = idx;
            end
            %temp_idx = mod(i, 7)+1;
            eb{idx} = errorbar(APbins, ones(1, length(APbins)),...
                .1*ones(1, length(APbins)), 'vertical',...
                'LineStyle', 'none');
            
            hold on
            if strcmp(lower(PlottingColors), 'gradient')
                set(eb{idx}, 'color', colors(temp_idx,:), 'LineWidth', 1);
            else
                set(eb{idx}, 'color', colors(temp_idx,:), 'LineWidth', 1);
            end
            set(get(get(eb{idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            if UseLines
                if strcmp(lower(PlottingColors), 'gradient')
                    prof{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(temp_idx,:));
                else
                    prof{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(temp_idx,:));
                end
            else
                if strcmp(lower(PlottingColors), 'gradient')
                    prof{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(temp_idx,:), 'linestyle', 'none');
                else
                    prof{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(temp_idx,:), 'linestyle', 'none');
                end
            end
            
            set(eb{idx},'Visible','off'); %'off' or 'on'
            set(prof{idx},'Visible','off'); %'off' or 'on'
        end
        
        
        
        hold off
        if ~UsePhysicalAPLength
            xlabel('Fraction Embryo Length')
            xlim([0, 1])
        else
            xlabel('Distance from the Anterior Pole (\mum)')
            xlim([0, max(this.APLengths)])
        end
        
        ylabel(ylab)
        ylim([GlobalPlotYmin, GlobalPlotYmax])
        
        %ylim([min(0, min(min(AllTimeOns+AllStdErrorTimeOns))*1.2), max(max(max(AllTimeOns+AllStdErrorTimeOns))*1.2, 1)])
        %
        
        title(FrameProfAx, {'',...
            ['Nuclear Cycle ',num2str(NC)]})
        
        
        
        %legend(FrameProfAx, legend_labels, 'Location', 'eastoutside')
        HasDataPlotted = false;
        for idx=1:NumTemperatures
            
            if UsePhysicalAPLength
                EmbryoIndex =  this.ProcessedExperiments(idx);
                APLength = this.APLengths(idx);
            else
                APLength = 1;
            end
            if ~all(isnan(AllNCParamSEs(idx, :)))
                use_idx = ~isnan(AllNCParams(idx,:)) &  (AllNCParams(idx, :)./AllNCParamSEs(idx, :) >= 1) & ...
                    (AllCounts(idx,:)>= this.MinimumBinCount);
            else
                use_idx = ~isnan(AllNCParams(idx,:)) & ...
                    (AllCounts(idx,:)>= this.MinimumBinCount);
            end
            if sum(use_idx) == 0 %| sum(DiffMeanFluoMat(i, use_idx, j) == 0)
                FrameProfAx.Children(end-(2*(idx-1)+1)).XData = APLength*APbins;
                FrameProfAx.Children(end-(2*(idx-1)+1)).YData = zeros(1, length(APbins));
                FrameProfAx.Children(end-(2*(idx-1))).XData = APLength*APbins;
                FrameProfAx.Children(end-(2*(idx-1))).YData = zeros(1, length(APbins));
                FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = .1*ones(1, length(APbins));
                FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta = .1*ones(1, length(APbins));
                set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','off'); %'off' or 'on'
                set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','off'); %'off' or 'on'
            else
                HasDataPlotted = true;
                FrameProfAx.Children(end-(2*(idx-1)+1)).YData = AllNCParams(idx,use_idx);
                FrameProfAx.Children(end-(2*(idx-1)+1)).XData = APbins(use_idx);
                FrameProfAx.Children(end-(2*(idx-1))).YData = AllNCParams(idx,use_idx);
                FrameProfAx.Children(end-(2*(idx-1))).XData = APbins(use_idx);
                FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = AllNCParamSEs(idx,use_idx);
                FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta  = AllNCParamSEs(idx,use_idx);
                set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','on'); %'off' or 'on'
                set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','on'); %'off' or 'on'
                if ~strcmp(lower(PlottingColors), 'gradient')
                    temp_idx = idx;
                    if plotted_temps(temp_idx) == 0
                        plotted_temps(temp_idx) = idx;
                    end
                end
                
            end
        end
        %try
        if exist('PlotTitle', 'var')
            
            
            title(FrameProfAx, {PlotTitle,...
                ['Nuclear Cycle ',num2str(NC)]})
            
        else
            
            title(FrameProfAx,  ['Nuclear Cycle ',num2str(NC)])
            
            
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
        else
            legend_labels = {};
            profs2 = {};
            for i= 1:length(plotted_temps)
                if plotted_temps(i) > 0
                    current_temp = temperatures(plotted_temps(i));
                    legend_labels{1, length(legend_labels)+1} =[ num2str(current_temp),'°C'];
                    profs2{1, length(profs2)+1} = prof{plotted_temps(i)};
                end
            end
            hlegend = legend([profs2{:}], legend_labels, 'Location', 'northeast',...
                'FontSize', 10);
        end
        
        saveas(FrameProfFig,[outdir3, filesep,...
            TraceType, PhysicalAPString, GradString,'_Binned', OutputString, '_NC',num2str(NC),'.png']);
        
    end
end


%% Calculate x  limits for 2nd subplot

AllTpoints = reshape(ParamTemperatures, 1,...
    size(ParamTemperatures, 1)*size(ParamTemperatures, 2)*size(ParamTemperatures, 3));
AllSETpoints = reshape(ParamSETemperatures, 1,...
    size(ParamTemperatures, 1)*size(ParamTemperatures, 2)*size(ParamTemperatures, 3));
AllCombinedTpoints = reshape(ParamTemperatures+ParamSETemperatures, 1,...
    size(ParamTemperatures, 1)*size(ParamTemperatures, 2)*size(ParamTemperatures, 3));
AllCombinedTpoints = AllCombinedTpoints(~isnan(AllCombinedTpoints));
AllTpoints = AllTpoints(~isnan(AllTpoints));
AllSETpoints = AllSETpoints(~isnan(AllSETpoints));
TemperatureVector = 1./(R*(AllTpoints + 273))+sqrt(((1./(R*(AllTpoints + 273).^2)).^2).*AllSETpoints.^2);
Subplot1DataXmin = min(AllCombinedTpoints)*0.95;
Subplot1DataXmax = max(AllCombinedTpoints)*1.05;
Subplot2DataXmin = min(TemperatureVector);
Subplot2DataXmax = max(TemperatureVector);
Subplot2Xspan = Subplot2DataXmax-Subplot2DataXmin;
Plot2Xmin = Subplot2DataXmin - Subplot2Xspan*.05;
Plot2Xmax = Subplot2DataXmax + Subplot2Xspan*.05;

if ~SkipBinnedParamsVsTemp
    for nc_idx=1:length(this.IncludedNCs)
        
        NC = this.IncludedNCs(nc_idx);
        % Prepare Traces for plotting
        
        NCMaxParams = NaN(1,NumTemperatures);
        AllNCParams = NaN(NumTemperatures, length(APbins));
        AllNCParamSEs = NaN(NumTemperatures, length(APbins));
        AllNCTemperatures = NaN(NumTemperatures, length(APbins));
        AllNCTemperatureSEs = NaN(NumTemperatures, length(APbins));
        AllCounts = NaN(NumTemperatures, length(APbins));
        for i=1:NumTemperatures
            SetParams = BinnedParams(i,:,NC-8).';
            SetSEParams = BinnedSEParams(i,:,NC-8).';
            SetTemps = ParamTemperatures(i,:,NC-8).';
            SetSETemps = ParamSETemperatures(i,:,NC-8).';
            SetCounts = Counts(i,:,NC-8).';
            
            IncludedBins = find(~isnan(SetParams) & (SetCounts >= this.MinimumBinCount) & (SetParams./SetSEParams >= 1)) ;
            if ~isempty(IncludedBins)
                AllNCParams(i,:) = SetParams;
                AllNCParamSEs(i,:) = SetSEParams;
                NCMaxParams(i,:) = max(SetParams(i,:) +SetSEParams(i,:));
                AllCounts(i,:) = SetCounts;
                AllNCTemperatures(i,:) =SetTemps;
                AllNCTemperatureSEs(i,:) =SetSETemps;
            end
        end
        
        if all(isnan(NCMaxParams))
            continue
        end
        outdir2 = [outdir,filesep,'Binned', OutputString];
        if ~exist(outdir2, 'dir')
            mkdir(outdir2)
        end
        
        outdir3 = [outdir2,filesep, datestr(now, 'yyyymmdd')];
        if ~exist(outdir3, 'dir')
            mkdir(outdir3)
        end
        for APindex = 1:length(APbins)
            APbin = APbins(APindex);
            clear eb prof eb2 profB
            plotted_temps = zeros(1, length(temperatures));
            %         plotted_temps = zeros(1, length(temperatures));
            
            close all
            
            
            
            
            FrameProfFig = figure(1);
            set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .4]);
            set(gcf,'color','w');
            FrameProfAx = subplot(1, 6, [1,2], gca);
            for idx =1:NumTemperatures
                
                
                eb{idx} = errorbar(APbins, ones(1, length(APbins)),...
                    .1*ones(1, length(APbins)), .1*ones(1, length(APbins)),...
                    .1*ones(1, length(APbins)), .1*ones(1, length(APbins)),...
                    'LineStyle', 'none');
                
                hold on
                
                set(eb{idx}, 'color', colors(idx,:), 'LineWidth', 1);
                
                set(get(get(eb{idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
                
                prof{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(idx,:), 'linestyle', 'none');
                
                
                
                set(eb{idx},'Visible','off'); %'off' or 'on'
                set(prof{idx},'Visible','off'); %'off' or 'on'
            end
            
            
            
            hold off
            
            xlabel('Temperature (°C)')
            xlim([Subplot1DataXmin, Subplot1DataXmax])
            
            ylabel(ylab)
            ylim([GlobalPlotYmin, GlobalPlotYmax])
            
            FrameProfAx2 = subplot(1, 6, [4,5,6]);
            
            for idx =1:NumTemperatures
                
                
                %temp_idx = mod(i, 7)+1;
                eb2{idx} = errorbar(APbins, ones(1, length(APbins)),...
                    .1*ones(1, length(APbins)), .1*ones(1, length(APbins)),...
                    .1*ones(1, length(APbins)), .1*ones(1, length(APbins)),...
                    'LineStyle', 'none');
                
                hold on
                set(eb2{idx}, 'color', colors(idx,:), 'LineWidth', 1);
                
                
                set(get(get(eb2{idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
                
                
                profB{idx} = plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(idx,:), 'linestyle', 'none');
                
                
                set(eb2{idx},'Visible','off'); %'off' or 'on'
                set(profB{idx},'Visible','off'); %'off' or 'on'
            end
            
            
            
            hold off
            
            xlabel('1/(RT) (mol/kJ)')
            xlim([Plot2Xmin, Plot2Xmax])
            
            
            
            ylabel(ylab)
            ylim([LogPlotYmin, GlobalPlotYmax])
            
            set(FrameProfAx2, 'YScale', 'log')
            
            %legend(FrameProfAx, legend_labels, 'Location', 'eastoutside')
            HasDataPlotted = false;
            for idx=1:NumTemperatures
                if ~isnan(AllNCParamSEs(idx, APindex))
                    NoDataCondition = isnan(AllNCParams(idx, APindex)) | (AllNCParams(idx,APindex)/AllNCParamSEs(idx, APindex) <  1) | ...
                        (AllCounts(idx,APindex)< this.MinimumBinCount);
                else
                    NoDataCondition = isnan(AllNCParams(idx, APindex)) | ...
                        (AllCounts(idx,APindex)< this.MinimumBinCount);
                end
                
                if  NoDataCondition
                    FrameProfAx.Children(end-(2*(idx-1)+1)).XData = 25;
                    FrameProfAx.Children(end-(2*(idx-1)+1)).YData = 1;
                    FrameProfAx.Children(end-(2*(idx-1))).XData = 25;
                    FrameProfAx.Children(end-(2*(idx-1))).YData = 1;
                    FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = 1;
                    FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta = 1;
                    FrameProfAx.Children(end-(2*(idx-1))).XPositiveDelta = 1;
                    FrameProfAx.Children(end-(2*(idx-1))).XNegativeDelta = 1;
                    set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','off'); %'off' or 'on'
                    set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','off'); %'off' or 'on'
                    
                    FrameProfAx2.Children(end-(2*(idx-1)+1)).XData = 1/(R*(25+273));
                    FrameProfAx2.Children(end-(2*(idx-1)+1)).YData = 1;
                    FrameProfAx2.Children(end-(2*(idx-1))).XData = 1/(R*(25+273));
                    FrameProfAx2.Children(end-(2*(idx-1))).YData = 1;
                    FrameProfAx2.Children(end-(2*(idx-1))).YPositiveDelta = 1;
                    FrameProfAx2.Children(end-(2*(idx-1))).YNegativeDelta = 1;
                    FrameProfAx2.Children(end-(2*(idx-1))).XPositiveDelta = 1;
                    FrameProfAx2.Children(end-(2*(idx-1))).XNegativeDelta = 1;
                    set(FrameProfAx2.Children(end-(2*(idx-1)+1)),'Visible','off'); %'off' or 'on'
                    set(FrameProfAx2.Children(end-(2*(idx-1))),'Visible','off'); %'off' or 'on'
                else
                    HasDataPlotted = true;
                    FrameProfAx.Children(end-(2*(idx-1)+1)).YData = AllNCParams(idx, APindex);
                    FrameProfAx.Children(end-(2*(idx-1)+1)).XData =AllNCTemperatures(idx, APindex);
                    FrameProfAx.Children(end-(2*(idx-1))).YData =  AllNCParams(idx, APindex);
                    FrameProfAx.Children(end-(2*(idx-1))).XData = AllNCTemperatures(idx, APindex);
                    FrameProfAx.Children(end-(2*(idx-1))).YPositiveDelta = AllNCParamSEs(idx,APindex);
                    FrameProfAx.Children(end-(2*(idx-1))).YNegativeDelta  = AllNCParamSEs(idx,APindex);
                    FrameProfAx.Children(end-(2*(idx-1))).XPositiveDelta = AllNCTemperatureSEs(idx,APindex);
                    FrameProfAx.Children(end-(2*(idx-1))).XNegativeDelta  = AllNCTemperatureSEs(idx,APindex);
                    
                    set(FrameProfAx.Children(end-(2*(idx-1)+1)),'Visible','on'); %'off' or 'on'
                    set(FrameProfAx.Children(end-(2*(idx-1))),'Visible','on'); %'off' or 'on'
                    
                    FrameProfAx2.Children(end-(2*(idx-1)+1)).YData = AllNCParams(idx, APindex);
                    FrameProfAx2.Children(end-(2*(idx-1)+1)).XData =1/(R*(AllNCTemperatures(idx, APindex)+273));
                    FrameProfAx2.Children(end-(2*(idx-1))).YData =  AllNCParams(idx, APindex);
                    FrameProfAx2.Children(end-(2*(idx-1))).XData = 1/(R*(AllNCTemperatures(idx, APindex)+273));
                    FrameProfAx2.Children(end-(2*(idx-1))).YPositiveDelta = AllNCParamSEs(idx,APindex);
                    FrameProfAx2.Children(end-(2*(idx-1))).YNegativeDelta  = AllNCParamSEs(idx,APindex);
                    FrameProfAx2.Children(end-(2*(idx-1))).XPositiveDelta = sqrt(AllNCTemperatureSEs(idx,APindex)^2*(1/(R*(AllNCTemperatures(idx, APindex)+273)^2))^2);
                    FrameProfAx2.Children(end-(2*(idx-1))).XNegativeDelta  = sqrt(AllNCTemperatureSEs(idx,APindex)^2*(1/(R*(AllNCTemperatures(idx, APindex)+273)^2))^2);
                    
                    set(FrameProfAx2.Children(end-(2*(idx-1)+1)),'Visible','on'); %'off' or 'on'
                    set(FrameProfAx2.Children(end-(2*(idx-1))),'Visible','on'); %'off' or 'on'
                    
                    
                    if plotted_temps(idx) == 0
                        plotted_temps(idx) = idx;
                    end
                    
                    
                end
            end
            
            
            if ~HasDataPlotted
                continue
            end
            legend_labels2 = {};
            profs2 = {};
            for i= 1:length(plotted_temps)
                if plotted_temps(i) > 0
                    current_temp = temperatures(i);
                    legend_labels2{1, length(legend_labels2)+1} =[ num2str(current_temp),'°C'];
                    profs2{1, length(profs2)+1} = profB{plotted_temps(i)};
                end
            end
            hlegend = legend(FrameProfAx2, [profs2{:}], legend_labels2, 'Location', 'eastoutside',...
                'FontSize', 10);
            %try
            if exist('PlotTitle', 'var')
                
                
                sgtitle({PlotTitle,...
                    ['Nuclear Cycle ',num2str(NC), ', Fraction Embryo Length: ', num2str(APbin)]})
                
            else
                
                sgtitle(['Nuclear Cycle ',num2str(NC), ', Fraction Embryo Length: ', num2str(APbin)])
                
            end
            
            
            saveas(FrameProfFig,[outdir3, filesep,...
                TraceType,'_Binned', OutputString,  '_NC',num2str(NC), 'AP', num2str(APindex), '.png']);
        end
        
    end
    close all
end
%%
if ~SkipAPSubplots
    
    if  strcmp(lower(PlottingColors), 'gradient')
        [~, colors] = getColorPalettes();
        GradString = '';
    end
    
    
    PlottedParamSEs((R2s < R2bound)| (PlottedParams > GlobalPlotYmax) | (PlottedParams < GlobalPlotYmin)) = NaN;
    PlottedParams((R2s < R2bound)| (PlottedParams > GlobalPlotYmax) | (PlottedParams < GlobalPlotYmin)) = NaN;
    
    WhereValid = squeeze(sum(~isnan(PlottedParams), 1));
    WhereValidAP = sum(WhereValid, 2).';
    ValidAPIndices = find(WhereValidAP);
    WhereValidNC = sum(WhereValid, 1);
    ValidNCIndices = find(WhereValidNC);
    
    
    SubplotDims = numSubplots(length(ValidAPIndices) + 1);
    
    SubFigDims = [0.9, 0.9*3072/1920*SubplotDims(1)/SubplotDims(2)*1.2];
    SubFigDims(2) = min([SubFigDims(2), 0.95]);
    SubFigDims = round(SubFigDims, 2);
    
    
    Nsets = size(PlottedParams, 1);
    for nc_idx=1:length(ValidNCIndices)
        
        NC = ValidNCIndices(nc_idx)+8;
        % Prepare Traces for plotting
        
        NCMaxParams = NaN(1,Nsets);
        AllNCParams = NaN(Nsets, length(ValidAPIndices));
        AllNCParamSEs = NaN(Nsets, length(ValidAPIndices));
        AllNCTemperatures = NaN(Nsets, length(ValidAPIndices));
        for i=1:Nsets
            SetParams = PlottedParams(i,ValidAPIndices,NC-8);
            SetSEParams = PlottedParamSEs(i,ValidAPIndices,NC-8);
            SetTemps = ones(size(SetParams))*this.Temp_sps(i);
            
            if ~all(isnan(SetSEParams))
                IncludedBins = find(~isnan(SetParams)  & (SetParams./SetSEParams >= 1)) ;
            else
                IncludedBins = find(~isnan(SetParams)) ;
            end
            if ~isempty(IncludedBins)
                AllNCParams(i,:) = SetParams;
                AllNCParamSEs(i,:) = SetSEParams;
                TempSEParams = SetSEParams;
                TempSEParams(isnan(TempSEParams )) = 0;
                NCMaxParams(i) = max(SetParams +TempSEParams);
                AllNCTemperatures(i,:) =SetTemps;
            end
        end
        
        if all(isnan(NCMaxParams))
            continue
        end
        outdir2 = [outdir,filesep, OutputString];
        if ~exist(outdir2, 'dir')
            mkdir(outdir2)
        end
        
        outdir3 = [outdir2,filesep, datestr(now, 'yyyymmdd')];
        if ~exist(outdir3, 'dir')
            mkdir(outdir3)
        end
        
        clear eb prof FrameProfAx
        eb = cell(Nsets, length(ValidAPIndices));
        prof = cell(Nsets, length(ValidAPIndices));
        FrameProfAx = cell(1,  length(ValidAPIndices));
        close all
        FrameProfFig = figure(1);
        set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.01, SubFigDims(1), SubFigDims(2)]);
        set(gcf,'color','w');
        
        for l = 1:length(ValidAPIndices)
            APindex = ValidAPIndices(l);
            APbin = APbins(APindex);
            
            
            
            
            
            
            
            
            
            if l == 1
                FrameProfAx{l} = subplot(SubplotDims(1), SubplotDims(2), l, gca);
            else
                FrameProfAx{l} = subplot(SubplotDims(1), SubplotDims(2), l);
            end
            
            
            
            for idx=1:Nsets
                temp_idx = find(temperatures == this.Temp_sps(idx));
                if ~isnan(AllNCParamSEs(idx, l))
                    NoDataCondition = isnan(AllNCParams(idx, l)) | (AllNCParams(idx,l)/AllNCParamSEs(idx, l) <  1);
                else
                    NoDataCondition = isnan(AllNCParams(idx, l));
                end
                if  NoDataCondition
                    
                    eb{idx, l} = errorbar(AllNCTemperatures(idx, l), AllNCParams(idx, l),...
                        AllNCParamSEs(idx, l),...
                        'LineStyle', 'none');
                    
                    hold on
                    
                    set(eb{idx,l}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                    
                    set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                    
                    
                    prof{idx,l} = plot(AllNCTemperatures(idx, l), AllNCParams(idx, l), '.-', 'Color', colors(temp_idx,:), 'linestyle', 'none');
                    
                    
                    set(eb{idx,l},'Visible','off'); %'off' or 'on'
                    set(prof{idx,l},'Visible','off'); %'off' or 'on'
                    
                    
                else
                    
                    eb{idx, l} = errorbar(AllNCTemperatures(idx, l), AllNCParams(idx, l),...
                        AllNCParamSEs(idx, l),...
                        'LineStyle', 'none');
                    
                    hold on
                    
                    set(eb{idx,l}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                    
                    set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                    
                    
                    prof{idx,l} = plot(AllNCTemperatures(idx, l), AllNCParams(idx, l), '.-', 'Color', colors(temp_idx,:), 'linestyle', 'none');
                    
                    
                    set(eb{idx,l},'Visible','on'); %'off' or 'on'
                    set(prof{idx,l},'Visible','on'); %'off' or 'on'
                    
                end
            end
            
            
            
            hold off
            
            xlabel('Temperature (°C)')
            xlim([Subplot1DataXmin, Subplot1DataXmax])
            
            ylabel(ylab)
            ylim([GlobalPlotYmin, GlobalPlotYmax])
            title(['AP: ', num2str(round(APbin, 3))])
        end
        if exist('PlotTitle', 'var')
            
            
            sgtitle({PlotTitle,...
                ['Nuclear Cycle ',num2str(NC)]})
            
        else
            
            sgtitle(['Nuclear Cycle ',num2str(NC)])
            
        end
        
        
        LegendAx = subplot(SubplotDims(1), SubplotDims(2), SubplotDims(1)*SubplotDims(2));
        for idx = 1:NumTemperatures
            plot([0], [0], '.-', 'Color', colors(idx,:), 'linestyle', 'none');
            hold on
        end
        hold off
        %xlim([-10, 10])
        %ylim([-10, 10])
        axis off
        
        legend_labels2 = {};
        profs2 = {};
        for i= 1:NumTemperatures
            current_temp = temperatures(i);
            legend_labels2{1, length(legend_labels2)+1} =[ num2str(current_temp),'°C'];
            
            %profs2{1, length(profs2)+1} = profB{plotted_temps(i)};
        end
        hlegend = legend(LegendAx, legend_labels2,...
            'FontSize', 12);%, 'Location', );
        hlegend.Position(1) = LegendAx.Position(1)+ LegendAx.Position(3)/2 - hlegend.Position(3)/2;
        hlegend.Position(2) = LegendAx.Position(2)+ LegendAx.Position(4)/2 - hlegend.Position(4)/2;
        
        %try
        
        saveas(FrameProfFig,[outdir3, filesep,...
            TraceType,'_Subplots', OutputString,  '_NC',num2str(NC),'.png']);
        
    end
    for nc_idx=1:length(ValidNCIndices)
        
        NC = ValidNCIndices(nc_idx)+8;
        % Prepare Traces for plotting
        
        NCMaxParams = NaN(1,Nsets);
        AllNCParams = NaN(Nsets, length(ValidAPIndices));
        AllNCParamSEs = NaN(Nsets, length(ValidAPIndices));
        AllNCTemperatures = NaN(Nsets, length(ValidAPIndices));
        for i=1:Nsets
            SetParams = PlottedParams(i,ValidAPIndices,NC-8);
            SetSEParams = PlottedParamSEs(i,ValidAPIndices,NC-8);
            SetTemps = ones(size(SetParams))*this.Temp_sps(i);
            
            
            IncludedBins = find(~isnan(SetParams)  & (SetParams./SetSEParams >= 1)) ;
            if ~isempty(IncludedBins)
                AllNCParams(i,:) = SetParams;
                AllNCParamSEs(i,:) = SetSEParams;
                NCMaxParams(i) = max(SetParams +SetSEParams);
                AllNCTemperatures(i,:) =SetTemps;
            end
        end
        
        if all(isnan(NCMaxParams))
            continue
        end
        outdir2 = [outdir,filesep,OutputString];
        if ~exist(outdir2, 'dir')
            mkdir(outdir2)
        end
        
        outdir3 = [outdir2,filesep, datestr(now, 'yyyymmdd')];
        if ~exist(outdir3, 'dir')
            mkdir(outdir3)
        end
        
        clear eb prof FrameProfAx
        eb = cell(NumTemperatures, length(ValidAPIndices));
        prof = cell(NumTemperatures, length(ValidAPIndices));
        FrameProfAx = cell(1,  length(ValidAPIndices));
        close all
        FrameProfFig = figure(1);
        set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.01, SubFigDims(1), SubFigDims(2)]);
        set(gcf,'color','w');
        
        for l = 1:length(ValidAPIndices)
            APindex = ValidAPIndices(l);
            APbin = APbins(APindex);
            
            
            
            
            
            
            
            
            
            if l == 1
                FrameProfAx{l} = subplot(SubplotDims(1), SubplotDims(2), l, gca);
            else
                FrameProfAx{l} = subplot(SubplotDims(1), SubplotDims(2), l);
            end
            
            
            
            for idx=1:Nsets
                temp_idx = find(temperatures == this.Temp_sps(idx));
                if ~isnan(AllNCParamSEs(idx, l))
                    NoDataCondition = isnan(AllNCParams(idx, l)) | (AllNCParams(idx,l)/AllNCParamSEs(idx, l) <  1);
                else
                    NoDataCondition = isnan(AllNCParams(idx, l));
                end
                if  NoDataCondition
                    
                    eb{idx, l} = errorbar(1/(R*(AllNCTemperatures(idx, l)+273)),...
                        AllNCParams(idx, l),...
                        AllNCParamSEs(idx, l),...
                        'LineStyle', 'none');
                    
                    hold on
                    
                    set(eb{idx,l}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                    
                    set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                    
                    
                    prof{idx,l} = plot(1/(R*(AllNCTemperatures(idx, l)+273)),...
                        AllNCParams(idx, l), '.-', 'Color', colors(temp_idx,:), 'linestyle', 'none');
                    
                    
                    set(eb{idx,l},'Visible','off'); %'off' or 'on'
                    set(prof{idx,l},'Visible','off'); %'off' or 'on'
                    
                    
                else
                    
                    eb{idx, l} = errorbar(1/(R*(AllNCTemperatures(idx, l)+273)),...
                        AllNCParams(idx, l),...
                        AllNCParamSEs(idx, l),...
                        'LineStyle', 'none');
                    
                    
                    hold on
                    
                    set(eb{idx,l}, 'color', colors(temp_idx,:), 'LineWidth', 1);
                    
                    set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                    
                    
                    prof{idx,l} = plot(1/(R*(AllNCTemperatures(idx, l)+273)),...
                        AllNCParams(idx, l), '.-', 'Color', colors(temp_idx,:), 'linestyle', 'none');
                    
                    set(eb{idx,l},'Visible','on'); %'off' or 'on'
                    set(prof{idx,l},'Visible','on'); %'off' or 'on'
                    
                end
            end
            
            
            
            hold off
            
            xlabel('1/(RT) (mol/kJ)')
            xlim([Plot2Xmin, Plot2Xmax])
            
            
            ylabel(ylab)
            ylim([LogPlotYmin, GlobalPlotYmax])
            
            set(FrameProfAx{l} , 'YScale', 'log')
            
            title(['AP: ', num2str(round(APbin, 3))])
        end
        if exist('PlotTitle', 'var')
            
            
            sgtitle({PlotTitle,...
                ['Nuclear Cycle ',num2str(NC)]})
            
        else
            
            sgtitle(['Nuclear Cycle ',num2str(NC)])
            
        end
        
        
        LegendAx = subplot(SubplotDims(1), SubplotDims(2), SubplotDims(1)*SubplotDims(2));
        for idx = 1:NumTemperatures
            plot([0], [0], '.-', 'Color', colors(idx,:), 'linestyle', 'none');
            hold on
        end
        hold off
        %xlim([-10, 10])
        %ylim([-10, 10])
        axis off
        
        legend_labels2 = {};
        profs2 = {};
        for i= 1:NumTemperatures
            current_temp = temperatures(i);
            legend_labels2{1, length(legend_labels2)+1} =[ num2str(current_temp),'°C'];
            
            %profs2{1, length(profs2)+1} = profB{plotted_temps(i)};
        end
        hlegend = legend(LegendAx, legend_labels2,...
            'FontSize', 12);%, 'Location', );
        hlegend.Position(1) = LegendAx.Position(1)+ LegendAx.Position(3)/2 - hlegend.Position(3)/2;
        hlegend.Position(2) = LegendAx.Position(2)+ LegendAx.Position(4)/2 - hlegend.Position(4)/2;
        
        %try
        
        
        
        saveas(FrameProfFig,[outdir3, filesep,...
            TraceType,'_LogSubplots', OutputString,  '_NC',num2str(NC),'.png']);
    end
    
end
%%
if ~SkipAPBinnedSubplots
    
    if  strcmp(lower(PlottingColors), 'gradient')
        [~, colors] = getColorPalettes();
        GradString = '';
    end
    
    PlottedParamSEs(R2s >= R2bound) = NaN;
    
    
    PlottedParamSEs((R2s < R2bound)| (PlottedParams > GlobalPlotYmax) | (PlottedParams < GlobalPlotYmin)) = NaN;
    PlottedParams((R2s < R2bound)| (PlottedParams > GlobalPlotYmax) | (PlottedParams < GlobalPlotYmin)) = NaN;
    
    WhereValid = squeeze(sum(~isnan(PlottedParams), 1));
    WhereValidAP = sum(WhereValid, 2).';
    ValidAPIndices = find(WhereValidAP);
    WhereValidNC = sum(WhereValid, 1);
    ValidNCIndices = find(WhereValidNC);
    
    
    SubplotDims = numSubplots(length(ValidAPIndices) + 1);
    
    SubFigDims = [0.9, 0.9*3072/1920*SubplotDims(1)/SubplotDims(2)*1.2];
    if SubFigDims(2) > 0.95
        SubFigDims = SubFigDims*0.95/SubFigDims(2);
        
    end
    SubFigDims = round(SubFigDims, 2);
    
    
    
    
    BinnedSEParams((Counts < this.MinimumBinCount)| (BinnedParams > GlobalPlotYmax) | (BinnedParams < GlobalPlotYmin)) = NaN;
    BinnedParams((Counts < this.MinimumBinCount)| (BinnedParams > GlobalPlotYmax) | (BinnedParams < GlobalPlotYmin)) = NaN;
    
    WhereBinnedValid = squeeze(sum(~isnan(BinnedParams), 1));
    WhereBinnedValidAP = sum(WhereBinnedValid, 2).';
    BinnedValidAPIndices = find(WhereBinnedValidAP);
    WhereBinnedValidNC = sum(WhereBinnedValid, 1);
    BinnedValidNCIndices = find(WhereBinnedValidNC);
    
    BinnedSubplotDims = numSubplots(length(BinnedValidAPIndices) + 1);
    % 3072×1920
    BinnedFigDims = [0.9, 0.9*3072/1920*BinnedSubplotDims(1)/BinnedSubplotDims(2)*1.2];
    BinnedFigDims(2) = min([BinnedFigDims(2), 0.95]);
    
    BinnedFigDims = round(BinnedFigDims, 2);
    
    
    for nc_idx=1:length(BinnedValidNCIndices)
        
        NC = BinnedValidNCIndices(nc_idx)+8;
        % Prepare Traces for plotting
        
        NCMaxParams = NaN(1,NumTemperatures);
        AllNCParams = NaN(NumTemperatures, length(BinnedValidAPIndices));
        AllNCParamSEs = NaN(NumTemperatures, length(BinnedValidAPIndices));
        AllNCTemperatures = NaN(NumTemperatures, length(BinnedValidAPIndices));
        AllNCTemperatureSEs = NaN(NumTemperatures, length(BinnedValidAPIndices));
        AllCounts = NaN(NumTemperatures, length(BinnedValidAPIndices));
        for i=1:NumTemperatures
            SetParams = BinnedParams(i,BinnedValidAPIndices,NC-8).';
            SetSEParams = BinnedSEParams(i,BinnedValidAPIndices,NC-8).';
            SetTemps = ParamTemperatures(i,BinnedValidAPIndices,NC-8).';
            SetSETemps = ParamSETemperatures(i,BinnedValidAPIndices,NC-8).';
            SetCounts = Counts(i,BinnedValidAPIndices,NC-8).';
            if ~all(isnan(SetSEParams))
                IncludedBins = find(~isnan(SetParams)  & (SetParams./SetSEParams >= 1)) ;
            else
                IncludedBins = find(~isnan(SetParams) );
            end
            if ~isempty(IncludedBins)
                AllNCParams(i,:) = SetParams;
                AllNCParamSEs(i,:) = SetSEParams;
                TempSEParams = SetSEParams;
                TempSEParams(isnan(TempSEParams)) = 0;
                NCMaxParams(i) = max(SetParams+TempSEParams);
                AllCounts(i,:) = SetCounts;
                AllNCTemperatures(i,:) =SetTemps;
                AllNCTemperatureSEs(i,:) =SetSETemps;
            end
        end
        
        if all(isnan(NCMaxParams))
            continue
        end
        outdir2 = [outdir,filesep,'Binned', OutputString];
        if ~exist(outdir2, 'dir')
            mkdir(outdir2)
        end
        
        outdir3 = [outdir2,filesep, datestr(now, 'yyyymmdd')];
        if ~exist(outdir3, 'dir')
            mkdir(outdir3)
        end
        
        clear eb prof FrameProfAx
        eb = cell(NumTemperatures, length(BinnedValidAPIndices));
        prof = cell(NumTemperatures, length(BinnedValidAPIndices));
        FrameProfAx = cell(1,  length(BinnedValidAPIndices));
        close all
        FrameProfFig = figure(1);
        set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.01, BinnedFigDims(1), BinnedFigDims(2)]);
        set(gcf,'color','w');
        
        for l = 1:length(BinnedValidAPIndices)
            APindex = BinnedValidAPIndices(l);
            APbin = APbins(APindex);
            
            
            
            
            
            
            
            
            
            if l == 1
                FrameProfAx{l} = subplot(BinnedSubplotDims(1), BinnedSubplotDims(2), l, gca);
            else
                FrameProfAx{l} = subplot(BinnedSubplotDims(1), BinnedSubplotDims(2), l);
            end
            
            
            
            for idx=1:NumTemperatures
                if ~isnan(AllNCParamSEs(idx, l) )
                    NoDataCondition = isnan(AllNCParams(idx, l)) | (AllNCParams(idx,l)/AllNCParamSEs(idx, l) <  1);
                else
                    NoDataCondition = isnan(AllNCParams(idx, l));
                end
                if  NoDataCondition
                    
                    eb{idx, l} = errorbar(AllNCTemperatures(idx, l), AllNCParams(idx, l),...
                        AllNCParamSEs(idx, l), AllNCParamSEs(idx, l),...
                        AllNCTemperatureSEs(idx, l), AllNCTemperatureSEs(idx, l),...
                        'LineStyle', 'none');
                    
                    hold on
                    
                    set(eb{idx,l}, 'color', colors(idx,:), 'LineWidth', 1);
                    
                    set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                    
                    
                    prof{idx,l} = plot(AllNCTemperatures(idx, l), AllNCParams(idx, l), '.-', 'Color', colors(idx,:), 'linestyle', 'none');
                    
                    
                    set(eb{idx,l},'Visible','off'); %'off' or 'on'
                    set(prof{idx,l},'Visible','off'); %'off' or 'on'
                    
                    
                else
                    
                    eb{idx, l} = errorbar(AllNCTemperatures(idx, l), AllNCParams(idx, l),...
                        AllNCParamSEs(idx, l), AllNCParamSEs(idx, l),...
                        AllNCTemperatureSEs(idx, l), AllNCTemperatureSEs(idx, l),...
                        'LineStyle', 'none');
                    
                    hold on
                    
                    set(eb{idx,l}, 'color', colors(idx,:), 'LineWidth', 1);
                    
                    set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                    
                    
                    prof{idx,l} = plot(AllNCTemperatures(idx, l), AllNCParams(idx, l), '.-', 'Color', colors(idx,:), 'linestyle', 'none');
                    
                    
                    set(eb{idx,l},'Visible','on'); %'off' or 'on'
                    set(prof{idx,l},'Visible','on'); %'off' or 'on'
                    
                end
            end
            
            
            
            hold off
            
            xlabel('Temperature (°C)')
            xlim([Subplot1DataXmin, Subplot1DataXmax])
            
            ylabel(ylab)
            ylim([GlobalPlotYmin, GlobalPlotYmax])
            title(['AP: ', num2str(round(APbin, 3))])
        end
        if exist('PlotTitle', 'var')
            
            
            sgtitle({PlotTitle,...
                ['Nuclear Cycle ',num2str(NC)]})
            
        else
            
            sgtitle(['Nuclear Cycle ',num2str(NC)])
            
        end
        
        
        LegendAx = subplot(BinnedSubplotDims(1), BinnedSubplotDims(2), BinnedSubplotDims(1)*BinnedSubplotDims(2));
        for idx = 1:NumTemperatures
            plot([0], [0], '.-', 'Color', colors(idx,:), 'linestyle', 'none');
            hold on
        end
        hold off
        %xlim([-10, 10])
        %ylim([-10, 10])
        axis off
        
        legend_labels2 = {};
        profs2 = {};
        for i= 1:NumTemperatures
            current_temp = temperatures(i);
            legend_labels2{1, length(legend_labels2)+1} =[ num2str(current_temp),'°C'];
            
            %profs2{1, length(profs2)+1} = profB{plotted_temps(i)};
        end
        hlegend = legend(LegendAx, legend_labels2,...
            'FontSize', 12);%, 'Location', );
        hlegend.Position(1) = LegendAx.Position(1)+ LegendAx.Position(3)/2 - hlegend.Position(3)/2;
        hlegend.Position(2) = LegendAx.Position(2)+ LegendAx.Position(4)/2 - hlegend.Position(4)/2;
        
        %try
        
        saveas(FrameProfFig,[outdir3, filesep,...
            TraceType,'_BinnedSubplots', OutputString,  '_NC',num2str(NC),'.png']);
        
    end
    for nc_idx=1:length(BinnedValidNCIndices)
        
        NC = BinnedValidNCIndices(nc_idx)+8;
        % Prepare Traces for plotting
        
        NCMaxParams = NaN(1,NumTemperatures);
        AllNCParams = NaN(NumTemperatures, length(BinnedValidAPIndices));
        AllNCParamSEs = NaN(NumTemperatures, length(BinnedValidAPIndices));
        AllNCTemperatures = NaN(NumTemperatures, length(BinnedValidAPIndices));
        AllNCTemperatureSEs = NaN(NumTemperatures, length(BinnedValidAPIndices));
        AllCounts = NaN(NumTemperatures, length(BinnedValidAPIndices));
        for i=1:NumTemperatures
            SetParams = BinnedParams(i,BinnedValidAPIndices,NC-8).';
            SetSEParams = BinnedSEParams(i,BinnedValidAPIndices,NC-8).';
            SetTemps = ParamTemperatures(i,BinnedValidAPIndices,NC-8).';
            SetSETemps = ParamSETemperatures(i,BinnedValidAPIndices,NC-8).';
            SetCounts = Counts(i,BinnedValidAPIndices,NC-8).';
            if ~all(isnan(SetSEParams))
                IncludedBins = find(~isnan(SetParams)  & (SetParams./SetSEParams >= 1)) ;
            else
                IncludedBins = find(~isnan(SetParams)  );
            end
            if ~isempty(IncludedBins)
                AllNCParams(i,:) = SetParams;
                AllNCParamSEs(i,:) = SetSEParams;
                TempSEParams = SetSEParams;
                TempSEParams(isnan(TempSEParams)) = 0;
                NCMaxParams(i) = max(SetParams +TempSEParams);
                AllCounts(i,:) = SetCounts;
                AllNCTemperatures(i,:) =SetTemps;
                AllNCTemperatureSEs(i,:) =SetSETemps;
            end
        end
        
        if all(isnan(NCMaxParams))
            continue
        end
        outdir2 = [outdir,filesep,'Binned', OutputString];
        if ~exist(outdir2, 'dir')
            mkdir(outdir2)
        end
        
        outdir3 = [outdir2,filesep, datestr(now, 'yyyymmdd')];
        if ~exist(outdir3, 'dir')
            mkdir(outdir3)
        end
        
        clear eb prof FrameProfAx
        eb = cell(NumTemperatures, length(BinnedValidAPIndices));
        prof = cell(NumTemperatures, length(BinnedValidAPIndices));
        FrameProfAx = cell(1,  length(BinnedValidAPIndices));
        close all
        FrameProfFig = figure(1);
        set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.01, BinnedFigDims(1), BinnedFigDims(2)]);
        set(gcf,'color','w');
        
        for l = 1:length(BinnedValidAPIndices)
            APindex = BinnedValidAPIndices(l);
            APbin = APbins(APindex);
            
            
            
            
            
            
            
            
            
            if l == 1
                FrameProfAx{l} = subplot(BinnedSubplotDims(1), BinnedSubplotDims(2), l, gca);
            else
                FrameProfAx{l} = subplot(BinnedSubplotDims(1), BinnedSubplotDims(2), l);
            end
            
            
            
            for idx=1:NumTemperatures
                if ~isnan(AllNCParamSEs(idx, l) )
                    NoDataCondition = isnan(AllNCParams(idx, l)) | (AllNCParams(idx,l)/AllNCParamSEs(idx, l) <  1);
                else
                    NoDataCondition = isnan(AllNCParams(idx, l));
                end
                
                if  NoDataCondition
                    
                    eb{idx, l} = errorbar(1/(R*(AllNCTemperatures(idx, l)+273)),...
                        AllNCParams(idx, l),...
                        AllNCParamSEs(idx, l), AllNCParamSEs(idx, l),...
                        sqrt(AllNCTemperatureSEs(idx,l)^2*(1/(R*(AllNCTemperatures(idx, l)+273)^2))^2),...
                        sqrt(AllNCTemperatureSEs(idx,l)^2*(1/(R*(AllNCTemperatures(idx, l)+273)^2))^2),...
                        'LineStyle', 'none');
                    
                    hold on
                    
                    set(eb{idx,l}, 'color', colors(idx,:), 'LineWidth', 1);
                    
                    set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                    
                    
                    prof{idx,l} = plot(1/(R*(AllNCTemperatures(idx, l)+273)),...
                        AllNCParams(idx, l), '.-', 'Color', colors(idx,:), 'linestyle', 'none');
                    
                    
                    set(eb{idx,l},'Visible','off'); %'off' or 'on'
                    set(prof{idx,l},'Visible','off'); %'off' or 'on'
                    
                    
                else
                    
                    eb{idx, l} = errorbar(1/(R*(AllNCTemperatures(idx, l)+273)),...
                        AllNCParams(idx, l),...
                        AllNCParamSEs(idx, l), AllNCParamSEs(idx, l),...
                        sqrt(AllNCTemperatureSEs(idx,l)^2*(1/(R*(AllNCTemperatures(idx, l)+273)^2))^2),...
                        sqrt(AllNCTemperatureSEs(idx,l)^2*(1/(R*(AllNCTemperatures(idx, l)+273)^2))^2),...
                        'LineStyle', 'none');
                    
                    
                    hold on
                    
                    set(eb{idx,l}, 'color', colors(idx,:), 'LineWidth', 1);
                    
                    set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                    
                    
                    prof{idx,l} = plot(1/(R*(AllNCTemperatures(idx, l)+273)),...
                        AllNCParams(idx, l), '.-', 'Color', colors(idx,:), 'linestyle', 'none');
                    
                    set(eb{idx,l},'Visible','on'); %'off' or 'on'
                    set(prof{idx,l},'Visible','on'); %'off' or 'on'
                    
                end
            end
            
            
            
            hold off
            
            xlabel('1/(RT) (mol/kJ)')
            xlim([Plot2Xmin, Plot2Xmax])
            
            
            ylabel(ylab)
            ylim([LogPlotYmin, GlobalPlotYmax])
            
            set(FrameProfAx{l} , 'YScale', 'log')
            
            title(['AP: ', num2str(round(APbin, 3))])
        end
        if exist('PlotTitle', 'var')
            
            
            sgtitle({PlotTitle,...
                ['Nuclear Cycle ',num2str(NC)]})
            
        else
            
            sgtitle(['Nuclear Cycle ',num2str(NC)])
            
        end
        
        
        LegendAx = subplot(BinnedSubplotDims(1), BinnedSubplotDims(2), BinnedSubplotDims(1)*BinnedSubplotDims(2));
        for idx = 1:NumTemperatures
            plot([0], [0], '.-', 'Color', colors(idx,:), 'linestyle', 'none');
            hold on
        end
        hold off
        %xlim([-10, 10])
        %ylim([-10, 10])
        axis off
        
        legend_labels2 = {};
        profs2 = {};
        for i= 1:NumTemperatures
            current_temp = temperatures(i);
            legend_labels2{1, length(legend_labels2)+1} =[ num2str(current_temp),'°C'];
            
            %profs2{1, length(profs2)+1} = profB{plotted_temps(i)};
        end
        hlegend = legend(LegendAx, legend_labels2,...
            'FontSize', 12);%, 'Location', );
        hlegend.Position(1) = LegendAx.Position(1)+ LegendAx.Position(3)/2 - hlegend.Position(3)/2;
        hlegend.Position(2) = LegendAx.Position(2)+ LegendAx.Position(4)/2 - hlegend.Position(4)/2;
        
        %try
        
        
        
        saveas(FrameProfFig,[outdir3, filesep,...
            TraceType,'_LogBinnedSubplots', OutputString,  '_NC',num2str(NC),'.png']);
    end
    
end


close all
end








