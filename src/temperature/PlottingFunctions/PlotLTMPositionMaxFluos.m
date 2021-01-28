function PlotLTMPositionMaxFluos(this, outdir, varargin)
%%

% PlotTitle, PlottingColors, UseDifferentColors,
% UseDiffProfiles, UsePhysicalAPLength
UsePhysicalAPLength = false;
UseDifferentColors = true;
UseLines = true;
R2bound = 0.9;

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
    elseif strcmp(lower(varargin{x}), 'tracetype')
        TraceType = lower(varargin{x+1});
        x = x+1;
    elseif strcmp(lower(varargin{x}), 'R2bound')
        R2bound = lower(varargin{x+1});
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
elseif ~strcmp(lower(TraceType), 'fluo3d') & ~strcmp(lower(TraceType), 'fluo')  & ~strcmp(lower(TraceType), 'anaphasealigned')&...
        ~strcmp(lower(TraceType), 'anaphasealigned3d') & ~strcmp(lower(TraceType), 'tbinned') & ~strcmp(lower(TraceType), 'tbinned3d')
    error('Invalid choice of trace type. Can use either "fluo", "fluo3d", "anaphasealigned", or "anaphasealigned3d".') % change to error
end

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

temperatures = flip(unique(this.Temp_sps(this.ProcessedExperiments)));

R = 8.314*10^(-3); % kJ * K^(-1)*mol^(-1)

APResolution = this.Experiments{1}.APResolution;
APbins = 0:APResolution:1;
NumSets = length(this.ProcessedExperiments);

legend_labels = this.LegendLabels(this.ProcessedExperiments);
MinimumTraceCount = this.MinimumTraceCount;

%% Load relevant parameters into memory
InitiationRates = getTrapezoidParameters(this, 'MeanInitiationRates', TraceType, false);
SEInitiationRates = getTrapezoidParameters(this, 'MeanInitiationRates', TraceType, true);
TimeOns = getTrapezoidParameters(this, 'TimeOns', TraceType, false);
SETimeOns = getTrapezoidParameters(this, 'TimeOns', TraceType, true);
TimeOffs = getTrapezoidParameters(this, 'TimeOffs', TraceType, false);
SETimeOffs = getTrapezoidParameters(this, 'TimeOffs', TraceType, true);
ElongationTimes = getTrapezoidParameters(this, 'ElongationTimes', TraceType, false);
SEElongationTimes = getTrapezoidParameters(this, 'ElongationTimes', TraceType, true);
UnloadingRates = getTrapezoidParameters(this, 'UnloadingRates', TraceType, false);
SEUnloadingRates = getTrapezoidParameters(this, 'UnloadingRates', TraceType, true);
R2s = getTrapezoidParameters(this, 'R2s', TraceType);
[MaxFluos, SEMaxFluos] = getMaxFluoMatForPlotting(this, TraceType);
%% function specific calculations
[~, APMaxPositions ] = max(MaxFluos, [], 2);
APMaxPositions = squeeze(APMaxPositions);
APMaxPositions(APMaxPositions == 1) = NaN;
PlottedParams = (APMaxPositions-1)*APResolution;

ylab = 'Fraction of Embryo Length for Max. Spot Fluo.';
OutputString = 'PositionMaxFluos';

GlobalPlotYmax = 1;
GlobalPlotYmin = 0;





for nc_idx=1:length(this.IncludedNCs)
    
    NC = this.IncludedNCs(nc_idx);
    % Prepare Traces for plotting
    
    

    AllNCParams = PlottedParams(this.ProcessedExperiments,NC-8).';

    if all(isnan(AllNCParams))
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

        clear  prof
        plotted_temps = zeros(1, length(temperatures));
        %         plotted_temps = zeros(1, length(temperatures));
        
        close all
        
        
        
        
        FrameProfFig = figure(1);
        set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .5, .4]);
        set(gcf,'color','w');
        FrameProfAx = gca;
        for idx =1:length(Temp_sp)
            
            temp_idx = find(temperatures == Temp_sp(idx));
            
        
            
            prof{idx} = plot(APbins, ones(1, length(APbins)), '.-',...
                'Color', colors(temp_idx,:), 'linestyle', 'none',...
                'MarkerSize', 20);
            
            hold on
     
            set(prof{idx},'Visible','off'); %'off' or 'on'
        end
        
        
        
        hold off
        
        xlabel('Temperature (°C)')
        xlim([15, 30])
        
        ylabel(ylab)
        ylim([GlobalPlotYmin, GlobalPlotYmax])%GlobalPlotYmax*1.05])
        
        
       
   
        for idx=1:length(Temp_sp)
            
            if  isnan(AllNCParams(idx)) 
                FrameProfAx.Children(end-idx+1).XData = 25;
                FrameProfAx.Children(end-idx+1).YData = 1;
            
                set(FrameProfAx.Children(end-idx+1),'Visible','off'); %'off' or 'on'
              
                
            else
 
                FrameProfAx.Children(end-idx+1).YData = AllNCParams(idx);
                FrameProfAx.Children(end-idx+1).XData = Temp_obs(idx);
               
                set(FrameProfAx.Children(end-idx+1),'Visible','on'); %'off' or 'on'
                
                temp_idx = find(temperatures == Temp_sp(idx));
                if plotted_temps(temp_idx) == 0
                    plotted_temps(temp_idx) = idx;
                end
                
                
            end
        end
        
     
        legend_labels2 = {};
        profs2 = {};
        for i= 1:length(plotted_temps)
            if plotted_temps(i) > 0
                current_temp = Temp_sp(plotted_temps(i));
                legend_labels2{1, length(legend_labels2)+1} =[ num2str(current_temp),'°C'];
                profs2{1, length(profs2)+1} = prof{plotted_temps(i)};
            end
        end
        hlegend = legend(FrameProfAx, [profs2{:}], legend_labels2, 'Location', 'northeast',...
            'FontSize', 10);
        %try
        if exist('PlotTitle', 'var')
            
            
            sgtitle({PlotTitle,...
                ['Nuclear Cycle ',num2str(NC)]})
            
        else
            
            sgtitle(['Nuclear Cycle ',num2str(NC)])
            
        end
        
        
        saveas(FrameProfFig,[outdir3, filesep,...
            'NC',num2str(NC), OutputString, '_',TraceType,'.png']);

end
end