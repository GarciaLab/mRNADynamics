function PlotLTMTrapParamsVsTemp(this, parameter, outdir, varargin)
%%

% PlotTitle, PlottingColors, UseDifferentColors,
% UseDiffProfiles, UsePhysicalAPLength
IncludeFits = true;

x = 1;
while x <= length(varargin)
    if strcmpi(varargin{x}, 'plottitle')
        PlotTitle = varargin{x+1};
        x = x+1;
    elseif strcmpi(varargin{x}, 'plottingcolors')
        PlottingColors = varargin{x+1};
        x = x+1;
    elseif strcmpi(varargin{x}, 'tracetype')
        TraceType = lower(varargin{x+1});
        x = x+1;
    elseif strcmpi(varargin{x}, 'ExcludeFits')
        IncludeFits = false;
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
    TraceType = 'Unaligned3D';
elseif strcmpi(TraceType, 'Fluo')
    TraceType = 'Unaligned';
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
        ~strcmpi(parameter, 'LoadingRates')& ~strcmpi(parameter, 'PostTranscriptionDuration')
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

R = this.R;
R2bound = this.R2bound;

temperatures = flip(unique(this.Temp_sps(this.ProcessedExperiments)));
NumTemperatures = length(temperatures);


APResolution = this.Experiments{1}.APResolution;
APbins = 0:APResolution:1;
NumAPbins = length(APbins);

NumSets = length(this.ExperimentPrefixes);
UseSet = ismember(1:NumSets, this.ProcessedExperiments);
Nsigfigs = 3;
legend_labels = this.LegendLabels;


MarkerStyles = {'o', 'd', 's', '>', '^','p', 'h', '*', 'x'};

TempMatches = cell(1, NumTemperatures);

for t_index = 1:NumTemperatures
    TempMatches{t_index} = find((this.Temp_sps == temperatures(t_index)) & UseSet);
end


%% Load relevant parameters into memory
[PlottedParams, PlottedParamSEs,R2s, ylab,OutputString,GlobalPlotYmax,GlobalPlotYmin,LogPlotYmin] = ...
    getPlottingVariables(this, parameter,  TraceType, R2bound);

% Calculate Plot Xlims

TemperatureVector = 1./(R*(temperatures + 273));
Subplot2DataXmin = min(TemperatureVector);
Subplot2DataXmax = max(TemperatureVector);
Subplot2Xspan = Subplot2DataXmax-Subplot2DataXmin;
Plot2Xmin = Subplot2DataXmin - Subplot2Xspan*.05;
Plot2Xmax = Subplot2DataXmax + Subplot2Xspan*.05;



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


close all
FrameProfFig = figure(1);
set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .9, .7]);
set(gcf,'color','w');


eb = cell(1, NumSets);
prof = cell(1, NumSets);
FrameProfAx = subplot(1, 3, 1, gca);
for SetIndex =1:NumSets
    if strcmp(lower(PlottingColors), 'gradient')
        ColIndex = find(abs(Temp_range-Temp_obs(SetIndex)) == min(abs(Temp_range-Temp_obs(SetIndex))));
    else
        ColIndex = find(temperatures == Temp_sp(SetIndex));
    end
    
    MarkerIndex = find(TempMatches{ColIndex} == SetIndex);
    if isempty(MarkerIndex)
        MarkerIndex = length(MarkerStyles);
    end
    
    eb{SetIndex} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)),...
        'vertical', 'LineStyle', 'none');
    hold on
    if strcmp(lower(PlottingColors), 'gradient')
        set(eb{SetIndex}, 'color', colors(ColIndex,:), 'LineWidth', 1, 'CapSize', 0);
    else
        set(eb{SetIndex}, 'color', colors(ColIndex,:), 'LineWidth', 1, 'CapSize', 0);
    end
    set(get(get(eb{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    if UseLines
        prof{SetIndex} = plot(APbins, ones(1, length(APbins)), [MarkerStyles{MarkerIndex}, '-'],...
            'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:), 'Color', colors(ColIndex,:),...
            'MarkerSize', 10);
        
    else
        prof{SetIndex} = plot(APbins, ones(1, length(APbins)), MarkerStyles{MarkerIndex},...
            'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:),...
            'MarkerSize', 10, 'linestyle', 'none');
        
    end
    
    set(eb{SetIndex},'Visible','off'); %'off' or 'on'
    set(prof{SetIndex},'Visible','off'); %'off' or 'on'
    if ~ismember(SetIndex, this.ProcessedExperiments)
        set(get(get(prof{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    end
    
end

grid on
hold off

xlabel('Temperature (°C)')
xlim([15, 30])

ylabel(ylab)
ylim([GlobalPlotYmin, GlobalPlotYmax])

FrameProfAx.YAxis.FontSize = 14;
FrameProfAx.XAxis.FontSize = 14;


%%
eb2 = cell(1, NumSets);
prof2 = cell(1, NumSets);
FrameProfAx2 = subplot(1, 3, 2);

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
for SetIndex =1:NumSets
    if strcmp(lower(PlottingColors), 'gradient')
        ColIndex = find(abs(Temp_range-Temp_obs(SetIndex)) == min(abs(Temp_range-Temp_obs(SetIndex))));
    else
        ColIndex = find(temperatures == Temp_sp(SetIndex));
    end
    
    MarkerIndex = find(TempMatches{ColIndex} == SetIndex);
    if isempty(MarkerIndex)
        MarkerIndex = length(MarkerStyles);
    end
    
    eb2{SetIndex} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)),...
        'vertical', 'LineStyle', 'none');
    hold on
    if strcmp(lower(PlottingColors), 'gradient')
        set(eb2{SetIndex}, 'color', colors(ColIndex,:), 'LineWidth', 1, 'CapSize', 0);
    else
        set(eb2{SetIndex}, 'color', colors(ColIndex,:), 'LineWidth', 1, 'CapSize', 0);
    end
    set(get(get(eb2{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    if UseLines
        prof2{SetIndex} = plot(APbins, ones(1, length(APbins)), [MarkerStyles{MarkerIndex}, '-'],...
            'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:), 'Color', colors(ColIndex,:),...
            'MarkerSize', 10);
        
    else
        prof2{SetIndex} = plot(APbins, ones(1, length(APbins)), MarkerStyles{MarkerIndex},...
            'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:),...
            'MarkerSize', 10, 'linestyle', 'none');
        
    end
    
    set(eb2{SetIndex},'Visible','off'); %'off' or 'on'
    set(prof2{SetIndex},'Visible','off'); %'off' or 'on'
    if ~ismember(SetIndex, this.ProcessedExperiments)
        set(get(get(prof2{SetIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    end
    
end

grid on
hold off
xlabel('1/(RT) (mol/kJ)')
xlim([Plot2Xmin, Plot2Xmax])

ylabel(ylab)
ylim([LogPlotYmin, GlobalPlotYmax])

set(FrameProfAx2, 'YScale', 'log')

FrameProfAx2.YAxis.FontSize = 14;
FrameProfAx2.XAxis.FontSize = 14;
sgt = sgtitle('Nuclear Cycle');

LegendAx = subplot(1, 3, 3);


pos1 = get(FrameProfAx, 'position');
pos1(1) = 0.05;
pos1(3) = .3;
set(FrameProfAx, 'position', pos1);

pos2 = get(FrameProfAx2, 'position');
pos2(1) = 0.4;
pos2(3) = .3;
set(FrameProfAx2, 'position', pos2);

pos3 = get(LegendAx, 'position');
pos3(1) = 0.75;
pos3(3) = .2;
set(LegendAx, 'position', pos3);




sgt.FontSize = 14;



%%


for NCIndex=1:length(this.IncludedNCs)
    NC = this.IncludedNCs(NCIndex);
    
    
    % Prepare Traces for plotting
    
    NCMaxParams = NaN(1, NumSets);
    AllNCParams = NaN(NumSets, NumAPbins);
    AllNCParamSEs = NaN(NumSets, NumAPbins);
    AllR2s = NaN(NumSets, NumAPbins);
    for SetIndex=1:NumSets
        SetParams = PlottedParams(SetIndex,:,NC-8).';
        SetSEParams = PlottedParamSEs(SetIndex,:,NC-8).';
        SetR2s = R2s(SetIndex,:,NC-8).';
        if ~all(isnan(SetSEParams))
            IncludedBins = find(~isnan(SetParams) & (SetR2s > R2bound) & (SetParams./SetSEParams >= 1)) ;
        else
            IncludedBins = find(~isnan(SetParams) & (SetR2s > R2bound)) ;
        end
        if ~isempty(IncludedBins)
            AllNCParams(SetIndex,:) = SetParams;
            AllNCParamSEs(SetIndex,:) = SetSEParams;
            TempSEParams = SetSEParams;
            TempSEParams(isnan(TempSEParams)) = 0;
            NCMaxParams(SetIndex) = max(SetParams+TempSEParams);
            AllR2s(SetIndex,:) = SetR2s;
        end
    end
    
    if all(isnan(NCMaxParams))
        continue
    end
    
    
    
    for APindex = 1:NumAPbins
        APbin = APbins(APindex);
        
        PlottedSets = zeros(1, NumSets, 'logical');
        for SetIndex=1:NumSets
            if ~ismember(SetIndex, this.ProcessedExperiments)
                continue
            end
            if ~isnan(AllNCParamSEs(SetIndex, APindex))
                NoDataCondition = isnan(AllNCParams(SetIndex, APindex)) | (AllNCParams(SetIndex, APindex)./AllNCParamSEs(SetIndex, APindex) <  1) | ...
                    (AllR2s(SetIndex,APindex)< R2bound);
            else
                NoDataCondition = isnan(AllNCParams(SetIndex, APindex)) | (AllR2s(SetIndex,APindex)< R2bound);
            end
            
            if  NoDataCondition
                FrameProfAx.Children(end-(2*(SetIndex-1)+1)).XData = 25;
                FrameProfAx.Children(end-(2*(SetIndex-1)+1)).YData = 1;
                FrameProfAx.Children(end-(2*(SetIndex-1))).XData = 25;
                FrameProfAx.Children(end-(2*(SetIndex-1))).YData = 1;
                FrameProfAx.Children(end-(2*(SetIndex-1))).YPositiveDelta = 1;
                FrameProfAx.Children(end-(2*(SetIndex-1))).YNegativeDelta = 1;
                set(FrameProfAx.Children(end-(2*(SetIndex-1)+1)),'Visible','off'); %'off' or 'on'
                set(FrameProfAx.Children(end-(2*(SetIndex-1))),'Visible','off'); %'off' or 'on'
                
                if IncludeFits
                    FrameProfAx2.Children(end-(2*(SetIndex-1)+1)-2).XData = 1/(R*(25+273));
                    FrameProfAx2.Children(end-(2*(SetIndex-1)+1)-2).YData = 1;
                    FrameProfAx2.Children(end-(2*(SetIndex-1))-2).XData = 1/(R*(25+273));
                    FrameProfAx2.Children(end-(2*(SetIndex-1))-2).YData = 1;
                    FrameProfAx2.Children(end-(2*(SetIndex-1))-2).YPositiveDelta = 1;
                    FrameProfAx2.Children(end-(2*(SetIndex-1))-2).YNegativeDelta = 1;
                    set(FrameProfAx2.Children(end-(2*(SetIndex-1)+1)-2),'Visible','off'); %'off' or 'on'
                    set(FrameProfAx2.Children(end-(2*(SetIndex-1))-2),'Visible','off'); %'off' or 'on'
                else
                    FrameProfAx2.Children(end-(2*(SetIndex-1)+1)).XData = 1/(R*(25+273));
                    FrameProfAx2.Children(end-(2*(SetIndex-1)+1)).YData = 1;
                    FrameProfAx2.Children(end-(2*(SetIndex-1))).XData = 1/(R*(25+273));
                    FrameProfAx2.Children(end-(2*(SetIndex-1))).YData = 1;
                    FrameProfAx2.Children(end-(2*(SetIndex-1))).YPositiveDelta = 1;
                    FrameProfAx2.Children(end-(2*(SetIndex-1))).YNegativeDelta = 1;
                    set(FrameProfAx2.Children(end-(2*(SetIndex-1)+1)),'Visible','off'); %'off' or 'on'
                    set(FrameProfAx2.Children(end-(2*(SetIndex-1))),'Visible','off'); %'off' or 'on'
                end
            else
                PlottedSets(SetIndex) = 1;
                FrameProfAx.Children(end-(2*(SetIndex-1)+1)).YData = AllNCParams(SetIndex, APindex);
                FrameProfAx.Children(end-(2*(SetIndex-1)+1)).XData =Temp_obs(SetIndex);
                FrameProfAx.Children(end-(2*(SetIndex-1))).YData =  AllNCParams(SetIndex, APindex);
                FrameProfAx.Children(end-(2*(SetIndex-1))).XData = Temp_obs(SetIndex);
                FrameProfAx.Children(end-(2*(SetIndex-1))).YPositiveDelta = AllNCParamSEs(SetIndex,APindex);
                FrameProfAx.Children(end-(2*(SetIndex-1))).YNegativeDelta  = AllNCParamSEs(SetIndex,APindex);
                
                set(FrameProfAx.Children(end-(2*(SetIndex-1)+1)),'Visible','on'); %'off' or 'on'
                set(FrameProfAx.Children(end-(2*(SetIndex-1))),'Visible','on'); %'off' or 'on'
                
                if IncludeFits
                    FrameProfAx2.Children(end-(2*(SetIndex-1)+1)-2).XData =1/(R*(Temp_obs(SetIndex)+273));
                    FrameProfAx2.Children(end-(2*(SetIndex-1)+1)-2).YData = AllNCParams(SetIndex, APindex);
                    FrameProfAx2.Children(end-(2*(SetIndex-1))-2).XData = 1/(R*(Temp_obs(SetIndex)+273));
                    FrameProfAx2.Children(end-(2*(SetIndex-1))-2).YData = AllNCParams(SetIndex, APindex);
                    FrameProfAx2.Children(end-(2*(SetIndex-1))-2).YPositiveDelta = AllNCParamSEs(SetIndex,APindex);
                    FrameProfAx2.Children(end-(2*(SetIndex-1))-2).YNegativeDelta = AllNCParamSEs(SetIndex,APindex);
                    set(FrameProfAx2.Children(end-(2*(SetIndex-1)+1)-2),'Visible','on'); %'off' or 'on'
                    set(FrameProfAx2.Children(end-(2*(SetIndex-1))-2),'Visible','on'); %'off' or 'on'
                else
                    FrameProfAx2.Children(end-(2*(SetIndex-1)+1)).YData = AllNCParams(SetIndex, APindex);
                    FrameProfAx2.Children(end-(2*(SetIndex-1)+1)).XData =1/(R*(Temp_obs(idx)+273));
                    FrameProfAx2.Children(end-(2*(SetIndex-1))).YData =  AllNCParams(SetIndex, APindex);
                    FrameProfAx2.Children(end-(2*(SetIndex-1))).XData = 1/(R*(Temp_obs(idx)+273));
                    FrameProfAx2.Children(end-(2*(SetIndex-1))).YPositiveDelta = AllNCParamSEs(SetIndex,APindex);
                    FrameProfAx2.Children(end-(2*(SetIndex-1))).YNegativeDelta  = AllNCParamSEs(SetIndex,APindex);
                    
                    set(FrameProfAx2.Children(end-(2*(SetIndex-1)+1)),'Visible','on'); %'off' or 'on'
                    set(FrameProfAx2.Children(end-(2*(SetIndex-1))),'Visible','on'); %'off' or 'on'
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
                    'FontSize', 14);
                
                
                
            else
                set(FrameProfAx2.Children(end),'Visible','off');
                set(FrameProfAx2.Children(end-1),'Visible','off');
                subplot2_labels{2} = '';
                subplot2_labels{1} = 'No Fit Info available';
                if exist('hlegend2', 'var')
                    set(hlegend2,'visible','off')
                end
            end
        else
            set(FrameProfAx2.Children(end),'Visible','off');
            set(FrameProfAx2.Children(end-1),'Visible','off');
            
            subplot2_labels{2} = '';
            subplot2_labels{1} = 'No Fit Info available';
            if exist('hlegend2', 'var')
                    set(hlegend2,'visible','off')
                end
            
            
        end
        
        if all(~PlottedSets)
            continue
        end
        
        %try
        if exist('PlotTitle', 'var')
            
            
            sgtitle({PlotTitle,...
                ['Nuclear Cycle ',num2str(NC), ', Fraction Embryo Length: ', num2str(APbin)]})
            
        else
            
            sgtitle(['Nuclear Cycle ',num2str(NC), ', Fraction Embryo Length: ', num2str(APbin)])
            
        end
        
        legend_profs = {};
        legend_labels3 = {};
        for SetIndex = 1:NumSets
            if ~PlottedSets(SetIndex)
               continue
            end
            ColIndex = find(temperatures == Temp_sp(SetIndex));
            MarkerIndex = find(TempMatches{ColIndex} == SetIndex);
            
            legend_profs{1, length(legend_profs)+1} = plot([0], [0], [MarkerStyles{MarkerIndex}, '-'],...
            'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:), 'Color', colors(ColIndex,:),...
            'MarkerSize', 10);
            hold on 
            
            legend_labels3{1, length(legend_labels3)+1} = legend_labels{SetIndex};
            
            
        end
        hold off
        axis off
        
        
        hlegend = legend(LegendAx,  legend_labels3,...
            'FontSize',14, 'Location', 'east');
        hlegend.Position(1) = LegendAx.Position(1)+ LegendAx.Position(3)/2 - hlegend.Position(3)/2;
        hlegend.Position(2) = LegendAx.Position(2)+ LegendAx.Position(4)/2 - hlegend.Position(4)/2;
        
       
        
        
        saveas(FrameProfFig,[outdir3, filesep,...
            TraceType,'_', OutputString,  '_NC',num2str(NC), 'AP', num2str(APindex), '.png']);
        
    end
end
close all