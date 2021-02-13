function PlotLTMActivationEnergiesVsAP(this, parameter, outdir, varargin)
%%

x = 1;
while x <= length(varargin)
    if strcmp(lower(varargin{x}), 'plottitle')
        PlotTitle = varargin{x+1};
        x = x+1;
    elseif strcmp(lower(varargin{x}), 'tracetype')
        TraceType = lower(varargin{x+1});
        x = x+1;
    end
    x = x+1;
end


if ~exist('TraceType', 'var')
    TraceType = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'Fluo3D')| strcmpi(TraceType, 'Unaligned3D')
    TraceType = 'Unaligned3D';
elseif strcmpi(TraceType, 'Fluo') | strcmpi(TraceType, 'Unaligned')
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


if strcmpi(parameter, 'TimeOns') | strcmpi(parameter, 'TimeOn')
    parameter = 'TimeOns';
elseif strcmpi(parameter, 'TimeOffs') | strcmpi(parameter, 'TimeOff')
    parameter = 'TimeOffs';
elseif strcmpi(parameter, 'TranscriptionWindows')| strcmpi(parameter, 'TranscriptionWindow')
    parameter = 'TranscriptionWindows';
elseif strcmpi(parameter, 'ElongationTimes') | strcmpi(parameter, 'ElongationTime')
    parameter = 'ElongationTimes';
elseif strcmpi(parameter, 'ElongationRates') | strcmpi(parameter, 'ElongationRate')
    parameter = 'ElongationRates';
elseif strcmpi(parameter, 'LoadingRates') | strcmpi(parameter, 'LoadingRate')
    parameter = 'LoadingRates';
elseif strcmpi(parameter, 'PostTranscriptionDurations') | strcmpi(parameter, 'PostTranscriptionDuration')
    parameter = 'PostTranscriptionDurations';
elseif strcmpi(parameter, 'PlateauHeights') | strcmpi(parameter, 'PlateauHeight')
    parameter = 'PlateauHeights';
elseif strcmpi(parameter, 'MaxFluos') | strcmpi(parameter, 'MaxFluo')
    parameter = 'MaxFluos';
end

%%

if ~exist(outdir, 'dir')
    mkdir(outdir)
end

Temp_obs = this.Temp_obs;
Temp_sp = this.Temp_sps;

%[~, colors] = getColorPalettes();


colors = customcolormap_preset('red-yellow-blue');
R_range = linspace(0, 1, size(colors, 1));
[bgcolors, ~] = getColorPalettes();


MarkerStyles = {'o', 'd', 's', '>', '^','p', 'h', '*', 'x'};


%%

R = this.R;




APResolution = this.Experiments{1}.APResolution;
APbins = 0:APResolution:1;



Nsigfigs = 3;



%% Load relevant parameters into memory
[ActivationEnergies, SEActivationEnergies, LogAs, SELogAs, R2s, PlotTitle] =...
    getActivationEnergyMatrices(this, parameter, TraceType);
% Calculate Plot Xlims

ObservedAPbins =  (sum(~isnan(ActivationEnergies),2) > 0).';
IncludedAPbins = APbins(ObservedAPbins);
NumAPbins = length(IncludedAPbins);
IncludedNCs = 9:14;
ObservedNCs =  (sum(~isnan(ActivationEnergies),1) > 0).';
IncludedNCs = IncludedNCs(ObservedNCs);

ActivationEnergies = ActivationEnergies(ObservedAPbins,:);
ActivationEnergies = ActivationEnergies(:,ObservedNCs);
SEActivationEnergies = SEActivationEnergies(ObservedAPbins,:);
SEActivationEnergies = SEActivationEnergies(:,ObservedNCs);
LogAs = LogAs(ObservedAPbins,:);
LogAs = LogAs(:,ObservedNCs);
SELogAs = SELogAs(ObservedAPbins,:);
SELogAs = SELogAs(:,ObservedNCs);
R2s = R2s(ObservedAPbins,:);
R2s = R2s(:,ObservedNCs);


PlotXmin = max([0, APResolution*(find(ObservedAPbins > 0, 1)-2)]);
PlotXmax = min([1, APResolution*(find(ObservedAPbins > 0, 1, 'last')+1)]);

GlobalPlotYmin = min(min(ActivationEnergies-SEActivationEnergies));
GlobalPlotYmax = max(max(ActivationEnergies+SEActivationEnergies));



%%
outdir2 = [outdir,filesep,'ActivationEnergies'];
if ~exist(outdir2, 'dir')
    mkdir(outdir2)
end
outdir3 = [outdir2,filesep,parameter];
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



%%
close all
for NCIndex=1:length(IncludedNCs)
    NC = IncludedNCs(NCIndex);
    
    
    eb = cell(1, NumAPbins);
    prof = cell(1, NumAPbins);
    close all
    FrameProfFig = figure(1);
    set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .9, .7]);
    set(gcf,'color','w');
    FrameProfAx = axes(FrameProfFig);
    for APIndex =1:NumAPbins
        ColIndex = find(abs(R_range-R2s(APIndex, NCIndex)) == min(abs(R_range-R2s(APIndex, NCIndex))));
        if isempty(ColIndex)
            continue
        end
        
        MarkerIndex = 1;
        if isempty(MarkerIndex)
            MarkerIndex = length(MarkerStyles);
        end
        
        eb{APIndex} = errorbar(IncludedAPbins(APIndex), ActivationEnergies(APIndex, NCIndex),  SEActivationEnergies(APIndex, NCIndex),...
            'vertical', 'LineStyle', 'none');
        hold on
        set(eb{APIndex}, 'color', colors(ColIndex,:), 'LineWidth', 1, 'CapSize', 0);
        
        set(get(get(eb{APIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        
        prof{APIndex} = plot(IncludedAPbins(APIndex), ActivationEnergies(APIndex, NCIndex), MarkerStyles{MarkerIndex},...
            'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:),...
            'MarkerSize', 10, 'linestyle', 'none');
        
        
        
    end
    
    grid on
    
    xlabel('Fraction Embryo Length')
    
    xlim([PlotXmin, PlotXmax])
    
    ylabel('Activation Energy (kj/mol)')
    ylim([GlobalPlotYmin, GlobalPlotYmax])
    
    FrameProfAx.YAxis.FontSize = 14;
    FrameProfAx.XAxis.FontSize = 14;
    
    if exist('PlotTitle', 'var')
        title(FrameProfAx, {PlotTitle,['Nuclear Cycle ',num2str(NC)]}, 'FontSize', 14)
    else
        title(FrameProfAx,  ['Nuclear Cycle ',num2str(NC)], 'FontSize', 14)
    end
    
    map = colormap(colors);
    h = colorbar;
    % %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
    hold off
    colorTitleHandle = get(h,'Title');
    titleString = 'R^2';
    set(colorTitleHandle ,'String',titleString, 'FontSize', 14);
    h.FontSize = 14;
    
    
    saveas(FrameProfFig,[outdir5, filesep,...
        parameter, '_ActivationEnergies_NC',num2str(NC),'.png']);
    
    %
end
%%

close all


eb = cell(length(IncludedNCs), NumAPbins);
prof = cell(length(IncludedNCs), NumAPbins);
plot_regions = cell(1, length(IncludedNCs));

FrameProfFig = figure(1);
set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .9, .7]);
set(gcf,'color','w');
FrameProfAx = axes(FrameProfFig);

for NCIndex=1:length(IncludedNCs)
    NC = IncludedNCs(NCIndex);
    bgXmin = PlotXmin+(NCIndex-1)*(PlotXmax-PlotXmin);
    bgXmax = PlotXmax+(NCIndex-1)*(PlotXmax-PlotXmin);
    plot_regions{NCIndex} = fill([bgXmin,bgXmin, bgXmax, bgXmax],...
        [GlobalPlotYmin,GlobalPlotYmax*1.1,GlobalPlotYmax*1.1, GlobalPlotYmin],...
        bgcolors(NCIndex+1,:));
    plot_regions{NCIndex}.FaceAlpha = 0.2;
    plot_regions{NCIndex}.EdgeColor = 'None';
    hold on
    
    for APIndex =1:NumAPbins
        ColIndex = find(abs(R_range-R2s(APIndex, NCIndex)) == min(abs(R_range-R2s(APIndex, NCIndex))));
        if isempty(ColIndex)
            continue
        end
        
        MarkerIndex = NCIndex;
        if isempty(MarkerIndex)
            MarkerIndex = length(MarkerStyles);
        end
        
        eb{NCIndex, APIndex} = errorbar(IncludedAPbins(APIndex)+(NCIndex-1)*(PlotXmax-PlotXmin), ActivationEnergies(APIndex, NCIndex),  SEActivationEnergies(APIndex, NCIndex),...
            'vertical', 'LineStyle', 'none');
        hold on
        set(eb{NCIndex,APIndex}, 'color', colors(ColIndex,:), 'LineWidth', 1, 'CapSize', 0);
        
        set(get(get(eb{NCIndex,APIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        
        prof{NCIndex,APIndex} = plot(IncludedAPbins(APIndex)+(NCIndex-1)*(PlotXmax-PlotXmin), ActivationEnergies(APIndex, NCIndex), MarkerStyles{MarkerIndex},...
            'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:),...
            'MarkerSize', 10, 'linestyle', 'none');
        
        
        
    end
    VertexPositions =  plot_regions{NCIndex}.Vertices;
    str = ['Nuclear Cycle ', num2str(NC)];
    text(VertexPositions(1,1)+0.025, GlobalPlotYmax*1.05,str,'FontSize',14);
    
    
    %
end

grid on

xlabel('Fraction Embryo Length')

xlim([PlotXmin, PlotXmax+(length(IncludedNCs)-1)*(PlotXmax-PlotXmin)])

ylabel('Activation Energy (kj/mol)')
ylim([GlobalPlotYmin, GlobalPlotYmax*1.1])

FrameProfAx.YAxis.FontSize = 14;
FrameProfAx.XAxis.FontSize = 14;


if exist('PlotTitle', 'var')
    title(FrameProfAx, PlotTitle, 'FontSize', 14)
    
end
map = colormap(colors);
h = colorbar;
% %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
hold off
colorTitleHandle = get(h,'Title');
titleString = 'R^2';
set(colorTitleHandle ,'String',titleString, 'FontSize', 14);
h.FontSize = 14;
%%
APticks = PlotXmin:APResolution:(PlotXmax+(length(IncludedNCs)-1)*(PlotXmax-PlotXmin));
xticks(APticks)
APlabels = cell(1,length(FrameProfAx.XTick));
for i = 1:length(APlabels)-1
    EvenDivisions = floor((APticks(i)-PlotXmin)/(PlotXmax-PlotXmin));
    APvalue = APticks(i)-EvenDivisions*(PlotXmax-PlotXmin);
    if ~mod(round(APvalue, 4), .1)
        APlabels{i} = num2str(APvalue);
    end
    %disp([num2str(APvalue), ', ', APlabels{i}])
end
xticklabels(APlabels)

%%


saveas(FrameProfFig,[outdir5, filesep,parameter, '_ActivationEnergies.png']);


%%
close all