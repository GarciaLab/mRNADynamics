function PlotHMMActivationEnergiesVsAP(this, parameter, outdir, varargin)
%%
UseRescaledFluo = false;
UseRescaledTiming = false;
x = 1;
while x <= length(varargin)
    if strcmpi(varargin{x}, 'userescaledfluo')
        UseRescaledFluo = true;
    elseif strcmpi(varargin{x}, 'userescaledtiming')
        UseRescaledTiming = true;
    end
    x = x+1;
end

if strcmpi(parameter, 'InitiationRates') | strcmpi(parameter, 'Loadingrates')
    parameter = 'InitiationRates';
elseif strcmpi(parameter, 'durations') | strcmpi(parameter, 'burstdurations')
    parameter = 'Durations';
elseif strcmpi(parameter, 'frequencies')| strcmpi(parameter, 'burstfrequencies')
    parameter = 'Frequencies';
end


%%

if ~exist(outdir, 'dir')
    mkdir(outdir)
end



%[~, colors] = getColorPalettes();


colors = customcolormap_preset('red-yellow-blue');
R_range = linspace(0, 1, size(colors, 1));
[bgcolors, ~] = getColorPalettes();
colors = brewermap(10,'Spectral');


MarkerStyles = {'o', 'd', 's', '>', '^','p', 'h', '*', 'x'};


%%

R =  8.314*10^(-3);




APResolution = 1/(size(this.InitiationRates, 3)-1);
APbins = 0:APResolution:1;



Nsigfigs = 3;



%% Load relevant parameters into memory
[ActivationEnergies, SEActivationEnergies, LogAs, SELogAs, R2s, PlotTitle] =...
    getHMMActivationEnergyMatrices(this, parameter, UseRescaledFluo, UseRescaledTiming);
% Calculate Plot Xlims



ObservedAPbins =  (sum(~isnan(ActivationEnergies),1) > 0);
IncludedAPbins = APbins(ObservedAPbins);
NumAPbins = length(IncludedAPbins);
IncludedTimeBins  = this.TimeVector;
ObservedTimeBins = (sum(~isnan(ActivationEnergies),2) > 0).';
IncludedTimeBins = IncludedTimeBins(ObservedTimeBins);
NumTimeBins= length(IncludedTimeBins);

ActivationEnergies = ActivationEnergies(:,ObservedAPbins);
ActivationEnergies = ActivationEnergies(ObservedTimeBins,:);
SEActivationEnergies = SEActivationEnergies(:,ObservedAPbins);
SEActivationEnergies = SEActivationEnergies(ObservedTimeBins,:);
LogAs = LogAs(:,ObservedAPbins);
LogAs = LogAs(ObservedTimeBins,:);
SELogAs = SELogAs(:,ObservedAPbins);
SELogAs = SELogAs(ObservedTimeBins,:);
R2s = R2s(:,ObservedAPbins);
R2s = R2s(ObservedTimeBins,:);


PlotXmin = max([0, APResolution*(find(ObservedAPbins > 0, 1)-2)]);
PlotXmax = min([1, APResolution*(find(ObservedAPbins > 0, 1, 'last')+1)]);


if strcmpi(parameter, 'durations')
    ActivationEnergies = -ActivationEnergies;
end
GlobalPlotYmin = min(min(ActivationEnergies-SEActivationEnergies));
GlobalPlotYmax = max(max(ActivationEnergies+SEActivationEnergies));
%GlobalPlotYmax = 20;


%%
outdir2 = [outdir,filesep,'ActivationEnergies'];
if ~exist(outdir2, 'dir')
    mkdir(outdir2)
end
if ~UseRescaledFluo & ~UseRescaledTiming
    outdir3 = [outdir2,filesep,parameter];
elseif ~UseRescaledTiming
    outdir3 = [outdir2,filesep,'RescaledFluo',parameter];
elseif ~UseRescaledFluo
    outdir3 = [outdir2,filesep,'RescaledTiming',parameter];
else
    outdir3 = [outdir2,filesep,'RescaledFluoRescaledTiming',parameter];
end
if ~exist(outdir3, 'dir')
    mkdir(outdir3)
end


outdir5 = [outdir3, filesep, datestr(now, 'yyyymmdd')];
if ~exist(outdir5, 'dir')
    mkdir(outdir5)
end



%%
close all
for TimeIndex=1:length(IncludedTimeBins)
    TimeBin = IncludedTimeBins(TimeIndex);
    
    
    eb = cell(1, NumAPbins);
    prof = cell(1, NumAPbins);
    close all
    FrameProfFig = figure(1);
    set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .9, .7]);
    set(gcf,'color','w');
    FrameProfAx = axes(FrameProfFig);
    for APIndex =1:NumAPbins
        ColIndex = TimeIndex;
        if isempty(ColIndex)
            continue
        end
        ColIndex = 1;
        MarkerIndex = 1;
        if isempty(MarkerIndex)
            MarkerIndex = length(MarkerStyles);
        end
        
        eb{APIndex} = errorbar(IncludedAPbins(APIndex), ActivationEnergies(TimeIndex, APIndex),  SEActivationEnergies(TimeIndex, APIndex),...
            'vertical', 'LineStyle', 'none');
        hold on
        set(eb{APIndex}, 'color', colors(ColIndex,:), 'LineWidth', 1, 'CapSize', 0);
        
        set(get(get(eb{APIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        
        prof{APIndex} = plot(IncludedAPbins(APIndex), ActivationEnergies(TimeIndex, APIndex), MarkerStyles{MarkerIndex},...
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
        title(FrameProfAx, {PlotTitle,[num2str(IncludedTimeBins(TimeIndex)),' min']}, 'FontSize', 14)
    else
        title(FrameProfAx,  [num2str(IncludedTimeBins(TimeIndex)),' min'], 'FontSize', 14)
    end
    
    %     map = colormap(colors);
    %     h = colorbar;
    %     % %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
    %     hold off
    %     colorTitleHandle = get(h,'Title');
    %     titleString = 'R^2';
    %     set(colorTitleHandle ,'String',titleString, 'FontSize', 14);
    %     h.FontSize = 14;
    
    
    saveas(FrameProfFig,[outdir5, filesep,...
        parameter, '_ActivationEnergies_TimeBin',num2str(TimeIndex),'.png']);
    
    %
end
%%
close all
eb = cell(NumTimeBins, NumAPbins);
prof = cell(NumTimeBins, NumAPbins);
close all
FrameProfFig = figure(1);
set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.05, .9, .7]);
set(gcf,'color','w');
TimeLabels = 1:length(ObservedTimeBins);
        TimeLabels = TimeLabels(ObservedTimeBins);
colors = brewermap(max(TimeLabels),'Spectral');
colormap(colors);
FrameProfAx = axes(FrameProfFig);
hold on
for TimeIndex=1:length(IncludedTimeBins)
    TimeBin = IncludedTimeBins(TimeIndex);
    
    
    
    for APIndex =1:NumAPbins
        
        ColIndex = TimeLabels(TimeIndex);
        if isempty(ColIndex)
            continue
        end
        %         ColIndex = 1;
        MarkerIndex = 1;
        if isempty(MarkerIndex)
            MarkerIndex = length(MarkerStyles);
        end
        
        eb{APIndex} = errorbar(IncludedAPbins(APIndex), ActivationEnergies(TimeIndex, APIndex),  SEActivationEnergies(TimeIndex, APIndex),...
            'vertical', 'LineStyle', 'none');
        hold on
        set(eb{APIndex}, 'color', colors(ColIndex,:), 'LineWidth', 1, 'CapSize', 0);
        
        set(get(get(eb{APIndex}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        
        prof{APIndex} = plot(IncludedAPbins(APIndex), ActivationEnergies(TimeIndex, APIndex), MarkerStyles{MarkerIndex},...
            'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:),...
            'MarkerSize', 10, 'linestyle', 'none');
        
        
        
    end
    
    
    %
end
grid on

xlabel('Fraction Embryo Length')

xlim([PlotXmin, PlotXmax])

ylabel('Activation Energy (kj/mol)')
ylim([GlobalPlotYmin, GlobalPlotYmax])
%ylim([-60 60])

FrameProfAx.YAxis.FontSize = 16;
FrameProfAx.XAxis.FontSize = 16;

if exist('PlotTitle', 'var')
    title(FrameProfAx, {PlotTitle,[num2str(IncludedTimeBins(TimeIndex)),' min']}, 'FontSize', 16)
else
    title(FrameProfAx,  [num2str(IncludedTimeBins(TimeIndex)),' min'], 'FontSize', 16)
end
set(FrameProfAx, 'Clim', [1, find(ObservedTimeBins, 1, 'last')+1])
h = colorbar;
h.Ticks = (1+0.5):(find(ObservedTimeBins, 1, 'last')+0.5);
h.TickLabels = this.TimeVector(1:max(TimeLabels));
%     map = colormap(colors);
%     h = colorbar;
%     % %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
%     hold off
%     colorTitleHandle = get(h,'Title');
%     titleString = 'R^2';
%     set(colorTitleHandle ,'String',titleString, 'FontSize', 14);
%     h.FontSize = 14;
if ~UseRescaledTiming
ylabel(h,'Time Cohort (min)')
else
    ylabel(h,'Normalized Time Cohort (min)')
end


hold off
set(FrameProfAx,'Fontsize',16)

saveas(FrameProfFig,[outdir5, filesep,...
    parameter, '_ActivationEnergies_AllBins.png']);


