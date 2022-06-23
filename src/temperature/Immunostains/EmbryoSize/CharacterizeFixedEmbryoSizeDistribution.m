function CharacterizeFixedEmbryoSizeDistribution(Prefixes, SetLabel)
 liveExperiments = cell(1, length(Prefixes));
 for i = 1:length(Prefixes)
     liveExperiments{i} = LiveExperiment(Prefixes{i});
 end
 FixedPixelSize_um = liveExperiments{1}.pixelSize_um;
 CT_CompiledEmbryos = CombineCompiledEmbryos(Prefixes);
 ApprovedEmbryos = CT_CompiledEmbryos.Approved;
 NEmbryos = length(ApprovedEmbryos);
 

FigFolder = ['S:/Gabriella/Dropbox/Figures/FixedEmbryoSizeMeasurements', filesep, SetLabel];
if ~isdir(FigFolder)
mkdir(FigFolder)
end

APLengths = NaN(1, NEmbryos);
APSlopes =  NaN(1, NEmbryos);
APIntercepts =  NaN(1, NEmbryos);
DorsalDistances =  NaN(1, NEmbryos);
VentralDistances =  NaN(1, NEmbryos);
DVLengths =  NaN(1, NEmbryos);

for i = 1:NEmbryos
    if CT_CompiledEmbryos.Approved(i) 
        APLengths(i) = FixedPixelSize_um*sqrt((CT_CompiledEmbryos.CoordAs(i,1)-CT_CompiledEmbryos.CoordPs(i,1))^2+...
            (CT_CompiledEmbryos.CoordAs(i,2)-CT_CompiledEmbryos.CoordPs(i,2))^2);
        APSlopes(i) = (CT_CompiledEmbryos.CoordAs(i,2)-CT_CompiledEmbryos.CoordPs(i,2))/(CT_CompiledEmbryos.CoordAs(i,1)-CT_CompiledEmbryos.CoordPs(i,1));
        APIntercepts(i) = CT_CompiledEmbryos.CoordAs(i,2)-APSlopes(i)*CT_CompiledEmbryos.CoordAs(i,1);
        DorsalDistances(i) = FixedPixelSize_um*abs(APIntercepts(i)+APSlopes(i)*CT_CompiledEmbryos.CoordDs(i,1)-CT_CompiledEmbryos.CoordDs(i,2))/sqrt(APSlopes(i)^2+1);
        VentralDistances(i) = FixedPixelSize_um*abs(APIntercepts(i)+APSlopes(i)*CT_CompiledEmbryos.CoordVs(i,1)-CT_CompiledEmbryos.CoordVs(i,2))/sqrt(APSlopes(i)^2+1);
        DVLengths(i) = DorsalDistances(i)+VentralDistances(i);
    end
end


AspectRatios = DVLengths./APLengths;
KeepEmbryos = CT_CompiledEmbryos.Approved & (CT_CompiledEmbryos.Flags == 0) & ~isnan(DVLengths);
APLengths = APLengths(KeepEmbryos);
APSlopes = APSlopes(KeepEmbryos);
APIntercepts = APIntercepts(KeepEmbryos);
DorsalDistances = DorsalDistances(KeepEmbryos);
VentralDistances = VentralDistances(KeepEmbryos);
DVLengths = DVLengths(KeepEmbryos);
AspectRatios = AspectRatios(KeepEmbryos);


MeanAspectRatio = mean(AspectRatios, 'omitnan');
StdAspectRatio = std(AspectRatios, 'omitnan');

MeanAPLength = mean(APLengths, 'omitnan');
StdAPLength = std(APLengths, 'omitnan');

MeanDVLength = mean(DVLengths, 'omitnan');
StdDVLength = std(DVLengths, 'omitnan');



%%
close all
APLengthMin = floor(min(APLengths)/5)*5;
APLengthMax = ceil(max(APLengths)/5)*5;
APLengthsHistFig = figure(1);
set(APLengthsHistFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
APLengthsHistAx = axes(APLengthsHistFig);
h = histogram(APLengths,APLengthMin:5:APLengthMax);
legend_label = ['AP Length = ', num2str(round(MeanAPLength, 1)), ' ± ', num2str(round(StdAPLength, 1)), ' microns'];
legend(legend_label, 'FontSize', 16);
grid on

hold off

xlabel('AP Length (microns)', 'FontSize', 16)
xlim([APLengthMin, APLengthMax+5])

ylabel('Counts', 'FontSize', 16)
ylim([0, ceil(max(h.Values)/5)*5])


APLengthsHistAx.YAxis.FontSize = 16;
APLengthsHistAx.XAxis.FontSize = 16;

outpath = [FigFolder,filesep, 'APLengthsHistogram.png'];
saveas(APLengthsHistFig,outpath);
        
%%
DVLengthMin = floor(min(DVLengths)/5)*5;
DVLengthMax = ceil(max(DVLengths)/5)*5;
close all
DVLengthsHistFig = figure(2);
set(DVLengthsHistFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
DVLengthsHistAx = axes(DVLengthsHistFig);
h = histogram(DVLengths,DVLengthMin:2.5:DVLengthMax);
legend_label = ['DV Length = ', num2str(round(MeanDVLength, 1)), ' ± ', num2str(round(StdDVLength, 1)), ' microns'];
legend(legend_label, 'FontSize', 16);
grid on

hold off

xlabel('DV Length (microns)', 'FontSize', 16)
xlim([DVLengthMin, DVLengthMax+5])

ylabel('Counts', 'FontSize', 16)
ylim([0, ceil(max(h.Values)/5)*5])



DVLengthsHistAx.YAxis.FontSize = 16;
DVLengthsHistAx.XAxis.FontSize = 16;

outpath = [FigFolder,filesep, 'DVLengthsHistogram.png'];
saveas(DVLengthsHistFig,outpath);
   
%%
ARBinSize = 0.005;
ARLengthMin = floor(min(AspectRatios)/ARBinSize)*ARBinSize;
ARLengthMax = ceil(max(AspectRatios)/ARBinSize)*ARBinSize;
close all
AspectRatiosHistFig = figure(3);
set(AspectRatiosHistFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
AspectRatiosHistAx = axes(AspectRatiosHistFig);
h = histogram(AspectRatios, ARLengthMin:ARBinSize:ARLengthMax);%,160:2.5:190)
legend_label = ['DV Length/AP Length = ', num2str(round(mean(AspectRatios), 2)), ' ± ', num2str(round(std(AspectRatios), 2))];
legend(legend_label, 'FontSize', 16);
grid on

hold off

xlabel('DV Length/AP Length', 'FontSize', 16)
xlim([ARLengthMin, ARLengthMax+ARBinSize])

ylabel('Counts', 'FontSize', 16)
ylim([0, ceil(max(h.Values)/5)*5+5])


AspectRatiosHistAx.YAxis.FontSize = 16;
AspectRatiosHistAx.XAxis.FontSize = 16;

outpath = [FigFolder,filesep, 'AspectRatiosHistogram.png'];
saveas(AspectRatiosHistFig,outpath);

%%


close all
APvsDVscatterFig = figure(4);
set(APvsDVscatterFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
APvsDVscatterAx = axes(APvsDVscatterFig);
scatter(APLengths, DVLengths, 50, 'r', 'filled')

grid on

hold off

xlabel('AP Length (microns)', 'FontSize', 16)
xlim([floor(min(APLengths)/5)*5, ceil(max(APLengths)/5)*5])

ylabel('DV Length (microns)', 'FontSize', 16)
ylim([floor(min(DVLengths)/5)*5, ceil(max(DVLengths)/5)*5])


APvsDVscatterAx.YAxis.FontSize = 16;
APvsDVscatterAx.XAxis.FontSize = 16;

outpath = [FigFolder,filesep, 'APvsDV_Scatter.png'];
saveas(APvsDVscatterFig,outpath);

%%
close all
DorsalvsVentralscatterFig = figure(5);
set(DorsalvsVentralscatterFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
DorsalvsVentralscatterAx = axes(DorsalvsVentralscatterFig);
scatter(DorsalDistances, VentralDistances, 50, 'r', 'filled')

grid on

hold off

xlabel('Dorsal Distance from Midline (microns)', 'FontSize', 16)
xlim([floor(min(DorsalDistances)/5)*5, ceil(max(DorsalDistances)/5)*5])

ylabel('Ventral Distance from Midline (microns)', 'FontSize', 16)
ylim([floor(min(VentralDistances)/5)*5, ceil(max(VentralDistances)/5)*5])


DorsalvsVentralscatterAx.YAxis.FontSize = 16;
DorsalvsVentralscatterAx.XAxis.FontSize = 16;

outpath = [FigFolder,filesep, 'DorsalvsVentral_Scatter.png'];
saveas(DorsalvsVentralscatterFig,outpath);
%%
DorsalBinSize = 2.5;
DorsalLengthMin = floor(min(DorsalDistances)/DorsalBinSize)*DorsalBinSize;
DorsalLengthMax = ceil(max(DorsalDistances)/DorsalBinSize)*DorsalBinSize;
close all
DorsalLengthsHistFig = figure(6);
set(DorsalLengthsHistFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
DorsalLengthsHistAx = axes(DorsalLengthsHistFig);
h = histogram(DorsalDistances,DorsalLengthMin:DorsalBinSize:DorsalLengthMax);
legend_label = ['Dorsal Distance = ', num2str(round(mean(DorsalDistances), 1)), ' ± ', num2str(round(std(DorsalDistances), 1)), ' microns'];
legend(legend_label, 'FontSize', 16);
grid on

hold off

xlabel('Dorsal Distance from Midline (microns)', 'FontSize', 16)
xlim([DorsalLengthMin, DorsalLengthMax+DorsalBinSize])

ylabel('Counts', 'FontSize', 16)
ylim([0, ceil(max(h.Values)/5)*5+5])


DorsalLengthsHistAx.YAxis.FontSize = 16;
DorsalLengthsHistAx.XAxis.FontSize = 16;

outpath = [FigFolder,filesep, 'DorsalDistancesHistogram.png'];
saveas(DorsalLengthsHistFig,outpath);

%%
DorsalBinSize = 2.5;
DorsalLengthMin = floor(min(VentralDistances)/DorsalBinSize)*DorsalBinSize;
DorsalLengthMax = ceil(max(VentralDistances)/DorsalBinSize)*DorsalBinSize;
close all
VentralLengthsHistFig = figure(7);
set(VentralLengthsHistFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
VentralLengthsHistAx = axes(VentralLengthsHistFig);
h = histogram(VentralDistances,DorsalLengthMin:DorsalBinSize:DorsalLengthMax);
legend_label = ['Ventral Distance = ', num2str(round(mean(VentralDistances), 1)), ' ± ', num2str(round(std(VentralDistances), 1)), ' microns'];
legend(legend_label, 'FontSize', 16);
grid on

hold off

xlabel('Ventral Distance from Midline (microns)', 'FontSize', 16)
xlim([DorsalLengthMin, DorsalLengthMax+DorsalBinSize])

ylabel('Counts', 'FontSize', 16)
ylim([0, ceil(max(h.Values)/5)*5+5])


VentralLengthsHistAx.YAxis.FontSize = 16;
VentralLengthsHistAx.XAxis.FontSize = 16;

outpath = [FigFolder,filesep, 'VentralDistancesHistogram.png'];
saveas(VentralLengthsHistFig,outpath);

