% EmbryoSizeMeasurements.m([Options])
% author: G. Martini
% date created: 5/11/22

% 
Prefix = '2022-04-07-FullSizeBrightfieldEmbryoImages';
rawDataPath = 'P:/Gabriella/LivemRNA\Data\RawDynamicsData';

rawDataFolder = 'P:/Gabriella/LivemRNA\Data\RawDynamicsData\2022-04-07\FullSizeBrightfieldEmbryoImages';
ProcPath = 'P:/Gabriella/LivemRNA\Data\ProcessedData';
DropboxFolder =  'S:/Gabriella/Dropbox\DynamicsResults';
PreProcPath =  'P:/Gabriella/LivemRNA\Data\PreProcessedData';
ExperimentType = {'IH'};
Channel1 = {'Brightfield'};
Channel2 = {char(0)};
Channel3 = {char(0)}; 
Channel4 = {char(0)}; 
Channel5 = {char(0)}; 

Channels = {Channel1, Channel2, Channel3, Channel4, Channel5};

[D, FileMode] = DetermineFileMode(rawDataFolder);
%   else
%     [D, FileMode] = DetermineFileMode(rootFolder);
%   end

%Create the output folder
OutputFolder = [PreProcPath, filesep, Prefix];
disp(['Creating folder: ', OutputFolder]);
mkdir(OutputFolder)
  
%Create the results folder
DropboxFolderName = [DropboxFolder, filesep, Prefix];
disp(['Creating folder: ', DropboxFolderName]);
mkdir(DropboxFolderName);


 %%
liveExperiment = LiveExperiment(Prefix);
disp('Exporting movie file...');

resultsFolder = liveExperiment.resultsFolder;

moviePrecision = 'uint16';
hisPrecision = 'uint16';


markandfind = true;

[XMLFolder, seriesPropertiesXML, seriesXML] = getSeriesFiles(rawDataFolder);
[LIFImages, LIFMeta] = loadLIFFile(rawDataFolder);
    if ~isempty(str2double(LIFMeta.getPixelsPhysicalSizeX(0))) &...
            ~isnan(str2double(LIFMeta.getPixelsPhysicalSizeX(0)))
        PixelSize_um = str2double(LIFMeta.getPixelsPhysicalSizeX(0));
    else
        try
            PixelSize_um = str2double(LIFMeta.getPixelsPhysicalSizeX(0).value);
        catch %no idea man
            PixelSize_um = str2double(LIFMeta.getPixelsPhysicalSizeX(1).value);
        end
    end
    MembraneZoomResultsFolder = [liveExperiment.resultsFolder, filesep, 'ZoomMembraneInfo'];
    if ~exist(MembraneZoomResultsFolder, 'dir')
        mkdir(MembraneZoomResultsFolder);
        
    end
    
    MembraneZoomPixelPath = [MembraneZoomResultsFolder, filesep, 'MembraneZoomPixelSize.mat'];
    save(MembraneZoomPixelPath,'PixelSize_um');
% 
% 
% try
%     if contains(seriesPropertiesXML(1).name, 'Mark_and_Find')
%         markandfind = true;
%     end
% catch % do nothing
% end
% 
% [LIFImages, LIFMeta] = loadLIFFile(rawDataFolder);
%  %Obtains frames information
% NEmbryos = size(LIFImages, 1);
% ySize = size(LIFImages{1,1}{1},1);
% xSize = size(LIFImages{1,1}{1},2);
% TIFStack = zeros(xSize, ySize, NEmbryos, 'double');
% for embryoIndex = 1:NEmbryos
%     TIFStack(:,:,embryoIndex) = LIFImages{embryoIndex,1}{1};
% end
% MaxPixel = max(max(max(TIFStack)));
% MinPixel = min(min(min(TIFStack)));
% TIFStack = (TIFStack-MinPixel)/(MaxPixel-MinPixel);
% TIFStack = uint8(round(255*TIFStack));
% 
% 
% saveNuclearProjection(TIFStack, [liveExperiment.preFolder, filesep, Prefix, '-Membrane.tif']);
% 
% DefineAPAxesForAllEmbryos(Prefix);
% CheckRotatedImages(Prefix);

%%
 load(MembraneZoomPixelPath,'PixelSize_um');
% 
liveExperiment = LiveExperiment(Prefix);
CompiledEmbryoPath = [liveExperiment.resultsFolder, filesep, 'CompiledEmbryos.Mat'];
if isfile(CompiledEmbryoPath)
    load(CompiledEmbryoPath,'CompiledEmbryos');
end

NEmbryos = size(CompiledEmbryos.CoordAs, 1);
CompiledEmbryos.DorsalDistances = NaN(size(CompiledEmbryos.APLengths));
CompiledEmbryos.VentralDistances = NaN(size(CompiledEmbryos.APLengths));
CompiledEmbryos.APLengths = NaN(size(CompiledEmbryos.APLengths));
CompiledEmbryos.DVLengths = NaN(size(CompiledEmbryos.APLengths));
for i = 1:NEmbryos
    if CompiledEmbryos.Approved(i)
        CompiledEmbryos.APLengths(i) = PixelSize_um*sqrt((CompiledEmbryos.CoordAs(i,1)-CompiledEmbryos.CoordPs(i,1))^2+...
            (CompiledEmbryos.CoordAs(i,2)-CompiledEmbryos.CoordPs(i,2))^2);
        CompiledEmbryos.DorsalDistances(i) = PixelSize_um*abs(CompiledEmbryos.APIntercepts(i)+CompiledEmbryos.APSlopes(i)*CompiledEmbryos.CoordDs(i,1)-CompiledEmbryos.CoordDs(i,2))/sqrt(CompiledEmbryos.APSlopes(i)^2+1);
        CompiledEmbryos.VentralDistances(i) = PixelSize_um*abs(CompiledEmbryos.APIntercepts(i)+CompiledEmbryos.APSlopes(i)*CompiledEmbryos.CoordVs(i,1)-CompiledEmbryos.CoordVs(i,2))/sqrt(CompiledEmbryos.APSlopes(i)^2+1);
        CompiledEmbryos.DVLengths(i) = CompiledEmbryos.DorsalDistances(i)+CompiledEmbryos.VentralDistances(i);
    end
end

save(CompiledEmbryoPath,'CompiledEmbryos');

KeepEmbryos = CompiledEmbryos.Approved & (CompiledEmbryos.Flags == 0);
FigFolder = 'S:/Gabriella/Dropbox/Figures/EmbryoSizeMeasurements';
mkdir(FigFolder)
APLengths = CompiledEmbryos.APLengths(KeepEmbryos);
DVLengths = CompiledEmbryos.DVLengths(KeepEmbryos);
VentralDistances = CompiledEmbryos.VentralDistances(KeepEmbryos);
DorsalDistances = CompiledEmbryos.DorsalDistances(KeepEmbryos);
NGoodEmbryos = length(APLengths);
MeanAPLength = mean(APLengths);
StdAPLength = std(APLengths);
MeanDVLength = mean(DVLengths);
StdDVLength = std(DVLengths);
scatter(APLengths, DVLengths)
AspectRatios = DVLengths./APLengths;
mkdir('S:/Gabriella/Dropbox/EmbryoSizeMeasurements');
SizeDataPath = 'S:/Gabriella/Dropbox/EmbryoSizeMeasurements/EmbryoSizeData.mat';
save(SizeDataPath,'APLengths','DVLengths','VentralDistances','DorsalDistances','NGoodEmbryos',...
    'MeanAPLength','StdAPLength','MeanDVLength','StdDVLength','AspectRatios');

%%
Prefixes = {'2022-04-11-BrightfieldMembraneFurrow-T25C-Embryo1',...
    '2022-04-11-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo1',...
    '2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo2',...
    '2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo3',...
    '2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo4',...
    '2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo5',...
    '2022-04-13-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo6',...
    '2022-04-18-BrightfieldMembraneFurrow-HisRFP-T27_5C-Embryo1',...
    '2022-04-18-BrightfieldMembraneFurrow-HisRFP-T27_5C-Embryo2',...
    '2022-04-18-BrightfieldMembraneFurrow-HisRFP-T27_5C-Embryo3',...
    '2022-04-19-BrightfieldMembraneFurrow-HisRFP-T22_5C-Embryo1',...
    '2022-04-19-BrightfieldMembraneFurrow-HisRFP-T22_5C-Embryo2',...
    '2022-04-17-BrightfieldMembraneFurrow-HisRFP-T20C-Embryo1',...
    '2022-04-17-BrightfieldMembraneFurrow-HisRFP-T20C-Embryo2',...
    '2022-04-17-BrightfieldMembraneFurrow-HisRFP-T20C-Embryo3',...
    '2022-04-15-BrightfieldMembraneFurrow-HisRFP-T17_5C-Embryo2',...
    '2022-04-15-BrightfieldMembraneFurrow-HisRFP-T17_5C-Embryo3',...
    '2022-04-16-BrightfieldMembraneFurrow-HisRFP-T17_5C-Embryo4'...
    };

T_sets = [25, 25, 25, 25, 25, 25, 27.5, 27.5, 27.5, 22.5, 22.5, 20, 20, 20, 17.5, 17.5, 17.5];
T_obs = [25.1, 25.1, 25.1, 25, 24.9, 24.9, 27.6, 27.5, 27.5, 22.5, 22.5, 20.3, 20.2, 20, 17.7, 17.8, 17.8];
Nlive = length(Prefixes);
LivePixelSizes = NaN(1, length(Prefixes));
Live_coordAs = NaN(length(Prefixes), 2);
Live_coordPs = NaN(length(Prefixes), 2);
Live_coordDs = NaN(length(Prefixes), 2);
Live_coordVs = NaN(length(Prefixes), 2);
Live_APLengths = NaN(1, length(Prefixes));
Live_DVLengths = NaN(1, length(Prefixes));
Live_VentralDistances = NaN(1, length(Prefixes));
Live_DorsalDistances = NaN(1, length(Prefixes));
Live_APSlopes = NaN(1, length(Prefixes));
Live_APIntercepts = NaN(1, length(Prefixes));
liveExperiments = cell(1, Nlive);
TestAPLengths =  NaN(1, length(Prefixes));
for embryo_index = 1:Nlive
    Pre = Prefixes{embryo_index};
    liveExperiments{embryo_index} = LiveExperiment(Pre);
    %GetAPAxisLength(Pre);
    load([liveExperiments{embryo_index}.resultsFolder,filesep,'APDetection.mat'],'coordA','coordP','coordD','coordV', 'APLength','DVLength', 'PixelSize');
    Live_coordAs(embryo_index,:) = coordA;
    Live_coordPs(embryo_index,:) = coordP;
    Live_coordDs(embryo_index,:) = coordD;
    Live_coordVs(embryo_index,:) = coordV;
    LivePixelSizes(embryo_index) = PixelSize;
    TestAPLengths(embryo_index) = APLength;
    Live_APSlopes(embryo_index) = (Live_coordAs(embryo_index,2)-Live_coordPs(embryo_index,2))/(Live_coordAs(embryo_index,1)-Live_coordPs(embryo_index,1));
    Live_APIntercepts(embryo_index) = Live_coordAs(embryo_index,2) - Live_APSlopes(embryo_index)*Live_coordAs(embryo_index,1);
    Live_APLengths(embryo_index) = LivePixelSizes(embryo_index)*sqrt((Live_coordAs(embryo_index,1)-Live_coordPs(embryo_index,1))^2+...
            (Live_coordAs(embryo_index,2)-Live_coordPs(embryo_index,2))^2);
 
    Live_DorsalDistances(embryo_index) = LivePixelSizes(embryo_index)*abs(Live_APIntercepts(embryo_index)+Live_APSlopes(embryo_index)*Live_coordDs(embryo_index,1)-Live_coordDs(embryo_index,2))/sqrt(Live_APSlopes(embryo_index)^2+1);
    Live_VentralDistances(embryo_index) = LivePixelSizes(embryo_index)*abs(Live_APIntercepts(embryo_index)+Live_APSlopes(embryo_index)*Live_coordVs(embryo_index,1)-Live_coordVs(embryo_index,2))/sqrt(Live_APSlopes(embryo_index)^2+1);
    Live_DVLengths(embryo_index) = Live_DorsalDistances(embryo_index)+Live_VentralDistances(embryo_index);
end
MeanLiveAPLength = mean(Live_APLengths);
StdLiveAPLength = std(Live_APLengths);
MeanLiveDVLength = mean(Live_DVLengths);
StdLiveDVLength = std(Live_DVLengths);
LiveAspectRatios = Live_DVLengths./Live_APLengths;
MembraneMovieSizeDataPath = 'S:/Gabriella/Dropbox/EmbryoSizeMeasurements/MembraneMovieEmbryoSizeData.mat';
save(MembraneMovieSizeDataPath,'Live_APLengths','Live_DVLengths','Live_VentralDistances','Live_DorsalDistances','Nlive',...
    'MeanLiveAPLength','StdLiveAPLength','MeanLiveDVLength','StdLiveDVLength','LiveAspectRatios');

%%
close all
APLengthsHistFig = figure(1);
set(APLengthsHistFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
APLengthsHistAx = axes(APLengthsHistFig);
histogram(APLengths,520:5:600)
legend_label = ['AP Length = ', num2str(round(MeanAPLength, 1)), ' ± ', num2str(round(StdAPLength, 1)), ' microns'];
legend(legend_label, 'FontSize', 16);
grid on

hold off

xlabel('AP Length (microns)', 'FontSize', 16)
xlim([520, 600])

ylabel('Counts', 'FontSize', 16)
ylim([0, 10])


APLengthsHistAx.YAxis.FontSize = 16;
APLengthsHistAx.XAxis.FontSize = 16;

outpath = [FigFolder,filesep, 'APLengthsHistogram.png'];
saveas(APLengthsHistFig,outpath);
        
%%
close all
DVLengthsHistFig = figure(2);
set(DVLengthsHistFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
DVLengthsHistAx = axes(DVLengthsHistFig);
histogram(DVLengths,160:2.5:190)
legend_label = ['DV Length = ', num2str(round(MeanDVLength, 1)), ' ± ', num2str(round(StdDVLength, 1)), ' microns'];
legend(legend_label, 'FontSize', 16);
grid on

hold off

xlabel('DV Length (microns)', 'FontSize', 16)
xlim([160, 190])

ylabel('Counts', 'FontSize', 16)
ylim([0, 20])


DVLengthsHistAx.YAxis.FontSize = 16;
DVLengthsHistAx.XAxis.FontSize = 16;

outpath = [FigFolder,filesep, 'DVLengthsHistogram.png'];
saveas(DVLengthsHistFig,outpath);
   
%%
close all
AspectRatiosHistFig = figure(3);
set(AspectRatiosHistFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
AspectRatiosHistAx = axes(AspectRatiosHistFig);
histogram(AspectRatios)%,160:2.5:190)
legend_label = ['DV Length/AP Length = ', num2str(round(mean(AspectRatios), 2)), ' ± ', num2str(round(std(AspectRatios), 2))];
legend(legend_label, 'FontSize', 16);
grid on

hold off

xlabel('DV Length/AP Length', 'FontSize', 16)
%xlim([160, 190])

ylabel('Counts', 'FontSize', 16)
%ylim([0, 20])


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
xlim([520, 600])

ylabel('DV Length (microns)', 'FontSize', 16)
ylim([160, 190])


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
%xlim([520, 600])

ylabel('Ventral Distance from Midline (microns)', 'FontSize', 16)
%ylim([160, 190])


DorsalvsVentralscatterAx.YAxis.FontSize = 16;
DorsalvsVentralscatterAx.XAxis.FontSize = 16;

outpath = [FigFolder,filesep, 'DorsalvsVentral_Scatter.png'];
saveas(DorsalvsVentralscatterFig,outpath);
%%
close all
DorsalLengthsHistFig = figure(6);
set(DorsalLengthsHistFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
DorsalLengthsHistAx = axes(DorsalLengthsHistFig);
histogram(DorsalDistances,50:2.5:80)
legend_label = ['Dorsal Distance = ', num2str(round(mean(DorsalDistances), 1)), ' ± ', num2str(round(std(DorsalDistances), 1)), ' microns'];
legend(legend_label, 'FontSize', 16);
grid on

hold off

xlabel('Dorsal Distance from Midline (microns)', 'FontSize', 16)
xlim([50, 85])

ylabel('Counts', 'FontSize', 16)
ylim([0, 12])


DorsalLengthsHistAx.YAxis.FontSize = 16;
DorsalLengthsHistAx.XAxis.FontSize = 16;

outpath = [FigFolder,filesep, 'DorsalDistancesHistogram.png'];
saveas(DorsalLengthsHistFig,outpath);

%%
close all
VentralLengthsHistFig = figure(7);
set(VentralLengthsHistFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
VentralLengthsHistAx = axes(VentralLengthsHistFig);
histogram(VentralDistances,85:2.5:125)
legend_label = ['Ventral Distance = ', num2str(round(mean(VentralDistances), 1)), ' ± ', num2str(round(std(VentralDistances), 1)), ' microns'];
legend(legend_label, 'FontSize', 16);
grid on

hold off

xlabel('Ventral Distance from Midline (microns)', 'FontSize', 16)
xlim([80, 130])

ylabel('Counts', 'FontSize', 16)
ylim([0, 12])


VentralLengthsHistAx.YAxis.FontSize = 16;
VentralLengthsHistAx.XAxis.FontSize = 16;

outpath = [FigFolder,filesep, 'VentralDistancesHistogram.png'];
saveas(VentralLengthsHistFig,outpath);

%%
close all
LiveAPLengthsHistFig = figure(1);
set(LiveAPLengthsHistFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
APLengthsHistAx = axes(LiveAPLengthsHistFig);
histogram(Live_APLengths,500:5:600)
legend_label = ['AP Length = ', num2str(round(MeanLiveAPLength, 1)), ' ± ', num2str(round(StdLiveAPLength, 1)), ' microns'];
legend(legend_label, 'FontSize', 16);
grid on

hold off

xlabel('AP Length (microns)', 'FontSize', 16)
xlim([500, 600])

ylabel('Counts', 'FontSize', 16)
ylim([0, 10])


APLengthsHistAx.YAxis.FontSize = 16;
APLengthsHistAx.XAxis.FontSize = 16;

outpath = [FigFolder,filesep, 'MembraneMovieAPLengthsHistogram.png'];
saveas(LiveAPLengthsHistFig,outpath);
        
%%
close all
LiveDVLengthsHistFig = figure(2);
set(LiveDVLengthsHistFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
DVLengthsHistAx = axes(LiveDVLengthsHistFig);
histogram(Live_DVLengths,170:2.5:190)
legend_label = ['DV Length = ', num2str(round(MeanLiveDVLength, 1)), ' ± ', num2str(round(StdLiveDVLength, 1)), ' microns'];
legend(legend_label, 'FontSize', 16);
grid on

hold off

xlabel('DV Length (microns)', 'FontSize', 16)
xlim([170, 190])

ylabel('Counts', 'FontSize', 16)
ylim([0, 8])


DVLengthsHistAx.YAxis.FontSize = 16;
DVLengthsHistAx.XAxis.FontSize = 16;

outpath = [FigFolder,filesep, 'MembraneMovieDVLengthsHistogram.png'];
saveas(LiveDVLengthsHistFig,outpath);
   
%%
close all
LiveAspectRatiosHistFig = figure(3);
set(LiveAspectRatiosHistFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
AspectRatiosHistAx = axes(LiveAspectRatiosHistFig);
histogram(LiveAspectRatios)%,160:2.5:190)
legend_label = ['DV Length/AP Length = ', num2str(round(mean(LiveAspectRatios), 2)), ' ± ', num2str(round(std(LiveAspectRatios), 2))];
legend(legend_label, 'FontSize', 16);
grid on

hold off

xlabel('DV Length/AP Length', 'FontSize', 16)
%xlim([160, 190])

ylabel('Counts', 'FontSize', 16)
ylim([0, 10])


AspectRatiosHistAx.YAxis.FontSize = 16;
AspectRatiosHistAx.XAxis.FontSize = 16;

outpath = [FigFolder,filesep, 'MembraneMovieAspectRatiosHistogram.png'];
saveas(LiveAspectRatiosHistFig,outpath);

%%
close all
LiveAPvsDVscatterFig = figure(4);
set(LiveAPvsDVscatterFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
APvsDVscatterAx = axes(LiveAPvsDVscatterFig);
scatter(Live_APLengths, Live_DVLengths, 50, 'r', 'filled')

grid on

hold off

xlabel('AP Length (microns)', 'FontSize', 16)
xlim([510, 570])

ylabel('DV Length (microns)', 'FontSize', 16)
ylim([160, 190])


APvsDVscatterAx.YAxis.FontSize = 16;
APvsDVscatterAx.XAxis.FontSize = 16;

outpath = [FigFolder,filesep, 'MembraneMovieAPvsDV_Scatter.png'];
saveas(LiveAPvsDVscatterFig,outpath);

%%
close all
LiveDorsalvsVentralscatterFig = figure(5);
set(LiveDorsalvsVentralscatterFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
DorsalvsVentralscatterAx = axes(LiveDorsalvsVentralscatterFig);
scatter(Live_DorsalDistances, Live_VentralDistances, 50, 'r', 'filled')

grid on

hold off

xlabel('Dorsal Distance from Midline (microns)', 'FontSize', 16)
xlim([55, 90])

ylabel('Ventral Distance from Midline (microns)', 'FontSize', 16)
%ylim([160, 190])


DorsalvsVentralscatterAx.YAxis.FontSize = 16;
DorsalvsVentralscatterAx.XAxis.FontSize = 16;

outpath = [FigFolder,filesep, 'MembraneMovieDorsalvsVentral_Scatter.png'];
saveas(LiveDorsalvsVentralscatterFig,outpath);
%%
close all
LiveDorsalLengthsHistFig = figure(6);
set(LiveDorsalLengthsHistFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
DorsalLengthsHistAx = axes(LiveDorsalLengthsHistFig);
histogram(Live_DorsalDistances,50:5:90)
legend_label = ['Dorsal Distance = ', num2str(round(mean(Live_DorsalDistances), 1)), ' ± ', num2str(round(std(Live_DorsalDistances), 1)), ' microns'];
legend(legend_label, 'FontSize', 16);
grid on

hold off

xlabel('Dorsal Distance from Midline (microns)', 'FontSize', 16)
xlim([55, 90])

ylabel('Counts', 'FontSize', 16)
ylim([0, 8])


DorsalLengthsHistAx.YAxis.FontSize = 16;
DorsalLengthsHistAx.XAxis.FontSize = 16;

outpath = [FigFolder,filesep, 'MembraneMovieDorsalDistancesHistogram.png'];
saveas(LiveDorsalLengthsHistFig,outpath);

%%
close all
LiveVentralLengthsHistFig = figure(7);
set(LiveVentralLengthsHistFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
set(gcf,'color','w');
VentralLengthsHistAx = axes(LiveVentralLengthsHistFig);
histogram(Live_VentralDistances,85:5:125)
legend_label = ['Ventral Distance = ', num2str(round(mean(Live_VentralDistances), 1)), ' ± ', num2str(round(std(Live_VentralDistances), 1)), ' microns'];
legend(legend_label, 'FontSize', 16);
grid on

hold off

xlabel('Ventral Distance from Midline (microns)', 'FontSize', 16)
xlim([85, 125])

ylabel('Counts', 'FontSize', 16)
ylim([0, 8])


VentralLengthsHistAx.YAxis.FontSize = 16;
VentralLengthsHistAx.XAxis.FontSize = 16;

outpath = [FigFolder,filesep, 'MembraneMovieVentralDistancesHistogram.png'];
saveas(LiveVentralLengthsHistFig,outpath);
   
   

   
   
