clear all
Prefixes={'2021-09-01-ywGM3_30m25Ccoll_80-100m25Cinc_Slide1',...
    '2021-09-01-ywGM5_30m25Ccoll_95-110m25Cinc_Slide1',...
    '2021-09-02-ywGM2_30m25Ccoll_110-130m25Cinc_Slide1',...
    '2021-09-03-ywGM4_30m25Ccoll_135-150m25Cinc_Slide1',...
    '2021-09-03-ywGM4_30m25Ccoll_135-150m25Cinc_Slide2',...
    '2021-08-31-yw_25C_Coll1hr_Mat2hr_ThreeStackStaging'};
outdir = 'S:/Gabriella/Dropbox/FixedTissueExperiments/';
load([outdir,filesep,'MembraneFurrowProfileInfo.mat']);
LiveAPLengths = APLengths;
LiveDVLengths = DVLengths;
NumSets = length(Prefixes);
CompiledEmbryos = {};
Ellipses ={};
liveExperiments = {};
for i = 1:NumSets
    liveExperiments{i} = LiveExperiment(Prefixes{i});
    Ellipses{i} = load([liveExperiments{i}.resultsFolder, filesep,'Ellipses.mat'], 'Ellipses');
    CompiledEmbryos{i} = load([liveExperiments{i}.resultsFolder, filesep,'CompiledEmbryos.mat'], 'CompiledEmbryos');
end

SetIDs = [];
TWindows = [];
EmbryoIDs = [];
MeanFurrows = [];
VarFurrows = [];
MeanCellWidths = [];
VarCellWidths = [];
NCs = [];
APLengths= [];
DVLengths = [];
CellCounts = [];
MedianCellSeparations = [];
for i = 1:NumSets
    NEmbryos = length(CompiledEmbryos{i}.CompiledEmbryos.nc);
    PixelSize = liveExperiments{i}.pixelSize_um;
    for j = 1:NEmbryos
        if CompiledEmbryos{i}.CompiledEmbryos.Approved(j)
            SetIDs(end+1) = i;
            if i < 5
                TWindows(end+1) = i;
            else
                TWindows(end+1) = i-1;
            end
            EmbryoIDs(end+1) = j;
            if ~isempty( CompiledEmbryos{i}.CompiledEmbryos.FurrowMeasurements{j})
                FurrowMeasurements = CompiledEmbryos{i}.CompiledEmbryos.FurrowMeasurements{j};
                FurrowDepths_pixels = sqrt((FurrowMeasurements(:,1)-FurrowMeasurements(:,3)).^2+...
                    (FurrowMeasurements(:,2)-FurrowMeasurements(:,4)).^2).';
                FurrowDepths = FurrowDepths_pixels*PixelSize;
                FurrowMean = mean(FurrowDepths);
                FurrowVar = var(FurrowDepths);
            else
                FurrowMean= NaN;
                FurrowVar = NaN;
            end
            MeanFurrows(end+1) = FurrowMean;
            VarFurrows(end+1) = FurrowVar;
            if ~isempty(CompiledEmbryos{i}.CompiledEmbryos.CellWidthMeasurements{j})
                WidthMeasurements = CompiledEmbryos{i}.CompiledEmbryos.CellWidthMeasurements{j};
                Widths_pixels = sqrt((WidthMeasurements(:,1)-WidthMeasurements(:,3)).^2+...
                    (WidthMeasurements(:,2)-WidthMeasurements(:,4)).^2).';
                Widths = Widths_pixels*PixelSize;
                WidthMean = mean(Widths);
                WidthVar = var(Widths);
            else
                WidthMean= NaN;
                WidthVar = NaN;
            end
            MeanCellWidths(end+1) = WidthMean;
            VarCellWidths(end+1) = WidthVar;
            NCs(end+1) = CompiledEmbryos{i}.CompiledEmbryos.nc(j);
            APLength = abs(CompiledEmbryos{i}.CompiledEmbryos.RotatedCoordAs(j,1)-...
                CompiledEmbryos{i}.CompiledEmbryos.RotatedCoordPs(j,1));
            APLengths(end+1) = CompiledEmbryos{i}.CompiledEmbryos.APLengths(j);
            DVLengths(end+1) = CompiledEmbryos{i}.CompiledEmbryos.DVLengths(j);
            if ~isempty(Ellipses{i}.Ellipses{j})
                CellData = Ellipses{i}.Ellipses{j};
                CellData = CellData(CellData(:,1) ~= 0,1:2);
                DVLength = abs(CompiledEmbryos{i}.CompiledEmbryos.RotatedCoordDs(j,2)-...
                    CompiledEmbryos{i}.CompiledEmbryos.RotatedCoordVs(j,2));
                APLength = abs(CompiledEmbryos{i}.CompiledEmbryos.RotatedCoordAs(j,1)-...
                    CompiledEmbryos{i}.CompiledEmbryos.RotatedCoordPs(j,1));
                yLims= [CompiledEmbryos{i}.CompiledEmbryos.RotatedCoordDs(j,2)-.1*DVLength,...
                    DVLength*1/3+CompiledEmbryos{i}.CompiledEmbryos.RotatedCoordDs(j,2)];
                xLims = [CompiledEmbryos{i}.CompiledEmbryos.RotatedCoordAs(j,1)+APLength*.25,...
                    CompiledEmbryos{i}.CompiledEmbryos.RotatedCoordAs(j,1)+APLength*.75];
                GoodCells = CellData(CellData(:,1) >=  xLims(1) &...
                    CellData(:,1) <=  xLims(2) & CellData(:,2) >=  yLims(1) & ...
                    CellData(:,2) <=  yLims(2),:);
                middley = median(GoodCells(:,2));
                ylims2 = [middley-10/PixelSize,middley+10/PixelSize];
                CellCounts(end+1) = size(GoodCells,1);
                GoodCells2 = GoodCells(GoodCells(:,2) >= ylims2(1) & GoodCells(:,2) <= ylims2(2),:); 
                [out,idx] = sort(GoodCells2(:,1));
                GoodCells2 = GoodCells2(idx,:);
                xDiff = diff(GoodCells2(:,1))*PixelSize;
                MedianCellSeparations(end+1) = median(xDiff);
            else
                CellCounts(end+1) = NaN;
                MedianCellSeparations(end+1) = NaN;
            end
        end
    end
end

DilationFactors = [mean(APLengths, 'omitnan')/mean(LiveAPLengths, 'omitnan'),...
    mean(DVLengths, 'omitnan')/mean(LiveDVLengths, 'omitnan')]; % Kind of a mess

InferredNCs = NaN(1,length(CellCounts));
InferredNCs(CellCounts >= 44) = 14;
InferredNCs(CellCounts >= 32 & CellCounts < 44) = 13;
InferredNCs(CellCounts < 32) = 12;
%%

Counts = NaN(5,3);
for t_index = 1:5
    for NC = 12:14
        Counts(t_index,NC-11) = sum(InferredNCs == NC & TWindows == t_index);
    end
end


%%
colors = brewermap(6,'Spectral');
close all
Fig3 = figure(3);

set(Fig3,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
    set(gcf,'color','w');
ax3 = axes(Fig3);

scatter(CellCounts, MeanFurrows, 50,'o','MarkerFaceColor',...
    colors(1,:),'MarkerEdgeColor','k')

grid on

hold off
% % 
% xlim([0, 50])
% ylim([0,1.05])
xlabel('Cell Count')
ylabel('Mean Furrow Depth (microns)')
% title('Gaussian Intensities for Frames 45 uW 25ºC')
ax3.FontSize = 18;
hold off
outpath = [plotdir, filesep, 'MeanFurrowDepthVsCellCount.png'];
saveas(Fig3,outpath)



