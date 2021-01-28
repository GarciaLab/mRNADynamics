function HealthSummary = SummarizeEmbryoHealth(Prefix, includePlots)
%% Load necessary info into memory
liveExperiment = LiveExperiment(Prefix);
FrameInfo = getFrameInfo(liveExperiment);
FrameTimes = [FrameInfo(:).Time]; % in seconds
nc_info = [liveExperiment.nc9, liveExperiment.nc10, liveExperiment.nc11,...
    liveExperiment.nc12, liveExperiment.nc13, liveExperiment.nc14, length(FrameInfo)];
PixelSize = liveExperiment.pixelSize_um;
nucleusDiameters = zeros(1, 6);
for nc=9:14
    nucleusDiameters(nc-8) = getDefaultParameters(FrameInfo,[['d', num2str(nc)]])/PixelSize; % in pixels
end
AddNuclearPosition(Prefix);
DetermineNucleiEndFrames(Prefix);
schnitzcells = CalculateNuclearMovement(Prefix);
schnitzcycles = [schnitzcells(:).cycle];
%% Initialize HealthSummaryInfo
HealthSummary = {};
HealthSummary.SchnitzCount = NaN(1, 5);
HealthSummary.FractionSickNuclei = NaN(1, 5);
HealthSummary.FractionRejectedNuclei = NaN(1,5);
HealthSummary.FractionCompleteNuclei = NaN(1,5);
HealthSummary.FractionFirstLastFrameNuclei = NaN(1,5);

HealthSummary.MeanTotalXDistanceTraveled = NaN(1,5);
HealthSummary.MeanTotalYDistanceTraveled = NaN(1,5);
HealthSummary.MeanTotalDistanceTraveled = NaN(1,5);
HealthSummary.MeanDistanceTraveledPerSecond = NaN(1,5);
HealthSummary.MeanTotalXDisplacement = NaN(1,5);
HealthSummary.MeanTotalYDisplacement = NaN(1,5);
HealthSummary.MeanTotalDisplacement = NaN(1,5);
HealthSummary.MeanDisplacementPerSecond = NaN(1,5);
HealthSummary.StdTotalXDistanceTraveled = NaN(1,5);
HealthSummary.StdTotalYDistanceTraveled = NaN(1,5);
HealthSummary.StdTotalDistanceTraveled = NaN(1,5);
HealthSummary.StdDistanceTraveledPerSecond = NaN(1,5);
HealthSummary.StdTotalXDisplacement = NaN(1,5);
HealthSummary.StdTotalYDisplacement = NaN(1,5);
HealthSummary.StdTotalDisplacement = NaN(1,5);
HealthSummary.StdDisplacementPerSecond = NaN(1,5);
[NCDivisionInfo,DivisionStdInfo] = CalculateSchnitzDivisionCycleTimes(Prefix);
HealthSummary.NCDivisionInfo = NCDivisionInfo;
HealthSummary.DivisionStdInfo = DivisionStdInfo;
HealthSummary.APboundaries = NaN(1,2);


%%
for NC = max(min(schnitzcycles), 10):max(schnitzcycles)
    NCSchnitzIndices = find(schnitzcycles == NC);
    NCschnitzcells = schnitzcells(NCSchnitzIndices);
    numSchnitzInCycle = length(NCSchnitzIndices);
    HealthSummary.SchnitzCount(NC-9) = numSchnitzInCycle;
    HealthSummary.FractionRejectedNuclei(NC-9) = sum([NCschnitzcells(:).Approved] == 0)/numSchnitzInCycle;
    HealthSummary.FractionSickNuclei(NC-9) = sum([NCschnitzcells(:).Flag] == 6)/numSchnitzInCycle;
    FractionCompleteNuclei = zeros(1, numSchnitzInCycle);
    FractionFirstLastFrameNuclei = zeros(1, numSchnitzInCycle);
    MeanTotalXDistanceTraveled = zeros(1, numSchnitzInCycle);
    MeanTotalYDistanceTraveled = zeros(1, numSchnitzInCycle);
    MeanTotalDistanceTraveled = zeros(1, numSchnitzInCycle);
    MeanDistanceTraveledPerSecond = zeros(1, numSchnitzInCycle);
    MeanTotalXDisplacement = zeros(1, numSchnitzInCycle);
    MeanTotalYDisplacement = zeros(1, numSchnitzInCycle);
    MeanTotalDisplacement = zeros(1, numSchnitzInCycle);
    MeanDisplacementPerSecond = zeros(1, numSchnitzInCycle);
    for schnitz_index=1:numSchnitzInCycle
        FractionCompleteNuclei(schnitz_index) = NCschnitzcells(schnitz_index).VelocityInfo.SchnitzHasAllFrames;
        FractionFirstLastFrameNuclei(schnitz_index) = NCschnitzcells(schnitz_index).VelocityInfo.SchnitzHasFirstAndLastCycleFrames;
        MeanTotalXDistanceTraveled(schnitz_index) = NCschnitzcells(schnitz_index).VelocityInfo.TotalXDistanceTraveled;
        MeanTotalYDistanceTraveled(schnitz_index) = NCschnitzcells(schnitz_index).VelocityInfo.TotalYDistanceTraveled;
        MeanTotalDistanceTraveled(schnitz_index) = NCschnitzcells(schnitz_index).VelocityInfo.TotalDistanceTraveled;
        MeanDistanceTraveledPerSecond(schnitz_index) = NCschnitzcells(schnitz_index).VelocityInfo.MeanDistanceTraveledPerSecond;
        MeanTotalXDisplacement(schnitz_index) = NCschnitzcells(schnitz_index).VelocityInfo.TotalXDisplacement;
        MeanTotalYDisplacement(schnitz_index) = NCschnitzcells(schnitz_index).VelocityInfo.TotalYDisplacement;
        MeanTotalDisplacement(schnitz_index) = NCschnitzcells(schnitz_index).VelocityInfo.TotalDisplacement;
        MeanDisplacementPerSecond(schnitz_index) = NCschnitzcells(schnitz_index).VelocityInfo.MeanDisplacementPerSecond;
    end
    HealthSummary.FractionCompleteNuclei(NC-9) = sum(FractionCompleteNuclei)/numSchnitzInCycle;
    HealthSummary.FractionFirstLastFrameNuclei(NC-9) = sum(FractionFirstLastFrameNuclei)/numSchnitzInCycle;
    HealthSummary.MeanTotalXDistanceTraveled(NC-9) = nanmean(MeanTotalXDistanceTraveled);
    HealthSummary.MeanTotalYDistanceTraveled(NC-9) = nanmean(MeanTotalYDistanceTraveled);
    HealthSummary.MeanTotalDistanceTraveled(NC-9) = nanmean(MeanTotalDistanceTraveled);
    HealthSummary.MeanDistanceTraveledPerSecond(NC-9) = nanmean(MeanDistanceTraveledPerSecond);
    HealthSummary.MeanTotalXDisplacement(NC-9) = nanmean(MeanTotalXDisplacement);
    HealthSummary.MeanTotalYDisplacement(NC-9) = nanmean(MeanTotalYDisplacement);
    HealthSummary.MeanTotalDisplacement(NC-9) = nanmean(MeanTotalDisplacement);
    HealthSummary.MeanDisplacementPerSecond(NC-9) = nanmean(MeanDisplacementPerSecond);
    HealthSummary.StdTotalXDistanceTraveled(NC-9) = nanstd(MeanTotalXDistanceTraveled)/length(MeanTotalXDistanceTraveled(~isnan(MeanTotalXDistanceTraveled)));;
    HealthSummary.StdTotalYDistanceTraveled(NC-9) = nanstd(MeanTotalYDistanceTraveled)/length(MeanTotalYDistanceTraveled(~isnan(MeanTotalYDistanceTraveled)));;
    HealthSummary.StdTotalDistanceTraveled(NC-9) = nanstd(MeanTotalDistanceTraveled)/length(MeanTotalDistanceTraveled(~isnan(MeanTotalDistanceTraveled)));;
    HealthSummary.StdDistanceTraveledPerSecond(NC-9) = nanstd(MeanDistanceTraveledPerSecond)/length(MeanDistanceTraveledPerSecond(~isnan(MeanDistanceTraveledPerSecond)));;
    HealthSummary.StdTotalXDisplacement(NC-9) = nanstd(MeanTotalXDisplacement)/length(MeanTotalXDisplacement(~isnan(MeanTotalXDisplacement)));
    HealthSummary.StdTotalYDisplacement(NC-9) = nanstd(MeanTotalYDisplacement)/length(MeanTotalYDisplacement(~isnan(MeanTotalYDisplacement)));
    HealthSummary.StdTotalDisplacement(NC-9) = nanstd(MeanTotalDisplacement)/length(MeanTotalDisplacement(~isnan(MeanTotalDisplacement)));
    HealthSummary.StdDisplacementPerSecond(NC-9) = nanstd(MeanDisplacementPerSecond)/length(MeanDisplacementPerSecond(~isnan(MeanDisplacementPerSecond)));
    
end

for schnitz_index=1:length(schnitzcells)
    if ~isempty(schnitzcells(schnitz_index).APpos)
        meanAPpos = nanmean(schnitzcells(schnitz_index).APpos);
        HealthSummary.APboundaries(1) = min(HealthSummary.APboundaries(1), meanAPpos);
        HealthSummary.APboundaries(2) = max(HealthSummary.APboundaries(2), meanAPpos);
    end
end
%%
HealthSummaryPath = [liveExperiment.resultsFolder, 'HealthSummary.mat'];
save(HealthSummaryPath, 'HealthSummary')



%% Plot Info
if includePlots
close all
subplot(2, 3, 1)
text(-0.55,0.9,['Embryo: '], 'FontWeight', 'bold'); axis off
text(-0.55,0.8,[Prefix]); axis off
text(-0.55,0.6,['AP range: ', num2str(round(HealthSummary.APboundaries(1), 2)), ' - ', num2str(round(HealthSummary.APboundaries(2), 2))]); 
subplot(2, 3, [2, 3])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.8, 0.8]);
errorbar(10:13, HealthSummary.NCDivisionInfo, HealthSummary.DivisionStdErrorInfo,'o',  'Color', 'r', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'r')
xlim([9, 14])
xticks(9:14)
try
ylim([0, max(HealthSummary.NCDivisionInfo+HealthSummary.DivisionStdErrorInfo)*1.1])
end
xlabel('division cycle')
ylabel('mean cycle duration')
% BarLabels = categorical({'Fraction of Sick Nuclei NC10','Fraction of Sick Nuclei NC11',...
%     'Fraction of Sick Nuclei NC12', 'Fraction of Sick Nuclei NC13',...
%     'Fraction of Sick Nuclei NC14','Fraction of Rejected Nuclei NC10',...
%     'Fraction of Rejected Nuclei NC11', 'Fraction of Rejected Nuclei NC12',...
%     'Fraction of Rejected Nuclei NC13', 'Fraction of Rejected Nuclei NC14',...
%     'Fraction of Complete Nuclei NC10',...
%     'Fraction of Complete Nuclei NC11', 'Fraction of Complete Nuclei NC12',...
%     'Fraction of Complete Nuclei NC13', 'Fraction of Complete Nuclei NC14',...
%     'Fraction of First/Last NC10',...
%     'Fraction of First/Last Nuclei NC11', 'Fraction of First/Last Nuclei NC12',...
%     'Fraction of First/Last Nuclei NC13', 'Fraction of First/Last Nuclei NC14'});
BarLabels = categorical({'Fraction of Sick Nuclei','Fraction of Rejected Nuclei',...
    'Fraction of Complete Nuclei',...
    'Fraction of First/Last'});

BarLabels = reordercats(BarLabels,{'Fraction of Sick Nuclei','Fraction of Rejected Nuclei',...
    'Fraction of Complete Nuclei',...
    'Fraction of First/Last'}); %specified order
FractionVector = [HealthSummary.FractionSickNuclei; HealthSummary.FractionRejectedNuclei;...
    HealthSummary.FractionCompleteNuclei; HealthSummary.FractionFirstLastFrameNuclei];
subplot(2,3, 4)
b = bar(BarLabels, FractionVector);
hold on 
ylabel('Fraction of Schnitz Cells')
ylim([0,1])
colors = get(gca,'ColorOrder');
nColors=5;                           % make variable so can change easily
labels={'Cycle 10';'Cycle 11';'Cycle 12';'Cycle 13';'Cycle 14'};
hBLG = bar(nan(2,nColors));         % the bar object array for legend
for i=1:nColors
    hBLG(i).FaceColor=colors(i,:);
end
%hLG=legend(hBLG,labels,'location','northeast');
hold off
subplot(2,3,5)
BarLabels2 = categorical({'Total X Distance Traveled','Total Y Distance Traveled',...
    'Total Distance Traveled',...
    'Final X Displacement', 'Final Y Displacement', 'Final Displacement'});

BarLabels2 = reordercats(BarLabels2,{'Total X Distance Traveled','Total Y Distance Traveled',...
    'Total Distance Traveled',...
    'Final X Displacement', 'Final Y Displacement', 'Final Displacement'}); %specified order
FractionVector2 = [HealthSummary.MeanTotalXDistanceTraveled; HealthSummary.MeanTotalYDistanceTraveled;...
    HealthSummary.MeanTotalDistanceTraveled; HealthSummary.MeanTotalXDisplacement;...
    HealthSummary.MeanTotalYDisplacement; HealthSummary.MeanTotalDisplacement];
StdVector2 = [HealthSummary.StdTotalXDistanceTraveled; HealthSummary.MeanTotalYDistanceTraveled;...
    HealthSummary.MeanTotalDistanceTraveled; HealthSummary.MeanTotalXDisplacement;...
    HealthSummary.MeanTotalYDisplacement; HealthSummary.MeanTotalDisplacement];
b2 = bar(BarLabels2, FractionVector2);
hold on
ylabel('microns')
nColors=5;                           % make variable so can change easily
labels={'Cycle 10';'Cycle 11';'Cycle 12';'Cycle 13';'Cycle 14'};
hBLG = bar(nan(2,nColors));         % the bar object array for legend
for i=1:nColors
    hBLG(i).FaceColor=colors(i,:);
end
%hLG=legend(hBLG,labels,'location','northeast');
hold off

subplot(2,3,6)
BarLabels3 = categorical({'Mean Distance Traveled per Second','Total Displacement per second'});

BarLabels3 = reordercats(BarLabels3,{'Mean Distance Traveled per Second','Total Displacement per second'}); %specified order
FractionVector3 = [HealthSummary.MeanDistanceTraveledPerSecond; HealthSummary.MeanDisplacementPerSecond];
b3 = bar(BarLabels3, FractionVector3);
hold on
ylabel('microns/s')
nColors=5;                           % make variable so can change easily
labels={'Cycle 10';'Cycle 11';'Cycle 12';'Cycle 13';'Cycle 14'};
hBLG = bar(nan(2,nColors));         % the bar object array for legend
for i=1:nColors
    hBLG(i).FaceColor=colors(i,:);
end
hLG=legend(hBLG,labels,'location','northeast');
hold off

HealthSummaryFigPath= [liveExperiment.resultsFolder, 'HealthSummary.png'];
saveas(gcf,HealthSummaryFigPath);
end
