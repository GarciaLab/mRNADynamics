% CalculateBackgroundAverages.m
% author: Gabriella Martini
% date added: 1/27/22
function [BackgroundZAverages, RawBackgroundZAverages, TotalBackgroundAverage] = ...
    CalculateBackgroundAverages(Prefix, Ellipses)
liveExperiment = LiveExperiment(Prefix);
FrameInfo = getFrameInfo(liveExperiment);
nuclear_diameter = getDefaultParameters(FrameInfo,'d14');
pixelSize = liveExperiment.pixelSize_um;
d_pixels = nuclear_diameter/pixelSize;
if ~exist('Ellipses', 'var')
outpath = [liveExperiment.resultsFolder, filesep, 'BackgroundEstimateEllipses.mat'];
load(outpath);
end
NumSlices = length(Ellipses);
RawBackgroundZAverages = NaN(1, NumSlices);
BackgroundZAverages = NaN(1, NumSlices);
AllBackgroundMeasures = [];
for j = 1:NumSlices
    if size(Ellipses{j},1) > 2
         RawBackgroundZAverages(j) = nanmean(Ellipses{j}(:,9));
    end
    AllBackgroundMeasures = [AllBackgroundMeasures Ellipses{j}(~isnan(Ellipses{j}(:,9)),9).']; 
end
TotalBackgroundAverage = nanmean(AllBackgroundMeasures);
for j = 2:NumSlices-1
    if ~isnan( RawBackgroundZAverages(j))
        BackgroundZAverages(j) =  RawBackgroundZAverages(j);
    else
        BackgroundZAverages(j) = TotalBackgroundAverage;
    end
end
