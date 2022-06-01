function PredictEmbryoNCs(Prefix,UseCustom, varargin)
%%
ShowPlots = false;
if ~exist('UseCustom', 'var')
    UseCustom = false;
end
liveExperiment= LiveExperiment(Prefix);

CompiledEmbryoPath = [liveExperiment.resultsFolder, filesep, 'CompiledEmbryos.Mat'];
load(CompiledEmbryoPath);
%%
ProcPath = liveExperiment.userProcFolder;
DropboxFolder = liveExperiment.userResultsFolder;

Channel1 = liveExperiment.Channel1;
Channel2 = liveExperiment.Channel2;
Channel3 = liveExperiment.Channel3;
Channel4 = liveExperiment.Channel4;
Channel5 = liveExperiment.Channel5;


xSize = liveExperiment.xDim;
ySize = liveExperiment.yDim;
PixelSize_um = liveExperiment.pixelSize_um;

%Get the nuclei segmentation data
Ellipses = getEllipses(liveExperiment);
try
    clrmp = single(hsv(20));
    clrmp = clrmp(randperm(length(clrmp)), :);
catch
    %in case the user doesn't have this colormap, just keep going.
end
Channels = {Channel1, Channel2, Channel3, Channel4, Channel5};

if UseCustom
    RotatedHisFile = [liveExperiment.preFolder, filesep, Prefix, '-CustomHis_Rotated.tif'];
else
    RotatedHisFile = [liveExperiment.preFolder, filesep, Prefix, '-His_Rotated.tif'];
end
RotatedMembraneFile =[liveExperiment.preFolder, filesep, Prefix, '-Membrane_Rotated.tif'];
hisMat = imreadStack2(RotatedHisFile, liveExperiment.yDim, liveExperiment.xDim, liveExperiment.nFrames);
memMat = imreadStack2(RotatedMembraneFile, liveExperiment.yDim, liveExperiment.xDim, liveExperiment.nFrames);

%%
hisMatCopy = hisMat;
xSize = size(hisMat,2);
ySize = size(hisMat,1);
close all
nEmbryos = size(hisMat, 3);
%Get information about the image size
% HisImage=imread([PreProcPath,filesep,Prefix,filesep,D(1).name]);
GoodEmbryos = 1:nEmbryos;
GoodEmbryos = GoodEmbryos(CompiledEmbryos.Approved);
for CurrentEmbryo = GoodEmbryos



HisImage = hisMat(:,:,CurrentEmbryo);
HisImage2 = hisMat(:,:,CurrentEmbryo);
MemImage = memMat(:,:,CurrentEmbryo);


DisplayRangeHis = [min(min(HisImage)), max(max(HisImage))];
DisplayRangeHis2 = [min(min(HisImage2)), max(max(HisImage2))];
if ShowPlots
close all

FullFigure=figure(1);
set(FullFigure,'units', 'normalized', 'position',[0.01, .3, .9, .5]);

fullAxes = axes(FullFigure,'Units', 'normalized', 'Position', [0 0 1 1]);
% Overlay=figure(2);
% set(Overlay,'units', 'normalized', 'position',[0.01, .1, .45, .65]);
%
% overlayAxes = axes(Overlay,'Units', 'normalized', 'Position', [0 0 1 1]);
%
%
% Original=figure(3);
% set(Original,'units', 'normalized', 'position',[0.51, .1, .45, .65]);
%
% originalAxes = axes(Original,'Units', 'normalized', 'Position', [0 0 1 1]);




tb = axtoolbar(fullAxes);
tb.Visible = 'off';
% tb2 = axtoolbar(originalAxes);
% tb2.Visible = 'off';
% tb3 = axtoolbar(originalAxes);
% tb3.Visible = 'off';

imFull = imshow(HisImage2,DisplayRangeHis2,'Border','Tight','Parent',fullAxes);

SwitchImageType = false;
hold(fullAxes,'on')

set(0, 'CurrentFigure', FullFigure)
end
%%
%Get the information about the centroids
[NCentroids,~]=size(Ellipses{CurrentEmbryo});

if ShowPlots
imFull.CData = HisImage2;
try
    caxis(fullAxes, DisplayRange);
    
end
%refresh ellipses plots by destroying and remaking
if exist('PlotHandle', 'var')
    cellfun(@delete, PlotHandle);
end

PlotHandle = cell(NCentroids, 1);
end
ellipseFrame = double(Ellipses{CurrentEmbryo});
CoordA = CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,:);
CoordP = CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,:);
CoordD = CompiledEmbryos.RotatedCoordDs(CurrentEmbryo,:);
CoordV = CompiledEmbryos.RotatedCoordVs(CurrentEmbryo,:);
%bIndex = boundary(ellipseFrame(:,1), ellipseFrame(:,2), 1);
APLength = CoordP(1)-CoordA(1);
DVLength = CoordV(2)-CoordD(2);
DorsalPoints = [];
for k = 1:NCentroids
    if (ellipseFrame(k, 1) > CoordA(1) + 0.1*APLength) & ...
            (ellipseFrame(k, 1) < CoordP(1) - 0.1*APLength)  & ...
            ellipseFrame(k, 2) < CoordA(2)
        DorsalPoints(end+1) = k;
    end
end

dorsalEllipses = ellipseFrame(DorsalPoints,:);
if ~isempty(dorsalEllipses)
dorsalEllipses = sortrows(dorsalEllipses, 1);
NuclearDiameter = 2*min(dorsalEllipses(:,3));
IncludedNuclei = [1];
NucleiCount = 1;
NuclearIndex = 1;
NextNuclearIndex = 2;
while NuclearIndex < size(dorsalEllipses, 1) & NextNuclearIndex <= size(dorsalEllipses, 1)
    if dorsalEllipses(NextNuclearIndex,1) <= dorsalEllipses(NuclearIndex,1)+0.5*NuclearDiameter
        NextNuclearIndex = NextNuclearIndex+1;
    else
        IncludedNuclei(end+1) = NextNuclearIndex;
        NucleiCount = NucleiCount+1;
        NuclearIndex = NextNuclearIndex;
        NextNuclearIndex = NuclearIndex+1;
    end
end
    
dorsalEllipses = dorsalEllipses(IncludedNuclei,:);
if ShowPlots
for k=1:size(dorsalEllipses, 1)
    n =k;
    %         PlotHandle{k} = drawellipse('Center',[ellipseFrame(n, 1) ellipseFrame(n, 2)],...
    %             'SemiAxes',[ellipseFrame(n, 3) ellipseFrame(n, 4)], ...
    %             'RotationAngle',ellipseFrame(n, 5) * (360/(2*pi)), 'FaceAlpha', 0,...
    %             'InteractionsAllowed', 'none', 'LabelVisible', 'hover', 'Label', num2str(ellipseFrame(n, 9)));
    colorhash = uint8(mod(round(dorsalEllipses(n, 1)+dorsalEllipses(n, 2)),20)+1);
    PlotHandle{k} = ellipse(2*dorsalEllipses(n, 3), 2*dorsalEllipses(n, 4),...
        dorsalEllipses(n, 5) * (360/(2*pi)), dorsalEllipses(n, 1),...
        dorsalEllipses(n, 2), clrmp(colorhash,:), 10, fullAxes, 0.5);
    
    
end

end
disp(['Current Embryo: ', num2str(CurrentEmbryo), ', Nuclear Count: ', num2str(NucleiCount)]);
if NucleiCount > 65
    CompiledEmbryos.nc(CurrentEmbryo) = 14;
elseif NucleiCount > 55
    CompiledEmbryos.nc(CurrentEmbryo) = 13;
end
else
    CompiledEmbryos.nc(CurrentEmbryo) = NaN;
end
end
if ShowPlots
    close all
end
save([DropboxFolder,filesep,Prefix,filesep,'CompiledEmbryos.mat'],'CompiledEmbryos', '-v6')
