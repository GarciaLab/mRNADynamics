function SelectMembraneZSlices(Prefix)
% Midpoints are in columns first (x) then rows (y) - y is smaller when
% higher up
liveExperiment = LiveExperiment(Prefix);
FrameInfo = getFrameInfo(liveExperiment);
load([liveExperiment.resultsFolder, filesep, 'MarkAndFindInfo.Mat'], 'MarkAndFindInfo');
PreProcFolder = liveExperiment.preFolder;


NEmbryos = MarkAndFindInfo.NSeries;
MembraneZoomResultsFolder = [liveExperiment.resultsFolder, filesep, 'ZoomMembraneInfo'];


MembraneZoomPixelPath = [MembraneZoomResultsFolder, filesep, 'MembraneZoomPixelSize.mat'];
load(MembraneZoomPixelPath,'PixelSize_um');
MembraneZoomImageSizePath = [MembraneZoomResultsFolder, filesep, 'MembraneZoomImageSize.mat'];
load(MembraneZoomImageSizePath);


CompiledEmbryoPath = [liveExperiment.resultsFolder, filesep, 'CompiledEmbryos.Mat'];
load(CompiledEmbryoPath);
CompiledEmbryos = UpdateAPAxisInfo(Prefix, CompiledEmbryos);
save(CompiledEmbryoPath,'CompiledEmbryos');

%%


RotatedMembraneMat = zeros(ySize/2, xSize, zSize, NEmbryos,'double');
for embryoIndex = 1:NEmbryos
    if CompiledEmbryos.Approved(embryoIndex) & CompiledEmbryos.MemCoordAs(embryoIndex,1) > 0 & CompiledEmbryos.MemCoordPs(embryoIndex,1) > 0 & CompiledEmbryos.MemCoordDs(embryoIndex,1) > 0
        
        NewName = [Prefix, '_Position',iIndex(embryoIndex,3),...
            '_RotatedZoomMembrane', '.tif'];
        RotatedMembraneMat(:,:,:,embryoIndex) = imreadStack2([PreProcFolder, filesep, NewName]);
    end
end

%%
DownsizedMembraneMat = zeros(ySize/2, xSize, zSize, NEmbryos,'uint8');
for embryoIndex = 1:NEmbryos
    if CompiledEmbryos.Approved(embryoIndex) & CompiledEmbryos.MemCoordAs(embryoIndex,1) > 0 & CompiledEmbryos.MemCoordPs(embryoIndex,1) > 0 & CompiledEmbryos.MemCoordDs(embryoIndex,1) > 0
        MaxPixelVal = max(max(max(RotatedMembraneMat(:,:,:,embryoIndex))));
        MinPixelVal = min(min(min(RotatedMembraneMat(:,:,:,embryoIndex))));
        
        DownsizedMembraneMat(:,:,:,embryoIndex) = uint8((RotatedMembraneMat(:,:,:,embryoIndex)-MinPixelVal)/(MaxPixelVal-MinPixelVal)*256-1);
    end
end

memImages = {};
IncludedImagesFile = [liveExperiment.resultsFolder, filesep, 'MembraneIncludeImagesInfo.mat'];
if isfile(IncludedImagesFile)
    load(IncludedImagesFile)
else
    includedImages = zeros(NEmbryos, zSize, 'logical');
end
MaxImages = zeros(ySize, xSize, NEmbryos, 'uint8');
MinImages = zeros(ySize, xSize, NEmbryos, 'uint8');
MedianImages = zeros(ySize, xSize, NEmbryos, 'uint8');


%%


embryoIndex=find(CompiledEmbryos.Approved,1);
DVLength = abs(CompiledEmbryos.MemRotatedCoordVs(embryoIndex, 2)-CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2));
ylims = [round(0.9*CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2)),...
    round(DVLength*1/3+CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2))];
APLength=abs(CompiledEmbryos.MemRotatedCoordPs(embryoIndex, 1)-CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1));
xlims = [round(CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)+APLength*.3),...
    round(CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)+APLength*.7)];

APGrid = 0.3:0.025:0.7;
APmarkers = APGrid*APLength+CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)-xlims(1);

xRange = xlims(1):xlims(2);
yRange = ylims(1):ylims(2);





zIndex = zSize;
cc=1;
EmbryoImage = DownsizedMembraneMat(yRange,xRange,zIndex,embryoIndex);
MaxImage = max(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)), [], 3);
MedImage = median(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)), 3);
MinImage = min(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)), [], 3);
if isempty(MaxImage)
    MaxImage = zeros(size(EmbryoImage), 'uint8');
    MedImage = zeros(size(EmbryoImage), 'uint8');
    MinImage = zeros(size(EmbryoImage), 'uint8');
end
DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
MaxDisplayRange=[min(min(MaxImage)),max(max(MaxImage))];
MinDisplayRange=[min(min(MinImage)),max(max(MinImage))];
MedDisplayRange=[min(min(MedImage)),max(max(MedImage))];

close all
EmbryoFigure=figure;
set(EmbryoFigure,'units', 'normalized', 'position',[0.1, .65, .7, .25]);

embryoAxes = axes(EmbryoFigure,'Units', 'normalized', 'Position', [0 0 1 1]);

% Show the first image
imEmbryo = imshow(EmbryoImage,DisplayRange,'Parent',embryoAxes);
hold on
GuidePlotHandles = cell(1,length(APmarkers));
for pl=1:length(APmarkers)
    GuidePlotHandles{pl} = plot(embryoAxes, [APmarkers(pl),APmarkers(pl)], [0 size(imEmbryo.CData,2)], 'b-+');
    GuidePlotHandles{pl}.Color(4) =0.3;
end



if (includedImages(embryoIndex, zIndex))
    set(EmbryoFigure,'Color','g')
else
    set(EmbryoFigure,'Color','r')
end
% axis image
axis off
hold off



FigureTitle=['Embryo Index: ',num2str(embryoIndex),'/',num2str(NEmbryos),...
    ', z Index: ',num2str(zIndex)];

set(EmbryoFigure,'Name',FigureTitle)


MaxFigure=figure;
set(MaxFigure,'units', 'normalized', 'position',[0.1, .35, .7, .2]);

maxAxes = axes(MaxFigure,'Units', 'normalized', 'Position', [0 0 1 1]);

% Show the first image
imMax = imshow(MaxImage,MaxDisplayRange,'Parent',maxAxes);



% axis image
axis off
hold off



MaxFigureTitle=['Max Projected Membrane Image'];

set(MaxFigure,'Name',MaxFigureTitle)


MinFigure=figure;
set(MinFigure,'units', 'normalized', 'position',[0.1, .35, .7, .2]);

minAxes = axes(MinFigure,'Units', 'normalized', 'Position', [0 0 1 1]);

% Show the first image
imMin = imshow(MinImage,MinDisplayRange,'Parent',minAxes);



% axis image
axis off
hold off



MinFigureTitle=['Min Projected Membrane Image'];

set(MinFigure,'Name',MinFigureTitle)

MedianFigure=figure;
set(MedianFigure,'units', 'normalized', 'position',[0.1, 0.05, .7, .2]);

medAxes = axes(MedianFigure,'Units', 'normalized', 'Position', [0 0 1 1]);

% Show the first image
imMed = imshow(MedImage,MedDisplayRange,'Parent',medAxes);

% axis image
axis off
hold off



MedFigureTitle=['Median Projected Membrane Image'];
set(MedianFigure,'Name',MedFigureTitle)

while (cc~='x')
    
    imEmbryo.CData = EmbryoImage;
    try
        caxis(embryoAxes, DisplayRange);
    end
    
    imMax.CData = MaxImage;
    try
        caxis(maxAxes, MaxDisplayRange);
    end
    
    imMin.CData = MinImage;
    try
        caxis(minAxes, MinDisplayRange);
    end
    
    imMed.CData = MedImage;
    try
        caxis(medAxes, MedDisplayRange);
    end
    
    for pl=1:length(APmarkers)
        GuidePlotHandles{pl}.YData =  [0 size(imEmbryo.CData,2)];
        GuidePlotHandles{pl}.XData  = [APmarkers(pl), APmarkers(pl)];
    end
    
    if includedImages(embryoIndex, zIndex)
        set(EmbryoFigure,'Color','g')
    else
        set(EmbryoFigure,'Color','r')
    end
    
    hold off
    
    
    
    FigureTitle=['Embryo Index: ',num2str(embryoIndex),'/',num2str(NEmbryos),...
        ', z Index: ',num2str(zIndex)];
    
    
    set(EmbryoFigure,'Name',FigureTitle)
    
    figure(EmbryoFigure)
    ct=waitforbuttonpress;
    cc=get(EmbryoFigure,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    
    if (ct~=0)&(cc=='q')        %Clear all AP information
        includedImages(embryoIndex, zIndex) = true;
        if zIndex > 1
            zIndex = zIndex-1;
        end
        DVLength = abs(CompiledEmbryos.MemRotatedCoordVs(embryoIndex, 2)-CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2));
        ylims = [round(0.9*CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2)),...
            round(DVLength*1/3+CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2))];
        APLength=abs(CompiledEmbryos.MemRotatedCoordPs(embryoIndex, 1)-CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1));
        xlims = [round(CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)+APLength*.3),...
            round(CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)+APLength*.7)];
        
        APGrid = 0.3:0.025:0.7;
        APmarkers = APGrid*APLength+CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)-xlims(1);
        
        xRange = xlims(1):xlims(2);
        yRange = ylims(1):ylims(2);
        
        
        
        
        
        EmbryoImage = DownsizedMembraneMat(yRange,xRange,zIndex,embryoIndex);
        MaxImage = max(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)), [], 3);
        MinImage = min(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)),[], 3);
        MedImage = median(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)), 3);
        if isempty(MaxImage)
            MaxImage = zeros(size(EmbryoImage), 'uint8');
            MinImage = zeros(size(EmbryoImage), 'uint8');
            MedImage = zeros(size(EmbryoImage), 'uint8');
        end
        DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
        MaxDisplayRange=[min(min(MaxImage)),max(max(MaxImage))];
        MedDisplayRange=[min(min(MedImage)),max(max(MedImage))];
        MinDisplayRange=[min(min(MinImage)),max(max(MinImage))];
    elseif (ct~=0)&(cc=='w')	%Select anterior end
        includedImages(embryoIndex, zIndex) = false;
        if zIndex > 1
            zIndex = zIndex-1;
        end
        DVLength = abs(CompiledEmbryos.MemRotatedCoordVs(embryoIndex, 2)-CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2));
        ylims = [round(0.9*CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2)),...
            round(DVLength*1/3+CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2))];
        APLength=abs(CompiledEmbryos.MemRotatedCoordPs(embryoIndex, 1)-CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1));
        xlims = [round(CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)+APLength*.3),...
            round(CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)+APLength*.7)];
        
        APGrid = 0.3:0.025:0.7;
        APmarkers = APGrid*APLength+CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)-xlims(1);
        
        xRange = xlims(1):xlims(2);
        yRange = ylims(1):ylims(2);
        
        
        
        
        
        EmbryoImage = DownsizedMembraneMat(yRange,xRange,zIndex,embryoIndex);
        MaxImage = max(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)), [], 3);
        MinImage = min(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)), [], 3);
        MedImage = median(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)), 3);
        if isempty(MaxImage)
            MaxImage = zeros(size(EmbryoImage), 'uint8');
            MinImage = zeros(size(EmbryoImage), 'uint8');
            MedImage = zeros(size(EmbryoImage), 'uint8');
        end
        DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
        MaxDisplayRange=[min(min(MaxImage)),max(max(MaxImage))];
        MedDisplayRange=[min(min(MedImage)),max(max(MedImage))];
        MinDisplayRange=[min(min(MinImage)),max(max(MinImage))];
    elseif (ct~=0)&(cc=='.')&((zIndex < zSize)|(embryoIndex < NEmbryos) )   %Increase contrast
        if zIndex > 1
            zIndex = zIndex-1;
        end
        DVLength = abs(CompiledEmbryos.MemRotatedCoordVs(embryoIndex, 2)-CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2));
        ylims = [round(0.9*CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2)),...
            round(DVLength*1/3+CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2))];
        APLength=abs(CompiledEmbryos.MemRotatedCoordPs(embryoIndex, 1)-CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1));
        xlims = [round(CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)+APLength*.3),...
            round(CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)+APLength*.7)];
        
        APGrid = 0.3:0.025:0.7;
        APmarkers = APGrid*APLength+CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)-xlims(1);
        
        xRange = xlims(1):xlims(2);
        yRange = ylims(1):ylims(2);
        
        
        
        
        
        EmbryoImage = DownsizedMembraneMat(yRange,xRange,zIndex,embryoIndex);
        MaxImage = max(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)), [], 3);
        MinImage = min(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)),[], 3);
        MedImage = median(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)), 3);
        if isempty(MaxImage)
            MaxImage = zeros(size(EmbryoImage), 'uint8');
            MedImage = zeros(size(EmbryoImage), 'uint8');
            MinImage = zeros(size(EmbryoImage), 'uint8');
        end
        DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
        MaxDisplayRange=[min(min(MaxImage)),max(max(MaxImage))];
        MedDisplayRange=[min(min(MedImage)),max(max(MedImage))];
        MinDisplayRange=[min(min(MinImage)),max(max(MinImage))];
    elseif (ct~=0)&(cc==',')&((embryoIndex > 1) | (zIndex > 1))   %Increase contrast
        if zIndex < zSize
            zIndex = zIndex + 1;
        end

        DVLength = abs(CompiledEmbryos.MemRotatedCoordVs(embryoIndex, 2)-CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2));
        ylims = [round(0.9*CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2)),...
            round(DVLength*1/3+CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2))];
        APLength=abs(CompiledEmbryos.MemRotatedCoordPs(embryoIndex, 1)-CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1));
        xlims = [round(CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)+APLength*.3),...
            round(CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)+APLength*.7)];
        
        APGrid = 0.3:0.025:0.7;
        APmarkers = APGrid*APLength+CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)-xlims(1);
        
        xRange = xlims(1):xlims(2);
        yRange = ylims(1):ylims(2);
        
        
        
        
        
        EmbryoImage = DownsizedMembraneMat(yRange,xRange,zIndex,embryoIndex);
        MaxImage = max(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)), [], 3);
        MinImage = min(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)),[], 3);
        MedImage = median(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)), 3);
        if isempty(MaxImage)
            MaxImage = zeros(size(EmbryoImage), 'uint8');
            MedImage = zeros(size(EmbryoImage), 'uint8');
            MinImage = zeros(size(EmbryoImage), 'uint8');
        end
        DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
        MaxDisplayRange=[min(min(MaxImage)),max(max(MaxImage))];
        MedDisplayRange=[min(min(MedImage)),max(max(MedImage))];
        MinDisplayRange=[min(min(MinImage)),max(max(MinImage))];
    elseif (ct~=0)&(cc=='>')&(embryoIndex < NEmbryos)
        if embryoIndex < max(find(CompiledEmbryos.Approved))
            embryoIndex = find(CompiledEmbryos.Approved(embryoIndex+1:end),1)+embryoIndex;
            zIndex = zSize;
        end
        DVLength = abs(CompiledEmbryos.MemRotatedCoordVs(embryoIndex, 2)-CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2));
        ylims = [round(0.9*CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2)),...
            round(DVLength*1/3+CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2))];
        APLength=abs(CompiledEmbryos.MemRotatedCoordPs(embryoIndex, 1)-CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1));
        xlims = [round(CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)+APLength*.3),...
            round(CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)+APLength*.7)];
        
        APGrid = 0.3:0.025:0.7;
        APmarkers = APGrid*APLength+CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)-xlims(1);
        
        xRange = xlims(1):xlims(2);
        yRange = ylims(1):ylims(2);
        
        
        
        
        
        EmbryoImage = DownsizedMembraneMat(yRange,xRange,zIndex,embryoIndex);
        MaxImage = max(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)), [], 3);
        MinImage = min(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)),[], 3);
        MedImage = median(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)), 3);
        if isempty(MaxImage)
            MaxImage = zeros(size(EmbryoImage), 'uint8');
            MedImage = zeros(size(EmbryoImage), 'uint8');
            MinImage = zeros(size(EmbryoImage), 'uint8');
        end
        DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
        MaxDisplayRange=[min(min(MaxImage)),max(max(MaxImage))];
        MedDisplayRange=[min(min(MedImage)),max(max(MedImage))];
        MinDisplayRange=[min(min(MinImage)),max(max(MinImage))];
    elseif (ct~=0)&(cc=='<')&(embryoIndex > 1)
        if embryoIndex > min(find(CompiledEmbryos.Approved))
            embryoIndex = find(CompiledEmbryos.Approved(1:embryoIndex-1),1, 'last');
            zIndex = zSize;
        end
        DVLength = abs(CompiledEmbryos.MemRotatedCoordVs(embryoIndex, 2)-CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2));
        ylims = [round(0.9*CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2)),...
            round(DVLength*1/3+CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2))];
        APLength=abs(CompiledEmbryos.MemRotatedCoordPs(embryoIndex, 1)-CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1));
        xlims = [round(CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)+APLength*.3),...
            round(CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)+APLength*.7)];
        
        APGrid = 0.3:0.025:0.7;
        APmarkers = APGrid*APLength+CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)-xlims(1);
        
        xRange = xlims(1):xlims(2);
        yRange = ylims(1):ylims(2);
        
        
        
        
        
        EmbryoImage = DownsizedMembraneMat(yRange,xRange,zIndex,embryoIndex);
        MaxImage = max(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)), [], 3);
        MinImage = min(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)),[], 3);
        MedImage = median(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)), 3);
        if isempty(MaxImage)
            MaxImage = zeros(size(EmbryoImage), 'uint8');
            MedImage = zeros(size(EmbryoImage), 'uint8');
            MinImage = zeros(size(EmbryoImage), 'uint8');
        end
        DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
        MaxDisplayRange=[min(min(MaxImage)),max(max(MaxImage))];
        MedDisplayRange=[min(min(MedImage)),max(max(MedImage))];
        MinDisplayRange=[min(min(MinImage)),max(max(MinImage))];
    elseif (ct~=0)&(cc=='m')    %Increase contrast
        DisplayRange(2)=DisplayRange(2)/2;
    elseif (ct~=0)&(cc=='n')    %Increase contrast
        DisplayRange(2)=DisplayRange(2)*2;
    elseif (ct~=0)&(cc=='j')
        try
            iJump = inputdlg('Embryo to jump to:', ...
                'Move to frame');
            iJump = str2double(iJump{1});
        catch
            iJump =embryoIndex;
        end
        if ismember(iJump, find(CompiledEmbryos.Approved))
            embryoIndex= iJump;
            zIndex = zSize;
        end
        
        DVLength = abs(CompiledEmbryos.MemRotatedCoordVs(embryoIndex, 2)-CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2));
        ylims = [round(0.9*CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2)),...
            round(DVLength*1/3+CompiledEmbryos.MemRotatedCoordDs(embryoIndex, 2))];
        APLength=abs(CompiledEmbryos.MemRotatedCoordPs(embryoIndex, 1)-CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1));
        xlims = [round(CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)+APLength*.3),...
            round(CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)+APLength*.7)];
        
        APGrid = 0.3:0.025:0.7;
        APmarkers = APGrid*APLength+CompiledEmbryos.MemRotatedCoordAs(embryoIndex, 1)-xlims(1);
        
        xRange = xlims(1):xlims(2);
        yRange = ylims(1):ylims(2);
        
        
        
        
        
        EmbryoImage = DownsizedMembraneMat(yRange,xRange,zIndex,embryoIndex);
        MaxImage = max(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)), [], 3);
        MinImage = min(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)),[], 3);
        MedImage = median(squeeze(DownsizedMembraneMat(yRange,xRange,includedImages(embryoIndex,:),embryoIndex)), 3);
        if isempty(MaxImage)
            MaxImage = zeros(size(EmbryoImage), 'uint8');
            MedImage = zeros(size(EmbryoImage), 'uint8');
            MinImage = zeros(size(EmbryoImage), 'uint8');
        end
        DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
        MaxDisplayRange=[min(min(MaxImage)),max(max(MaxImage))];
        MedDisplayRange=[min(min(MedImage)),max(max(MedImage))];
        MinDisplayRange=[min(min(MinImage)),max(max(MinImage))];
    elseif (ct~=0)&(cc=='r')    %Reset the contrast
        DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
        
    end
end


%%
close all
MaxMat = zeros(ySize/2, xSize, NEmbryos, 'uint8');
MedMat = zeros(ySize/2, xSize, NEmbryos, 'uint8');
MinMat = zeros(ySize/2, xSize, NEmbryos, 'uint8');
for embryoIndex = 1:NEmbryos
    if CompiledEmbryos.Approved(embryoIndex)
        try
        MaxMat(:,:,embryoIndex) = max(squeeze(DownsizedMembraneMat(:,:,includedImages(embryoIndex,:),embryoIndex)), [], 3);
        MedMat(:,:,embryoIndex) = median(squeeze(DownsizedMembraneMat(:,:,includedImages(embryoIndex,:),embryoIndex)), 3);
        MinMat(:,:,embryoIndex) = min(squeeze(DownsizedMembraneMat(:,:,includedImages(embryoIndex,:),embryoIndex)),[], 3);
        end
    end
end


RotatedMaxMembraneCustomFile = [liveExperiment.preFolder, filesep,Prefix, '-MaxZoomMembraneRotated.tif'];
    saveNuclearProjection(MaxMat,RotatedMaxMembraneCustomFile);
    
    
RotatedMedMembraneCustomFile = [liveExperiment.preFolder, filesep,Prefix, '-MedianZoomMembraneRotated.tif'];
    saveNuclearProjection(MedMat,RotatedMedMembraneCustomFile);
    
    
RotatedMinMembraneCustomFile = [liveExperiment.preFolder, filesep,Prefix, '-MinZoomMembraneRotated.tif'];
    saveNuclearProjection(MinMat,RotatedMinMembraneCustomFile);
   
IncludedImagesFile = [liveExperiment.resultsFolder, filesep, 'MembraneIncludeImagesInfo.mat'];
save(IncludedImagesFile,'includedImages');



