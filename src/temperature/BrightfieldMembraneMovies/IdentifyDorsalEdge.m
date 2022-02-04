function IdentifyDorsalEdge(Prefix, varargin)


close all
cleanupObj = onCleanup(@myCleanupFun);

liveExperiment = LiveExperiment(Prefix);

FrameInfo = getFrameInfo(liveExperiment);
FrameTimes_seconds = [FrameInfo(:).Time];
FrameTimes_minutes = FrameTimes_seconds/60;

ProcPath = liveExperiment.userProcFolder;
DropboxFolder = liveExperiment.userResultsFolder;

Channel1 = liveExperiment.Channel1;
anaphaseFrames = liveExperiment.anaphaseFrames;
nc9 = anaphaseFrames(1);
nc10 = anaphaseFrames(2);
nc11 = anaphaseFrames(3);
nc12 = anaphaseFrames(4);
nc13 = anaphaseFrames(5);
nc14 = anaphaseFrames(6);

xSize = liveExperiment.xDim;
ySize = liveExperiment.yDim;
PixelSize_um = liveExperiment.pixelSize_um;


membraneMat = getMembraneMat(liveExperiment);

nFrames = size(membraneMat, 3);
%Get information about the image size
% HisImage=imread([PreProcPath,filesep,Prefix,filesep,D(1).name]);


%%
counter= 1;
CurrentFrame=nFrames;
MemImage = membraneMat(:,:,CurrentFrame);
DisplayRange=[min(min(MemImage)),max(max(MemImage))];

Overlay=figure;
set(Overlay,'units', 'normalized', 'position',[0.01, .1, .8, .8]);

overlayAxes = axes(Overlay,'Units', 'normalized', 'Position', [0 0 1 1]);

tb = axtoolbar(overlayAxes);
tb.Visible = 'off';


currentCharacter=1;



% Show the first image
imOverlay = imshow(MemImage,DisplayRange,'Parent',overlayAxes);

FigureTitle='Click points along dorsal embryonic boundary.';
%%
title(overlayAxes, FigureTitle);


projFlag = false;
set(0, 'CurrentFigure', Overlay)

if isfile([DropboxFolder,filesep,Prefix,filesep,'EmbryoOrientationInfo.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'EmbryoOrientationInfo.mat']);
end
if ~exist('EdgePoints', 'var')
    EdgePoints = NaN(0, 2, 'double');
end

if ~exist('p', 'var')
    p = [];
end
ShowFit =false;

while (currentCharacter~='x')
    [NEdgePoints,~]=size(EdgePoints);
    
    imOverlay.CData = MemImage;
    hold(overlayAxes, 'on')
    if exist('FitPlotHandle', 'var')
        delete(FitPlotHandle);
    end
    if NEdgePoints > 2
        p = polyfit(EdgePoints(:,1), EdgePoints(:,2),1);
        if ShowFit
            FitPlotHandle = plot(overlayAxes, 0:xSize, polyval(p, 0:xSize), 'g-');
        end
    end
    
    
    
    
    
    %refresh ellipses plots by destroying and remaking
    if exist('PlotHandle', 'var')
        cellfun(@delete, PlotHandle);
    end
    PlotHandle = cell(NEdgePoints, 1);
    for k=1:NEdgePoints
        
        PlotHandle{k} = plot(overlayAxes, EdgePoints(k,1),EdgePoints(k,2) ,  'r.', 'MarkerSize',20);
        
    end
    
    
    
    
    %%
    
    tb = axtoolbar(overlayAxes);
    tb.Visible = 'off';
    
    
    ct=waitforbuttonpress;
    currentCharacter=get(Overlay,'currentcharacter');
    currentMouse=get(overlayAxes,'CurrentPoint');
    
    
    
    if (ct~=0)&(currentCharacter=='p')
        [x,y] = ginput;
        for point_index = 1:length(x)
            EdgePoints(end+1,:)= [x(point_index),y(point_index)];
        end
        
    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'normal'))
        currentCharacter=1;
        if (currentMouse(1,2)>0)&(currentMouse(1,1)>0)&(currentMouse(1,2)<=ySize)&(currentMouse(1,1)<=xSize)
            
            %Add a circle to this location with the median radius of the
            %ellipses found in this frame
            
            %(x, y, a, b, theta, maxcontourvalue, time,
            
            
            EdgePoints(end+1,:)= [currentMouse(1,1),currentMouse(1,2)];
            
            
        end
        
        
        
        
    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'alt'))
        currentCharacter=1;
        if (currentMouse(1,2)>0)&(currentMouse(1,1)>0)&(currentMouse(1,2)<=ySize)&(currentMouse(1,1)<=xSize)
            %Find out which ellipses we clicked on so we can delete it
            
            %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
            Distances=sqrt((EdgePoints(:,1)-currentMouse(1,1)).^2+...
                (EdgePoints(:,2)-currentMouse(1,2)).^2);
            [~,MinIndex]=min(Distances);
            
            EdgePoints=[EdgePoints(1:MinIndex-1,:);...
                EdgePoints(MinIndex+1:end,:)];
        end
        
        
        
    elseif (ct~=0)&(currentCharacter=='m')    %Increase contrast
        DisplayRange(2)=DisplayRange(2)/1.5;
        
    elseif (ct~=0)&(currentCharacter=='n')    %Decrease contrast
        DisplayRange(2)=DisplayRange(2)*1.5;
        
    elseif (ct~=0)&(currentCharacter=='r')    %Reset the contrast
        DisplayRange=[min(min(MemImage)),max(max(MemImage))];
    elseif (ct~=0)&(currentCharacter=='f')    %Reset the contrast
        
        ShowFit = ~ShowFit;
        
    elseif (ct~=0)&(currentCharacter=='d')    %Delete all ellipses in the current frame
        EdgePoints=[];
        
        
        
        
        
        
        
        
    elseif (ct~=0)&(currentCharacter=='0')    %Debug mode
        keyboard
        
    end
end
%%
BoundaryAngle = atan(-p(1));
BoundaryAngleDegrees = BoundaryAngle*180/(pi);


%%
close all

SpotFig=figure;
set(SpotFig,'units', 'normalized', 'position',[0.01, .1, .8, .8]);


spotAxes = axes(SpotFig,'Units', 'normalized', 'Position', [0 0 1 1]);

tb = axtoolbar(spotAxes);
tb.Visible = 'off';

BlankImage = 255*ones(size(MemImage), 'uint8');

% Show the first image
imSpot = imshow(BlankImage,[0 255],'Parent',spotAxes);
hold(spotAxes, 'on')





%refresh ellipses plots by destroying and remaking
if exist('PlotHandle', 'var')
    cellfun(@delete, PlotHandle);
end
PlotHandle = cell(NEdgePoints, 1);
for k=1:NEdgePoints
    
    PlotHandle{k} = plot(spotAxes, EdgePoints(k,1),EdgePoints(k,2) ,  'k.', 'MarkerSize',20);
    
end



export_fig(spotAxes,[DropboxFolder, filesep,Prefix, filesep, 'BoundarySpotsImage.png'],...
    '-a1', '-png','-preserve_size');
close all

spotImage =imread([DropboxFolder, filesep,Prefix, filesep, 'BoundarySpotsImage.png']);
if size(spotImage,1) ~= size(MemImage,1)
    spotImage = imresize(spotImage,size(MemImage,1)/size(spotImage,1));
end






% Show the first image

%%
imRotatedSpots = imrotate(spotImage, -BoundaryAngleDegrees, 'bilinear', 'crop');
%imSpot2 = imshow(imRotatedSpots,'Parent',spotAxes2);
resize_scale_factor = 5;;
RotatedImage = double(imRotatedSpots);
%  se = strel('disk', 2);
% thresholdedResizedImage = imdilate(resizedRotatedImage, se);
RotatedImage(RotatedImage < 0) = 0;
RotatedImage(RotatedImage >255) = 255;
thresholdedImage = RotatedImage < 1;
[im_label, n_spots] = bwlabel(thresholdedImage);
%%
spot_stats = regionprops('table', im_label, 'Centroid', 'Area','Circularity');
spot_areas = cat(1,spot_stats.Area).';
spot_centroids =cat(1, spot_stats.Centroid);
approx_spot_area= median(spot_areas);
approx_spot_row = median(spot_centroids(:,2));
spot_row_range = [approx_spot_row-.05*ySize, approx_spot_row+.05*ySize];
spot_area_range = [approx_spot_area/5,approx_spot_area*5];
good_spots = spot_areas >= spot_area_range(1) &spot_areas <= spot_area_range(2) & ...
    spot_centroids(:,2).' >= spot_row_range(1)   & spot_centroids(:,2).' <= spot_row_range(2) ;
close all
figure(1)
imshow(thresholdedImage)
hold on
plot(spot_centroids(good_spots,1),spot_centroids(good_spots,2), 'r.')
hold off


RotatedEdgePoints = spot_centroids(good_spots,:);
NRotatedEdgePoints = size(RotatedEdgePoints,1);
BoundaryPosition = mean(RotatedEdgePoints(:,2));

save([DropboxFolder,filesep,Prefix,filesep,'EmbryoOrientationInfo.mat'],'EdgePoints', 'p',...
    'BoundaryAngle','BoundaryAngleDegrees', 'RotatedEdgePoints','BoundaryPosition', '-v6')

