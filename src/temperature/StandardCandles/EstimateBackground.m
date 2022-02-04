function EstimateBackground(Prefix, MaxPixelValue, LoadEllipses)
if ~exist('MaxPixelValue', 'var')
    MaxPixelValue = 20;
end
if ~exist('LoadEllipses', 'var')
    LoadEllipses = true;
end
liveExperiment = LiveExperiment(Prefix);
FrameInfo = getFrameInfo(liveExperiment);
nuclear_diameter = getDefaultParameters(FrameInfo,'d14');
pixelSize = liveExperiment.pixelSize_um;
d_pixels = nuclear_diameter/pixelSize;
outpath = [liveExperiment.resultsFolder, filesep, 'BackgroundEstimateEllipses.mat'];
xDim = liveExperiment.xDim;
yDim = liveExperiment.yDim;
zDim = FrameInfo(1).NumberSlices + 2;
NumFrames = length(FrameInfo);
TifList = dir([liveExperiment.preFolder, filesep, '*ch01.tif']);
MovieCells = cell(1, NumFrames);
NumSlices = 0;
DisplayRange = [0, MaxPixelValue];
FilterRadius = uint16(round(d_pixels/4));
FilterSize = 2*FilterRadius+1;
[y,x] =meshgrid(1:FilterSize, 1:FilterSize);% end
dist_mat = sqrt(((double(y)-double(FilterRadius+1)).^2 +(double(x)-double(FilterRadius+1)).^2));
bg_filter =dist_mat <= FilterRadius;

for i = 1:NumFrames
    MovieCells{i} = imreadStack([TifList(i).folder, filesep,  TifList(i).name]);
    NumSlices = NumSlices + size(MovieCells{i}, 3);
end
if LoadEllipses
    try
        load(outpath)
    catch
        Ellipses = cell(NumSlices, 1);
        nEllipses = 0;
        for zIndex = 1:NumSlices
            Ellipses{zIndex} = zeros(nEllipses,9, 'double');
        end
        
    end
else
    Ellipses = cell(NumSlices, 1);
    nEllipses = 0;
    for zIndex = 1:NumSlices
        Ellipses{zIndex} = zeros(nEllipses,9, 'double');
    end
end
%%
Overlay=figure;
set(Overlay,'units', 'normalized', 'position',[0.01, .1, .8, .8]);

overlayAxes = axes(Overlay,'Units', 'normalized', 'Position', [0 0 1 1]);
tb = axtoolbar(overlayAxes);
tb.Visible = 'off';

CurrentZ=2;

zImage = MovieCells{1}(:,:,CurrentZ);
currentCharacter=1;

% Show the first image
imOverlay = imshow(zImage,DisplayRange,'Border','Tight','Parent',overlayAxes);


set(0, 'CurrentFigure', Overlay)

while (currentCharacter~='x')
    
    
    zImage = MovieCells{1}(:,:,CurrentZ);
    
    
    %Get the information about the centroids
    [NCentroids,~]=size(Ellipses{CurrentZ});
    
    
    imOverlay.CData = zImage;
    try
        caxis(overlayAxes, DisplayRange);
    end
    
    %refresh ellipses plots by destroying and remaking
    if exist('PlotHandle', 'var')
        cellfun(@delete, PlotHandle);
    end
    
    PlotHandle = cell(NCentroids, 1);
    ellipseFrame = double(Ellipses{CurrentZ});
    for k=1:NCentroids
        n = k;
        %         PlotHandle{k} = drawellipse('Center',[ellipseFrame(n, 1) ellipseFrame(n, 2)],...
        %             'SemiAxes',[ellipseFrame(n, 3) ellipseFrame(n, 4)], ...
        %             'RotationAngle',ellipseFrame(n, 5) * (360/(2*pi)), 'FaceAlpha', 0,...
        %             'InteractionsAllowed', 'none', 'LabelVisible', 'hover', 'Label', num2str(ellipseFrame(n, 9)));
        
        PlotHandle{k} = ellipse(ellipseFrame(n, 3), ellipseFrame(n, 4),...
            ellipseFrame(n, 5) * (360/(2*pi)), ellipseFrame(n, 1),...
            ellipseFrame(n, 2), 'k', 10, overlayAxes);
        
        
        
        
        
        %                 set(PlotHandle{k}, 'StripeColor', 'w', 'Color', 'w','Linewidth', 1);
        set(PlotHandle{k}, 'Color', 'r','Linewidth', 1);
        
        
        
        
    end
    
    
    FigureTitle=['Z: ',num2str(CurrentZ),'/',num2str(NumSlices)];
    
    
    set(Overlay,'Name',FigureTitle)
    
    
    
    
    
    
    
    %%
    
    tb = axtoolbar(overlayAxes);
    tb.Visible = 'off';
    
    ct=waitforbuttonpress;
    currentCharacter=get(Overlay,'currentcharacter');
    currentMouse=get(overlayAxes,'CurrentPoint');
    
    
    
    
    if (ct~=0)&(currentCharacter=='.')&(CurrentZ<NumSlices)
        CurrentZ=CurrentZ+1;
    elseif (ct~=0)&(currentCharacter==',')&(CurrentZ>1)
        CurrentZ=CurrentZ-1;
    elseif (ct~=0)&(currentCharacter=='s')
        save(outpath,'Ellipses', '-v6')
        disp('Ellipses saved.')
    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'normal'))
        currentCharacter=1;
        if (currentMouse(1,2)>0)&(currentMouse(1,1)>0)&(currentMouse(1,2)<=yDim)&(currentMouse(1,1)<=xDim)
            
            %Add a circle to this location with the median radius of the
            %ellipses found in this frame
            
            %(x, y, a, b, theta, maxcontourvalue, time,
            %particle_id)
            
            
            centroid = [uint16(round(currentMouse(1,1))),uint16(round(currentMouse(1,2)))];
            snip = double(zImage(centroid(2)-FilterRadius:centroid(2)+FilterRadius,...
                centroid(1)-FilterRadius:centroid(1)+FilterRadius));
            bgd = sum(sum(snip.*double(bg_filter)))/double(sum(sum(bg_filter)));
            
            Ellipses{CurrentZ}(end+1,:)=...
                [double(currentMouse(1,1)),double(currentMouse(1,2)),double(FilterRadius),double(FilterRadius),0,0,0,0,bgd];
            
            
        end
        
        
        
        
    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'alt'))
        currentCharacter=1;
        if ~isempty(Ellipses{CurrentZ})&(currentMouse(1,2)>0)&(currentMouse(1,1)>0)&(currentMouse(1,2)<=yDim)&(currentMouse(1,1)<=xDim)
            %Find out which ellipses we clicked on so we can delete it
            
            %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
            Distances=sqrt((Ellipses{CurrentZ}(:,1)-currentMouse(1,1)).^2+...
                (Ellipses{CurrentZ}(:,2)-currentMouse(1,2)).^2);
            [~,MinIndex]=min(Distances);
            
            Ellipses{CurrentZ}=[Ellipses{CurrentZ}(1:MinIndex-1,:);...
                Ellipses{CurrentZ}(MinIndex+1:end,:)];
        end
        
    elseif (ct~=0)&(currentCharacter=='j')
        iJump=input('Z to jump to: ');
        if (floor(iJump)>0)&(iJump<=NumSlices)
            CurrentZ=iJump;
        else
            disp('Z out of range.');
        end
        
    elseif (ct~=0)&(currentCharacter=='m')    %Increase contrast
        DisplayRange(2)=DisplayRange(2)/1.5;
        
    elseif (ct~=0)&(currentCharacter=='n')    %Decrease contrast
        DisplayRange(2)=DisplayRange(2)*1.5;
        
    elseif (ct~=0)&(currentCharacter=='r')    %Reset the contrast
        DisplayRange=[0, MaxPixelValue];
        
    elseif (ct~=0)&(currentCharacter=='d')    %Delete all ellipses in the current frame
        Ellipses{CurrentZ}=[];
        
        
        
        
        
    elseif (ct~=0)&(currentCharacter=='0')    %Debug mode
        keyboard
        
    end
end
close all

%%
for CurrentZ = 1:NumSlices
    zImage = MovieCells{1}(:,:,CurrentZ);
    for k = 1:size((Ellipses{CurrentZ}),1)
        centroid = [uint16(round(Ellipses{CurrentZ}(k,1))),uint16(round(Ellipses{CurrentZ}(k,2)))];
        snip = double(zImage(centroid(2)-FilterRadius:centroid(2)+FilterRadius,...
            centroid(1)-FilterRadius:centroid(1)+FilterRadius));
        bgd = sum(sum(snip.*double(bg_filter)))/double(sum(sum(bg_filter)));
        
        Ellipses{CurrentZ}(k,:)=...
            [Ellipses{CurrentZ}(k,1),Ellipses{CurrentZ}(k,2),double(FilterRadius),double(FilterRadius),0,0,0,0,bgd];
    end
end
%%
[BackgroundZAverages, RawBackgroundZAverages, TotalBackgroundAverage] = ...
    CalculateBackgroundAverages(Prefix, Ellipses);
save(outpath,'Ellipses','BackgroundZAverages','RawBackgroundZAverages',...
    'TotalBackgroundAverage','-v6')


