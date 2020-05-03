function movieMat = CheckNucleiSegmentation(Prefix, varargin)
%
%
%
%Usage:
%
% .  - Move a frame forward
% ,  - Move a frame backwards
% >  - Move 5 frames forward
% <  - Move 5 frames backwards
% j  - Jump to a frame
% q  - Move a cycle forward
% w  - Move a cycle backwards
% d  - Delete all ellipses in the current frame
% c  - Copy all ellipses from previous frame
% v  - Copy all ellipses from next frame
% s  - Save current analysis
% ~  - Create a different projection for the nuclear image
% m  - Increase contrast
% n  - Decrease contrast
% r  - Reset contrast setting
% / - Adjust ellipse centroids
% x  - Exit and save
% 9  - Debug mode
%
%
%right click  - delete region
%left click - add region with default nc radius and angle
%

cleanupObj = onCleanup(@myCleanupFun);


noAdd = false;
nWorkers = 1;
fish = false;
preMovie = false;
chooseHis = false;
yToRetrackPrompt = true;

for k = 1:length(varargin)
    if strcmpi(varargin{k}, 'noAdd') | strcmpi(varargin{k}, 'fish') | strcmpi(varargin{k}, 'markandfind')
        noAdd = true;
        fish = true;
    elseif strcmpi(varargin{k}, 'nWorkers')
        nWorkers = varargin{k+1};
    elseif strcmpi(varargin{k}, 'chooseHis')
        chooseHis = varargin{k+1};
    elseif strcmpi(varargin{k}, 'colormap')
        cmap = varargin{k+1};
    elseif strcmpi(varargin{k}, 'premovie')
        preMovie = true;
        movieMat = varargin{k+1};
    elseif strcmpi(varargin{k}, 'yToRetrackPrompt')
        yToRetrackPrompt = true;
    end
end

liveExperiment = LiveExperiment(Prefix);

ProcPath = liveExperiment.userProcFolder;
DropboxFolder = liveExperiment.userResultsFolder;

Channel1 = liveExperiment.Channel1;
Channel2 = liveExperiment.Channel2;
Channel3 = liveExperiment.Channel3;
anaphaseFrames = liveExperiment.anaphaseFrames;
nc9 = anaphaseFrames(1);
nc10 = anaphaseFrames(2);
nc11 = anaphaseFrames(3);
nc12 = anaphaseFrames(4);
nc13 = anaphaseFrames(5);
nc14 = anaphaseFrames(6);

xSize = liveExperiment.xDim;
ySize = liveExperiment.yDim;
pixelSize_um = liveExperiment.pixelSize_um;

%Get the nuclei segmentation data
Ellipses = getEllipses(liveExperiment);
schnitzcells = getSchnitzcells(liveExperiment);
%Load the reference histogram for the fake histone channel
load('ReferenceHist.mat', 'ReferenceHist')


schnitzcellsFile = [liveExperiment.resultsFolder, filesep, Prefix, '_lin.mat'];

[Ellipses, schnitzcells] = addSchnitzIndexToEllipses(Ellipses, schnitzcells);
save2([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'], Ellipses);
save2(schnitzcellsFile, schnitzcells)

Channels = {Channel1, Channel2, Channel3};

if chooseHis
    uiopen([ProcPath, filesep, Prefix,'_',filesep,'*.mat']);
    if exist('probHis_fiji', 'var')
        hisMat = probHis_fiji;
        clear probHis_fiji;
    elseif exist('probHis_matlab', 'var')
        hisMat = probHis_matlab;
        clear probHis_matlab;
    elseif exist('probHis', 'var')
        hisMat = probHis;
        clear probHis;
    end
else
    hisMat = getHisMat(liveExperiment);
end

nFrames = size(hisMat, 3);
%Get information about the image size
% HisImage=imread([PreProcPath,filesep,Prefix,filesep,D(1).name]);
HisImage = hisMat(:,:,1);
DisplayRange=[min(min(HisImage)),max(max(HisImage))];

nc = [];

%Make a vector containing the nc corresponding to each frame
for k=1:nFrames
    if k<nc9
        nc(k)=8;
    elseif (k>=nc9)&(k<nc10)
        nc(k)=9;
    elseif (k>=nc10)&(k<nc11)
        nc(k)=10;
    elseif (k>=nc11)&(k<nc12)
        nc(k)=11;
    elseif (k>=nc12)&(k<nc13)
        nc(k)=12;
    elseif (k>=nc13)&(k<nc14) %#ok<*AND2>
        nc(k)=13;
    elseif k>=nc14
        nc(k)=14;
    end
end



%%
Overlay=figure;
set(Overlay,'units', 'normalized', 'position',[0.01, .2, .5, .5]);

overlayAxes = axes(Overlay,'Units', 'normalized', 'Position', [0 0 1 1]);

OriginalImage=figure;
set(OriginalImage,'units', 'normalized', 'position',[0.55, .2, .5, .5]);

originalAxes = axes(OriginalImage,'Units', 'normalized', 'Position', [0 0 1 1]);

tb = axtoolbar(overlayAxes);
tb.Visible = 'off';
tb2 = axtoolbar(originalAxes);
tb2.Visible = 'off';

try
    clrmp = single(hsv(length(schnitzcells)));
    clrmp = clrmp(randperm(length(clrmp)), :);
catch
%in case the user doesn't have this colormap, just keep going.
end

CurrentFrame=1;
cc=1;

% Show the first image
imOverlay = imshow(HisImage,DisplayRange,'Border','Tight','Parent',overlayAxes);
imOriginal = imshow(HisImage,DisplayRange,'Border','Tight','Parent',originalAxes);

projFlag = false;
set(0, 'CurrentFigure', Overlay)

while (cc~='x')
    
    %Load subsequent images
    if ~projFlag
        HisImage = hisMat(:, :, CurrentFrame);
    else
        HisImage = Projection(:, :,CurrentFrame);
    end
    
    
    %Get the information about the centroids
    [NCentroids,~]=size(Ellipses{CurrentFrame});
    
        
    imOverlay.CData = HisImage;
    try
        caxis(overlayAxes, DisplayRange);
        caxis(originalAxes, DisplayRange);
    end
    
    %refresh ellipses plots by destroying and remaking
    if exist('PlotHandle', 'var')
        cellfun(@delete, PlotHandle);
    end
    
    PlotHandle = cell(NCentroids, 1);
    ellipseFrame = Ellipses{CurrentFrame};
    if ~fish
        for k=1:NCentroids
            n = k;
            PlotHandle{k} = drawellipse('Center',[ellipseFrame(n, 1) ellipseFrame(n, 2)],...
                'SemiAxes',[ellipseFrame(n, 3) ellipseFrame(n, 4)], ...
                'RotationAngle',ellipseFrame(n, 5) * (360/(2*pi)), 'FaceAlpha', 0,...
                'InteractionsAllowed', 'none');
            
            if size(Ellipses{CurrentFrame}, 2) > 8
                schnitzInd = Ellipses{CurrentFrame}(k, 9);
            else
                schnitzInd = getSchnitz(Ellipses{CurrentFrame}(k,:), schnitzcells, CurrentFrame);
                if ~isempty(schnitzInd)
                    Ellipses{CurrentFrame}(k, 9) = schnitzInd;
                else
                    Ellipses{CurrentFrame}(k, 9) = 0;
                end
            end
            
            if schnitzInd ~= 0
                set(PlotHandle{k}, 'StripeColor', clrmp(schnitzInd, :),...
                    'Color', clrmp(schnitzInd, :),'Linewidth', 2);
            else
                set(PlotHandle{k}, 'StripeColor', 'w', 'Color', 'w','Linewidth', 2);
            end
            
        end
    else
        for k=1:NCentroids
            n = k;
            
            PlotHandle{k} = drawellipse('Center',[ellipseFrame(n, 1) ellipseFrame(n, 2)],...
                'SemiAxes',[ellipseFrame(n, 3) ellipseFrame(n, 4)], ...
                'RotationAngle',ellipseFrame(n, 5) * (360/(2*pi)), 'FaceAlpha', 0,...
                'InteractionsAllowed', 'none');
            
        end
    end
    
    try
    FigureTitle=['Frame: ',num2str(CurrentFrame),'/',num2str(nFrames),...
        ', nc: ',num2str(nc(CurrentFrame))];
    catch
          FigureTitle=['Frame: ',num2str(CurrentFrame),'/',num2str(nFrames)];
    end
    
    set(Overlay,'Name',FigureTitle)
    
    
    imOriginal.CData = HisImage;
    
    
    tb = axtoolbar(overlayAxes);
    tb.Visible = 'off';
    tb2 = axtoolbar(originalAxes);
    tb2.Visible = 'off';
    
    ct=waitforbuttonpress;
    cc=get(Overlay,'currentcharacter');
    cm=get(overlayAxes,'CurrentPoint');
    
    
    
    
    if (ct~=0)&(cc=='.')&(CurrentFrame<nFrames)
        CurrentFrame=CurrentFrame+1;
    elseif (ct~=0)&(cc==',')&(CurrentFrame>1)
        CurrentFrame=CurrentFrame-1;
    elseif (ct~=0)&(cc=='>')&(CurrentFrame+5<nFrames)
        CurrentFrame=CurrentFrame+5;
    elseif (ct~=0)&(cc=='<')&(CurrentFrame-4>1)
        CurrentFrame=CurrentFrame-5;
    elseif (ct~=0)&(cc=='s')
        save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses', '-v6')
        disp('Ellipses saved.')
    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'normal'))
        cc=1;
        if (cm(1,2)>0)&(cm(1,1)>0)&(cm(1,2)<=ySize)&(cm(1,1)<=xSize)
            
            %Add a circle to this location with the mean radius of the
            %ellipses found in this frame
            
            %(x, y, a, b, theta, maxcontourvalue, time,
            %particle_id)
            
            MeanRadius = computeMeanRadius(Ellipses, CurrentFrame, nFrames);
           
            try
                Ellipses{CurrentFrame}(end+1,:)=...
                    [cm(1,1),cm(1,2),MeanRadius,MeanRadius,0,0,0,0,0];
            catch
                Ellipses{CurrentFrame}(end+1,:)=...
                    [cm(1,1),cm(1,2),MeanRadius,MeanRadius,0,0,0,0];
            end
        end
        
        
        
        
    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'alt'))
        cc=1;
        if (cm(1,2)>0)&(cm(1,1)>0)&(cm(1,2)<=ySize)&(cm(1,1)<=xSize)
            %Find out which ellipses we clicked on so we can delete it
            
            %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
            Distances=sqrt((Ellipses{CurrentFrame}(:,1)-cm(1,1)).^2+...
                (Ellipses{CurrentFrame}(:,2)-cm(1,2)).^2);
            [~,MinIndex]=min(Distances);
            
            Ellipses{CurrentFrame}=[Ellipses{CurrentFrame}(1:MinIndex-1,:);...
                Ellipses{CurrentFrame}(MinIndex+1:end,:)];
        end
        
    elseif (ct~=0)&(cc=='j')
        iJump=input('Frame to jump to: ');
        if (floor(iJump)>0)&(iJump<=nFrames)
            CurrentFrame=iJump;
        end
        
    elseif (ct~=0)&(cc=='m')    %Increase contrast
        DisplayRange(2)=DisplayRange(2)/1.5;
        
    elseif (ct~=0)&(cc=='n')    %Decrease contrast
        DisplayRange(2)=DisplayRange(2)*1.5;
        
    elseif (ct~=0)&(cc=='r')    %Reset the contrast
        DisplayRange=[min(min(HisImage)),max(max(HisImage))];
        
    elseif (ct~=0)&(cc=='d')    %Delete all ellipses in the current frame
        Ellipses{CurrentFrame}=[];
    elseif (ct~=0)&(cc=='D')    %Delete all ellipses in hand-drawn ROI
        roi = drawrectangle(overlayAxes);
        EllipsesCopy = Ellipses;
        EllipsesCopy{CurrentFrame} = [];
        for c = 1:NCentroids
            r = [Ellipses{CurrentFrame}(c, 1), Ellipses{CurrentFrame}(c, 2)];
            if ~inROI(roi, r(1), r(2))
                EllipsesCopy{CurrentFrame}(c, :) = Ellipses{CurrentFrame}(c, :);
            end
        end
        Ellipses = EllipsesCopy;
        delete(roi);
        clear EllipsesCopy;
        
    elseif (ct~=0)&(cc=='c') & CurrentFrame > 1
        %copy nuclear information from previous frame
        
        Ellipses{CurrentFrame} = Ellipses{CurrentFrame-1};
        Ellipses{CurrentFrame} =...
            registerEllipses(Ellipses{CurrentFrame},...
            HisImage, hisMat(:, :, CurrentFrame-1));
        
    elseif (ct~=0)&(cc=='v') & CurrentFrame < nFrames
        %copy nuclear information from next frame
        
        Ellipses{CurrentFrame} = Ellipses{CurrentFrame+1};
        Ellipses{CurrentFrame} =...
            registerEllipses(Ellipses{CurrentFrame},...
            HisImage, hisMat(:, :, CurrentFrame+1));
        
        
    elseif (ct~=0)&(cc=='{')
        %resegment from scratch
        
        Ellipses{CurrentFrame}=[];
        [centers, radii, ~] =...
            findEllipsesByKMeans(HisImage, 'displayFigures', false);
        
        for k = 1:length(radii)
            Ellipses{CurrentFrame}(k, :) = [centers(k,1),centers(k,2),radii(k),radii(k),...
                0,0,0,0];
        end
        
    elseif (ct~=0)&(cc=='~')
        
        ProjectionType = 'midsumprojection';
        
        [~, ~, Projection] = chooseNuclearChannels2(...
            movieMat, 'ProjectionType', ProjectionType,'Channels',...
            Channels,'ReferenceHist', ReferenceHist);
        
        DisplayRange = [mean(mean(squeeze(Projection(:, :, CurrentFrame)))),...
            max(max(squeeze(Projection(:, :, CurrentFrame)))) ];
        
        disp('changed projection');
        
    elseif (ct~=0)&(cc=='g')  %copy nuclear information from next frame
        mitDuration = 10; % ~10 frames before and after anaphase
        for frame = CurrentFrame - mitDuration:CurrentFrame
            Ellipses{frame} = Ellipses{CurrentFrame-mitDuration-1};
        end
        for frame = CurrentFrame + 1:CurrentFrame + mitDuration
            Ellipses{frame} = Ellipses{CurrentFrame+mitDuration+1};
        end
    elseif (ct~=0)&(cc=='q') %go to next nc
        nextncframes = find(nc == (nc(CurrentFrame)+1));
        if ~isempty(nextncframes)
            CurrentFrame = nextncframes(1);
        end
    elseif (ct~=0)&(cc=='w') %go to previous nc
        previousncframes = find(nc == (nc(CurrentFrame)-1));
        if ~isempty(previousncframes)
            CurrentFrame = previousncframes(1);
        end
    elseif (ct~=0)&(cc=='\')  %resegment with ksnakecircles
        
        [~, circles] = kSnakeCircles(HisImage, pixelSize_um);
        circles(:, 6:9) = zeros(size(circles, 1), 4);
        Ellipses{CurrentFrame} = circles;
        
    elseif (ct~=0)&(cc=='0')    %Debug mode
        keyboard
        
    end
end

save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses', '-v6')


%Decide whether we need to re-track
if yToRetrackPrompt
    reTrackAnswer = 'y';
else
    userPrompt = 'Did you make changes to nuclei and thus require re-tracking? (y/n)';
    reTrackAnswer = inputdlg(userPrompt);
end

if contains(reTrackAnswer,'n')
    disp('Ellipses saved. Per user input, not re-tracking. Exiting.')
else
    opts = {};  if fish opts = [opts, 'markandfind']; end
    disp('Ellipses saved. Running TrackNuclei to incorporate changes.')
    TrackNuclei(Prefix,'retrack', 'nWorkers', 1, opts{:});
end

end


 function MeanRadius = computeMeanRadius(Ellipses, CurrentFrame, nFrames)
 
 radius = @(x,f) nanmean( (1/2)*(x{f}(:, 3) + x{f}(:, 4)) );
 
 if ~isempty(Ellipses{CurrentFrame})
     for k = 1:size(Ellipses{CurrentFrame}, 1)
         if Ellipses{CurrentFrame}(k, 3) == 0
             Ellipses{CurrentFrame}(k, :) = nan;
         end
     end
     MeanRadius = radius(Ellipses, CurrentFrame);
 elseif CurrentFrame+1 < nFrames && ~isempty(Ellipses{CurrentFrame+1})
     for k = 1:size(Ellipses{CurrentFrame+1}, 1)
         if Ellipses{CurrentFrame+1}(k, 3) == 0
             Ellipses{CurrentFrame+1}(k, :) = nan;
         end
     end
     MeanRadius = radius(Ellipses, CurrentFrame+1);
     
 elseif CurrentFrame-1 >1 && ~isempty(Ellipses{CurrentFrame-1})
     for k = 1:size(Ellipses{CurrentFrame-1}, 1)
         if Ellipses{CurrentFrame-1}(k, 3) == 0
             Ellipses{CurrentFrame-1}(k, :) = nan;
         end
     end
     MeanRadius = radius(Ellipses, CurrentFrame-1);
 else
     MeanRadius = 20; %magic number just to avoid errors in weird situations (units of pixels)
 end
                
end