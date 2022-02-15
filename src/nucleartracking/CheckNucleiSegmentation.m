function movieMat = CheckNucleiSegmentation(Prefix, varargin)
%%
% DESCRIPTION
% Opens a user interface that allows for manual curation of nuclear
% segmentation and tracking results
%
% PARAMETERS
% Prefix: Prefix of the dataset being analyzed
%
% OPTIONS
% 'yToRetrackPrompt', true/false: If followed by true, will automatically 
%                                 rerun TrackNuclei upon exiting the GUI. 
%                                 If followed by false, will open a user
%                                 dialog that asks whether or not you need
%                                 to rerun TrackNuclei. By default, the 
%                                 former happens.
% 'chooseHis': If you used Weka for classifying nuclei, this option will
%              plot the probHis.tif files instead of the raw His.tif images
%              for easier manual curation
% 'nWorkers': set the number of workers for a parallel pool (as of
%             2020-07-27, this option does nothing)
% 'noAdd', 'fish', or 'markandfind': changes some things for compatibility
%                                    with mark-and-find Leica data (e.g.
%                                    FISH experiments)
% 'drawTraces': no idea what this does
% 'premovie', movieMat: as of 2020-7-27, this option seems to do nothing
%
%
% GUI COMMANDS
% .  - Move a frame forward
% ,  - Move a frame backwards
% >  - Move 5 frames forward
% <  - Move 5 frames backwards
% j  - Jump to a frame
% q  - Move a cycle forward
% w  - Move a cycle backwards
% right click  - delete region
% left click - add region with default nc radius and angle
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
% OUTPUT
% Ellipses.mat: saved to the folder 'Dropbox\Prefix\'
%
%
% Author (contact): uknown (hggarcia@berkeley.edu)
% Created: XXXX-XX-XX
% Last Updated: XXXX-XX-XX
% Documented by: Meghan Turner (meghan_turner@berkeley.edu)
%

cleanupObj = onCleanup(@myCleanupFun);


noAdd = false;
nWorkers = 1;
fish = false;
preMovie = false;
chooseHis = false;
yToRetrackPrompt = false;
drawTraces = false;

for k = 1:length(varargin)
    if strcmpi(varargin{k}, 'noAdd') | strcmpi(varargin{k}, 'fish') | strcmpi(varargin{k}, 'markandfind')
        noAdd = true;
        fish = true;
    elseif strcmpi(varargin{k}, 'nWorkers')
        nWorkers = varargin{k+1};
    elseif strcmpi(varargin{k}, 'chooseHis')
        chooseHis = true;
%     elseif strcmpi(varargin{k}, 'colormap')
%         cmap = varargin{k+1};
    elseif strcmpi(varargin{k}, 'premovie')
        preMovie = true;
        movieMat = varargin{k+1};
    elseif strcmpi(varargin{k}, 'yToRetrackPrompt')
        yToRetrackPrompt = varargin{k+1};
    elseif strcmpi(varargin{k}, 'drawTraces')
        drawTraces = true;
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
PixelSize_um = liveExperiment.pixelSize_um;

%Get the nuclei segmentation data
Ellipses = getEllipses(liveExperiment);
schnitzcells = getSchnitzcells(liveExperiment);
%Load the reference histogram for the fake histone channel
load('ReferenceHist.mat', 'ReferenceHist')


%i may bring this back later -AR
if false
    schnitzcellsFile = [liveExperiment.resultsFolder, filesep, Prefix, '_lin.mat'];

    [Ellipses, schnitzcells] = addSchnitzIndexToEllipses(Ellipses, schnitzcells);
    save2([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'], Ellipses);
    save2(schnitzcellsFile, schnitzcells)
end

Channels = {Channel1, Channel2, Channel3};

if chooseHis
    hisMat = imreadStack([liveExperiment.procFolder, filesep, 'probHis.tif']);
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
    elseif (k>=nc9)&(k<nc10 || isnan(nc10))
        nc(k)=9;
    elseif (k>=nc10)&(k<nc11 || isnan(nc11))
        nc(k)=10;
    elseif (k>=nc11)&(k<nc12 || isnan(nc12))
        nc(k)=11;
    elseif (k>=nc12)& (k<nc13 || isnan(nc13))
        nc(k)=12;
    elseif (k>=nc13)&( k<nc14 || isnan(nc14) ) %#ok<*AND2>
        nc(k)=13;
    elseif k>=nc14
        nc(k)=14;
    end
end



%%
Overlay=figure;
set(Overlay,'units', 'normalized', 'position',[0.01, .5, .4, .4]);

overlayAxes = axes(Overlay,'Units', 'normalized', 'Position', [0 0 1 1]);

%%
OriginalImage=figure;
set(OriginalImage,'units', 'normalized', 'position',[0.01, 0.05, .4, .4]);
originalAxes = axes(OriginalImage,'Units', 'normalized', 'Position', [0 0 1 1]);
set(OriginalImage,'menubar','none')
set(OriginalImage,'NumberTitle','off');
%%
if drawTraces
    schnitzTrackingFigure = figure;
    t = tiledlayout(schnitzTrackingFigure, 1, 2);
    schnitzXTrackingAxes = nexttile(t);
    schnitzYTrackingAxes = nexttile(t);
    set(schnitzTrackingFigure,'units', 'normalized', 'position',[0.6, .2, .3, .5]);
    title(schnitzXTrackingAxes, 'X over time')
    xlabel(t, 'frame')
    title(schnitzYTrackingAxes, 'Y over time')
    ylabel(t, 'centroid (pixels)')
end
%%

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
currentCharacter=1;

% Show the first image
imOverlay = imshow(HisImage,DisplayRange,'Border','Tight','Parent',overlayAxes);
imOriginal = imshow(HisImage,DisplayRange,'Border','Tight','Parent',originalAxes);

projFlag = false;
set(0, 'CurrentFigure', Overlay)

while (currentCharacter~='x')
    
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
    for k=1:NCentroids
        n = k;
%         PlotHandle{k} = drawellipse('Center',[ellipseFrame(n, 1) ellipseFrame(n, 2)],...
%             'SemiAxes',[ellipseFrame(n, 3) ellipseFrame(n, 4)], ...
%             'RotationAngle',ellipseFrame(n, 5) * (360/(2*pi)), 'FaceAlpha', 0,...
%             'InteractionsAllowed', 'none', 'LabelVisible', 'hover', 'Label', num2str(ellipseFrame(n, 9)));
        
        PlotHandle{k} = ellipse(ellipseFrame(n, 3), ellipseFrame(n, 4),...
            ellipseFrame(n, 5) * (360/(2*pi)), ellipseFrame(n, 1),...
            ellipseFrame(n, 2), 'k', 10, overlayAxes);
        
        if ~fish
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
%                 set(PlotHandle{k}, 'StripeColor', clrmp(schnitzInd, :),...
%                     'Color', clrmp(schnitzInd, :),'Linewidth', 1);
                 set(PlotHandle{k},...
                    'Color', clrmp(schnitzInd, :),'Linewidth', 1);
            else
%                 set(PlotHandle{k}, 'StripeColor', 'w', 'Color', 'w','Linewidth', 1);
                set(PlotHandle{k}, 'Color', 'w','Linewidth', 1);
            end
        end
        
        if drawTraces
            %plot tracking information in the third figure
            plot(schnitzXTrackingAxes, ...
                schnitzcells(schnitzInd).frames, schnitzcells(schnitzInd).cenx,...
            'Color',  clrmp(schnitzInd, :), 'Linewidth', 3)
            hold(schnitzXTrackingAxes, 'on');
            plot(schnitzYTrackingAxes, ...
                schnitzcells(schnitzInd).frames, schnitzcells(schnitzInd).ceny,...
                'Color',  clrmp(schnitzInd, :), 'Linewidth', 3)
            hold(schnitzYTrackingAxes, 'on');
        end
        
        
    end
    
    if drawTraces
        hold(schnitzXTrackingAxes, 'off');
        hold(schnitzYTrackingAxes, 'off');
    end
    
    try
        FigureTitle=['Frame: ',num2str(CurrentFrame),'/',num2str(nFrames),...
            ', nc: ',num2str(nc(CurrentFrame))];
    catch
        FigureTitle=['Frame: ',num2str(CurrentFrame),'/',num2str(nFrames)];
    end
    
    set(Overlay,'Name',FigureTitle)
    
    
    imOriginal.CData = HisImage;
    
     
    
    
    
    %%
    
    tb = axtoolbar(overlayAxes);
    tb.Visible = 'off';
    tb2 = axtoolbar(originalAxes);
    tb2.Visible = 'off';
    
    ct=waitforbuttonpress;
    currentCharacter=get(Overlay,'currentcharacter');
    currentMouse=get(overlayAxes,'CurrentPoint');
    
    
    
    
    if (ct~=0)&(currentCharacter=='.')&(CurrentFrame<nFrames)
        CurrentFrame=CurrentFrame+1;
    elseif (ct~=0)&(currentCharacter==',')&(CurrentFrame>1)
        CurrentFrame=CurrentFrame-1;
    elseif (ct~=0)&(currentCharacter=='>')&(CurrentFrame+5<nFrames)
        CurrentFrame=CurrentFrame+5;
    elseif (ct~=0)&(currentCharacter=='<')&(CurrentFrame-4>1)
        CurrentFrame=CurrentFrame-5;
    elseif (ct~=0)&(currentCharacter=='s')
        save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses', '-v6')
        disp('Ellipses saved.')
    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'normal'))
        currentCharacter=1;
        if (currentMouse(1,2)>0)&(currentMouse(1,1)>0)&(currentMouse(1,2)<=ySize)&(currentMouse(1,1)<=xSize)
            
            %Add a circle to this location with the median radius of the
            %ellipses found in this frame
            
            %(x, y, a, b, theta, maxcontourvalue, time,
            %particle_id)
            
            averageRadius = computeAverageRadius(Ellipses, CurrentFrame, nFrames);
            
            try
                Ellipses{CurrentFrame}(end+1,:)=...
                    [currentMouse(1,1),currentMouse(1,2),averageRadius,averageRadius,0,0,0,0,0];
            catch
                Ellipses{CurrentFrame}(end+1,:)=...
                    [currentMouse(1,1),currentMouse(1,2),averageRadius,averageRadius,0,0,0,0];
            end
            
        end
        
        
        
        
    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'alt'))
        currentCharacter=1;
        if (currentMouse(1,2)>0)&(currentMouse(1,1)>0)&(currentMouse(1,2)<=ySize)&(currentMouse(1,1)<=xSize)
            %Find out which ellipses we clicked on so we can delete it
            
            %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
            Distances=sqrt((Ellipses{CurrentFrame}(:,1)-currentMouse(1,1)).^2+...
                (Ellipses{CurrentFrame}(:,2)-currentMouse(1,2)).^2);
            [~,MinIndex]=min(Distances);
            
            Ellipses{CurrentFrame}=[Ellipses{CurrentFrame}(1:MinIndex-1,:);...
                Ellipses{CurrentFrame}(MinIndex+1:end,:)];
        end
        
    elseif (ct~=0)&(currentCharacter=='j')
        iJump=input('Frame to jump to: ');
        if (floor(iJump)>0)&(iJump<=nFrames)
            CurrentFrame=iJump;
        else
            disp('Frame out of range.');
        end
        
    elseif (ct~=0)&(currentCharacter=='m')    %Increase contrast
        DisplayRange(2)=DisplayRange(2)/1.5;
        
    elseif (ct~=0)&(currentCharacter=='n')    %Decrease contrast
        DisplayRange(2)=DisplayRange(2)*1.5;
        
    elseif (ct~=0)&(currentCharacter=='r')    %Reset the contrast
        DisplayRange=[min(min(HisImage)),max(max(HisImage))];
        
    elseif (ct~=0)&(currentCharacter=='d')    %Delete all ellipses in the current frame
        Ellipses{CurrentFrame}=[];
    elseif (ct~=0)&(currentCharacter=='D')    %Delete all ellipses in hand-drawn ROI
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
        
    elseif (ct~=0)&(currentCharacter=='c') & CurrentFrame > 1
        %copy nuclear information from previous frame
        
        Ellipses{CurrentFrame} = Ellipses{CurrentFrame-1};
        Ellipses{CurrentFrame} =...
            registerEllipses(Ellipses{CurrentFrame},...
            HisImage, hisMat(:, :, CurrentFrame-1));
        
    elseif (ct~=0)&(currentCharacter=='v') & CurrentFrame < nFrames
        %copy nuclear information from next frame
        
        Ellipses{CurrentFrame} = Ellipses{CurrentFrame+1};
        Ellipses{CurrentFrame} =...
            registerEllipses(Ellipses{CurrentFrame},...
            HisImage, hisMat(:, :, CurrentFrame+1));
        
        
    elseif (ct~=0)&(currentCharacter=='{')
        %resegment from scratch
        
        Ellipses{CurrentFrame}=[];
        [centers, radii, ~] =...
            findEllipsesByKMeans(HisImage, 'displayFigures', false);
        
        for k = 1:length(radii)
            Ellipses{CurrentFrame}(k, :) = [centers(k,1),centers(k,2),radii(k),radii(k),...
                0,0,0,0];
        end
        
    elseif (ct~=0)&(currentCharacter=='~')
        
        ProjectionType = 'midsumprojection';
        
        movieMat = getMovieMat(liveExperiment); 
        [~, ~, Projection] = chooseNuclearChannels2(...
            movieMat, 'ProjectionType', ProjectionType,'Channels',...
            Channels,'ReferenceHist', ReferenceHist);
        
        DisplayRange = [mean(mean(Projection(:, :, CurrentFrame))),...
            max(max(Projection(:, :, CurrentFrame))) ];
        
        disp('changed projection');
        
    elseif (ct~=0)&(currentCharacter=='g')  %copy nuclear information from next frame
        mitDuration = 10; % ~10 frames before and after anaphase
        for frame = CurrentFrame - mitDuration:CurrentFrame
            Ellipses{frame} = Ellipses{CurrentFrame-mitDuration-1};
        end
        for frame = CurrentFrame + 1:CurrentFrame + mitDuration
            Ellipses{frame} = Ellipses{CurrentFrame+mitDuration+1};
        end
    elseif (ct~=0)&(currentCharacter=='q') %go to next nc
        nextncframes = find(nc == (nc(CurrentFrame)+1));
        if ~isempty(nextncframes)
            CurrentFrame = nextncframes(1);
        end
    elseif (ct~=0)&(currentCharacter=='w') %go to previous nc
        previousncframes = find(nc == (nc(CurrentFrame)-1));
        if ~isempty(previousncframes)
            CurrentFrame = previousncframes(1);
        end
    elseif (ct~=0)&(currentCharacter=='\')  %resegment with ksnakecircles
        
        [~, circles] = kSnakeCircles(HisImage, PixelSize_um);
        circles(:, 6:9) = zeros(size(circles, 1), 4);
        Ellipses{CurrentFrame} = circles;
        
    elseif (ct~=0)&(currentCharacter=='`')  %perform active contouring
        
    Ellipses{CurrentFrame} = adjustNuclearContours(...
        Ellipses{CurrentFrame}, HisImage, liveExperiment.pixelSize_um);


    elseif (ct~=0)&(currentCharacter=='0')    %Debug mode
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


function averageRadius = computeAverageRadius(Ellipses, CurrentFrame, nFrames)

radius = @(x,f) nanmedian( (1/2)*(x{f}(:, 3) + x{f}(:, 4)) );

if ~isempty(Ellipses{CurrentFrame})
    for k = 1:size(Ellipses{CurrentFrame}, 1)
        if Ellipses{CurrentFrame}(k, 3) == 0
            Ellipses{CurrentFrame}(k, :) = nan;
        end
    end
    
    averageRadius = radius(Ellipses, CurrentFrame);
    
elseif CurrentFrame+1 < nFrames && ~isempty(Ellipses{CurrentFrame+1})
    for k = 1:size(Ellipses{CurrentFrame+1}, 1)
        if Ellipses{CurrentFrame+1}(k, 3) == 0
            Ellipses{CurrentFrame+1}(k, :) = nan;
        end
    end
    
    averageRadius = radius(Ellipses, CurrentFrame+1);
    
elseif CurrentFrame-1 >1 && ~isempty(Ellipses{CurrentFrame-1})
    for k = 1:size(Ellipses{CurrentFrame-1}, 1)
        if Ellipses{CurrentFrame-1}(k, 3) == 0
            Ellipses{CurrentFrame-1}(k, :) = nan;
        end
    end
    
    averageRadius = radius(Ellipses, CurrentFrame-1);
    
else
    averageRadius = 20; %magic number just to avoid errors in weird situations (units of pixels)
end

end