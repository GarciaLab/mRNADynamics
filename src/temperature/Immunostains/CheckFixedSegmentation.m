function CheckFixedSegmentation(Prefix, UseCustom, varargin)
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
% Documented by: Gabriella Martini (martini@berkeley.edu)
%

cleanupObj = onCleanup(@myCleanupFun);

%%

if ~exist('UseCustom', 'var')
    UseCustom = false;
end
liveExperiment= LiveExperiment(Prefix);


CompiledEmbryoPath = [liveExperiment.resultsFolder, filesep, 'CompiledEmbryos.Mat'];
load(CompiledEmbryoPath);

ProcPath = liveExperiment.userProcFolder;
DropboxFolder = liveExperiment.userResultsFolder;

Channel1 = liveExperiment.Channel1;
Channel2 = liveExperiment.Channel2;
Channel3 = liveExperiment.Channel3;

xSize = liveExperiment.xDim;
ySize = liveExperiment.yDim;
PixelSize_um = liveExperiment.pixelSize_um;

%Get the nuclei segmentation data
Ellipses = getEllipses(liveExperiment);
%Load the reference histogram for the fake histone channel
% load('ReferenceHist.mat', 'ReferenceHist')
try
    clrmp = single(hsv(20));
    clrmp = clrmp(randperm(length(clrmp)), :);
catch
    %in case the user doesn't have this colormap, just keep going.
end
Channels = {Channel1, Channel2, Channel3};

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
CurrentEmbryo = 1;
CurrentOctant = 1;
OctantXs = [1:xSize/4:xSize 1:xSize/4:xSize];
OctantYs = [1 1 1 1 1+ySize/2 1+ySize/2 1+ySize/2 1+ySize/2];
xRangeOctant = [OctantXs(CurrentOctant), OctantXs(CurrentOctant)+xSize/4-1];
yRangeOctant = [OctantYs(CurrentOctant), OctantYs(CurrentOctant)+ySize/2-1];

HisImage = hisMat(:,:,CurrentEmbryo);
HisImage2 = hisMat(:,:,CurrentEmbryo);
MemImage = memMat(:,:,CurrentEmbryo);


DisplayRangeHis = [min(min(HisImage)), max(max(HisImage))];
DisplayRangeHis2 = [min(min(HisImage2)), max(max(HisImage2))];
DisplayRangeMem = [min(min(MemImage)), max(max(MemImage))];

close all

FullFigure=figure(1);
set(FullFigure,'units', 'normalized', 'position',[0.01, .3, .9, .5]);

fullAxes = axes(FullFigure,'Units', 'normalized', 'Position', [0 0 1 1]);
Overlay=figure(2);
set(Overlay,'units', 'normalized', 'position',[0.01, .1, .45, .65]);

overlayAxes = axes(Overlay,'Units', 'normalized', 'Position', [0 0 1 1]);


Original=figure(3);
set(Original,'units', 'normalized', 'position',[0.51, .1, .45, .65]);

originalAxes = axes(Original,'Units', 'normalized', 'Position', [0 0 1 1]);




tb = axtoolbar(overlayAxes);
tb.Visible = 'off';
tb2 = axtoolbar(originalAxes);
tb2.Visible = 'off';
tb3 = axtoolbar(originalAxes);
tb3.Visible = 'off';

imFull = imshow(HisImage2,DisplayRangeHis2,'Border','Tight','Parent',fullAxes);
imOriginal = imshow(MemImage,DisplayRangeMem,'Border','Tight','Parent',originalAxes);

imOverlay = imshow(HisImage,DisplayRangeHis,'Border','Tight','Parent',overlayAxes);


SwitchImageType = false;
hold(overlayAxes,'on')

xlim(overlayAxes,xRangeOctant)
ylim(overlayAxes,yRangeOctant)
xlim(originalAxes,xRangeOctant)
ylim(originalAxes,yRangeOctant)
set(0, 'CurrentFigure', Original)
set(0, 'CurrentFigure', Overlay)
currentCharacter = 1;
while (currentCharacter~='x')
    
    %Load subsequent images
    
    
    %Get the information about the centroids
    [NCentroids,~]=size(Ellipses{CurrentEmbryo});
    
    imFull.CData = HisImage2;
    imOverlay.CData = HisImage;
    imOriginal.CData = MemImage;
    try
        caxis(overlayAxes, DisplayRange);
        
    end
    
    try
        caxis(originalAxes, DisplayRangeMem);
        
    end
    try
        caxis(fullAxes, DisplayRange);
        
    end
    %refresh ellipses plots by destroying and remaking
    if exist('PlotHandle', 'var')
        cellfun(@delete, PlotHandle);
    end
    
    PlotHandle = cell(NCentroids, 1);
    ellipseFrame = double(Ellipses{CurrentEmbryo});
    for k=1:NCentroids
        n = k;
        %         PlotHandle{k} = drawellipse('Center',[ellipseFrame(n, 1) ellipseFrame(n, 2)],...
        %             'SemiAxes',[ellipseFrame(n, 3) ellipseFrame(n, 4)], ...
        %             'RotationAngle',ellipseFrame(n, 5) * (360/(2*pi)), 'FaceAlpha', 0,...
        %             'InteractionsAllowed', 'none', 'LabelVisible', 'hover', 'Label', num2str(ellipseFrame(n, 9)));
        colorhash = uint8(mod(round(ellipseFrame(n, 1)+ellipseFrame(n, 2)),20)+1);
        PlotHandle{k} = ellipse(2*ellipseFrame(n, 3), 2*ellipseFrame(n, 4),...
            ellipseFrame(n, 5) * (360/(2*pi)), ellipseFrame(n, 1),...
            ellipseFrame(n, 2), clrmp(colorhash,:), 10, overlayAxes, 0.5);
        
    end
    
    
    
    
    
    
    
    
    FigureTitle=['Frame: ',num2str(CurrentEmbryo),'/',num2str(nEmbryos)];
    
    
    set(Overlay,'Name',FigureTitle)
    
    
    
    
    
    
    if (CompiledEmbryos.Approved(CurrentEmbryo))
        set(Overlay,'Color','g')
    else
        set(Overlay,'Color','r')
    end
    
    %%
    
    tb = axtoolbar(overlayAxes);
    tb.Visible = 'off';
    
    ct=waitforbuttonpress;
    currentCharacter=get(Overlay,'currentcharacter');
    currentMouse=get(overlayAxes,'CurrentPoint');
    
    
    
    
    if (ct~=0)&(currentCharacter=='.')&(CurrentOctant<8)
        CurrentOctant=CurrentOctant+1;
        xRangeOctant = [OctantXs(CurrentOctant), OctantXs(CurrentOctant)+xSize/4-1];
        yRangeOctant = [OctantYs(CurrentOctant), OctantYs(CurrentOctant)+ySize/2-1];
        
        xlim(overlayAxes,xRangeOctant)
        ylim(overlayAxes,yRangeOctant)
        xlim(originalAxes,xRangeOctant)
        ylim(originalAxes,yRangeOctant)
        HisImage = hisMat(:, :, CurrentEmbryo);
        DisplayRange = [min(min(HisImage)), max(max(HisImage))];
        HisImage2 = hisMatCopy(:, :, CurrentEmbryo);
        
        DisplayRange2 = [min(min(HisImage2)), max(max(HisImage2))];
        
        MemImage = memMat(:, :, CurrentEmbryo);
        DisplayRangeMem = [min(min(MemImage)), max(max(MemImage))];
    elseif (ct~=0)&(currentCharacter=='.')&(CurrentOctant==8)&(CurrentEmbryo < nEmbryos)
        NewIndex = find(1:nEmbryos > CurrentEmbryo & CompiledEmbryos.Approved, 1);
        if isempty(NewIndex)
            CurrentEmbryo = nEmbryos;
        else
            CurrentEmbryo = NewIndex;
        end
        CurrentOctant = 1;
        xRangeOctant = [OctantXs(CurrentOctant), OctantXs(CurrentOctant)+xSize/4-1];
        yRangeOctant = [OctantYs(CurrentOctant), OctantYs(CurrentOctant)+ySize/2-1];
        
        xlim(overlayAxes,xRangeOctant)
        ylim(overlayAxes,yRangeOctant)
        xlim(originalAxes,xRangeOctant)
        ylim(originalAxes,yRangeOctant)
        HisImage = hisMat(:, :, CurrentEmbryo);
        DisplayRange = [min(min(HisImage)), max(max(HisImage))];
        HisImage2 = hisMatCopy(:, :, CurrentEmbryo);
        DisplayRange2 = [min(min(HisImage2)), max(max(HisImage2))];
        MemImage = memMat(:, :, CurrentEmbryo);
        DisplayRangeMem = [min(min(MemImage)), max(max(MemImage))];
    elseif (ct~=0)&(currentCharacter==',')&(CurrentOctant>1)
        CurrentOctant=CurrentOctant-1;
        xRangeOctant = [OctantXs(CurrentOctant), OctantXs(CurrentOctant)+xSize/4-1];
        yRangeOctant = [OctantYs(CurrentOctant), OctantYs(CurrentOctant)+ySize/2-1];
        
        xlim(overlayAxes,xRangeOctant)
        ylim(overlayAxes,yRangeOctant)
        xlim(originalAxes,xRangeOctant)
        ylim(originalAxes,yRangeOctant)
    elseif (ct~=0)&(currentCharacter==',')&(CurrentOctant==1)&(CurrentEmbryo > 1)
        NewIndex = find(1:nEmbryos < CurrentEmbryo & CompiledEmbryos.Approved, 1,'last');
        if isempty(NewIndex)
            CurrentEmbryo = 1;
        else
            CurrentEmbryo = NewIndex;
        end
        CurrentOctant = 8;
        xRangeOctant = [OctantXs(CurrentOctant), OctantXs(CurrentOctant)+xSize/4-1];
        yRangeOctant = [OctantYs(CurrentOctant), OctantYs(CurrentOctant)+ySize/2-1];
        
        xlim(overlayAxes,xRangeOctant)
        ylim(overlayAxes,yRangeOctant)
        xlim(originalAxes,xRangeOctant)
        ylim(originalAxes,yRangeOctant)
        HisImage = hisMat(:, :, CurrentEmbryo);
        DisplayRange = [min(min(HisImage)), max(max(HisImage))];
        HisImage2 = hisMatCopy(:, :, CurrentEmbryo);
        DisplayRange2 = [min(min(HisImage2)), max(max(HisImage2))];
        MemImage = memMat(:, :, CurrentEmbryo);
        DisplayRangeMem = [min(min(MemImage)), max(max(MemImage))];
    elseif (ct~=0)&(currentCharacter=='>')&(CurrentEmbryo < nEmbryos)
        NewIndex = find(1:nEmbryos > CurrentEmbryo & CompiledEmbryos.Approved, 1);
        if isempty(NewIndex)
            CurrentEmbryo = nEmbryos;
        else
            CurrentEmbryo = NewIndex;
        end
        CurrentOctant = 1;
        xRangeOctant = [OctantXs(CurrentOctant), OctantXs(CurrentOctant)+xSize/4-1];
        yRangeOctant = [OctantYs(CurrentOctant), OctantYs(CurrentOctant)+ySize/2-1];
        
        xlim(overlayAxes,xRangeOctant)
        ylim(overlayAxes,yRangeOctant)
        xlim(originalAxes,xRangeOctant)
        ylim(originalAxes,yRangeOctant)
        HisImage = hisMat(:, :, CurrentEmbryo);
        DisplayRange = [min(min(HisImage)), max(max(HisImage))];
        HisImage2 = hisMatCopy(:, :, CurrentEmbryo);
        DisplayRange2 = [min(min(HisImage2)), max(max(HisImage2))];
        MemImage = memMat(:, :, CurrentEmbryo);
        DisplayRangeMem = [min(min(MemImage)), max(max(MemImage))];
    elseif (ct~=0)&(currentCharacter=='<')&(CurrentEmbryo > 1)
        NewIndex = find(1:nEmbryos < CurrentEmbryo & CompiledEmbryos.Approved, 1,'last');
        if isempty(NewIndex)
            CurrentEmbryo = 1;
        else
            CurrentEmbryo = NewIndex;
        end
        CurrentOctant = 1;
        xRangeOctant = [OctantXs(CurrentOctant), OctantXs(CurrentOctant)+xSize/4-1];
        yRangeOctant = [OctantYs(CurrentOctant), OctantYs(CurrentOctant)+ySize/2-1];
        
        xlim(overlayAxes,xRangeOctant)
        ylim(overlayAxes,yRangeOctant)
        xlim(originalAxes,xRangeOctant)
        ylim(originalAxes,yRangeOctant)
        HisImage = hisMat(:, :, CurrentEmbryo);
        DisplayRange = [min(min(HisImage)), max(max(HisImage))];
        HisImage2 = hisMatCopy(:, :, CurrentEmbryo);
        DisplayRange2 = [min(min(HisImage2)), max(max(HisImage2))];
        MemImage = memMat(:, :, CurrentEmbryo);
        DisplayRangeMem = [min(min(MemImage)), max(max(MemImage))];
    elseif (ct~=0)&(currentCharacter=='j')
        iJump=input('Frame to jump to: ');
        if (floor(iJump)>0)&(iJump<=nEmbryos)
            CurrentEmbryo=iJump;
        else
            disp('Frame out of range.');
        end
        CurrentOctant = 1;
        xRangeOctant = [OctantXs(CurrentOctant), OctantXs(CurrentOctant)+xSize/4-1];
        yRangeOctant = [OctantYs(CurrentOctant), OctantYs(CurrentOctant)+ySize/2-1];
        
        xlim(overlayAxes,xRangeOctant)
        ylim(overlayAxes,yRangeOctant)
        xlim(originalAxes,xRangeOctant)
        ylim(originalAxes,yRangeOctant)
        HisImage = hisMat(:, :, CurrentEmbryo);
        DisplayRange = [min(min(HisImage)), max(max(HisImage))];
        HisImage2 = hisMatCopy(:, :, CurrentEmbryo);
        DisplayRange2 = [min(min(HisImage2)), max(max(HisImage2))];
        MemImage = memMat(:, :, CurrentEmbryo);
        DisplayRangeMem = [min(min(MemImage)), max(max(MemImage))];
    elseif (ct~=0)&(currentCharacter=='t')
        if ~SwitchImageType
            SwitchImageType = true;
            hisMat = memMat;
            HisImage = hisMat(:, :, CurrentEmbryo);
            DisplayRange = [min(min(HisImage)), max(max(HisImage))];
            HisImage2 = hisMatCopy(:, :, CurrentEmbryo);
            DisplayRange2 = [min(min(HisImage2)), max(max(HisImage2))];
            MemImage = memMat(:, :, CurrentEmbryo);
            DisplayRangeMem = [min(min(MemImage)), max(max(MemImage))];
        else
            SwitchImageType = false;
            hisMat = hisMatCopy;
            HisImage = hisMat(:, :, CurrentEmbryo);
            HisImage2 = hisMatCopy(:, :, CurrentEmbryo);
            DisplayRange2 = [min(min(HisImage2)), max(max(HisImage2))];
            DisplayRange = [min(min(HisImage)), max(max(HisImage))];
            MemImage = memMat(:, :, CurrentEmbryo);
            DisplayRangeMem = [min(min(MemImage)), max(max(MemImage))];
        end
    elseif (ct~=0)&(currentCharacter=='q')
        CompiledEmbryos.Approved(CurrentEmbryo) = true;
    elseif (ct~=0)&(currentCharacter=='w')
        CompiledEmbryos.Approved(CurrentEmbryo) = false;
    elseif (ct~=0)&(currentCharacter=='s')
        Ellipses = UpdateEllipsesOriginalCoordinates(Ellipses, Prefix, CompiledEmbryos,liveExperiment);
        save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses', '-v6')
        save([DropboxFolder,filesep,Prefix,filesep,'CompiledEmbryos.mat'],'CompiledEmbryos', '-v6')
        disp('Ellipses saved.')
    elseif (ct~=0)&(currentCharacter == 'f')
        try
            [flag, flag_string]  = chooseFixedEmbryoFlag;
            disp(flag_string)
            CompiledEmbryos.Flags(CurrentEmbryo) =flag;
        catch
            disp('No Flag Selected')
        end
    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'normal'))
        currentCharacter=1;
        if (currentMouse(1,2)>0)&(currentMouse(1,1)>0)&(currentMouse(1,2)<=ySize)&(currentMouse(1,1)<=xSize)
            
            %Add a circle to this location with the median radius of the
            %ellipses found in this frame
            
            %(x, y, a, b, theta, maxcontourvalue, time,
            %particle_id)
            
            averageRadius = max(Ellipses{CurrentEmbryo}(:,3));
            
            
            
            Ellipses{CurrentEmbryo}(end+1,:)=...
                [currentMouse(1,1),currentMouse(1,2),averageRadius,averageRadius,0,0,0,0];
            
            
        end
        
        
        
        
    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'alt'))
        currentCharacter=1;
        if (currentMouse(1,2)>0)&(currentMouse(1,1)>0)&(currentMouse(1,2)<=ySize)&(currentMouse(1,1)<=xSize)
            %Find out which ellipses we clicked on so we can delete it
            
            %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
            Distances=sqrt((Ellipses{CurrentEmbryo}(:,1)-currentMouse(1,1)).^2+...
                (Ellipses{CurrentEmbryo}(:,2)-currentMouse(1,2)).^2);
            [~,MinIndex]=min(Distances);
            
            Ellipses{CurrentEmbryo}=[Ellipses{CurrentEmbryo}(1:MinIndex-1,:);...
                Ellipses{CurrentEmbryo}(MinIndex+1:end,:)];
        end
        
        
        
    elseif (ct~=0)&(currentCharacter=='m')    %Increase contrast
        DisplayRange(2)=DisplayRange(2)/1.5;
        
    elseif (ct~=0)&(currentCharacter=='n')    %Decrease contrast
        DisplayRange(2)=DisplayRange(2)*1.5;
        
    elseif (ct~=0)&(currentCharacter=='r')    %Reset the contrast
        DisplayRange=[min(min(HisImage)),max(max(HisImage))];
        
        %     elseif (ct~=0)&(currentCharacter=='d')    %Delete all ellipses in the current frame
        %         Ellipses{CurrentFrame}=[];
    elseif (ct~=0)&(currentCharacter=='D')    %Delete all ellipses in hand-drawn ROI
        roi = drawrectangle(overlayAxes);
        EllipsesCopy = Ellipses;
        EllipsesCopy{CurrentEmbryo} = [];
        for c = 1:NCentroids
            r = [Ellipses{CurrentEmbryo}(c, 1), Ellipses{CurrentEmbryo}(c, 2)];
            if ~inROI(roi, r(1), r(2))
                EllipsesCopy{CurrentEmbryo}(c, :) = Ellipses{CurrentEmbryo}(c, :);
            end
        end
        Ellipses = EllipsesCopy;
        delete(roi);
        clear EllipsesCopy;
        
        
        
        
        
        
    elseif (ct~=0)&(currentCharacter=='0')    %Debug mode
        keyboard
        
    end
end
close all
Ellipses = UpdateEllipsesOriginalCoordinates(Ellipses,  Prefix, CompiledEmbryos,liveExperiment);
save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses', '-v6')
save([DropboxFolder,filesep,Prefix,filesep,'CompiledEmbryos.mat'],'CompiledEmbryos', '-v6')
