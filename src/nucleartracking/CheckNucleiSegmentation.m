function movieMat = CheckNucleiSegmentation(Prefix, varargin)
%
%
%To do:
%1) Allow to edit the size and angle of an ellipse
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




close all

%Load the folder information
[SourcePath,ProcPath,DefaultDropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

noAdd = false;
nWorkers = 1;
fish = false;
preMovie = false;
chooseHis = false;
map = gray;

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'noAdd') | strcmpi(varargin{i}, 'fish') | strcmpi(varargin{i}, 'markandfind')
        noAdd = true;
        fish = true;
    elseif strcmpi(varargin{i}, 'nWorkers')
        nWorkers = varargin{i+1};
    elseif strcmpi(varargin{i}, 'chooseHis')
        chooseHis = varargin{i+1};
    elseif strcmpi(varargin{i}, 'colormap')
        map = varargin{i+1};  
    elseif strcmpi(varargin{i}, 'premovie')
        preMovie = true;
        movieMat = varargin{i+1};
    end
end

startParallelPool(nWorkers, 0, 1);


[~,ProcPath,DropboxFolder,~,PreProcPath]=...
    DetermineLocalFolders(Prefix);


%Set the source folders
Folder=[ProcPath,filesep,Prefix,'_',filesep,'preanalysis',filesep];
FileName=['CompactResults_',Prefix,'_.mat'];

%Set the destination folders
OutputFolder=[DropboxFolder,filesep,Prefix];
FilePrefix=FileName(16:end-4);
DataFolder=[Folder,'..',filesep,'..',filesep,'..',filesep,'Data',filesep,FilePrefix(1:end-1)];


%Find out how many frames we have
% D=dir([PreProcPath,filesep,Prefix,filesep,Prefix,'-His_*.tif']);



[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
    Channel1, Channel2, Objective, Power, DataFolderFromDataColumn, DropboxFolderName, Comments,...
    nc9, nc10, nc11, nc12, nc13, nc14, CF, Channel3] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);

load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo');

[xSize, ySize, pixelSize, ~, ~,...
    nFrames, ~, ~] = getFrameInfoParams(FrameInfo);

%Get the nuclei segmentation data
load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'], 'Ellipses');
load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'], 'schnitzcells');
%Load the reference histogram for the fake histone channel
load('ReferenceHist.mat', 'ReferenceHist')

hasSchnitzInd =size(Ellipses{1},2) == 9;

if ~hasSchnitzInd
    Ellipses = addSchnitzIndexToEllipses(Ellipses, schnitzcells);
    save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'], 'Ellipses');
end

Channels = {Channel1{1}, Channel2{1}, Channel3{1}};
nCh = sum(~cellfun(@isempty, Channels));

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
    [~,hisMat, maxMat, medMat, midMat]...
                = makeMovieMats(Prefix, PreProcPath, nWorkers, FrameInfo, 'loadMovie', false);
end

nFrames = size(hisMat, 3);
%Get information about the image size
% HisImage=imread([PreProcPath,filesep,Prefix,filesep,D(1).name]);
HisImage = squeeze(hisMat(:,:,1));
DisplayRange=[min(min(HisImage)),max(max(HisImage))];



%Make a vector containing the nc corresponding to each frame
for i=1:nFrames
    if i<nc9
        nc(i)=8;
    elseif (i>=nc9)&(i<nc10)
        nc(i)=9;
    elseif (i>=nc10)&(i<nc11)
        nc(i)=10;
    elseif (i>=nc11)&(i<nc12)
        nc(i)=11;
    elseif (i>=nc12)&(i<nc13)
        nc(i)=12;
    elseif (i>=nc13)&(i<nc14) %#ok<*AND2>
        nc(i)=13;
    elseif i>=nc14
        nc(i)=14;
    end
end



%%


%%
Overlay=figure;
% set(Overlay,'units', 'normalized', 'position',[0.01, .55, .75, .33]);
set(Overlay,'units', 'normalized', 'position',[0.01, .2, .5, .5]);

overlayAxes = axes(Overlay,'Units', 'normalized', 'Position', [0 0 1 1]);

OriginalImage=figure;
% set(OriginalImage,'units', 'normalized', 'position',[0.01, .1, .75, .33]);
set(OriginalImage,'units', 'normalized', 'position',[0.55, .2, .5, .5]);

originalAxes = axes(OriginalImage,'Units', 'normalized', 'Position', [0 0 1 1]);

tb = axtoolbar(overlayAxes);
tb.Visible = 'off';
tb2 = axtoolbar(originalAxes);
tb2.Visible = 'off';

try
    clrmp = single(hsv(length(schnitzcells)));
    clrmp = clrmp(randperm(length(clrmp)), :);
end

CurrentFrame=1;
cc=1;

% Show the first image
imOverlay = imshow(HisImage,DisplayRange,'Border','Tight','Parent',overlayAxes);
% colormap(overlayAxes,map);
imOriginal = imshow(HisImage,DisplayRange,'Border','Tight','Parent',originalAxes);
% set(overlayAxes,'Units', 'normalized', 'Position', [0 0 1 1]);
% % set(originalAxes,'Units', 'normalized', 'Position', [0 0 1 1]);
% imOverlay = imagescUpdate(overlayAxes, HisImage, []);
% set(overlayAxes,'Units', 'normalized', 'Position', [0 0 1 1]);

projFlag = false;
set(0, 'CurrentFigure', Overlay)

while (cc~='x')
    
    %Load subsequent images
    if ~projFlag 

            HisImage = squeeze(hisMat(:, :, CurrentFrame));
    else
        HisImage = squeeze(Projection(:, :,CurrentFrame));
    end
    
    
    %Get the information about the centroids
    [NCentroids,~]=size(Ellipses{CurrentFrame});
    
    %     imshow(HisImage,DisplayRange,'Border','Tight','Parent',overlayAxes)
    imOverlay.CData = HisImage;
    try
        caxis(overlayAxes, DisplayRange);
        caxis(originalAxes, DisplayRange);
    end
    axesHandlesToChildObjects = findobj(overlayAxes, 'Type', 'line');
    if ~isempty(axesHandlesToChildObjects)
        delete(axesHandlesToChildObjects);
    end
    %     hold(overlayAxes, 'on')
    PlotHandle=zeros(NCentroids, 1);
    if ~fish
        for i=1:NCentroids
            PlotHandle(i)=ellipse(Ellipses{CurrentFrame}(i,3),...
                Ellipses{CurrentFrame}(i,4),...
                Ellipses{CurrentFrame}(i,5),Ellipses{CurrentFrame}(i,1)+1,...
                Ellipses{CurrentFrame}(i,2)+1,[],20,overlayAxes);
            if size(Ellipses{CurrentFrame}, 2) > 8
                schnitzInd = Ellipses{CurrentFrame}(i, 9);
            else
                schnitzInd = getSchnitz(Ellipses{CurrentFrame}(i,:), schnitzcells, CurrentFrame);
                if ~isempty(schnitzInd)
                    Ellipses{CurrentFrame}(i, 9) = schnitzInd;
                else
                    Ellipses{CurrentFrame}(i, 9) = 0;
                end
            end
            if schnitzInd ~=0
                set(PlotHandle(i), 'Color', clrmp(schnitzInd, :),'Linewidth', 2);
            else
                set(PlotHandle(i), 'Color', 'w','Linewidth', 3);
                new_handle = copyobj(PlotHandle(i),overlayAxes);
                set(new_handle, 'Color', 'k','Linewidth', 2);
            end
        end
    else
        for i=1:NCentroids
            PlotHandle(i)=ellipse(Ellipses{CurrentFrame}(i,3),...
                Ellipses{CurrentFrame}(i,4),...
                Ellipses{CurrentFrame}(i,5),Ellipses{CurrentFrame}(i,1)+1,...
                Ellipses{CurrentFrame}(i,2)+1, 'g', 4,overlayAxes, .05);
            %              set(PlotHandle(i), 'Color', 'g','Linewidth', .5);
        end
        %         for i=1:NCentroids
        %             set(PlotHandle(i), 'Color', 'w','Linewidth', .5);
        %         end
    end
    %     hold(overlayAxes, 'off')
    %     set(PlotHandle,'Color','r', 'Linewidth', 3)
    
    
    
    FigureTitle=['Frame: ',num2str(CurrentFrame),'/',num2str(nFrames),...
        ', nc: ',num2str(nc(CurrentFrame))];
    set(Overlay,'Name',FigureTitle)
    
    
    %     imshow(HisImage,DisplayRange,'Border','Tight''Parent',originalAxes)
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
        save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses')
        disp('Ellipses saved.')
    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'normal'))
        cc=1;
        if (cm(1,2)>0)&(cm(1,1)>0)&(cm(1,2)<=ySize)&(cm(1,1)<=xSize)
            
            %Add a circle to this location with the mean radius of the
            %ellipses found in this frame
            
            %(x, y, a, b, theta, maxcontourvalue, time,
            %particle_id)
            if ~isempty(Ellipses{CurrentFrame})
                MeanRadius=mean((Ellipses{CurrentFrame}(:,3)+Ellipses{CurrentFrame}(:,4))/2);
            elseif ~isempty(Ellipses{CurrentFrame+1})
                MeanRadius=mean((Ellipses{CurrentFrame+1}(:,3)+Ellipses{CurrentFrame+1}(:,4))/2);
            elseif ~isempty(Ellipses{CurrentFrame-1})
                MeanRadius=mean((Ellipses{CurrentFrame-1}(:,3)+Ellipses{CurrentFrame-1}(:,4))/2);
            end
            
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
    elseif (ct~=0)&(cc=='c') & CurrentFrame > 1 %copy nuclear information from previous frame
        Ellipses{CurrentFrame} = Ellipses{CurrentFrame-1};
    elseif (ct~=0)&(cc=='v') & CurrentFrame < nFrames %copy nuclear information from next frame
        Ellipses{CurrentFrame} = Ellipses{CurrentFrame+1};
    elseif (ct~=0)&(cc=='{') %resegment from scratch 
        
        Ellipses{CurrentFrame}=[];
%         [centers, radii, mask] = maskNuclei2(HisImage);
         [centers, radii, mask] = findEllipsesByKMeans(HisImage, 'displayFigures', false);

        for i = 1:length(radii)
            Ellipses{CurrentFrame}(i, :) = [centers(i,1),centers(i,2),radii(i),radii(i),0,0,0,0];
        end
        
    elseif (ct~=0)&(cc=='~')
 
         [ProjectionType, nonInverted, inverted] = makeNuclearProjection_CNT(nCh);
        disp('calculating projection...')
        nuclearMovie = nan(xSize, ySize, nFrames, nCh, 'double'); % ch z t x y
        %ch z t x y
        for ch = 1:nCh
            if inverted(ch)
                if strcmpi(ProjectionType, 'maxprojection')
                    nuclearMovie(:, :, :, ch) = imcomplement(maxMat(:, :, :, ch));
                elseif strcmpi(ProjectionType, 'medianprojection')
                    nuclearMovie(:, :, :, ch) = imcomplement(medMat(:, :, :, ch));
                elseif strcmpi(ProjectionType, 'midprojection')
                    nuclearMovie(:, :, :, ch) = imcomplement(midMat(:, :, :, ch));
                end
            end
            if nonInverted(ch) & ~inverted(ch)
                if strcmpi(ProjectionType, 'maxprojection')
                    nuclearMovie(:, :, :, ch) = maxMat(:, :, :, ch);
                elseif strcmpi(ProjectionType, 'medprojection')
                    nuclearMovie(:, :, :, ch) = medMat(:, :, :, ch);
                elseif strcmpi(ProjectionType, 'midprojection')
                    nuclearMovie(:, :, :, ch) = midMat(:, :, :, ch);
                end
            end
            % Use the reference histogram to scale the Projection (This part
            % might need some more optimization later-YJK)
            nuclearMovie(:, :, :, ch) = histeq(mat2gray(nuclearMovie(:, :, :, ch)), ReferenceHist);
        end
        
        % Get average of all Projections
        Projection = squeeze(nanmean(nuclearMovie, 1));
        projFlag = true;
        
        DisplayRange = [mean(mean(squeeze(Projection(:, :, CurrentFrame)))), max(max(squeeze(Projection(:, :, CurrentFrame)))) ];
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
    elseif (ct~=0)&(cc=='/')  %adjust ellipse centroids
       
        Ellipses{CurrentFrame} = adjustEllipseCentroidsFrame(Ellipses{CurrentFrame}, squeeze(hisMat(:, :, CurrentFrame)), 'pixelSize', pixelSize);

    elseif (ct~=0)&(cc=='0')    %Debug mode
        keyboard
        
    end
end



save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses', '-v7.3', '-nocompression')
close all;

%Decide whether we need to re-track
userPrompt = 'Did you make changes to nuclei and thus require re-tracking? (y/n)';
reTrackAnswer = inputdlg(userPrompt);
if contains(reTrackAnswer,'n')
    disp('Ellipses saved. Per user input, not re-tracking. Exiting.')
else
    opts = {};  if fish opts = [opts, 'markandfind']; end
    disp('Ellipses saved. Running TrackNuclei to incorporate changes.')
    TrackNuclei(Prefix,'NoBulkShift','ExpandedSpaceTolerance', 1.5, 'retrack', 'nWorkers', 1, opts{:});
end

end
