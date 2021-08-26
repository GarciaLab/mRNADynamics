function [Particles, SpotFilter] = AddParticlePosition(Prefix, varargin)
%
% DESCRIPTION
% Locates particles from a zoomed-in movie within full embryo images using
% spatial cross-correlation.
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
%
% OPTIONS
% 'SkipAlignment': If you want to skip alignment
% 'ManualAlignment': If you want to manually align the zoomed in and full
%                       embyro images
% 'NoAP': Just adds X and Y information
% 'SelectChannel': Prompts user to select the channel to use for alignment
% 'optionalResults':
% 'yToManualAlignmentPrompt':
% 'correctDV': If you want the DV shift calculated
%
% MANUAL ALIGNMENT CONTROLS
% . - Move to the right
% > - Move to the right further
% , - Move to the left
% < - Move to the left further
% a - Move up
% A - Move up further
% z - Move down
% Z - Move down further
% x - Save and exit
%
% V2: Changed this function to use a correlation in order to center the images.

% Default set of variables to save
saveVars={'coordA', 'coordP', 'coordAZoom', 'coordPZoom'};

[SkipAlignment, ManualAlignment, NoAP, SelectChannel, InvertHis,...
    optionalResults, yToManualAlignmentPrompt, correctDV] = determineAddParticlePositionOptions(varargin);


liveExperiment = LiveExperiment(Prefix);

DropboxFolder = liveExperiment.userResultsFolder;
PreProcPath = liveExperiment.userPreFolder;
RawDynamicsPath = liveExperiment.userRawFolder;

Channel1 = liveExperiment.Channel1;
Channel2 = liveExperiment.Channel2;
Channel3 = liveExperiment.Channel3;

APResolution = liveExperiment.APResolution;

[Particles, Spots, SpotFilter, NChannels] = loadParticlesIfExists(DropboxFolder, Prefix);

% See if we had any lineage/nuclear information
hisDir=dir([PreProcPath,filesep,Prefix,filesep,'*his*']);
if ~isempty(hisDir)
    histoneChannelPresent = true;
else
    histoneChannelPresent = false;
end

% isn't getMicroscope the same as DetermineFileMode?
[FileMode, EmbryoName, projectDate] = getMicroscope(Prefix);


if ~NoAP
    rawPrefixPath = [RawDynamicsPath,filesep,projectDate,filesep,EmbryoName,filesep];
    fullEmbryoPath = [rawPrefixPath, 'FullEmbryo', filesep];

    %If you want to select which channel to load as alignment
    ChannelToLoad = determineChannelToLoad(SelectChannel, liveExperiment.Channels);
    
    %Get information about all images. This depends on the microscope used.
    
    %Get the information about the zoom
    if strcmp(FileMode,'TIF')
        [ZoomImage, SurfImage, ZoomRatio, SurfInfo, SurfColumns, SurfRows, Rows, Columns] =...
            APP_getTifZoomInfo(DropboxFolder, Prefix, rawPrefixPath, fullEmbryoPath, ChannelToLoad);
    else
        %Figure out the different channels
        nuclearChannels = contains({Channel1,Channel2,Channel3},'nuclear','IgnoreCase',true);
        nuclearChannelIndices = find(nuclearChannels);
        invertedChannels = contains({Channel1,Channel2,Channel3},'inverted','IgnoreCase',true);
        if ~any(nuclearChannels)
            error('You should check the MovieDatabase.csv to see whether ":Nuclear" or "invertedNuclear" is plugged into your channels.')
        end
        
        %Find the zoomed movie pixel size
        D=dir([rawPrefixPath,'*.',FileMode(1:3)]);
        
        if strcmp(FileMode, 'DSPIN')            %CS20170912
            D=dir([rawPrefixPath,'*.nd']);
        end
        
        %Load only the metadata from the zoomed images
        MetaReader=bfGetReader([rawPrefixPath,D(end).name]);
        MetaZoom=MetaReader.getMetadataStore();
        try
            PixelSizeZoom=str2double(MetaZoom.getPixelsPhysicalSizeX(0).value);
        catch
            PixelSizeZoom=str2double(MetaZoom.getPixelsPhysicalSizeX(0));
        end
        
        %Find the full embryo pixel size and load the image
        D=dir([fullEmbryoPath,filesep,'*surf*.',FileMode(1:3)]);
        if strcmp(FileMode, 'DSPIN')            %CS20170912
            D=dir([RawDynamicsPath,filesep,projectDate,filesep,...
                EmbryoName, filesep,'FullEmbryo',filesep,'*surf*']);
        end
        
        surfFile = [fullEmbryoPath,filesep,D(end).name];
        ImageTempRaw=bfopen(surfFile);
        MetaFullEmbryo= ImageTempRaw{:, 4};
        PixelSizeFullEmbryoSurf=str2double(MetaFullEmbryo.getPixelsPhysicalSizeX(0) );
        try
            PixelSizeFullEmbryoSurf=str2double(MetaFullEmbryo.getPixelsPhysicalSizeX(0).value);
        catch
            PixelSizeFullEmbryoSurf=str2double(MetaFullEmbryo.getPixelsPhysicalSizeX(0));
        end
        
        %Check that the surface and midsaggital images have the same zoom
        midDirectory=dir([fullEmbryoPath,filesep,'*mid*.',FileMode(1:3)]);

        ImageTempMid=bfopen([fullEmbryoPath,filesep,midDirectory(end).name]);
        MetaFullEmbryoMid= ImageTempMid{:, 4};
        
        %This if for BioFormats backwards compatibility
        % Note 7/27/20: This is returning NaN so the check below doesn't
        % actually do anything. 
        if ~isempty(str2double(MetaFullEmbryoMid.getPixelsPhysicalSizeX(0)))
            PixelSizeFullEmbryoMid=str2double(MetaFullEmbryoMid.getPixelsPhysicalSizeX(0).value);
        else
            PixelSizeFullEmbryoMid=str2double(MetaFullEmbryoMid.getPixelsPhysicalSizeX(0));
        end
        
        
        %In principle, we would be comparing PixelSizeFullEmbryo==PixelSizeFullEmbryoMid
        %However, some issues of machine precision made this not work
        %sometimes.
        if abs(PixelSizeFullEmbryoSurf/PixelSizeFullEmbryoMid-1)>0.01
            error('The surface and midsaggital images were not taken with the same pixel size')
        end
        
        
        %How many channels and slices do we have?
        NChannelsMeta=MetaFullEmbryo.getChannelCount(0);

        %Look for the image with the largest size. In this way, we avoid
        %loading individual tiles in the case of a tile scan.
        
        ImageTemp = ImageTempRaw{1,1};

        for i=1:size(ImageTemp,1)
            ImageSizes(i)=size(ImageTemp{i,1},1);
        end
        
        [MaxSize,~]=max(ImageSizes);
        ElligibleImages = find(ImageSizes==MaxSize); %NL 2021-08-07: we want to use all images that are elligible 
        
        projectionTemp = [];
        % generate projections for each relevant channel
        for n = 1:length(nuclearChannelIndices)
            TempChannel = nuclearChannelIndices(n);
            InvertFlag = invertedChannels(TempChannel)==1;
            iter = 1;
            MaxTemp = [];
        
            % iterate through and store slices
            for i=TempChannel:NChannelsMeta:length(ImageSizes)
                if ismember(i,ElligibleImages)
                    MaxTemp(:,:,iter)=ImageTemp{i,1};
                    iter = iter + 1;
                end
            end
            i_slice = mat2gray(mean(MaxTemp,3))*255;
            if InvertFlag
                i_mask = getEmbryoMaskLive(i_slice,liveExperiment.pixelSize_um);
                i_slice(~i_mask) = 255;
                i_slice = imcomplement(i_slice);
            end
            projectionTemp(:,:,n) = i_slice;
        end
        SurfImage = mean(double(projectionTemp),3);

        %For Nikon spinnind disk data load the stitched surface image that
        %was made by FindAPAxisFullEmbryo
        if strcmp(FileMode,'DSPIN')
            SurfImage = bfopen([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryoSurf.tif']);
            SurfImage = SurfImage{1}{1};
        end
        
        %Rotates the full embryo image to match the rotation of the zoomed
        %time series
        zoom_angle = 0;
        full_embryo_angle = 0;
        
        if strcmp(FileMode,'LIFExport')
            
            zoom_angle = getZoomAngle(Prefix, rawPrefixPath);
            
            full_embryo_angle = getFullEmbryoAngle(fullEmbryoPath,...
                surfFile, Prefix);
            
            
        elseif strcmp(FileMode,'LSM')|strcmp(FileMode,'CZI')|strcmp(FileMode,'DSPIN') %CS20170912
            LSMSurf=bfopen([fullEmbryoPath,D(1).name]);
            LSMMeta=LSMSurf{:,4};
            LSMMeta2=LSMSurf{:,2};
            
            %Get the metadata about rotation. This will depend on whether
            %we have a LSM or CZI file
            if strcmp(FileMode,'LSM')
                %Figure out the rotation of the full embryo image
                full_embryo_angle = LSMMeta2.get('Recording Rotation #1');
                
                %Figure out the rotation of the zoomed-in image
                DLSMZoom=dir([rawPrefixPath,'*.lsm']);
                LSMZoom=bfopen([rawPrefixPath, DLSMZoom(1).name]);
                LSMMetaZoom2=LSMZoom{:,2};
                zoom_angle=LSMMetaZoom2.get('Recording Rotation #1');
            elseif strcmp(FileMode,'CZI')
                %Figure out the rotation of the full embryo image
                full_embryo_angle=str2num(LSMMeta2.get('Global HardwareSetting|ParameterCollection|RoiRotation #1'));
                
                %Figure out the rotation of the zoomed-in image
                DLSMZoom=dir([rawPrefixPath,'*.czi']);
                LSMZoom=bfopen([rawPrefixPath,...
                    DLSMZoom(1).name]);
                LSMMetaZoom2=LSMZoom{:,2};
                zoom_angle=str2double(LSMMetaZoom2.get('Global HardwareSetting|ParameterCollection|RoiRotation #1'));
            end
        end

        SurfImage = imrotate(SurfImage, -zoom_angle + full_embryo_angle);
        
        mkdir([DropboxFolder,filesep,Prefix, filesep, 'DV']);
        maxSurfSavePath = [DropboxFolder,filesep,Prefix, filesep, 'DV', filesep, 'surf_max.tif'];
        imwrite(SurfImage,maxSurfSavePath);
        
        
        clear ImageTemp
        
        %Zoom factor
        MovieZoom=PixelSizeFullEmbryoSurf(1)/PixelSizeZoom(1);
        SurfZoom=1;     %We'll call the zoom of the full embryo image 1
        
        %Get the size of the zoom image
        Columns = str2double(MetaZoom.getPixelsSizeX(0));
        Rows = str2double(MetaZoom.getPixelsSizeY(0));
        surf_size = size(SurfImage);
        SurfColumns= surf_size(2);
        SurfRows=surf_size(1);
        ZoomRatio = MovieZoom / SurfZoom;
    end
    
    %Get the information about the AP axis as well as the image shifts
    %used for the stitching of the two halves of the embryo
    loadMatFile([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'], true);
    
    
    %Make a folder to store the images
    mkdir([DropboxFolder,filesep,Prefix,filesep,'APDetection']);
    
    
    hisMat = getHisMat(liveExperiment);

    ZoomImage = hisMat(:, :, end-1);
    
    %Do a correlation between the zoomed in and zoomed out surface images
    %to figure out the shift.
    FullEmbryo=imread([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo.tif']);
    
    if ~SkipAlignment && histoneChannelPresent
        
        ShiftColumn=0;
        ShiftRow=0;
        NucMaskZoomIn = false(size(ZoomImage));
        NucMaskZoomOut = false(size(SurfImage));
        
        if ZoomRatio < 24  %ZoomRatio > 1 && ZoomRatio < 24. AR 12/4/17- where did this number come from
            
            %Enlarge the zoomed out image so we can do the cross-correlation
            SurfImageResized=imresize(SurfImage, ZoomRatio);
            
            %Calculate the correlation matrix and find the maximum
            im1 = ZoomImage;
            im2 = SurfImageResized;
            
            C_raw = normxcorr2(imgaussfilt(im1,2), imgaussfilt(im2,2));            

            [i_mask, ~] = getEmbryoMaskLive(im2,liveExperiment.pixelSize_um/ZoomRatio);
            i_mask2 = false(size(C_raw));
            
            x_ref = ceil(size(im1,2)/2);
            y_ref = ceil(size(im1,1)/2);
            
            i_mask2(y_ref+1:size(i_mask,1)+y_ref,x_ref+1:size(i_mask,2)+x_ref) = i_mask;
            se = strel('disk',x_ref);
            i_mask_dil = imdilate(i_mask2,se);
            
            C = C_raw;
            C(~i_mask_dil) = 0;
            
            [Max2,MaxRows]=max(C);
            [~,MaxColumn]=max(Max2);
            MaxRow=MaxRows(MaxColumn);
            [CRows,CColumns]=size(C);
            
            
            %This shift is now converted to the zoom out distances. If we
            %want to translate to the zoomed in coordinates we need to
            %multiply again by ZoomRatio.
            ShiftRow=round((MaxRow-(CRows/2+1))/ZoomRatio);
            ShiftColumn=round((MaxColumn-(CColumns/2+1))/ZoomRatio);
           
            %How well did we do with the alignment?
            [RowsResized,ColumnsResized]=size(SurfImageResized);
            [RowsZoom,ColumnsZoom]=size(ZoomImage);
            
            %If manual alignment was done before then load the results
            if exist('ManualAlignmentDone', 'var')
                if ManualAlignmentDone
                    disp('Manual alignment results saved.')
                    if ~yToManualAlignmentPrompt
                        Answer=input('Would you like to use them (y/n)? ','s');
                    else
                        Answer = 'y';
                        ManualAlignment = false;
                    end
                    if strcmpi(Answer,'y')
                        load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'],'ShiftRow','ShiftColumn')
                    elseif strcmpi(Answer,'n')
                        disp('Deleting manual alignment results')
                        clear ManualAlignmentDone
                    else
                        msgbox('Error: Answer not recognized');
                        error('Answer not recognized')
                    end
                end
            end
            
            %If we can, we'll crop the surface image to center the overlay
            %on the area of acquisition. If not, we'll just show both images
            %without cropping.
            
            if ((round(RowsResized/2-RowsZoom/2+ShiftRow*ZoomRatio))>0)&...
                    (round(RowsResized/2+RowsZoom/2-1+ShiftRow*ZoomRatio)<=RowsResized)&...
                    (round(ColumnsResized/2-ColumnsZoom/2+ShiftColumn*ZoomRatio)>0)&...
                    (round(ColumnsResized/2+ColumnsZoom/2-1+ShiftColumn*ZoomRatio)<=ColumnsResized)
                RowsResizedRange=...
                    round(RowsResized/2-RowsZoom/2+ShiftRow*ZoomRatio):...
                    round(RowsResized/2+RowsZoom/2-1+ShiftRow*ZoomRatio);
                ColumnsResizedRange=...
                    round(ColumnsResized/2-ColumnsZoom/2+ShiftColumn*ZoomRatio):...
                    round(ColumnsResized/2+ColumnsZoom/2-1+ShiftColumn*ZoomRatio);
            else
                RowsResizedRange=1:size(ZoomImage,1);%RowsResized;
                ColumnsResizedRange=1:size(ZoomImage,2);%ColumnsResized;
            end
            
            
            %Make an overlay of the zoomed in and zoomed out real
            %images as well as of a quickly segmented nuclear mask. If this
            %fails, we'll switch to ManualAlignment.
            
            try
                %Real image overlay
                %Crop the zoomed out image to match the zoomed in one
                SurfImageResizeZoom=...
                    SurfImageResized(RowsResizedRange,ColumnsResizedRange);
                ImOverlay=cat(3,mat2gray(SurfImageResizeZoom),...
                    +mat2gray(ZoomImage),zeros(size(SurfImageResizeZoom)));
                
                %Show the correlation image, but crop it a little bit
                contFig = figure(2);
                contAxes = axes(contFig);
                
                %HG: Note that I changed the ranges here by two pixels at least.
                lX = uint16((CRows+1)/2-RowsZoom+1);
                uX = uint16((CRows-1)/2+RowsZoom);
                lY = uint16((CColumns+1)/2-ColumnsZoom+1);
                uY = uint16((CColumns-1)/2+ColumnsZoom);
                
                contourf(contAxes,imresize(abs(C(lX:uX,lY:uY)),.5));
                saveas(contFig, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'AlignmentCorrelation.tif']);
            catch
                warning('Could not generate correlation image. Switching to manual alignment')
                ManualAlignment=true;
            end
            
        else
            warning('Not able to do the cross correlation. Switching to manual alignment mode.')
            ManualAlignment=true;
        end
        
        
        %See if we need the manual alignment
        if ManualAlignment
            %See if we need to load the previous manual alignment results
            [ShiftColumn,ShiftRow]=ManualAPCorrection(SurfImage,ZoomImage,C,ZoomRatio,ShiftRow,ShiftColumn,...
                FullEmbryo, ZoomRatio, SurfRows,Rows, Columns, coordA, coordP, SurfColumns);
            ManualAlignmentDone = true;
        end
        
    else
        warning('Not able to do the cross correlation. Assuming no shift between surface-level and movie-level images.')
    end
    
    
    ImageCenter=[SurfRows/2,SurfColumns/2];
    
    %This is for the acquisition/zoom image
    TopLeft=[ImageCenter(1)-Rows/ZoomRatio/2+ShiftRow,...
        ImageCenter(2)-Columns/ZoomRatio/2+ShiftColumn];
    BottomRight=[ImageCenter(1)+Rows/ZoomRatio/2+ShiftRow,...
        ImageCenter(2)+Columns/ZoomRatio/2+ShiftColumn];
    coordAZoom=(coordA-[TopLeft(2),TopLeft(1)])*ZoomRatio;
    coordPZoom=(coordP-[TopLeft(2),TopLeft(1)])*ZoomRatio;
    if exist('coordD', 'var')
        coordDZoom=(coordD-[TopLeft(2),TopLeft(1)])*ZoomRatio;
        coordVZoom=(coordV-[TopLeft(2),TopLeft(1)])*ZoomRatio;
    end
    
    fullFigure = figure(7);
    fullAxes = axes(fullFigure);
    imshow(imadjust(mat2gray(FullEmbryo)),'DisplayRange',[],'InitialMagnification',100,'Parent', fullAxes)
    hold(fullAxes, 'on')
    rectangle('Position',[TopLeft([2,1]),BottomRight([2,1])-TopLeft([2,1])],'EdgeColor','r')
    plot(coordA(1),coordA(2),'.g','MarkerSize',30)
    plot(coordP(1),coordP(2),'.r','MarkerSize',30)
    plot([coordA(1),coordP(1)],[coordA(2),coordP(2)],'-b')
    hold(fullAxes, 'off')
    saveas(fullFigure, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryoArea.tif']);
    
    zoomImageFigure = figure(9);
    zoomImageAxes = axes(zoomImageFigure);
    imshow(imadjust(ZoomImage),[], 'Parent', zoomImageAxes)
    %imshow(NucMaskZoomIn)
    %imshow(NucMaskZoomOutResizedCropped)
    hold(zoomImageAxes,'on')
    plot(zoomImageAxes,[coordAZoom(1),coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'-b')
    plot(zoomImageAxes,[coordAZoom(1)+1,coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'--r')
    plot(zoomImageAxes,[coordAZoom(1)-1,coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'-.g')
    hold(zoomImageAxes,'off')
    saveas(zoomImageFigure, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'ZoomedEmbryoAP.tif']);
    
    
    %With AP coordinates in hand we can now determine the AP position of
    %all particles. Look I my notes in "Calculating AP positions" in Notability
    %for details of the calculation.
    
    %Angle between the x-axis and the AP-axis
    APAngle=atan2((coordPZoom(2)-coordAZoom(2)),(coordPZoom(1)-coordAZoom(1)));
    APLength = norm(coordAZoom-coordPZoom);
    if exist('coordD', 'var')
        DVLength = norm(coordDZoom-coordVZoom);
    else
        DVLength = APLength/2; %rough estimate but surprisingly accurate
    end
    saveVars = [saveVars, 'APLength', 'DVLength'];
    
    APPosImage=zeros(size(ZoomImage));
    [Rows,Columns]=size(ZoomImage);
    
    for y=1:Rows
        for x=1:Columns
            Angle=atan2((y-coordAZoom(2)),(x-coordAZoom(1)));
            
            Distance=sqrt((coordAZoom(2)-y).^2+(coordAZoom(1)-x).^2);
            APPosition=Distance.*cos(Angle-APAngle);
            APPosImage(y,x)=APPosition/APLength;
        end
    end
    
    
    APbinID = 0:APResolution:1;
    
    
    APPosBinImage=zeros(size(APPosImage));
    for i=1:(length(APbinID)-1)
        FilteredMask=(APbinID(i)<=APPosImage)&(APbinID(i+1)>APPosImage);
        
        APPosBinImage=APPosBinImage+FilteredMask*i;
    end
    
    
    zoomOverlayFigure = figure(10);
    zoomOverlayAxes = axes(zoomOverlayFigure);
    ZoomOverlay=cat(3,mat2gray(ZoomImage)/2+mat2gray(APPosBinImage)/2,...
        mat2gray(ZoomImage)/2,mat2gray(ZoomImage)/2);
    imshow(ZoomOverlay, 'Parent', zoomOverlayAxes)
    hold(zoomOverlayAxes, 'on')
    plot(zoomOverlayAxes,[coordAZoom(1),coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'-b')
    plot(zoomOverlayAxes,[coordAZoom(1)+1,coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'--r')
    plot(zoomOverlayAxes,[coordAZoom(1)-1,coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'-.g')
    hold(zoomOverlayAxes,'off')
    
    DV_correction = 0;
    if exist([DropboxFolder,filesep,Prefix,filesep,'DV',filesep,'Classified_image.tif'], 'file')
        correctDV = 1;
    end
    if correctDV
        if exist([DropboxFolder,filesep,Prefix,filesep,'DV',filesep,'DV_correction.mat'], 'file')
            load([DropboxFolder,filesep,Prefix,filesep,'DV',filesep,'DV_correction.mat'],'DV_correction');
        else
            try
                DV_correction = FindDVShift_full(Prefix);
            catch
                disp('failed to apply dv correction');
            end
            save([DropboxFolder,filesep,Prefix,filesep,'DV',filesep,'DV_correction.mat'],'DV_correction');
        end
        saveVars = [saveVars, 'DV_correction'];
    end
    
    
    DVLength = APLength/2;
    if exist([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'], 'file')
        for ChN=1:NChannels
            for i=1:length(Particles{ChN})
                %Angle between the x-axis and the particle using the A position as a
                %zero
                Angles=atan2((Particles{ChN}(i).yPos-coordAZoom(2)),...
                    (Particles{ChN}(i).xPos-coordAZoom(1)));
                
                %Distance between the points and the A point
                Distances=sqrt((coordAZoom(2)-Particles{ChN}(i).yPos).^2+(coordAZoom(1)-Particles{ChN}(i).xPos).^2);
                APPositions=Distances.*cos(Angles-APAngle);
                Particles{ChN}(i).APpos=APPositions/APLength;
                
                %Determine the distance perpendicular to the AP axis. This is a
                %proxy for a DV axis.
                DVPositions = Distances.*sin(Angles-APAngle); %units of pixels on the surface of blastoderm.
                if correctDV
                    %ventral midline is dv_correction pixels away from the
                    %AP axis.
                    Particles{ChN}(i).DVpos=abs(DVPositions-DV_correction)/DVLength;
                    %so DVpos is pixels away from the ventral midline
                    %across the blastoderm.
                else
                    Particles{ChN}(i).DVpos=DVPositions/DVLength;
                end
            end
        end
    end
    
    %Save AP detection information
    
    %Information about shifts
    if exist('xShift', 'var')
        saveVars=[saveVars,'xShift','yShift'];
    elseif  exist('xShift1', 'var')
        saveVars=[saveVars,'xShift1','yShift1',...
            'xShift2','yShift2'];
    end
    %Rotation information
    if exist('zoom_angle', 'var')
        ImageRotation=zoom_angle;
        saveVars=[saveVars,'ImageRotation'];
    end
    
    save([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'],saveVars{:})
    
    if exist('ManualAlignmentDone', 'var')
        if ManualAlignmentDone
            save([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'],'ManualAlignmentDone',...
                'ShiftColumn','ShiftRow','-append')
        end
    end
    
end




%Save particle information
if exist([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'], 'file')
    
    %Bring the one channel case back to the legacy setting
    if NChannels==1
        Particles=Particles{1};
    end
    try
        save([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'],...
            'Particles','SpotFilter', '-v6');
    catch
        save([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'],...
            'Particles','SpotFilter', '-v7.3', '-nocompression');
    end
end

if correctDV
    CheckDivisionTimes(Prefix, 'lazy');
    %for convenience.
end

try
    Ellipses = getEllipses(liveExperiment);
    schnitzcells = getSchnitzcells(liveExperiment);
    
    ellipsesOld = Ellipses;
    schnitzcellsOld = schnitzcells; 
    
    [EllipsePos, APAngle, APLength]...
        = convertToFractionalEmbryoLength(Prefix);
    
    for s = 1:length(schnitzcells)
        for f = 1:length(schnitzcells(s).frames)
            ellipseInd = schnitzcells(s).cellno(f);
            schnitzcells(s).APpos(f) = EllipsePos{f}(ellipseInd);
        end
    end
    
    ellipsesSizeUnchanged(ellipsesOld, Ellipses);
    schnitzcellsSizeUnchanged(schnitzcellsOld, schnitzcells);
    
    
    
    save2([liveExperiment.resultsFolder, filesep,'Ellipses.mat'], Ellipses);
    save2([liveExperiment.resultsFolder, filesep,Prefix,'_lin.mat'], schnitzcells);
    
catch
    warning('failed to add AP positions to nuclear structures.')
end


end

