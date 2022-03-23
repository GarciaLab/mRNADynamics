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
%V2: Changed this function to use a correlation in order to center the
%images.

cleanupObj = onCleanup(@myCleanupFun);


%Default set of variables to save
saveVars={'coordA','coordP','coordAZoom','coordPZoom'};


SkipAlignment=false;
ManualAlignment=true;
NoAP=false;
SelectChannel=0;
InvertHis=false;
optionalResults = '';
yToManualAlignmentPrompt = false;
correctDV = false;


for i=1:length(varargin)
    switch varargin{i}
        case {'SkipAlignment'}
            disp('Skipping alignment step')
            SkipAlignment=1;
        case {'ManualAlignment'}
            ManualAlignment=true;
        case {'NoAP'}
            NoAP=1;
        case {'SelectChannel'}
            SelectChannel=1;
        case {'optionalResults'}
            optionalResults = varargin{i+1};
        case {'yToManualAlignmentPrompt'}
            yToManualAlignmentPrompt = 1;
        case {'correctDV'}
            correctDV = true;
    end
end

liveExperiment = LiveExperiment(Prefix);

DropboxFolder = liveExperiment.userResultsFolder;
PreProcPath = liveExperiment.userPreFolder;
RawDynamicsPath = liveExperiment.userRawFolder;

Channel1 = liveExperiment.Channel1;
Channel2 = liveExperiment.Channel2;
Channel3 = liveExperiment.Channel3;

APResolution = liveExperiment.APResolution;

if exist([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'], 'file')
    load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'], 'Particles', 'SpotFilter')
    load([DropboxFolder,filesep,Prefix,filesep,'Spots.mat'], 'Spots')
    
    %Create the particle array. This is done so that we can support multiple
    %channels. Also figure out the number of channels
    if iscell(Particles)
        NChannels=length(Particles);
    else
        Particles={Particles};
        Spots={Spots};
        NChannels=1;
    end
    
    %Now, get the particle positions (if they're not there already).
    for ChN=1:NChannels
        Particles = addPositionsToParticles(Particles, Spots, ChN);
    end
    
    if isfield(Particles{ChN},'DVpos')
        warning('Particles.mat already has DV positions stored. They will be rewritten')
    end
    if isfield(Particles{ChN},'APpos')
        warning('Particles.mat already has AP positions stored. They will be rewritten')
    end
    
else
    warning('No Particles.mat found. Just updating APDetection.mat')
end


%See if we had any lineage/nuclear information
hisDir=dir([PreProcPath,filesep,Prefix,filesep,'*his*']);
if ~isempty(hisDir)
    histoneChannelPresent = true;
else
    histoneChannelPresent = false;
end



[FileMode, EmbryoName, projectDate] = getMicroscope(Prefix);
rawPrefixPath = [RawDynamicsPath,filesep,projectDate,filesep,EmbryoName,filesep];
fullEmbryoPath = [rawPrefixPath, 'FullEmbryo', filesep];

if ~NoAP
    %If you want to select which channel to load as alignment
    ChannelToLoad = determineChannelToLoad(SelectChannel, liveExperiment.Channels);
    
    %Get information about all images. This depends on the microscope used.
    
    %Get the information about the zoom
    if strcmp(FileMode,'TIF')
        D=dir([rawPrefixPath,'*.tif']);
        ImageInfo = imfinfo([rawPrefixPath,D(1).name]);
        
        %Figure out the zoom factor
        MovieZoom=ExtractInformationField(ImageInfo(1),'state.acq.zoomFactor=');
        MovieZoom=str2double(MovieZoom);
        
        
        %Get the zoomed out surface image and its dimensions from the FullEmbryo folder
        D=dir([fullEmbryoPath,filesep,'*.tif']);
        SurfName=D(find(~cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
        SurfImage=imread([rawPrefixPath,...
            'FullEmbryo',filesep,SurfName],ChannelToLoad);
        
        %Get the size of the zoom image
        Rows = str2double(ExtractInformationField(ImageInfo(1), 'state.acq.linesPerFrame='));
        Columns = str2double(ExtractInformationField(ImageInfo(1), 'state.acq.pixelsPerLine='));
        
        SurfInfo = imfinfo([RawDynamicsPath, filesep, projectDate, filesep, EmbryoName, filesep, 'FullEmbryo', filesep, SurfName]);
        SurfZoom = ExtractInformationField(SurfInfo(1), 'state.acq.zoomFactor=');
        SurfZoom = str2double(SurfZoom);
        
        SurfRows = str2double(ExtractInformationField(SurfInfo(1), 'state.acq.linesPerFrame='));
        SurfColumns = str2double(ExtractInformationField(SurfInfo(1), 'state.acq.pixelsPerLine='));
        
        %HG: I had to add this for some surface images that were edited
        %with ImageJ and lost their metadata. Note that I'm hardcoding the
        %zoom of the low magnification images.
        if isnan(SurfRows)
            SurfRows=SurfInfo(1).Height;
            SurfColumns=SurfInfo(1).Width;
            SurfZoom=1;
        end
        
        
        %Get the full embryo image
        FullEmbryo=imread([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo.tif']);
        
        %Ratio between the two zoom levels
        ZoomRatio = MovieZoom / SurfZoom;
        ResizeFactor = max([Rows/SurfRows*ZoomRatio, Columns/SurfColumns*ZoomRatio]);
        % ES 2013-10-30: the reason I have to define ResizeFactor differently
        % from ZoomRatio is because you can't necessarily infer the
        % microns-per-pixel resolution from the zoom alone: it also depends on
        % the dimensions of the image. This may not work for all possible
        % resolutions, though...
        ZoomRatio = ResizeFactor;
        
    else
        SurfName=[];
        
        %Figure out the different channels
        %NuclearChannel=contains([Channel1,Channel2,Channel3],'nuclear','IgnoreCase',true);
        
        
        if any(contains([Channel1,Channel2,Channel3],'nuclear','IgnoreCase',true))
            if SelectChannel
                HisChannel = ChannelToLoad;
                InvertHis=0;
            elseif ~any(contains([Channel1,Channel2,Channel3],'inverted','IgnoreCase',true))
                HisChannel = ChannelToLoad;
                InvertHis=0;
            elseif any(contains([Channel1,Channel2,Channel3],'inverted','IgnoreCase',true))
                HisChannel= ChannelToLoad;
                InvertHis=1;
            else
                error('You should check the MovieDatabase.csv to see whether ":Nuclear" or "invertedNuclear" is plugged into your channels.')
            end
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
        ImageTemp=bfopen(surfFile);
        MetaFullEmbryo= ImageTemp{:, 4};
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
        NSlices=str2double(MetaFullEmbryo.getPixelsSizeZ(0));
        clear MaxTemp
        
        %Do a maximum projection
        
        %Look for the image with the largest size. In this way, we avoid
        %loading individual tiles in the case of a tile scan.
        for i=1:size(ImageTemp,1)
            SizesImages(i)=size(ImageTemp{i,1}{1,1},1);
        end
        [~,ImageCellToUse]=max(SizesImages);
        
        
        %By looking at the last image we make sure we're avoiding the
        %individual tiles if we're dealing with tile scan
        for i=HisChannel:NChannelsMeta:size(ImageTemp{ImageCellToUse,1},1)
            MaxTemp(:,:,i)=ImageTemp{ImageCellToUse,1}{i,1};
        end
        
        if InvertHis
            SurfImage=MaxTemp(:,:,HisChannel+round(NSlices/2)*NChannelsMeta);
        else
            SurfImage=max(MaxTemp,[],3);
        end
        
        %For Nikon spinnind disk data load the stitched surface image that
        %was made by FindAPAxisFullEmbryo
        if strcmp(FileMode,'DSPIN')
            SurfImage = bfopen([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryoSurf.tif']);
            SurfImage = SurfImage{1}{1};
        end
        %%
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
            
            %Get the surface image
            %SurfImage=LSMSurf{1}{HisChannel,1};
            
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
        %%
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
    %%
    
    
    %Get the information about the AP axis as well as the image shifts
    %used for the stitching of the two halves of the embryo
    load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])
    
    
    %Make a folder to store the images
    mkdir([DropboxFolder,filesep,Prefix,filesep,'APDetection'])
    
    
    
    % Load the last frame of the zoomed-in image, for calculating the
    % cross-correlation. This part now only requires the same
    % ChannelToLoad, which is defined above using ":Nuclear" or
    % ":invertedNuclear"
    
    %Get the surface image in the zoomed case by looking at the last
    %frame of our movie
    % From "Nuclear" or "invertedNuclear", it should've made His images
    % when you run ExportDataforLivemRNA. If you edit the channels in the
    % MovieDatabase after you ran the ExportDataForLivemRNA, then run this
    % script, the code might freak out. I'll put a warning message about
    % that.
    %     DHis=dir([PreProcPath,filesep,Prefix,filesep,Prefix,'-His*.tif']);
    hisMat = getHisMat(liveExperiment);
    %     if ~isempty(DHis)
    %         try
    %3D stack
    %             hisMat = imreadStack([PreProcPath,filesep,Prefix,filesep,Prefix,'-His.tif']);
    ZoomImage = hisMat(:, :, end-1);
    %         catch
    %             %single planes
    %             ZoomImage=imread([PreProcPath,filesep,Prefix,filesep,DHis.name]);
    %         end
    %     else
    %         disp('Did you run ExportDataForLivemRNA again, after editing the MovieDatabase.csv with ":Nuclear" ("invertedNuclear")?')
    % This might be the case, for instance, if you're just trying
    % to find AP information about an image without using FISH
    % code. In that case, just extract the nuclei from the last
    % raw image.
    %         DGFP = dir([rawPrefixPath, '*.tif']);
    %         ImageInfo = imfinfo([rawPrefixPath, DGFP(end).name]);
    %         NumFramesAndSlices = length(ImageInfo)/2;
    %         RawImage3M = NaN(Rows, Columns, NumFramesAndSlices);
    %         for lImageIndex = 1:NumFramesAndSlices
    %             RawImage3M(:, :, lImageIndex) = imread([rawPrefixPath,DGFP(end).name],'Index', 2*(lImageIndex-1) + ChannelToLoad);
    %         end
    %         ZoomImage = max(RawImage3M, [], 3) / 255;
    %     end
    
    
    
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
            
            %             if InvertHis
            %                 im1 = imcomplement(im1);
            %                 im2 = imcomplement(im2);
            %             end
            
            %             try
            %                 C = gather(normxcorr2(gpuArray(im1), gpuArray(im2)));
            %             catch
            C = normxcorr2(im1, im2);
            %             end
            
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
                RowsResizedRange=1:RowsResized;
                ColumnsResizedRange=1:ColumnsResized;
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
                
%                 %Nuclear mask overlay
%                 NucMaskZoomOut= GetNuclearMask(SurfImage,2.5,0);
%                 NucMaskZoomOutResized=imresize(NucMaskZoomOut, ZoomRatio);
%                 NucMaskZoomOutResizedCropped=...
%                     NucMaskZoomOutResized(RowsResizedRange,ColumnsResizedRange);
%                 NucMaskZoomIn=GetNuclearMask(ZoomImage,8,2);
%                 ImOverlayMask=cat(3,mat2gray(NucMaskZoomOutResizedCropped),...
%                     +mat2gray(NucMaskZoomIn),zeros(size(NucMaskZoomOutResizedCropped)));
%                 
%                 alOvFig = figure(1);
%                 imOv = subplot(2,1,1);
%                 imshow(ImOverlay, 'Parent', imOv)
%                 imOvMask = subplot(2,1,2);
%                 imshow(ImOverlayMask, 'Parent', imOvMask)
%                 saveas(alOvFig, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'AlignmentOverlay.tif']);
%                 
                
                
                
                %Show the correlation image, but crop it a little bit
                contFig = figure(2);
                contAxes = axes(contFig);
                %HG: Note that I changed the ranges here by two pixels at
                %least.
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

function ChannelToLoad = determineChannelToLoad(SelectChannel, Channels)

if SelectChannel
    list = string(Channels);
    [indx,tf] = listdlg('PromptString','Select the channel to use for alignment:','ListString',list);
    ChannelToLoad = indx;
else
    % From now, we will use a better way to define the channel for
    % alignment (used for cross-correlation).
    % Find channels with ":Nuclear"
    ChannelToLoadTemp= contains([Channels(1),Channels(2),Channels(3)],'nuclear','IgnoreCase',true);
    
    % Define the Channel to load, for calculating the cross-correlation
    % In future, we can think about combining multiple channels for
    % calculating the cross-correlation to get more accurate estimates.
    % For now, let's pick only one channel for this. For multiple
    % channels, let's first pick the first channel. This can be fine in
    % most cases, since we normally use lower wavelength for sth we
    % care more, or we get better signal from those.
    if sum(ChannelToLoadTemp) && sum(ChannelToLoadTemp)==1
        ChannelToLoad=find(ChannelToLoadTemp);
    elseif sum(ChannelToLoadTemp) && length(ChannelToLoadTemp)>=2
        ChannelToLoad=find(ChannelToLoadTemp);
        ChannelToLoad = ChannelToLoad(1);
    else
        error('No histone channel found. Was it defined in MovieDatabase as :Nuclear or :InvertedNuclear?')
    end
    
end



end