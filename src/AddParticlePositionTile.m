function [Particles, SpotFilter] = AddParticlePositionTile(Prefix, varargin)
%% author: Gabriella Martini
% date created: 1/29/20
% date last modified: 8/3/20
%Default set of variables to save

% DESCRIPTION
% Aligns movie with full embryo image created 
% by stitching tiled images using spatial cross-correlation.
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
% 'UseFullEmbryoBoundaries': restrict zoomed in movie to be completely
% contained by the boundaries of the full embryo image.
% 'ImposeEmbryoMask': restruct zoomed in movie to be completely contained
% within the boundaries of the embryo with some tolerance. 
% 'embryo_mask_tolerance': number between 0 and 1 which corresponds to the
% percentage of the movie that must be contained within the embryo when
% using the 'ImposeEmbryoMask' flag



%% 


%Default set of variables to save
saveVars={'coordA','coordP','coordAZoom','coordPZoom'};
 

NoAP=false;
SelectChannel=0;
InvertHis=false;
optionalResults = '';
UseFullEmbryoBoundaries = false;
ImposeEmbryoMask = false; % Note that this won't work if a lot of zoomImage is outside of embryo. Use embryo_mask_tolerance to adjust how rigid this feature is.
embryo_mask_tolerance = .9; %Change value using varargin
ImposeManualBounds = false;

close all
%% 
i = 1;
while i <= length(varargin)
    switch varargin{i}
        case {'NoAP'}
            NoAP=1;
            i = i + 1;
        case {'SelectChannel'}
            SelectChannel=1;
            i = i + 1;
        case {'optionalResults'}
            optionalResults = varargin{i+1};
            i = i + 2;
        case {'UseFullEmbryoBoundaries'}
            UseFullEmbryoBoundaries = true;
            i = i + 1;
        case {'ImposeEmbryoMask'}
            ImposeEmbryoMask = true;
            i = i + 1;
        case {'ImposeManualBounds'}
            ImposeManualBounds = true;
            i = i + 1;
        case {'embryo_mask_tolerance'}
            embryo_mask_tolerance = varargin{i+1};
            i = i + 2;
        otherwise
            error(['"',varargin{i}, '" not a known flag.']);
            
    end
end
%% 



%Get the relevant folders for this data set
[RawDynamicsPath, ~, DefaultDropboxFolder, DropboxFolder, ~, PreProcPath,...
    configValues,...
    movieDatabasePath, movieDatabase] = DetermineAllLocalFolders(Prefix, optionalResults);

[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoopEnd, APResolution,...
    Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
    nc9, nc10, nc11, nc12, nc13, nc14, CF,Channel3,prophase,metaphase, anaphase] = getExperimentDataFromMovieDatabase(Prefix, movieDatabase);

%% 

if exist([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'], 'file')
    load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
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
        for i=1:length(Particles{ChN})
            for j=1:length(Particles{ChN}(i).Frame)
                [x,y,z]=SpotsXYZ(Spots{ChN}(Particles{ChN}(i).Frame(j)));
                if ~isempty(x)
                    Particles{ChN}(i).xPos(j)=x(Particles{ChN}(i).Index(j));
                    Particles{ChN}(i).yPos(j)=y(Particles{ChN}(i).Index(j));
                    Particles{ChN}(i).zPos(j)=z(Particles{ChN}(i).Index(j));
                end
            end
        end

        if isfield(Particles{ChN},'APpos')
            warning('Particles.mat already has AP positions stored. They will be rewritten')
        end
    end
    
else
    warning('No Particles.mat found. Just updating APDetection.mat')
end

%% 

%See if we had any lineage/nuclear information
hisDir=dir([PreProcPath,filesep,Prefix,filesep,'*-His_*']);
if ~isempty(hisDir)
    histoneChannelPresent = true;
else
    histoneChannelPresent = false;
end



[FileMode, EmbryoName, projectDate] = getMicroscope(Prefix);
rawPrefixPath = [RawDynamicsPath,filesep,projectDate,filesep,EmbryoName,filesep];
fullEmbryoPath = [rawPrefixPath, 'FullEmbryo', filesep];

%% 


if ~strcmp(FileMode,'LIFExport')
    msg = 'Non_LIFExport filetypes not supported.';
    error(msg)
end
    

if ~NoAP
    %If you want to select which channel to load as alignment.
    if SelectChannel
        list = string({Channel1,Channel2,Channel3});
        [indx,tf] = listdlg('PromptString','Select the channel to use for alignment:','ListString',list);
        ChannelToLoad = indx;
    else
        % From now, we will use a better way to define the channel for
        % alignment (used for cross-correlation).
        % Find channels with ":Nuclear"
        ChannelToLoadTemp=contains([Channel1,Channel2,Channel3],'nuclear','IgnoreCase',true);
        
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
    
    %Get information about all images. This depends on the microscope used.
    


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

    ImageTemp=bfopen([fullEmbryoPath,filesep,D(end).name]);
    MetaFullEmbryo= ImageTemp{:, 4};
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



    %MidImage =imread([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo.tif']);
    SurfImage = imread([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryoSurf.tif']);


    %%
    %Rotates the full embryo image to match the rotation of the zoomed
    %time series
    zoom_angle = 0;
    full_embryo_angle = 0;


    if isfolder([rawPrefixPath, 'MetaData'])
        xml_file_path = dir([rawPrefixPath,'MetaData', filesep, '*.xml']);
        xml_file = xml_file_path(1).name;
        xDoc = searchXML([rawPrefixPath,'MetaData', filesep, xml_file]);
        zoom_angle = str2double(evalin('base','rot'));
    else
        warning('No time series metadata found.')
    end
    if isfolder([fullEmbryoPath,'MetaData'])
        xml_file_path2 = dir([fullEmbryoPath,'MetaData', filesep,'*Surf*.xml']);
        xml_file2 = xml_file_path2(1).name;
        xDoc2 = searchXML([fullEmbryoPath,'MetaData', filesep, xml_file2]);
        full_embryo_angle = str2double(evalin('base','rot'));
    else
        warning('No full embryo metadata found.')
    end

    evalin('base','clear rot')

    %%
    SurfImage = imrotate(SurfImage, -zoom_angle + full_embryo_angle);


    % This is the original z stack full embryo surface image file from bfopen 
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
    DHis=dir([PreProcPath,filesep,Prefix,filesep,Prefix,'-His*.tif']);
    if ~isempty(DHis)
        ZoomImage=imread([PreProcPath,filesep,Prefix,filesep,DHis(end-1).name]);
    else
        error('Did you run ExportDataForLivemRNA again, after editing the MovieDatabase.csv with ":Nuclear" ("invertedNuclear")?')
    end
    

    
    %Do a correlation between the zoomed in and zoomed out surface images
    %to figure out the shift.
    
    % This is the mid-saggital image
    FullEmbryo=imread([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo.tif']);
    %% 
    
    if histoneChannelPresent
        
        ShiftColumn=0;
        ShiftRow=0;
        NucMaskZoomIn = false(size(ZoomImage));
        NucMaskZoomOut = false(size(SurfImage));
        
        
        if ZoomRatio > 24
            warning('Not able to do the cross correlation. Assuming no shift between surface-level and movie-level images.')
        else
            if ZoomRatio <= 1
                ZoomImageResized = imresize(ZoomImage, 1/ZoomRatio);
                im1 = ZoomImageResized;
                im2 = SurfImage;
                % Here "Resized" refers to the shared coordinate system 
                [RowsSurfResized,ColumnsSurfResized]=size(SurfImage);
                [RowsZoomResized,ColumnsZoomResized]=size(ZoomImageResized);
            else
                %Enlarge the zoomed out image so we can do the cross-correlation
                SurfImageResized=imresize(SurfImage, ZoomRatio);

                %Calculate the correlation matrix and find the maximum
                im1 = ZoomImage;
                im2 = SurfImageResized;

                % Here "Resized" refers to the shared coordinate system 
                [RowsSurfResized,ColumnsSurfResized]=size(SurfImageResized);
                [RowsZoomResized,ColumnsZoomResized]=size(ZoomImage);
            end

            C = normxcorr2(im1, im2);
            [CRows,CColumns]=size(C);

            if UseFullEmbryoBoundaries
                row_subset = RowsZoomResized:(CRows-RowsZoomResized);
                col_subset = ColumnsZoomResized:(CColumns-ColumnsZoomResized);
                CMask = zeros(size(C));
                CMask(row_subset, col_subset) = 1;
                C(CMask == 0) = min(min(C));
            elseif ImposeEmbryoMask
                load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo')
                % change to be flexible nucleus diameter 
                nucleusDiameter = getDefaultParameters(FrameInfo,['d14']);
                im2_padded = zeros(RowsSurfResized+10-rem(RowsSurfResized, 10), ColumnsSurfResized+10-rem(ColumnsSurfResized, 10), 'uint16');
                im2_padded(1:RowsSurfResized, 1:ColumnsSurfResized) = im2;
                I_resize = imresize(im2_padded, .1);
                f_sigma_resize = round(nucleusDiameter /( min([PixelSizeFullEmbryoSurf, PixelSizeZoom])*10));
                I_blurred_resize = imfilter(I_resize,...
                     fspecial('gaussian',2*f_sigma_resize,f_sigma_resize),'symmetric','conv');
                embryoMask_resized = imbinarize(I_blurred_resize);

                SubRowSize = (RowsZoomResized-rem(RowsZoomResized, 10))/10;
                SubColumnSize = (ColumnsZoomResized-rem(ColumnsZoomResized, 10))/10;
                temp_mask = zeros(size(embryoMask_resized, 1)-SubRowSize+1, size(embryoMask_resized, 2)-SubColumnSize+1);

                ZoomImageSize = SubRowSize*SubColumnSize;
                for r=1:size(temp_mask, 1)
                    for c=1:size(temp_mask, 2)
                        ZoomMask = zeros(size(embryoMask_resized));
                        ZoomMask(r:(r+SubRowSize-1), c:(c+SubColumnSize-1)) = 1;
                        if sum(sum(and(ZoomMask,embryoMask_resized))) >= embryo_mask_tolerance*ZoomImageSize
                            temp_mask(r, c) = 1;
                        end

                    end
                end
                CMask_sub = imresize(temp_mask, 10);
                CMask_sub(CMask_sub > 0) = 1;
                CMask_sub_ch = bwconvhull(CMask_sub);
                CMask = zeros(size(C));
                CMask(RowsZoomResized:(RowsZoomResized+size(CMask_sub_ch, 1)-1),ColumnsZoomResized:(ColumnsZoomResized+size(CMask_sub_ch, 2)-1))=CMask_sub_ch;
                C(CMask == 0) = min(min(C));
            elseif ImposeManualBounds
                [minRow, maxRow, minColumn, maxColumn] = RestrictOverlapWindow(C);
                CMask = zeros(size(C));
                CMask(minRow:maxRow, minColumn:maxColumn) = 1;
                C(CMask == 0) = min(min(C));
            end

            [Max2,MaxRows]=max(C);
            [~,MaxColumn]=max(Max2);
            MaxRow=MaxRows(MaxColumn);

            if ZoomRatio < 1
                %This shift kept in the zoom in coordinates (full embryo image coordinates).
                %If we want to translate to the zoomed out coordinates we need to
                %divide (I think!) again by ZoomRatio.
                ShiftRow=round((MaxRow-(CRows-1)/2));
                ShiftColumn=round((MaxColumn-(CColumns-1)/2));
            else
                ShiftRow=round((MaxRow-(CRows-1)/2)/ZoomRatio);
                ShiftColumn=round((MaxColumn-(CColumns-1)/2)/ZoomRatio);
            end

            %How well did we do with the alignment?



            %If we can, we'll crop the surface image to center the overlay
            %on the area of acquisition. If not, we'll just show both images
            %without cropping.

            if ((round(RowsSurfResized/2-RowsZoomResized/2+ShiftRow))>0)&...
                    (round(RowsSurfResized/2+RowsZoomResized/2-1+ShiftRow)<=RowsSurfResized)&...
                    (round(ColumnsSurfResized/2-ColumnsZoomResized/2+ShiftColumn)>0)&...
                    (round(ColumnsSurfResized/2+ColumnsZoomResized/2-1+ShiftColumn)<=ColumnsSurfResized)
                RowsSurfRange=...
                    round(RowsSurfResized/2-RowsZoomResized/2+ShiftRow):...
                    round(RowsSurfResized/2+RowsZoomResized/2-1+ShiftRow);
                ColumnsSurfRange=...
                    round(ColumnsSurfResized/2-ColumnsZoomResized/2+ShiftColumn):...
                    round(ColumnsSurfResized/2+ColumnsZoomResized/2-1+ShiftColumn);
            else
                RowsSurfRange=1:RowsSurfResized;
                ColumnsSurfRange=1:ColumnsSurfResized;
            end


            
            %Make an overlay of the zoomed in and zoomed out real
            %images as well as of a quickly segmented nuclear mask. 
            
            if ZoomRatio < 1
                SurfImageCrop=...
                    SurfImage(RowsSurfRange,ColumnsSurfRange);
            else 
                SurfImageCrop=...
                    SurfImageResized(RowsSurfRange,ColumnsSurfRange);
            end

            try
                %Real image overlay
                %Crop the zoomed out image to match the zoomed in one
                %im1 = ZoomImageResized or ZoomImage;
                %im2 = SurfImageResized or SurfImage;
                SurfImageCrop=...
                    im2(RowsSurfRange,ColumnsSurfRange);
                
                
                ImOverlay=cat(3,mat2gray(SurfImageCrop),...
                    +mat2gray(im1),zeros(size(SurfImageCrop)));
                
                load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo')
                nucleusDiameter = getDefaultParameters(FrameInfo,['d14']);
                NucMaskPixelSize = min([PixelSizeFullEmbryoSurf, PixelSizeZoom]);
                %Nuclear mask overlay
                NucMaskSurf= GetNuclearMaskTile(SurfImageCrop,nucleusDiameter,NucMaskPixelSize, ImposeEmbryoMask);
                NucMaskZoom=GetNuclearMaskTile(im1,nucleusDiameter, NucMaskPixelSize,ImposeEmbryoMask);
                ImOverlayMask=cat(3,mat2gray(NucMaskSurf),...
                    +mat2gray(NucMaskZoom),zeros(size(NucMaskSurf)));

                alOvFig = figure(1);
                imOv = subplot(2,1,1);
                imshow(ImOverlay, 'Parent', imOv)
                imOvMask = subplot(2,1,2);
                imshow(ImOverlayMask, 'Parent', imOvMask)
                saveas(alOvFig, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'AlignmentOverlay.tif']);
                
                
                
       
                %Show the correlation image, but crop it a little bit
                contFig = figure(2);
                contAxes = axes(contFig);
                %HG: Note that I changed the ranges here by two pixels at
                %least.
                lX = uint16((CRows+1)/2-RowsZoomResized+1);
                uX = uint16((CRows-1)/2+RowsZoomResized);
                lY = uint16((CColumns+1)/2-ColumnsZoomResized+1);
                uY = uint16((CColumns-1)/2+ColumnsZoomResized);
                
                contourf(contAxes,imresize(abs(C(lX:uX,lY:uY)),.5));
                saveas(contFig, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'AlignmentCorrelation.tif']);
            catch
                warning('Could not generate correlation image.')
            end
       
    
        end
        
        

       
    end
    
    % Need to figure out unit conversion here 
    ImageCenter=[SurfRows/2,SurfColumns/2];
    
    %This is for the acquisition/zoom image
    TopLeft=[ImageCenter(1)-Rows/ZoomRatio/2+ShiftRow,...
        ImageCenter(2)-Columns/ZoomRatio/2+ShiftColumn];
    BottomRight=[ImageCenter(1)+Rows/ZoomRatio/2+ShiftRow,...
        ImageCenter(2)+Columns/ZoomRatio/2+ShiftColumn];
    coordAZoom=(coordA-[TopLeft(2),TopLeft(1)])*ZoomRatio;
	coordPZoom=(coordP-[TopLeft(2),TopLeft(1)])*ZoomRatio;

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

    saveVars = [saveVars, 'APLength'];

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
    
    save([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'],'Particles','SpotFilter');
end


end
%% 



function [coordT, coordB, coordL, coordR] = RestrictOverlapWindow(C)
close all
coordL = 1;
coordR = size(C, 2);
coordT = 1;
coordB = size(C, 1);
APImageFig=figure;
apAx = axes(APImageFig);
%Now, do the correction
cc=1;

while (cc~='x')
    
    
    imagesc(C, 'Parent', apAx)
    %imshow(APImage,DisplayRange)
    
    axis image
    %axis off
    title({'left boundary (green), right boundary (red),', 'top boundary (blue), bottom boundary (yellow); original'})
    hold on
    
    if exist('coordL', 'var')
        xline(coordL,'g')
    end
    
    if exist('coordR', 'var')
        xline(coordR,'r')
    end
    
    if exist('coordT', 'var')
        yline(coordT,'b')
    end
    
    if exist('coordB', 'var')
        yline(coordB,'y')
    end
%     
%     try
%         plot(coordA(1),coordA(2),'g.','MarkerSize',20);
%     catch
%         %not sure what happened here.
%     end
%     
%     try
%         plot(coordP(1),coordP(2),'r.','MarkerSize',20);
%     catch
%         %not sure what happened here.
%     end
    

    
    hold off
    
    figure(APImageFig)
    ct=waitforbuttonpress;
    cc=get(APImageFig,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    
    if (ct~=0)&(cc=='c')        %Clear all AP information
        leftboundary=[];
        rightboundary=[];
        topboundary=[];
        bottomboundary=[];
    elseif (ct~=0)&(cc=='l')	%Select anterior end
        [coordL,coordLy]=ginputc(1,'Color',[1,1,1]);
    elseif (ct~=0)&(cc=='r')    %Select posterior end
        [coordR,coordRy]=ginputc(1,'Color',[1,1,1]);
    elseif (ct~=0)&(cc=='t')    %Select posterior end
        [coordTx,coordT]=ginputc(1,'Color',[1,1,1]);
    elseif (ct~=0)&(cc=='b')    %Select posterior end
        [coordBx,coordB]=ginputc(1,'Color',[1,1,1]);
    end
end

close all

end
