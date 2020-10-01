function FindAPAxisTile_deprecated(Prefix, varargin)
% author: Gabriella Martini
% date created: 12/28/19
% date last modified: 7/28/20

%% SOMETHING SEEMS TO HAVE GONE WRONG WITH THIS FILE.

% DESCRIPTION
% Locates particles from a zoomed-in movie within full embryo images using
% spatial cross-correlation.
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
%
% OPTIONS
% 'NoAP': Just adds X and Y information % Add support for this
% 'SelectChannel': Prompts user to select the channel to use for alignment
% 'optionalResults': 
% 'yToManualAlignmentPrompt': 
% 'correctDV': If you want the DV shift calculated
% UseFullEmbryoBoundaries = false;
% ImposeEmbryoMask = false; % Note that this won't work if a lot of zoomImage is outside of embryo. Use embryo_mask_tolerance to adjust how rigid this feature is.
% embryo_mask_tolerance = .9; %Change value using varargin
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

%Default set of variables to save
saveVars={'coordA','coordP','coordAZoom','coordPZoom'};


SkipAlignment=false;
ManualAlignment=false;
NoAP=false;
SelectChannel=0;
InvertHis=false;
optionalResults = '';
yToManualAlignmentPrompt = false;
correctDV = false;
UseFullEmbryoBoundaries = false;
ImposeEmbryoMask = false; % Note that this won't work if a lot of zoomImage is outside of embryo. Use embryo_mask_tolerance to adjust how rigid this feature is.
embryo_mask_tolerance = .9; %Change value using varargin
%ImposeManualBounds = false;

close all



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
        if isfield(Particles{ChN},'DVpos')
            warning('Particles.mat already has DV positions stored. They will be rewritten')
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
    
    
    %Get the information about the zoom
    if strcmp(FileMode,'TIF')
       msg = '.tif files not supported.';
       error(msg)
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
             msg = 'DSPIN files not supported.';
             error(msg)
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
            msg = 'DSPIN files not supported.';
            error(msg)
        end
        
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
        
        if strcmp(FileMode,'LIFExport')
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
        elseif strcmp(FileMode,'LSM')|strcmp(FileMode,'CZI')|strcmp(FileMode,'DSPIN') %CS20170912
           msg='No Support for these files types';
           error(msg)
        end
        %%
        SurfImage = imrotate(SurfImage, -zoom_angle + full_embryo_angle);
        
        mkdir([DropboxFolder,filesep,Prefix, filesep, 'DV']);
        maxSurfSavePath = [DropboxFolder,filesep,Prefix, filesep, 'DV', filesep, 'surf_max.tif'];
        imwrite(SurfImage,maxSurfSavePath);
        
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
    end
    
    
    
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
        % This might be the case, for instance, if you're just trying
        % to find AP information about an image without using FISH
        % code. In that case, just extract the nuclei from the last
        % raw image.
        % DGFP = dir([rawPrefixPath, '*.tif']);
        % ImageInfo = imfinfo([rawPrefixPath, DGFP(end).name]);
        % NumFramesAndSlices = length(ImageInfo)/2;
        % RawImage3M = NaN(Rows, Columns, NumFramesAndSlices);
        % for lImageIndex = 1:NumFramesAndSlices
        %     RawImage3M(:, :, lImageIndex) = imread([rawPrefixPath,DGFP(end).name],'Index', 2*(lImageIndex-1) + ChannelToLoad);
        % end
        % ZoomImage = max(RawImage3M, [], 3) / 256;
    end
    

    
    %Do a correlation between the zoomed in and zoomed out surface images
    %to figure out the shift.
    
    % This is the mid-saggital image
    FullEmbryo=imread([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo.tif']);
    
    if ~SkipAlignment && histoneChannelPresent
        
        ShiftColumn=0;
        ShiftRow=0;
        NucMaskZoomIn = false(size(ZoomImage));
        NucMaskZoomOut = false(size(SurfImage));
        
        
        if ZoomRatio <= 1
            ZoomImageResized = imresize(ZoomImage, 1/ZoomRatio);
            im1 = ZoomImageResized;
            im2 = SurfImage;
            C = normxcorr2(im1, im2);
            [CRows,CColumns]=size(C);
            [RowsSurf,ColumnsSurf]=size(SurfImage);
            [RowsZoomResized,ColumnsZoomResized]=size(ZoomImageResized);
            %             end
            if UseFullEmbryoBoundaries
                row_subset = RowsZoomResized:(CRows-RowsZoomResized);
                col_subset = ColumnsZoomResized:(CColumns-ColumnsZoomResized);
                Csub = C(row_subset, col_subset);
                [Max2sub,MaxRowsSub]=max(Csub);
                [~,MaxColumnSub]=max(Max2sub);
                MaxRowSub=MaxRowsSub(MaxColumnSub);
            elseif ImposeEmbryoMask
                load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo')
                % change to be flexible nucleus diameter 
                nucleusDiameter = getDefaultParameters(FrameInfo,['d14']);
                im2_padded = zeros(SurfRows+10-rem(SurfRows, 10), SurfColumns+10-rem(SurfColumns, 10), 'uint16');
                im2_padded(1:SurfRows, 1:SurfColumns) = im2;
                I_resize = imresize(im2_padded, .1);
                f_sigma_resize = round(nucleusDiameter /( PixelSizeFullEmbryoSurf*10));
                I_blurred_resize = imfilter(I_resize,...
                     fspecial('gaussian',2*f_sigma_resize,f_sigma_resize),'symmetric','conv');
                embryoMask_resized = imbinarize(I_blurred_resize);

                SubRowSize = (RowsZoomResized-rem(RowsZoomResized, 10))/10;
                SubColumnSize = (ColumnsZoomResized-rem(ColumnsZoomResized, 10))/10;
                temp_mask = zeros(size(embryoMask_resized, 1)-SubRowSize+1, size(embryoMask_resized, 2)-SubColumnSize+1);

                ZoomImageSize = SubRowSize*SubColumnSize;
                tolerance = .95;
                for r=1:size(temp_mask, 1)
                    %display(['r: ', num2str(r)]);
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
            end
            
            
            [Max2,MaxRows]=max(C);
            [~,MaxColumn]=max(Max2);
            MaxRow=MaxRows(MaxColumn);
            
            
            %This shift kept in the zoom in coordinates (full embryo image coordinates).
            %If we want to translate to the zoomed out coordinates we need to
            %divide (I think!) again by ZoomRatio.
            ShiftRow=round((MaxRow-(CRows-1)/2));
            ShiftColumn=round((MaxColumn-(CColumns-1)/2));
            
            %How well did we do with the alignment?
            
            
            
            %If we can, we'll crop the surface image to center the overlay
            %on the area of acquisition. If not, we'll just show both images
            %without cropping.
            
            if ((round(RowsSurf/2-RowsZoomResized/2+ShiftRow))>0)&...
                    (round(RowsSurf/2+RowsZoomResized/2-1+ShiftRow)<=RowsSurf)&...
                    (round(ColumnsSurf/2-ColumnsZoomResized/2+ShiftColumn)>0)&...
                    (round(ColumnsSurf/2+ColumnsZoomResized/2-1+ShiftColumn)<=ColumnsSurf)
                RowsSurfRange=...
                    round(RowsSurf/2-RowsZoomResized/2+ShiftRow):...
                    round(RowsSurf/2+RowsZoomResized/2-1+ShiftRow);
                ColumnsSurfRange=...
                    round(ColumnsSurf/2-ColumnsZoomResized/2+ShiftColumn):...
                    round(ColumnsSurf/2+ColumnsZoomResized/2-1+ShiftColumn);
            else
                RowsSurfRange=1:RowsSurf;
                ColumnsSurfRange=1:ColumnsSurf;
            end
            
            
            %Make an overlay of the zoomed in and zoomed out real
            %images as well as of a quickly segmented nuclear mask. If this
            %fails, we'll swtich to ManualAlignment.
            
            try
                %Real image overlay
                %Crop the zoomed out image to match the zoomed in one
                SurfImageCrop=...
                    SurfImage(RowsSurfRange,ColumnsSurfRange);
                ImOverlay=cat(3,mat2gray(SurfImageCrop),...
                    +mat2gray(ZoomImageResized),zeros(size(SurfImageCrop)));
                
                load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo')
                nucleusDiameter = getDefaultParameters(FrameInfo,['d14']);
                NucMaskPixelSize = min([PixelSizeFullEmbryoMid, PixelSizeZoom]);
                %Nuclear mask overlay
                NucMaskSurf= GetNuclearMaskTile(SurfImageCrop,nucleusDiameter,NucMaskPixelSize, ImposeEmbryoMask);
                NucMaskZoomResized=GetNuclearMaskTile(ZoomImageResized,nucleusDiameter, NucMaskPixelSize,ImposeEmbryoMask);
                ImOverlayMask=cat(3,mat2gray(NucMaskSurf),...
                    +mat2gray(NucMaskZoomResized),zeros(size(NucMaskSurf)));

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
                warning('Could not generate correlation image. Switching to manual alignment')
                ManualAlignment=true;
            end
        elseif ZoomRatio < 24  %ZoomRatio > 1 && ZoomRatio < 24. AR 12/4/17- where did this number come from
            
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
            
            % Smaller image is first variable
            C = normxcorr2(im1, im2);
            %             end
            
            [Max2,MaxRows]=max(C);
            [~,MaxColumn]=max(Max2);
            MaxRow=MaxRows(MaxColumn);
            [CRows,CColumns]=size(C);
            
            
            %This shift is now converted to the zoom out distances. If we
            %want to translate to the zoomed in coordinates we need to
            %multiply again by ZoomRatio.
            ShiftRow=round((MaxRow-(CRows-1)/2)/ZoomRatio);
            ShiftColumn=round((MaxColumn-(CColumns-1)/2)/ZoomRatio);
            
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
            %fails, we'll swtich to ManualAlignment.
            
            try
                %Real image overlay
                %Crop the zoomed out image to match the zoomed in one
                SurfImageResizeCrop=...
                    SurfImageResized(RowsResizedRange,ColumnsResizedRange);
                ImOverlay=cat(3,mat2gray(SurfImageResizeZoom),...
                    +mat2gray(ZoomImage),zeros(size(SurfImageResizeZoom)));
                load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo')
                nucleusDiameter = getDefaultParameters(FrameInfo,['d14']);
                NucMaskPixelSize = min([PixelSizeFullEmbryoMid, PixelSizeZoom]);
                %Nuclear mask overlay
                NucMaskSurfResized= GetNuclearMaskTile(SurfImageResizeCrop,nucleusDiameter,NucMaskPixelSize, ImposeEmbryoMask);
                NucMaskZoom=GetNuclearMaskTile(ZoomImage,nucleusDiameter, NucMaskPixelSize,ImposeEmbryoMask);
                ImOverlayMask=cat(3,mat2gray(NucMaskZoomOutResizedCropped),...
                    +mat2gray(NucMaskZoomIn),zeros(size(NucMaskZoomOutResizedCropped)));

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
    
    if ZoomRatio < 1
        % Need to figure out unit conversion here 
        ImageCenter=[SurfRows/2,SurfColumns/2];

        %This is for the acquisition/zoom image
        TopLeft=[ImageCenter(1)-RowsZoomResized/2+ShiftRow,...
            ImageCenter(2)-ColumnsZoomResized/2+ShiftColumn];
        BottomRight=[ImageCenter(1)+RowsZoomResized/2+ShiftRow,...
            ImageCenter(2)+ColumnsZoomResized/2+ShiftColumn];
        coordAZoom=(coordA-[TopLeft(2),TopLeft(1)])*ZoomRatio;
        coordPZoom=(coordP-[TopLeft(2),TopLeft(1)])*ZoomRatio;
        if exist('coordD', 'var')
            coordDZoom=(coordD-[TopLeft(2),TopLeft(1)])*ZoomRatio;
            coordVZoom=(coordV-[TopLeft(2),TopLeft(1)])*ZoomRatio;
        end
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
    
    save([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'],'Particles','SpotFilter');
end

if correctDV
    CheckDivisionTimes(Prefix, 'lazy');
    %for convenience. 
end
%end
