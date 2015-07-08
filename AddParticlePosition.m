function AddParticlePosition(varargin)

%First parameter should be the prefix. The other parameters can be:
%SkipAlignment
%ManualAlignment
%NoAP: Just add X and Y information

%V2: Changed this function to use a correlation in order to center the
%images.


%Get the relevant folders for this data set
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(varargin{1});


SkipAlignment=0;
ManualAlignment=0;
NoAP=0;
InvertHis=0;

close all

if ~isempty(varargin)
    Prefix=varargin{1};
    for i=2:length(varargin)
        switch varargin{i}
            case {'SkipAlignment'}
                display('Skipping alignment step')
                SkipAlignment=1;
            case {'ManualAlignment'}
                ManualAlignment=1;
            case {'NoAP'}
                NoAP=1;
        end
    end
else
    FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
    Dashes=strfind(FolderTemp,filesep);
    Prefix=FolderTemp((Dashes(end)+1):end);
end


if exist([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
    
    %Create the particle array. This is done so that we can support multiple
    %channels. Also figure out the number of channels
    if iscell('Particles')
        NChannels=length(Particles);
    else
        Particles={Particles};
        NChannels=1;
    end
    
    %Now, get the particle positions (if they're not there already). Notice
    %that the code pulls out the position information from fad. This is because
    %of historical reasons mostly.
    for ChN=1:NChannels
        for i=1:length(Particles{ChN})
            for j=1:length(Particles{ChN}(i).Frame)
                [x,y]=fad2xyzFit(Particles{ChN}(i).Frame(j),fad(ChN), 'addMargin'); 
                Particles{ChN}(i).xPos(j)=x(Particles{ChN}(i).Index(j));
                Particles{ChN}(i).yPos(j)=y(Particles{ChN}(i).Index(j));
            end
        end
        if isfield(Particles{ChN},'APpos')
            warning('Particles.mat already has AP positions stored. They will be rewritten')
        end
    end

else
    warning('No Particles.mat found. Just updating APDetection.mat')
end


%See if we had any lineage/nuclear information
D=dir([PreProcPath,filesep,Prefix,filesep,'*-His_*']);
if length(D)>0
    HistoneChannel=1;
else
    HistoneChannel=0;
end


%Determine whether we're dealing with 2-photon data from Princeton or LSM
%data. 2-photon data uses TIF files. In LSM mode multiple files will be
%combined into one.
%Find out the date it was taken
Dashes=findstr(Prefix,'-');
Date=Prefix(1:Dashes(3)-1);
EmbryoName=Prefix(Dashes(3)+1:end);

DTIF=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'*.tif']);
DLSM=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'*.lsm']);
DLIF=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'*.lif']);
DLAT=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'IsLatticeData.txt']);

if (length(DTIF)>0)&(length(DLSM)==0)
    if length(DLIF)==0
        if length(DLAT)==0
            display('2-photon @ Princeton data mode')
            D=DTIF;
            FileMode='TIF';
        else
            display('Lattice Light Sheet data mode')
            D=DTIF;
            FileMode='LAT';
        end
    else
        display('LIF export mode')
        D=DTIF;
        FileMode='LIFExport';
    end
elseif (length(DTIF)==0)&(length(DLSM)>0)
    display('LSM mode')
    D=DLSM;
    FileMode='LSM';
else
    error('File type not recognized')
end




%Figure out what type of experiment we have
[XLSNum,XLSTxt]=xlsread([DropboxFolder,filesep,'MovieDatabase.xlsx']);
DataFolderColumn=find(strcmp(XLSTxt(1,:),'DataFolder'));
ExperimentTypeColumn=find(strcmp(XLSTxt(1,:),'ExperimentType'));
Channel1Column=find(strcmp(XLSTxt(1,:),'Channel1'));
Channel2Column=find(strcmp(XLSTxt(1,:),'Channel2'));

% Convert the prefix into the string used in the XLS file
Dashes = strfind(Prefix, '-');
PrefixRow = find(strcmp(XLSTxt(:, DataFolderColumn),...
    [Prefix(1:Dashes(3)-1), '\', Prefix(Dashes(3)+1:end)]));
if isempty(PrefixRow)
    PrefixRow = find(strcmp(XLSTxt(:, DataFolderColumn),...
        [Prefix(1:Dashes(3)-1), '/', Prefix(Dashes(3)+1:end)]));
    if isempty(PrefixRow)
        error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
    end
end

ExperimentType=XLSTxt(PrefixRow,ExperimentTypeColumn);
Channel1=XLSTxt(PrefixRow,Channel1Column);
Channel2=XLSTxt(PrefixRow,Channel2Column);


ChannelToLoadTemp=(~cellfun(@isempty,strfind({lower(Channel1{1}),lower(Channel2{1})},'mcherry'))|...
    ~cellfun(@isempty,strfind({lower(Channel1{1}),lower(Channel2{1})},'his')));


if sum(ChannelToLoadTemp)
    ChannelToLoad=find(ChannelToLoadTemp);
else
    error('No histone channel found. Was it defined in MovieDatabase.XLS?')
end



%Figure out how our field of view maps to the AP coordinates.

%First, figure out how our field of view maps to the stitched embryo
%image

%Were the images taken on the left or right half of the embryo?
if ~NoAP

    %Get information about all images. This depends on the microsocpe used.
    
    %Get the information about the zoom
    if strcmp(FileMode,'TIF')
        D=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'*.tif']);
        ImageInfo = imfinfo([SourcePath,filesep,Date,filesep,EmbryoName,filesep,D(1).name]);
        
        %Figure out the zoom factor
        MovieZoom=ExtractInformationField(ImageInfo(1),'state.acq.zoomFactor=');
        MovieZoom=str2num(MovieZoom);
    
        
        %Get the zoomed out surface image and its dimensions from the FullEmbryo folder
        D=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*.tif']);
        SurfName=D(find(~cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
        SurfImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,SurfName],ChannelToLoad); 
        
        %Get the size of the zoom image
        Rows = str2double(ExtractInformationField(ImageInfo(1), 'state.acq.linesPerFrame='));
        Columns = str2double(ExtractInformationField(ImageInfo(1), 'state.acq.pixelsPerLine='));
        
        SurfInfo = imfinfo([SourcePath, filesep, Date, filesep, EmbryoName, filesep, 'FullEmbryo', filesep, SurfName]);
        SurfZoom = ExtractInformationField(SurfInfo(1), 'state.acq.zoomFactor=');
        SurfZoom = str2double(SurfZoom);

        SurfRows = str2double(ExtractInformationField(SurfInfo(1), 'state.acq.linesPerFrame='));
        SurfColumns = str2double(ExtractInformationField(SurfInfo(1), 'state.acq.pixelsPerLine='));
        
        %Get the full emrbyo image
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
        
%     elseif strcmp(FileMode,'LSM') 
%         
%         D=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'*.lsm']);
%         ImageInfo = imfinfo([SourcePath,filesep,Date,filesep,EmbryoName,filesep,D(1).name]);
%         
%         %Figure out the zoom factor
%         MovieZoom=ExtractInformationField(ImageInfo(1),'state.acq.zoomFactor=');
%         MovieZoom=str2num(MovieZoom);
%     
%         
%         %Get the zoomed out surface image and its dimensions from the FullEmbryo folder
%         D=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*.tif']);
%         SurfName=D(find(~cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
%         SurfImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,SurfName],ChannelToLoad); 
%         
%         %Get the size of the zoom image
%         Rows = str2double(ExtractInformationField(ImageInfo(1), 'state.acq.linesPerFrame='));
%         Columns = str2double(ExtractInformationField(ImageInfo(1), 'state.acq.pixelsPerLine='));
%         
%         SurfInfo = imfinfo([SourcePath, filesep, Date, filesep, EmbryoName, filesep, 'FullEmbryo', filesep, SurfName]);
%         SurfZoom = ExtractInformationField(SurfInfo(1), 'state.acq.zoomFactor=');
%         SurfZoom = str2double(SurfZoom);
% 
%         SurfRows = str2double(ExtractInformationField(SurfInfo(1), 'state.acq.linesPerFrame='));
%         SurfColumns = str2double(ExtractInformationField(SurfInfo(1), 'state.acq.pixelsPerLine='));
%         
%         %Get the full embryo image
%         FullEmbryo=imread([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo.tif']);
%     
        
    elseif strcmp(FileMode,'LSM')|strcmp(FileMode,'LIFExport')
%     elseif strcmp(FileMode,'LIFExport')

        %This is so that the code doesn't freak out later
        SurfName=[];
        
        %Figure out the different channels
        if  ~isempty(strfind(lower(Channel1{1}),'his'))
            HisChannel=1;
            InvertHis=0;
        elseif ~isempty(strfind(lower(Channel1{1}),'mcherry'))
            HisChannel=1;
            InvertHis=1;
        elseif ~isempty(strfind(lower(Channel2{1}),'his'))
            HisChannel=2;
            InvertHis=0;
        elseif ~isempty(strfind(lower(Channel2{1}),'mcherry'))
            HisChannel=2;
            InvertHis=1;
        else
            error('LIF Mode error: Channel name not recognized. Check MovieDatabase.XLSX')
        end
                
        %Find the zoomed movie pixel size
        D=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'*.',FileMode(1:3)]);
        %Load only the metadata from the zoomed images
        MetaReader=bfGetReader([SourcePath,filesep,Date,filesep,EmbryoName,filesep,D(end).name]);
        MetaZoom=MetaReader.getMetadataStore();
        PixelSizeZoom=str2num(MetaZoom.getPixelsPhysicalSizeX(0));

        %Find the full embryo pixel size and load the image
        D=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*surf*.',FileMode(1:3)]);
        ImageTemp=bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,D(end).name]);
        MetaFullEmbryo= ImageTemp{:, 4};
        PixelSizeFullEmbryo=str2num(MetaFullEmbryo.getPixelsPhysicalSizeX(0));
        
        %How many channels and slices do we have?
        NChannelsMeta=MetaFullEmbryo.getChannelCount(0);
        NSlices=str2num(MetaFullEmbryo.getPixelsSizeZ(0));
        clear MaxTemp
  
        %Do a maximum projections

        %By looking at the last image we make sure we're avoiding the
        %individual tiles if we're dealing with tile scan
        for i=HisChannel:NChannelsMeta:size(ImageTemp{end,1},1)
                MaxTemp(:,:,i)=ImageTemp{end,1}{i,1};
        end
        if InvertHis
            SurfImage=MaxTemp(:,:,HisChannel-1+round(NSlices/2)*NChannelsMeta-1);
        else
            SurfImage=max(MaxTemp,[],3);
        end
        Dashes=findstr(Prefix,'-');
        Date=Prefix(1:Dashes(3)-1);
        EmbryoName=Prefix(Dashes(3)+1:end);
        
        %Rotates the full embryo image to match the rotation of the zoomed
        %time series
        zoom_angle = 0;
        full_embryo_angle = 0;

        if isdir([SourcePath, filesep, Date, filesep, EmbryoName, filesep, 'MetaData']) 
            xml_file_path = dir([SourcePath, filesep, Date, filesep, EmbryoName, filesep, 'MetaData', filesep, '*.xml']);
            xml_file = xml_file_path(1).name;
            xDoc = searchXML([SourcePath, filesep, Date, filesep, EmbryoName, filesep, 'MetaData', filesep, xml_file]);
            zoom_angle = str2double(evalin('base','rot'));
        else 
            warning('No time series metadata found.')
        end
        if isdir([SourcePath, filesep, Date, filesep, EmbryoName, filesep, 'FullEmbryo', filesep...
                'MetaData'])     
            xml_file_path2 = dir([SourcePath, filesep, Date, filesep, EmbryoName, filesep, 'FullEmbryo',...
                filesep, 'MetaData', filesep,'*Surf*.xml']);
            xml_file2 = xml_file_path2(1).name;
            xDoc2 = searchXML([SourcePath, filesep, Date, filesep, EmbryoName, filesep,'FullEmbryo', filesep,...
                    'MetaData', filesep, xml_file2]);
            full_embryo_angle = str2double(evalin('base','rot'));
        else 
            warning('No full embryo metadata found.')
        end
        
        evalin('base','clear rot')
        SurfImage = imrotate(SurfImage, -zoom_angle + full_embryo_angle);
        clear ImageTemp

        %Zoom factor
        MovieZoom=PixelSizeFullEmbryo(1)/PixelSizeZoom(1);
        SurfZoom=1;     %We'll call the zoom of the full embryo image 1
        
        %Get the size of the zoom image
        Columns = str2num(MetaZoom.getPixelsSizeX(0));
        Rows = str2num(MetaZoom.getPixelsSizeY(0));
        %Get the size of the full embryo image
%         SurfColumns=str2num(MetaFullEmbryo.getPixelsSizeX(0));
%         SurfRows=str2num(MetaFullEmbryo.getPixelsSizeY(0));
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

    if HistoneChannel
        ChannelToLoad=2;

        %Get the surface image in the zoomed case by looking at the last
        %frame of our movie
        DHis=dir([PreProcPath,filesep,Prefix,filesep,Prefix,'-His*.tif']);
        ZoomImage=imread([PreProcPath,filesep,Prefix,filesep,DHis(end).name]);
    else
        ChannelToLoad=1;

        %Get the surface image in the zoomed case
        DGFP=dir([PreProcPath,filesep,Prefix,filesep,Prefix,'*_z*.tif']);
        if ~isempty(DGFP)
            ZoomImage=imread([PreProcPath,filesep,Prefix,filesep,DGFP(end).name],ChannelToLoad);
        else
            % This might be the case, for instance, if you're just trying
            % to find AP information about an image without using FISH
            % code. In that case, just extract the nuclei from the last
            % raw image.
            DGFP = dir([SourcePath, filesep, Date, filesep, EmbryoName, filesep, '*.tif']);
            ImageInfo = imfinfo([SourcePath, filesep, Date, filesep, EmbryoName, filesep, DGFP(end).name]);
            NumFramesAndSlices = length(ImageInfo)/2;
            RawImage3M = NaN(Rows, Columns, NumFramesAndSlices);
            for lImageIndex = 1:NumFramesAndSlices
                RawImage3M(:, :, lImageIndex) = imread([SourcePath, filesep, Date, filesep, EmbryoName, filesep, DGFP(end).name],...
                    'Index', 2*(lImageIndex-1) + ChannelToLoad);
            end
            ZoomImage = max(RawImage3M, [], 3) / 255;
        end
    end

      

   
    
    %If there is no zoom information on the surface image then look into
    %the temp folder. This is because sometimes we edit images in ImageJ
    %which leads to losing the zoom information.
    if isnan(SurfZoom)
        Dtemp=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'temp',filesep,'*.tif']);
        LeftFileIndex=find(~cellfun('isempty',strfind(lower({Dtemp.name}),'left'))&...
            cellfun('isempty',strfind(lower({Dtemp.name}),'surf')));
        ImageInfo = imfinfo([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,...
            'temp',filesep,Dtemp(LeftFileIndex).name]);
        SurfZoom=str2double(ExtractInformationField(ImageInfo(1),'state.acq.zoomFactor='));
         
        SurfRows=SurfInfo.Height;
        SurfColumns=SurfInfo.Width;
    end
   
    
    
    %Do a correlation between the zoomed in and zoomed out surface images
    %to figure out the shift.
    
    FullEmbryo=imread([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo.tif']);
    if ~SkipAlignment && HistoneChannel
        if ZoomRatio > 1 && ZoomRatio < 24 
            
            %Enlarge the zoomed out image so we can do the cross-correlation
            SurfImageResized=imresize(SurfImage, ZoomRatio);

            %Calculate the correlation matrix and find the maximum
            if InvertHis
                warning('I still need to fix this part')
                C = normxcorr2(imcomplement(ZoomImage), SurfImageResized);
            else            
                C = normxcorr2(ZoomImage, SurfImageResized);
            end

            [Max2,MaxRows]=max(C);
            [Dummy,MaxColumn]=max(Max2);
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
            if exist('ManualAlignmentDone')
                if ManualAlignmentDone
                    display('Manual alignment results saved. Using them.')
                    load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'],'ShiftRow','ShiftColumn')
                end
            end
            
            
            %See if we need the manual alignment
            if ManualAlignment
                %See if we need to load the previous manual alignment results
                [ShiftColumn,ShiftRow]=ManualAPCorrection(SurfImage,ZoomImage,C,ZoomRatio,ShiftRow,ShiftColumn,...
                     FullEmbryo, ZoomRatio, SurfRows,Rows, Columns, coordA, coordP, SurfColumns);
                ManualAlignmentDone=1;
            end
            
            try
                %Make an overlay of the zoomed in and zoomed out real
                %images as well as of a quickly segmented nuclear mask
                
                %Real image overlay
                %Crop the zoomed out image to match the zoomed in one
                SurfImageResizeZoom=...
                    SurfImageResized(round(RowsResized/2-RowsZoom/2+ShiftRow*ZoomRatio):round(RowsResized/2+RowsZoom/2-1+ShiftRow*ZoomRatio),...
                    round(ColumnsResized/2-ColumnsZoom/2+ShiftColumn*ZoomRatio):round(ColumnsResized/2+ColumnsZoom/2-1+ShiftColumn*ZoomRatio));
                ImOverlay=cat(3,mat2gray(SurfImageResizeZoom),...
                    +mat2gray(ZoomImage),zeros(size(SurfImageResizeZoom)));

                %Nuclear mask overlay
                NucMaskZoomOut=GetNuclearMask(SurfImage,2.5,0);
                NucMaskZoomOutResized=imresize(NucMaskZoomOut, ZoomRatio);
                NucMaskZoomOutResizedCropped=...
                    NucMaskZoomOutResized(round(RowsResized/2-RowsZoom/2+ShiftRow*ZoomRatio):round(RowsResized/2+RowsZoom/2-1+ShiftRow*ZoomRatio),...
                    round(ColumnsResized/2-ColumnsZoom/2+ShiftColumn*ZoomRatio):round(ColumnsResized/2+ColumnsZoom/2-1+ShiftColumn*ZoomRatio));
                NucMaskZoomIn=GetNuclearMask(ZoomImage,8,2);
                ImOverlayMask=cat(3,mat2gray(NucMaskZoomOutResizedCropped),...
                    +mat2gray(NucMaskZoomIn),zeros(size(NucMaskZoomOutResizedCropped)));


                figure(1)
                subplot(2,1,1)
                imshow(ImOverlay)
                subplot(2,1,2)
                imshow(ImOverlayMask)
 
                
                saveas(gcf, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'AlignmentOverlay.tif']);

                %Show the correlation image, but crop it a little bit
                figure(2)
                contourf(abs(C((CRows-1)/2-RowsZoom:(CRows-1)/2+RowsZoom,...
                    (CColumns-1)/2-ColumnsZoom:(CColumns-1)/2+ColumnsZoom)))
                saveas(gcf, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'AlignmentCorrelation.tif']);
            catch
                warning('Could not generate correlation image. Switching to manual alignment')
                
                ManualAlignment=1;
                ShiftColumn=0;
                ShiftRow=0;
                NucMaskZoomIn = false(size(ZoomImage));
                NucMaskZoomOut = false(size(SurfImage));
            end

        else

            warning('Not able to do the cross correlation. Switching to manual alignment mode.')

            ManualAlignment=1;
            ShiftColumn=0;
            ShiftRow=0;
            NucMaskZoomIn = false(size(ZoomImage));
            NucMaskZoomOut = false(size(SurfImage));
        end
    else
        warning('Not able to do the cross correlation. Assuming no shift between surface-level and movie-level images.')
        
        ManualAlignment = 0;
        ShiftColumn=0;
        ShiftRow=0;
        NucMaskZoomIn = false(size(ZoomImage));
        NucMaskZoomOut = false(size(SurfImage));
    end
   
    %Now figure out how the shift of the zoomed in and zoomed out surface
    %images translates to the whole embryo image (which is most of the
    %times stitched). This is necessary to figure out the AP position.

    
    %For full embryo images stitched out of two separate images, we need to
    %look at each case: the zoom in version being on the right or on the
    %left.
    %
    %We'll overlay the zoomed out surface and mid saggital images to check we got
    %things right.


    %Patch the surface image to fit the full embryo image.
    PatchedSurfaceImage=zeros(size(FullEmbryo));
    

    %Are we dealing with three images?
    if ~isempty(findstr(lower(SurfName),'center'))
        %Is this a left-right orientation or top-down?
        if sum(~cellfun('isempty',strfind(lower({D.name}),'left'))) && sum(~cellfun('isempty',strfind(lower({D.name}),'right')))
        
            HalfName=D(find(sum(~cellfun('isempty',strfind(lower({D.name}),'center'))&~cellfun('isempty',strfind(lower({D.name}),'surf'))))).name;
            HalfImageSurf=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,HalfName],ChannelToLoad);
            [Rows1x,Columns1x]=size(HalfImageSurf);


            ImageCenter=[Rows1x/2,Columns1x/2];


            %This is for the half image
            TopLeftHalf = round([ImageCenter(1)-Rows/ZoomRatio/2+ShiftRow+1,...
                ImageCenter(2)-Columns/ZoomRatio/2+ShiftColumn+1]);
            BottomRightHalf = round([ImageCenter(1)+Rows/ZoomRatio/2+ShiftRow,...
                ImageCenter(2)+Columns/ZoomRatio/2+ShiftColumn]);
           
            %This is for the full image
            TopLeft=TopLeftHalf+[0,Columns1x-xShift1];
            BottomRight=BottomRightHalf+[0,Columns1x-xShift1];
            
            %AP position mapped onto the zoomed out image
            coordAHalf=coordA+[-Columns1x+xShift1,0];
            coordPHalf=coordP+[-Columns1x+xShift1,0];
            
            
        elseif sum(~cellfun('isempty',strfind(lower({D.name}),'top'))) && sum(~cellfun('isempty',strfind(lower({D.name}),'bottom')))
        
            HalfName=D(find(sum(~cellfun('isempty',strfind(lower({D.name}),'center'))&~cellfun('isempty',strfind(lower({D.name}),'surf'))))).name;
            HalfImageSurf=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,HalfName],ChannelToLoad);
            [Rows1x,Columns1x]=size(HalfImageSurf);

            %Get the imaging region
            ImageCenter=[Rows1x/2,Columns1x/2];
            %Imaged region mapped onto the zoomed out image
            TopLeftHalf = round([ImageCenter(1)-Rows/ZoomRatio/2+ShiftRow+1,...
                ImageCenter(2)-Columns/ZoomRatio/2+ShiftColumn+1]);
            BottomRightHalf = round([ImageCenter(1)+Rows/ZoomRatio/2+ShiftRow,...
                ImageCenter(2)+Columns/ZoomRatio/2+ShiftColumn]);
            
            %AP position mapped onto the zoomed out image. Careful, AP
            %coordinates are defined (x,y).
            coordAHalf=coordA+[0,-(Rows1x-yShift1)];
            coordPHalf=coordP+[0,-(Rows1x-yShift1)];
            
            %This is for the full image
            TopLeft=TopLeftHalf+[Rows1x-yShift1,0];
            BottomRight=BottomRightHalf+[Rows1x-yShift1,0];
            
        end
    
    %Only two images:
    %Are we dealing with a left-right orientation or with top-bottom?
    elseif ~isempty(findstr(lower(SurfName),'top'))|~isempty(findstr(lower(SurfName),'bottom'))
        %The information from the top-bottom stitching of the two images is as follows:
        %xShift and yShift are the shifts used to stitch the images.
        %yShift is the displacement of the bottom image with respect to the
        %top image. Positive yShift moves the bottom image up.
        %xShift is the displacement of the top image with respect to the right
        %image. Positive xShift moves the top image to the right. Note that if we're
        %aligning images on the bottom we don't need to worry about this in
        %terms of the overlap of the surface and mid images.
        
        %If the zoomed in image coincides with the bottom zoomed out image
        if sum(~cellfun('isempty',strfind(lower({D.name}),'bottom'))&~cellfun('isempty',strfind(lower({D.name}),'surf')))
        
       
            %Load the half image at the surface
            HalfNameSurf=D(find(~cellfun('isempty',strfind(lower({D.name}),'bottom'))&~cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
            HalfImageSurf=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,HalfNameSurf],ChannelToLoad);
            [Rows1x,Columns1x]=size(HalfImageSurf);

            TopMidImageName=D(find(~cellfun('isempty',strfind(lower({D.name}),'top'))&cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
            BottomMidImageName=D(find(~cellfun('isempty',strfind(lower({D.name}),'bottom'))&cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
            TopMidImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,TopMidImageName],ChannelToLoad);
            BottomMidImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,BottomMidImageName],ChannelToLoad);


            %This is for the half image
          
            %Get the imaging region
            ImageCenter=[Rows1x/2,Columns1x/2];
            %Imaged region mapped onto the zoomed out image
            TopLeftHalf = round([ImageCenter(1)-Rows/ZoomRatio/2+ShiftRow+1,...
                ImageCenter(2)-Columns/ZoomRatio/2+ShiftColumn+1]);
            BottomRightHalf = round([ImageCenter(1)+Rows/ZoomRatio/2+ShiftRow,...
                ImageCenter(2)+Columns/ZoomRatio/2+ShiftColumn]);
            
            %AP position mapped onto the zoomed out image. Careful, AP
            %coordinates are defined (x,y).
            coordAHalf=coordA+[0,-(Rows1x-yShift)];
            coordPHalf=coordP+[0,-(Rows1x-yShift)];
            
            %This is for the full image
            TopLeft=TopLeftHalf+[Rows1x-yShift,0];
            BottomRight=BottomRightHalf+[Rows1x-yShift,0];

    
        elseif sum(~cellfun('isempty',strfind(lower({D.name}),'top'))&...
                ~cellfun('isempty',strfind(lower({D.name}),'surf')))

            %Load the half image at the midsaggital plane
            %HalfName=D(find(sum(cellfun('isempty',strfind(lower({D.name}),'right'))&cellfun('isempty',strfind(lower({D.name}),'surf'))))).name;

            %Load the half image at the surface
            HalfNameSurf=D(find(~cellfun('isempty',strfind(lower({D.name}),'top'))&~cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
            HalfImageSurf=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,HalfNameSurf],ChannelToLoad);
            [Rows1x,Columns1x]=size(HalfImageSurf);


            TopMidImageName=D(find(~cellfun('isempty',strfind(lower({D.name}),'top'))&cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
            BottomMidImageName=D(find(~cellfun('isempty',strfind(lower({D.name}),'bottom'))&cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
            TopMidImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,TopMidImageName],ChannelToLoad);
            BottomMidImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,BottomMidImageName],ChannelToLoad);

            
            %This is for the half image

            %Get the imaging region
            ImageCenter=[Rows1x/2,Columns1x/2];
            %Imaged region mapped onto the zoomed out image
            TopLeftHalf = round([ImageCenter(1)-Rows/ZoomRatio/2+ShiftRow+1,...
                ImageCenter(2)-Columns/ZoomRatio/2+ShiftColumn+1]);
            BottomRightHalf = round([ImageCenter(1)+Rows/ZoomRatio/2+ShiftRow,...
                ImageCenter(2)+Columns/ZoomRatio/2+ShiftColumn]);


            %AP position mapped onto the zoomed out image
            coordAHalf=coordA-[xShift,0];
            coordPHalf=coordP-[xShift,0];


            %This is for the full image

            TopLeft=TopLeftHalf+[0,+xShift];
            BottomRight=BottomRightHalf+[0,+xShift];

            
            
            
            
        
        end    
            
        
    
    else
        %The information from the left-right stitching of the two images is as follows:
        %xShift and yShift are the shifts used to stitch the images.
        %xShift is the displacement of the right image with respect to the left
        %image. Positive xShift moves the right image towards the left.
        %yShift is the displacement of the left image with respect to the right
        %image. Positive yShift moves the left image up. Note that if we're
        %aligning images on the right we don't need to worry about this in
        %terms of the overlap of the surface and mid images.


        %If the zoomed in image coincides with the right zoomed out image
        if sum(~cellfun('isempty',strfind(lower({D.name}),'right'))&~cellfun('isempty',strfind(lower({D.name}),'surf')))

            %Load the half image at the midsaggital plane
            %HalfName=D(find(sum(cellfun('isempty',strfind(lower({D.name}),'right'))&cellfun('isempty',strfind(lower({D.name}),'surf'))))).name;

            %Load the half image at the surface
            HalfNameSurf=D(find(~cellfun('isempty',strfind(lower({D.name}),'right'))&~cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
            HalfImageSurf=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,HalfNameSurf],ChannelToLoad);
            [Rows1x,Columns1x]=size(HalfImageSurf);



            LeftMidImageName=D(find(~cellfun('isempty',strfind(lower({D.name}),'left'))&cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
            RightMidImageName=D(find(~cellfun('isempty',strfind(lower({D.name}),'right'))&cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
            LeftMidImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,LeftMidImageName],ChannelToLoad);
            RightMidImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,RightMidImageName],ChannelToLoad);


            %The code below is meant for troubleshooting. I basically used it
            %to figure out whether the alignment between the two half images
            %was correct.
            %What I found out is that the margin removal was causing issues. I
            %need to get back to it later.
            %close all

    %         %Stitch the full embryo and overlay with the surface image - This
    %         %seems to work!
    %         FullEmbryoStitch = mat2gray(imstitch(LeftMidImage,RightMidImage, xShift, yShift,[1 2]),[0,100]);
    %         PatchedSurfImage=zeros(size(FullEmbryoStitch));
    %         
    %         PatchedSurfImage(:,...
    %             Columns1x-xShift+1:2*Columns1x-xShift)=...
    %             mat2gray(HalfImageSurf,[0,100]);
    %        
    %         OverlayFullEmbryoStitch=cat(3,PatchedSurfImage==1,FullEmbryoStitch==1,zeros(size(FullEmbryoStitch)));
    %         
    %         
    %         
    %         
    %         %This one seems to coincide with what I get if I do the overlay in
    %         %ImageJ
    %         figure(3)
    %         imshow(OverlayFullEmbryoStitch)
    %         
    %         
    %         
    %         %Grab the FullEmbryo image and overlay it with the surface one
    %         
    %         PatchedMidImage=zeros(Rows1x*2,Columns1x*2);
    %         PatchedSurfImage=zeros(Rows1x*2,Columns1x*2);
    %         
    % 
    %         
    %         PatchedMidImage(Rows1x/2+1:size(FullEmbryo,1)+Rows1x/2,...
    %             Columns1x/2+1:size(FullEmbryo,2)+Columns1x/2)=mat2gray(FullEmbryo,[0,100]);
    %         PatchedSurfImage(Rows1x/2+1:Rows1x+Rows1x/2,...
    %             Columns1x*3/2+1-xShift:Columns1x*3/2-xShift+Columns1x)=...
    %             mat2gray(HalfImageSurf,[0,100]);
    %         
    %         OverlayFullEmbryo=cat(3,PatchedSurfImage==1,PatchedMidImage==1,zeros(Rows1x*2,Columns1x*2));
    %         
    %         figure(4)
    %         imshow(OverlayFullEmbryo)


            %This is for the half image

            %Get the imaging region
            ImageCenter=[Rows1x/2,Columns1x/2];
            %Imaged region mapped onto the zoomed out image
            TopLeftHalf = round([ImageCenter(1)-Rows/ZoomRatio/2+ShiftRow+1,...
                ImageCenter(2)-Columns/ZoomRatio/2+ShiftColumn+1]);
            BottomRightHalf = round([ImageCenter(1)+Rows/ZoomRatio/2+ShiftRow,...
                ImageCenter(2)+Columns/ZoomRatio/2+ShiftColumn]);
            %AP position mapped onto the zoomed out image
            coordAHalf=coordA+[-Columns1x+xShift,0];
            coordPHalf=coordP+[-Columns1x+xShift,0];


            %Start by overlaying the zoom in figure on top of the zoom out
            %figure. Also add the rectangle.

            %Create an overlay of the mask in zoom in and zoom out
            NucMaskZoomInOverlay=zeros(size(SurfImage));

            NucMaskZoomInOverlay(TopLeftHalf(1):BottomRightHalf(1),TopLeftHalf(2):BottomRightHalf(2))=...
                        imresize(NucMaskZoomIn,...
                        [length(TopLeftHalf(1):BottomRightHalf(1)), length(TopLeftHalf(2):BottomRightHalf(2))]);

            SurfOutMaskInOverlay=cat(3,imadjust(mat2gray(NucMaskZoomOut)),imadjust(mat2gray(NucMaskZoomInOverlay)),...
                zeros(size(SurfImage)));


            figure(5)
            imshow(SurfOutMaskInOverlay)
            hold on
            rectangle('Position',[TopLeftHalf([2,1]),BottomRightHalf([2,1])-TopLeftHalf([2,1])],'EdgeColor','r')
            plot(coordAHalf(1),coordAHalf(2),'.g','MarkerSize',30)
            plot(coordPHalf(1),coordPHalf(2),'.r','MarkerSize',30)
            plot([coordAHalf(1),coordPHalf(1)],[coordAHalf(2),coordPHalf(2)],'-b')
            hold off


            %This is for the full image

            TopLeft=TopLeftHalf+[0,Columns1x-xShift];
            BottomRight=BottomRightHalf+[0,Columns1x-xShift];

            figure(6)
            imshow(imadjust(mat2gray(FullEmbryo)),'DisplayRange',[],'InitialMagnification',100)
            hold on
            rectangle('Position',[TopLeft([2,1]),BottomRight([2,1])-TopLeft([2,1])],'EdgeColor','r')
            plot(coordA(1),coordA(2),'.g','MarkerSize',30)
            plot(coordP(1),coordP(2),'.r','MarkerSize',30)
            plot([coordA(1),coordP(1)],[coordA(2),coordP(2)],'-b')
            hold off


        elseif sum(~cellfun('isempty',strfind(lower({D.name}),'left'))&...
                ~cellfun('isempty',strfind(lower({D.name}),'surf')))

            %Load the half image at the midsaggital plane
            %HalfName=D(find(sum(cellfun('isempty',strfind(lower({D.name}),'right'))&cellfun('isempty',strfind(lower({D.name}),'surf'))))).name;

            %Load the half image at the surface
            HalfNameSurf=D(find(~cellfun('isempty',strfind(lower({D.name}),'left'))&~cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
            HalfImageSurf=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,HalfNameSurf],ChannelToLoad);
            [Rows1x,Columns1x]=size(HalfImageSurf);


            LeftMidImageName=D(find(~cellfun('isempty',strfind(lower({D.name}),'left'))&cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
            RightMidImageName=D(find(~cellfun('isempty',strfind(lower({D.name}),'right'))&cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
            LeftMidImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,LeftMidImageName],ChannelToLoad);
            RightMidImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,RightMidImageName],ChannelToLoad);

            %Stitch the full embryo and overlay with the surface image - This
            %seems to work!
            FullEmbryoStitch = mat2gray(imstitch(LeftMidImage,RightMidImage, xShift, yShift,[1 2]),[0,50]);
            PatchedSurfImage=zeros(size(FullEmbryoStitch));

            PatchedSurfImage(:,...
               1:Columns1x)=...
                mat2gray(circshift(HalfImageSurf,[-yShift,0]),[0,50]);



            OverlayFullEmbryoStitch=cat(3,PatchedSurfImage==1,FullEmbryoStitch==1,zeros(size(FullEmbryoStitch)));


            %The code below is meant for troubleshooting. I basically used it
            %to figure out whether the alignment between the two half images
            %was correct.
            %What I found out is that the margin removal was causing issues. I
            %need to get back to it later.
            %close all

            %This one seems to coincide with what I get if I do the overlay in
            %ImageJ
            figure(3)
            imshow(OverlayFullEmbryoStitch)



            %Grab the FullEmbryo image and overlay it with the surface one

            PatchedMidImage=zeros(Rows1x*2,Columns1x*2);
            PatchedSurfImage=zeros(Rows1x*2,Columns1x*2);



            PatchedMidImage(Rows1x/2+1:size(FullEmbryo,1)+Rows1x/2,...
                Columns1x/2+1:size(FullEmbryo,2)+Columns1x/2)=mat2gray(FullEmbryo,[0,50]);
            PatchedSurfImage(Rows1x/2+1-yShift:Rows1x+Rows1x/2-yShift,...
                Columns1x/2+1:Columns1x/2+Columns1x)=...
                mat2gray(HalfImageSurf,[0,50]);

            OverlayFullEmbryo=cat(3,PatchedSurfImage==1,PatchedMidImage==1,zeros(Rows1x*2,Columns1x*2));

            figure(4)
            imshow(OverlayFullEmbryo)


            %This is for the half image

            %Get the imaging region
            ImageCenter=[Rows1x/2,Columns1x/2];
            %Imaged region mapped onto the zoomed out image
            TopLeftHalf = round([ImageCenter(1)-Rows/ZoomRatio/2+ShiftRow+1,...
                ImageCenter(2)-Columns/ZoomRatio/2+ShiftColumn+1]);
            BottomRightHalf = round([ImageCenter(1)+Rows/ZoomRatio/2+ShiftRow,...
                ImageCenter(2)+Columns/ZoomRatio/2+ShiftColumn]);


            %AP position mapped onto the zoomed out image
            coordAHalf=coordA+[0,yShift];
            coordPHalf=coordP+[0,yShift];


            %Start by overlaying the zoom in figure on top of the zoom out
            %figure. Also add the rectangle.

            %Create an overlay of the mask in zoom in and zoom out
            NucMaskZoomInOverlay=zeros(size(SurfImage));

            NucMaskZoomInOverlay(TopLeftHalf(1):BottomRightHalf(1),TopLeftHalf(2):BottomRightHalf(2))=...
                        imresize(NucMaskZoomIn,...
                        [length(TopLeftHalf(1):BottomRightHalf(1)), length(TopLeftHalf(2):BottomRightHalf(2))]);

            SurfOutMaskInOverlay=cat(3,imadjust(mat2gray(NucMaskZoomOut)),imadjust(mat2gray(NucMaskZoomInOverlay)),...
                zeros(size(SurfImage)));


            figure(5)
            imshow(SurfOutMaskInOverlay)
            hold on
            rectangle('Position',[TopLeftHalf([2,1]),BottomRightHalf([2,1])-TopLeftHalf([2,1])],'EdgeColor','r')
            plot(coordAHalf(1),coordAHalf(2),'.g','MarkerSize',30)
            plot(coordPHalf(1),coordPHalf(2),'.r','MarkerSize',30)
            plot([coordAHalf(1),coordPHalf(1)],[coordAHalf(2),coordPHalf(2)],'-b')
            hold off


            %This is for the full image

            TopLeft=TopLeftHalf+[-yShift,0];
            BottomRight=BottomRightHalf+[-yShift,0];

            figure(6)
            imshow(imadjust(mat2gray(FullEmbryo)),'DisplayRange',[],'InitialMagnification',100)
            hold on
            rectangle('Position',[TopLeft([2,1]),BottomRight([2,1])-TopLeft([2,1])],'EdgeColor','r')
            plot(coordA(1),coordA(2),'.g','MarkerSize',30)
            plot(coordP(1),coordP(2),'.r','MarkerSize',30)
            plot([coordA(1),coordP(1)],[coordA(2),coordP(2)],'-b')
            hold off

        else
            %error('Problem with the surface file (or its naming) in the source data folder "FullEmbryo"')
        end
    end



    %Plot the area where we imaged on top of the embryo

    %Check if the embryo could actually fit in one of the images. If that's the
    %case we need to shift the the AP poisitions and boxes.
    if strcmp(FileMode,'TIF')
        
        
        
        if sum(size(HalfImageSurf)==size(FullEmbryo))==2

            warning('Have HG check this part of the code')

            ImageCenter=[Rows1x/2,Columns1x/2];

            %This is for the full image
            TopLeft=[ImageCenter(1)-Rows/ZoomRatio/2,ImageCenter(2)-Columns/ZoomRatio/2]...
                -[yShift,xShift];
            BottomRight=[ImageCenter(1)+Rows/ZoomRatio/2,ImageCenter(2)+Columns/ZoomRatio/2]...
                -[yShift,xShift];

            %This is for the acquisition image
            TopLeftHalf=[ImageCenter(1)-Rows/ZoomRatio/2+ShiftRow,...
                ImageCenter(2)-Columns/ZoomRatio/2+ShiftColumn];
            BottomRightHalf=[ImageCenter(1)+Rows/ZoomRatio/2+ShiftRow,...
                ImageCenter(2)+Columns/ZoomRatio/2+ShiftColumn];

            coordAHalf=coordA+[xShift,yShift];
            coordPHalf=coordP+[xShift,yShift];


        end
    else
        ImageCenter=[SurfRows/2,SurfColumns/2];

        %This is for the acquisition image
        TopLeftHalf=[ImageCenter(1)-Rows/ZoomRatio/2+ShiftRow,...
            ImageCenter(2)-Columns/ZoomRatio/2+ShiftColumn];
        BottomRightHalf=[ImageCenter(1)+Rows/ZoomRatio/2+ShiftRow,...
            ImageCenter(2)+Columns/ZoomRatio/2+ShiftColumn];
        
        %HG: It's silly to be defining things using the Half coordinates. I
        %need to do better here
        TopLeft=TopLeftHalf;
        BottomRight=BottomRightHalf;
        

        coordAHalf=coordA;
        coordPHalf=coordP;
    end


    figure(7)
    imshow(imadjust(mat2gray(FullEmbryo)),'DisplayRange',[],'InitialMagnification',100)
    hold on
    rectangle('Position',[TopLeft([2,1]),BottomRight([2,1])-TopLeft([2,1])],'EdgeColor','r')
    plot(coordA(1),coordA(2),'.g','MarkerSize',30)
    plot(coordP(1),coordP(2),'.r','MarkerSize',30)
    plot([coordA(1),coordP(1)],[coordA(2),coordP(2)],'-b')
    hold off
    saveas(gcf, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryoArea.tif']);

%     figure(8)
%     imshow(imadjust(SurfImage),'DisplayRange',[],'InitialMagnification',100)
%     hold on
%     rectangle('Position',[TopLeftHalf([2,1]),BottomRightHalf([2,1])-TopLeftHalf([2,1])],'EdgeColor','r')
%     plot(coordAHalf(1),coordAHalf(2),'.g','MarkerSize',30)
%     plot(coordPHalf(1),coordPHalf(2),'.r','MarkerSize',30)
%     plot([coordAHalf(1),coordPHalf(1)],[coordAHalf(2),coordPHalf(2)],'-b')
%     plot([1],[1],'.y','MarkerSize',50)
%     hold off
%     saveas(gcf, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'HalfEmbryoArea.tif']);
%     

    


    if strcmp(FileMode,'TIF')
        figure(8)
        imshow(imadjust(SurfImage),'DisplayRange',[],'InitialMagnification',100)
        hold on
        rectangle('Position',[TopLeftHalf([2,1]),BottomRightHalf([2,1])-TopLeftHalf([2,1])],'EdgeColor','r')
        plot(coordAHalf(1),coordAHalf(2),'.g','MarkerSize',30)
        plot(coordPHalf(1),coordPHalf(2),'.r','MarkerSize',30)
        plot([coordAHalf(1),coordPHalf(1)],[coordAHalf(2),coordPHalf(2)],'-b')
        plot([1],[1],'.y','MarkerSize',50)
        hold off
        saveas(gcf, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'HalfEmbryoArea.tif']);
        
        %We have the position of the anterior and posterior in the coordinates of the 
        %field of view we took. Do the mapping of the imaging region with respect to the AP axis.

        %Convert them to the corresponding zoom factor. In order to do this I
        %need to measure A and P with respect to the center of the rectangle.
        %Notice that there might be a +/-1 issue in the positioning given that
        %Matlab starts indexing from 1. However, as shown in the images below,
        %this doesn't make any real difference.   
        coordAZoom=(coordAHalf-[TopLeftHalf(2),TopLeftHalf(1)])*ZoomRatio;
        coordPZoom=(coordPHalf-[TopLeftHalf(2),TopLeftHalf(1)])*ZoomRatio;
    else
        figure(8)
        imshow(imadjust(SurfImage),'DisplayRange',[],'InitialMagnification',100)
        hold on
        rectangle('Position',[TopLeft([2,1]),BottomRight([2,1])-TopLeft([2,1])],'EdgeColor','r')
        plot(coordA(1),coordA(2),'.g','MarkerSize',30)
        plot(coordP(1),coordP(2),'.r','MarkerSize',30)
        plot([coordA(1),coordP(1)],[coordA(2),coordP(2)],'-b')
        plot([1],[1],'.y','MarkerSize',50)
        hold off
        saveas(gcf, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'HalfEmbryoArea.tif'])
        
        coordAZoom=(coordA-[TopLeft(2),TopLeft(1)])*ZoomRatio;
        coordPZoom=(coordP-[TopLeft(2),TopLeft(1)])*ZoomRatio;
    end
    

 
    

    
    
    figure(9)
    imshow(imadjust(ZoomImage),[])
    %imshow(NucMaskZoomIn)
    %imshow(NucMaskZoomOutResizedCropped)
    hold on
    plot([coordAZoom(1),coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'-b')
    plot([coordAZoom(1)+1,coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'--r')
    plot([coordAZoom(1)-1,coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'-.g')
    hold off
    saveas(gcf, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'ZoomedEmbryoAP.tif']);


    %With AP coordinates in hand we can now determine the AP position of
    %all particles. Look I my notes in "Calculating AP positions" in Notability
    %for details of the calculation.



    %Angle between the x-axis and the AP-axis
    APAngle=atan((coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1)));
    %Correction for if APAngle is in quadrants II or III
    if coordPZoom(1)-coordAZoom(1) < 0
        APAngle = APAngle + pi;
    end
    
    
    APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);


    APPosImage=zeros(size(ZoomImage));
    [Rows,Columns]=size(ZoomImage);

    for i=1:Rows
        for j=1:Columns
            Angle=atan((i-coordAZoom(2))./(j-coordAZoom(1)));
            if j-coordAZoom(1) < 0
                Angle = Angle + pi;
            end
            % Correction for if Angle is in quadrants II or III
            
            Distance=sqrt((coordAZoom(2)-i).^2+(coordAZoom(1)-j).^2);
            APPosition=Distance.*cos(Angle-APAngle);
            APPosImage(i,j)=APPosition/APLength;
        end
    end
    [XLSNum,XLSTxt,XLSRaw]=xlsread([DropboxFolder,filesep,'MovieDatabase.xlsx']);
    APResolutionColumn = find(strcmp(XLSRaw(1,:),'APResolution'));
    APResolution = XLSRaw{PrefixRow,APResolutionColumn};

    
    APbinID=0:APResolution:1;


    APPosBinImage=zeros(size(APPosImage));
    for i=1:(length(APbinID)-1)
        FilteredMask=(APbinID(i)<=APPosImage)&(APbinID(i+1)>APPosImage);

        APPosBinImage=APPosBinImage+FilteredMask*i;
    end

    ZoomOverlay=cat(3,mat2gray(ZoomImage)/2+mat2gray(APPosBinImage)/2,...
        mat2gray(ZoomImage)/2,mat2gray(ZoomImage)/2);

    figure(10)
    imshow(ZoomOverlay)
    hold on
    plot([coordAZoom(1),coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'-b')
    plot([coordAZoom(1)+1,coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'--r')
    plot([coordAZoom(1)-1,coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'-.g')
    hold off


    if exist([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
        for ChN=1:NChannels
            for i=1:length(Particles{ChN})
                %Angle between the x-axis and the particle using the A position as a
                %zero
                Angles=atan((Particles{ChN}(i).yPos-coordAZoom(2))./(Particles{ChN}(i).xPos-coordAZoom(1)));
                if Particles{ChN}(i).xPos-coordAZoom(1) < 0
                    Angles = Angles + pi;
                end
                % Correction for if Angles is in quadrants II or III

                %Distance between the points and the A point
                Distances=sqrt((coordAZoom(2)-Particles{ChN}(i).yPos).^2+(coordAZoom(1)-Particles{ChN}(i).xPos).^2);
                APPositions=Distances.*cos(Angles-APAngle);
                Particles{ChN}(i).APpos=APPositions/APLength;

                %Determine the distance perpendicular to the AP axis. This is a
                %proxy for a DV axis.

                Particles{ChN}(i).DVpos=Distances.*sin(Angles-APAngle);

            end
        end
    end




    
    
    %Save AP detection information
    if exist('xShift')
        save([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'],'coordA','coordP',...
            'xShift','yShift','coordAZoom','coordPZoom') 
    elseif  exist('xShift1')
        save([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'],'coordA','coordP',...
            'xShift1','yShift1','xShift2','yShift2','coordAZoom','coordPZoom') 
    else
        save([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'],'coordA','coordP',...
            'coordAZoom','coordPZoom') 
    end

    if ManualAlignment
        if ManualAlignmentDone
            save([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'],'ManualAlignmentDone',...
                'ShiftColumn','ShiftRow','-append')
        end
    end
    

end




%Save particle information
if exist([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
    
    %Bring the one channel case back to the legacy setting
    if NChannels==1
        Particles=Particles{1};
    end
    
    if exist('Threshold1')
        save([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'],'Particles','fad','fad2',...
            'Threshold1','Threshold2');
    else
        save([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'],'Particles','fad','fad2');
    end
end
    
 


