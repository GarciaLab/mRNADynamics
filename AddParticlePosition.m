function AddParticlePosition(varargin)

%First parameter should be the prefix. The other parameters can be:
%SkipAlignment
%ManualAlignment
%NoAP: Just add X and Y information

%V2: Changed this function to use a correlation in order to center the
%images.


%Get the relevant folders for this data set
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders(varargin{1});


SkipAlignment=0;
ManualAlignment=0;
NoAP=0;

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
    Dashes=strfind(FolderTemp,'\');
    Prefix=FolderTemp((Dashes(end)+1):end);
end


if exist([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
    %Now, get the particle positions (if they're not there already). Notice
    %that the code pulls out the position information from fad. This is because
    %of historical reasons mostly.
    for i=1:length(Particles)
        for j=1:length(Particles(i).Frame)
            [x,y]=fad2xyzFit(Particles(i).Frame(j),fad, 'addMargin'); 
            Particles(i).xPos(j)=x(Particles(i).Index(j));
            Particles(i).yPos(j)=y(Particles(i).Index(j));
        end
    end
    if isfield(Particles,'APpos')
        warning('Particles.mat already has AP positions stored. They will be rewritten')
    end

else
    warning('No Particles.mat found. Just updating APDetection.mat')
end


%See if we had any lineage/nuclear information
D=dir([FISHPath,filesep,'Data',filesep,Prefix,filesep,'*-His_*']);
if length(D)>0
    HistoneChannel=1;
else
    HistoneChannel=0;
end




%Figure out how our field of view maps to the AP coordinates.

%First, figure out how our field of view maps to the stitched embryo
%image

%Were the images taken on the left or right half of the embryo?
if ~NoAP

    %Get information about all images
    
    %Find out the date it was taken
    Dashes=findstr(Prefix,'-');
    Date=Prefix(1:Dashes(3)-1);
    EmbryoName=Prefix(Dashes(3)+1:end);

    D=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'*.tif']);

    %Get the information about the zoom
    ImageInfo = imfinfo([SourcePath,filesep,Date,filesep,EmbryoName,filesep,D(1).name]);

    %Get the information about the AP axis as well as the image shifts
    %used for the stitching of the two halves of the embryo
    load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])

    %See if manual alignment was performed on this set. If so we'll skip the
    %automated alignment
    if exist('ManualAlignment')
        display('Manual alignment results saved. Using them.')
        ManualAlignmentDone=1;
    else
        ManualAlignmentDone=0;
    end

    %Figure out the zoom factor
    MovieZoom=ExtractInformationField(ImageInfo(1),'state.acq.zoomFactor=');
    MovieZoom=str2num(MovieZoom);
    
    %Get the size of the zoom image
    Rows = str2double(ExtractInformationField(ImageInfo(1), 'state.acq.linesPerFrame='));
    Columns = str2double(ExtractInformationField(ImageInfo(1), 'state.acq.pixelsPerLine='));

    %Make a folder to store the images
    mkdir([DropboxFolder,filesep,Prefix,filesep,'APDetection'])

    if HistoneChannel
        ChannelToLoad=2;

        %Get the surface image in the zoomed case by looking at the last
        %frame of our movie
        D=dir([FISHPath,filesep,'Data',filesep,Prefix,filesep,Prefix,'-His*.tif']);
        ZoomImage=imread([FISHPath,filesep,'Data',filesep,Prefix,filesep,D(end).name]);
    else
        ChannelToLoad=1;

        %Get the surface image in the zoomed case
        D=dir([FISHPath,filesep,'Data',filesep,Prefix,filesep,Prefix,'*_z*.tif']);
        ZoomImage=imread([FISHPath,filesep,'Data',filesep,Prefix,filesep,D(end-10).name],ChannelToLoad);
    end

    %Get the zoomed out surface image and its dimensions from the FullEmbryo folder
    D=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*.tif']);
    SurfName=D(find(~cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
    SurfImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,SurfName],ChannelToLoad);   

    SurfInfo = imfinfo([SourcePath, filesep, Date, filesep, EmbryoName, filesep, '\FullEmbryo\', SurfName]);
    SurfZoom = ExtractInformationField(SurfInfo(1), 'state.acq.zoomFactor=');
    SurfZoom = str2double(SurfZoom);
    
    SurfRows = str2double(ExtractInformationField(SurfInfo(1), 'state.acq.linesPerFrame='));
    SurfColumns = str2double(ExtractInformationField(SurfInfo(1), 'state.acq.pixelsPerLine='));
    
    
    %If there is no zoom information on the surface image then look into
    %the temp folder. This is because sometimes we edit images in ImageJ
    %which leads to losing the zoom information.
    if isnan(SurfZoom)
        Dtemp=dir([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo\temp\*.tif']);
        LeftFileIndex=find(~cellfun('isempty',strfind(lower({Dtemp.name}),'left'))&...
            cellfun('isempty',strfind(lower({Dtemp.name}),'surf')));
        ImageInfo = imfinfo([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo\temp',filesep,Dtemp(LeftFileIndex).name]);
        SurfZoom=str2double(ExtractInformationField(ImageInfo(1),'state.acq.zoomFactor='));
         
        SurfRows=SurfInfo.Height;
        SurfColumns=SurfInfo.Width;
    end
   
    %Ratio between the two zoom levels
    ZoomRatio = MovieZoom / SurfZoom;

    
    %Do a correlation between the zoomed in and zoomed out surface images
    %to figure out the shift.
    
    
    if ~SkipAlignment&HistoneChannel
        if ZoomRatio > 1 && ZoomRatio < 24 
            
            %Enlarge the zoomed out image so we can do the cross-correlation
            ResizeFactor = max([Rows/SurfRows*ZoomRatio, Columns/SurfColumns*ZoomRatio]);
            % ES 2013-10-30: the reason I have to define ResizeFactor
            % differently from ZoomRatio is because you can't necessarily
            % infer the microns-per-pixel resolution from the zoom alone:
            % it also depends on the dimensions of the image. This may not
            % work for all possible resolutions, though...
            ZoomRatio = ResizeFactor;
            SurfImageResized=imresize(SurfImage, ResizeFactor);
            
            %Calculate the correlation matrix and find the maximum
            C = normxcorr2(ZoomImage, SurfImageResized);
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

            try
                %Make an overlay of the zoomed in and zoomed out real
                %images as well as of a quickly segmented nuclear mask
                
                %Real image overlay
                %Crop the zoomed out image to match the zoomed in one
                SurfImageResizeZoom=...
                    SurfImageResized(RowsResized/2-RowsZoom/2+ShiftRow*ZoomRatio:RowsResized/2+RowsZoom/2-1+ShiftRow*ZoomRatio,...
                    ColumnsResized/2-ColumnsZoom/2+ShiftColumn*ZoomRatio:ColumnsResized/2+ColumnsZoom/2-1+ShiftColumn*ZoomRatio);
                ImOverlay=cat(3,mat2gray(SurfImageResizeZoom),...
                    +mat2gray(ZoomImage),zeros(size(SurfImageResizeZoom)));

                %Nuclear mask overlay
                NucMaskZoomOut=GetNuclearMask(SurfImage,2.5,0);
                NucMaskZoomOutResized=imresize(NucMaskZoomOut, ResizeFactor);
                NucMaskZoomOutResizedCropped=...
                    NucMaskZoomOutResized(RowsResized/2-RowsZoom/2+ShiftRow*ZoomRatio:RowsResized/2+RowsZoom/2-1+ShiftRow*ZoomRatio,...
                    ColumnsResized/2-ColumnsZoom/2+ShiftColumn*ZoomRatio:ColumnsResized/2+ColumnsZoom/2-1+ShiftColumn*ZoomRatio);
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
            end

        else

            warning('Not doing cross-correlation for AP finding. Instead just looking at the center')

            %Add this flag to do the alignment manually later on
            if (~ManualAlignmentDone)|(ManualAlignment)
                ManualAlignment=1;
                ShiftColumn=0;
                ShiftRow=0;
            end
            
            ShiftColumn=0;
            ShiftRow=0;

        end


    elseif ~ManualAlignmentDone
        ShiftRow=0;
        ShiftColumn=0;
    end

%ManualAPCorrection(SurfImage,ZoomImage,C,ResizeFactor,ShiftRow,ShiftColumn)
    
    
    
    %Now figure out how the shift of the zoomed in and zoomed out surface
    %images translates to the whole embryo image (which is most of the
    %times stitched). This is necessary to figure out the AP position.

    
    %For full embryo images stitched out of two separate images, we need to
    %look at each case: the zoom in version being on the right or on the
    %left.
    
    %Load the full embryo image
    FullEmbryo=imread([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo.tif']);
    
    %We'll overlay the zoomed out surface and mid saggital images to check we got
    %things right.

    %Patch the surface image to fit the full embryo image.
    PatchedSurfaceImage=zeros(size(FullEmbryo));
    
    %The information from the stitching of the two images is as follows:
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
        TopLeftHalf=[ImageCenter(1)-Rows/ZoomRatio/2+ShiftRow+1,...
            ImageCenter(2)-Columns/ZoomRatio/2+ShiftColumn+1];
        BottomRightHalf=[ImageCenter(1)+Rows/ZoomRatio/2+ShiftRow,...
            ImageCenter(2)+Columns/ZoomRatio/2+ShiftColumn];
        %AP position mapped onto the zoomed out image
        coordAHalf=coordA+[-Columns+xShift,0];
        coordPHalf=coordP+[-Columns+xShift,0];


        %Start by overlaying the zoom in figure on top of the zoom out
        %figure. Also add the rectangle.
        
        %Create an overlay of the mask in zoom in and zoom out
        NucMaskZoomInOverlay=zeros(size(SurfImage));
        
        NucMaskZoomInOverlay(TopLeftHalf(1):BottomRightHalf(1),TopLeftHalf(2):BottomRightHalf(2))=...
                    imresize(NucMaskZoomIn,1/ZoomRatio);
        
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
        
%         TopLeft=[size(FullEmbryo,1)-Rows1x/2-Rows/ZoomRatio/2+ShiftRow,...
%             size(FullEmbryo,2)-Columns1x/2+Columns/ZoomRatio/2+xShift/ZoomRatio+ShiftColumn];
%         BottomRight=[size(FullEmbryo,1)-Rows1x/2+Rows/ZoomRatio/2+ShiftRow,...
%             size(FullEmbryo,2)-Columns1x/2+Columns/ZoomRatio*1.5+xShift/ZoomRatio+ShiftColumn];
        
        
        figure(6)
        imshow(imadjust(mat2gray(FullEmbryo)),'DisplayRange',[],'InitialMagnification',100)
        hold on
        rectangle('Position',[TopLeft([2,1]),BottomRight([2,1])-TopLeft([2,1])],'EdgeColor','r')
        plot(coordA(1),coordA(2),'.g','MarkerSize',30)
        plot(coordP(1),coordP(2),'.r','MarkerSize',30)
        plot([coordA(1),coordP(1)],[coordA(2),coordP(2)],'-b')
        hold off
            

    elseif sum(~cellfun('isempty',strfind(lower({D.name}),'left'))&~cellfun('isempty',strfind(lower({D.name}),'surf')))

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
        TopLeftHalf=[ImageCenter(1)-Rows/ZoomRatio/2+ShiftRow+1,...
            ImageCenter(2)-Columns/ZoomRatio/2+ShiftColumn+1];
        BottomRightHalf=[ImageCenter(1)+Rows/ZoomRatio/2+ShiftRow,...
            ImageCenter(2)+Columns/ZoomRatio/2+ShiftColumn];
        
%         TopLeftHalf=[ImageCenter(1)-Rows/ZoomRatio/2+1,...
%             ImageCenter(2)-Columns/ZoomRatio/2+1];
%         BottomRightHalf=[ImageCenter(1)+Rows/ZoomRatio/2,...
%             ImageCenter(2)+Columns/ZoomRatio/2];
        
        %AP position mapped onto the zoomed out image
        coordAHalf=coordA+[0,yShift];
        coordPHalf=coordP+[0,yShift];


        %Start by overlaying the zoom in figure on top of the zoom out
        %figure. Also add the rectangle.
        
        %Create an overlay of the mask in zoom in and zoom out
        NucMaskZoomInOverlay=zeros(size(SurfImage));
        
        NucMaskZoomInOverlay(TopLeftHalf(1):BottomRightHalf(1),TopLeftHalf(2):BottomRightHalf(2))=...
                    imresize(NucMaskZoomIn,1/ZoomRatio);
        
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
        
%         TopLeft=[size(FullEmbryo,1)-Rows1x/2-Rows/ZoomRatio/2+ShiftRow,...
%             size(FullEmbryo,2)-Columns1x/2+Columns/ZoomRatio/2+xShift/ZoomRatio+ShiftColumn];
%         BottomRight=[size(FullEmbryo,1)-Rows1x/2+Rows/ZoomRatio/2+ShiftRow,...
%             size(FullEmbryo,2)-Columns1x/2+Columns/ZoomRatio*1.5+xShift/ZoomRatio+ShiftColumn];
        
        
        figure(6)
        imshow(imadjust(mat2gray(FullEmbryo)),'DisplayRange',[],'InitialMagnification',100)
        hold on
        rectangle('Position',[TopLeft([2,1]),BottomRight([2,1])-TopLeft([2,1])],'EdgeColor','r')
        plot(coordA(1),coordA(2),'.g','MarkerSize',30)
        plot(coordP(1),coordP(2),'.r','MarkerSize',30)
        plot([coordA(1),coordP(1)],[coordA(2),coordP(2)],'-b')
        hold off
            
        
        

    %For the case where we needed to stitch three images
    elseif sum(~cellfun('isempty',strfind(lower({D.name}),'center'))&~cellfun('isempty',strfind(lower({D.name}),'surf')))

        
        warning('Have HG check this case in the code')
        
        
        HalfName=D(find(sum(~cellfun('isempty',strfind(lower({D.name}),'center'))&~cellfun('isempty',strfind(lower({D.name}),'surf'))))).name;
        HalfImageSurf=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,HalfName],ChannelToLoad);
        [Rows1x,Columns1x]=size(HalfImageSurf);


        ImageCenter=[Rows1x/2,Columns1x/2];



        %This is for the full image
        TopLeft=[ImageCenter(1)-Rows/ZoomRatio/2+ShiftRow,...
            Columns1x*3/2-xShift1-Columns/ZoomRatio/2+ShiftColumn];

        BottomRight=[ImageCenter(1)+Rows/ZoomRatio/2+ShiftRow,...
            Columns1x*3/2-xShift1+Columns/ZoomRatio/2+ShiftColumn];

        %This is for the half image
        TopLeftHalf=[ImageCenter(1)-Rows/ZoomRatio/2+ShiftRow,...
            ImageCenter(2)-Columns/ZoomRatio/2+ShiftColumn];
        BottomRightHalf=[ImageCenter(1)+Rows/ZoomRatio/2+ShiftRow,...
            ImageCenter(2)+Columns/ZoomRatio/2+ShiftColumn];

        coordAHalf=coordA+[-Columns1x+xShift1,0];
        coordPHalf=coordP+[-Columns1x+xShift1,0];





    else
        error('Problem with the surface file (or its naming) in the source data folder "FullEmbryo"')
    end



    %Plot the area where we imaged on top of the embryo

    %Check if the embryo could actually fit in one of the images. If that's the
    %case we need to shift the the AP poisitions and boxes.
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




    figure(7)
    imshow(imadjust(mat2gray(FullEmbryo)),'DisplayRange',[],'InitialMagnification',100)
    hold on
    rectangle('Position',[TopLeft([2,1]),BottomRight([2,1])-TopLeft([2,1])],'EdgeColor','r')
    plot(coordA(1),coordA(2),'.g','MarkerSize',30)
    plot(coordP(1),coordP(2),'.r','MarkerSize',30)
    plot([coordA(1),coordP(1)],[coordA(2),coordP(2)],'-b')
    hold off
    saveas(gcf, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryoArea.tif']);

    

    


    %Now, compare it to the surface picture
    SurfName=D(find(~cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
    SurfImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,SurfName],ChannelToLoad);

    %Did we already have information from manual alignment?
%     if ManualAlignmentDone
%         TopLeftHalf=[ImageCenter(1)-Rows/ZoomRatio/2+ShiftRow,...
%                 ImageCenter(2)-Columns1x/ZoomRatio/2+ShiftColumn];
%         BottomRightHalf=[ImageCenter(1)+Rows/ZoomRatio/2+ShiftRow,...
%             ImageCenter(2)+Columns1x/ZoomRatio/2+ShiftColumn];
%     end

    
    
    
    %Did the user want to do manual alignment?
    if ManualAlignment
       [ShiftColumn,ShiftRow,TopLeftHalf,BottomRightHalf]=ManualAPCorrection(SurfImage,TopLeftHalf,BottomRightHalf,...
           ShiftColumn,ShiftRow,coordAHalf,coordPHalf,Rows1x,Columns1x,ZoomRatio);
    end



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

%     coordAZoom=(coordAHalf-ImageCenter)*ZoomRatio+[Columns/2,Rows/2]-[ShiftColumn,ShiftRow]*ZoomRatio;
%     coordPZoom=(coordPHalf-ImageCenter)*ZoomRatio+[Columns/2,Rows/2]-[ShiftColumn,ShiftRow]*ZoomRatio;

    coordAZoom=(coordAHalf-[TopLeftHalf(2),TopLeftHalf(1)])*ZoomRatio;
    coordPZoom=(coordPHalf-[TopLeftHalf(2),TopLeftHalf(1)])*ZoomRatio;
    
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
    if coordPZoom(1)-coordAZoom(1) < 0
        APAngle = APAngle + pi;
    end
    % Correction for if APAngle is in quadrants II or III
    
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


    %Divide the image into AP bins. The size of the bin will depend on the
    %experiment
    if strfind(lower(Prefix),'eve')     %Eve2 experiments
        APResolution=0.01;
    %hb or kni BAC experiments
    elseif ~isempty(strfind(lower(Prefix),'hbbac'))|...
            ~isempty(strfind(lower(Prefix),'knibac'))     
        APResolution=0.015;
    else                                %All other experiments
        APResolution=0.025;
    end
    
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
        for i=1:length(Particles)
            %Angle between the x-axis and the particle using the A position as a
            %zero
            Angles=atan((Particles(i).yPos-coordAZoom(2))./(Particles(i).xPos-coordAZoom(1)));
            if Particles(i).xPos-coordAZoom(1) < 0
                Angles = Angles + pi;
            end
            % Correction for if Angles is in quadrants II or III
            
            %Distance between the points and the A point
            Distances=sqrt((coordAZoom(2)-Particles(i).yPos).^2+(coordAZoom(1)-Particles(i).xPos).^2);
            APPositions=Distances.*cos(Angles-APAngle);
            Particles(i).APpos=APPositions/APLength;
            
            %Determine the distance perpendicular to the AP axis. This is a
            %proxy for a DV axis.
            
            Particles(i).DVpos=Distances.*sin(Angles-APAngle);
            
        end
    end




    
    
    %Save AP detection information
    if exist('xShift')
        save([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'],'coordA','coordP',...
            'xShift','yShift','coordAZoom','coordPZoom') 
    else
        save([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'],'coordA','coordP',...
            'xShift1','yShift1','xShift2','yShift2','coordAZoom','coordPZoom') 
    end

    if ManualAlignment|ManualAlignmentDone
        save([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'],'ManualAlignmentDone',...
            'ShiftColumn','ShiftRow','-append')
    end
    

end




%Save particle information
if exist([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
    if exist('Threshold1')
        save([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'],'Particles','fad','fad2',...
            'Threshold1','Threshold2');
    else
        save([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'],'Particles','fad','fad2');
    end
end
    
 


