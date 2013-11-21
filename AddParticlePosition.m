function AddParticlePositionV2(varargin)

%First parameter should be the prefix. The other parameters can be:
%SkipAlignment
%ManualAlignment
%NoAP: Just add X and Y information

%V2: Changed this function to use a correlation in order to center the
%images.


% ES 2013-10-29: Required for multiple users to be able to analyze data on
% one computer
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
    %NOTE: I switched this because of the retracking with the new flat field
    %for some data sets
    %if ~isfield(Particles,'xPos')
        for i=1:length(Particles)
            for j=1:length(Particles(i).Frame)
                [x,y]=fad2xyzFit(Particles(i).Frame(j),fad, 'addMargin'); 
                Particles(i).xPos(j)=x(Particles(i).Index(j));
                Particles(i).yPos(j)=y(Particles(i).Index(j));
            end
        end
    %end



    if isfield(Particles,'APpos')
        warning('Particles.mat already has AP positions stored.')
        %pause
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
    Zoom=ExtractInformationField(ImageInfo(1),'state.acq.zoomFactor=');
    Zoom=str2num(Zoom);

    if Zoom==4
        Rows=256;
        Columns=512;
    elseif Zoom==8
        Rows=256;
        Columns=256;
    elseif Zoom==16
        Rows=128;
        Columns=128;
    else
        error('Include zoom information in AddParticlePosition.m')
    end







    %Make a folder to store the images
    mkdir([DropboxFolder,filesep,Prefix,filesep,'APDetection'])

    if HistoneChannel
        ChannelToLoad=2;

        %Get the surface image in the zoomed case
        D=dir([FISHPath,filesep,'Data',filesep,Prefix,filesep,Prefix,'-His*.tif']);
        ZoomImage=imread([FISHPath,filesep,'Data',filesep,Prefix,filesep,D(end-10).name]);


    else
        ChannelToLoad=1;

        %Get the surface image in the zoomed case
        D=dir([FISHPath,filesep,'Data',filesep,Prefix,filesep,Prefix,'*_z*.tif']);
        ZoomImage=imread([FISHPath,filesep,'Data',filesep,Prefix,filesep,D(end-10).name],ChannelToLoad);

    end

    %Get the surface image in the 1x
    D=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*.tif']);
    SurfName=D(find(~cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
    SurfImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,SurfName],ChannelToLoad);   




    if ~SkipAlignment&HistoneChannel
        if Zoom==4
            %Enlarge the 1x image so we can do the cross-correlation
            SurfImageResized=imresize(SurfImage,Zoom);

            %Calculate the correlation marrix and find the maximum
            C = normxcorr2(ZoomImage, SurfImageResized);
            [Max2,MaxRows]=max(C);
            [Dummy,MaxColumn]=max(Max2);
            MaxRow=MaxRows(MaxColumn);
            [CRows,CColumns]=size(C);

            ShiftRow=round((MaxRow-(CRows/2+1))/Zoom);
            ShiftColumn=round((MaxColumn-(CColumns/2+1))/Zoom);


            %How well did we do with the alignment?
            [RowsResized,ColumnsResized]=size(SurfImageResized);
            [RowsZoom,ColumnsZoom]=size(ZoomImage);

            SurfImageResizeZoom=SurfImageResized(RowsResized/2-RowsZoom/2+ShiftRow*Zoom:RowsResized/2+RowsZoom/2-1+ShiftRow*Zoom,...
                ColumnsResized/2-ColumnsZoom/2+ShiftColumn*Zoom:ColumnsResized/2+ColumnsZoom/2-1+ShiftColumn*Zoom);

            ImOverlay=cat(3,mat2gray(SurfImageResizeZoom),...
                +mat2gray(ZoomImage),zeros(size(SurfImageResizeZoom)));


            figure(4)
            imshow(ImOverlay)
            saveas(gcf, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'AlignmentOverlay.tif']);




            figure(5)
            contourf(abs(C((CRows-1)/2-RowsZoom/2:(CRows-1)/2+RowsZoom/2,...
                (CColumns-1)/2-ColumnsZoom/2:(CColumns-1)/2+ColumnsZoom/2)))
            saveas(gcf, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'AlignmentCorrelation.tif']);

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




    if sum(~cellfun('isempty',strfind(lower({D.name}),'right'))&~cellfun('isempty',strfind(lower({D.name}),'surf')))

        %Load the half image at the midsaggital plane
        HalfName=D(find(sum(cellfun('isempty',strfind(lower({D.name}),'right'))&cellfun('isempty',strfind(lower({D.name}),'surf'))))).name;
        %HalfImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo\',HalfName]);
        %[Rows1x,Columns1x]=size(HalfImage);

        %Load the half image at the surface
        HalfNameSurf=D(find(sum(cellfun('isempty',strfind(lower({D.name}),'right'))&cellfun('isempty',strfind(lower({D.name}),'surf'))))).name;
        HalfImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,HalfName],ChannelToLoad);
        [Rows1x,Columns1x]=size(HalfImage);


    %     %Load the zoomed-in image at the surface
    %     DZoom=dir([SourcePath,filesep,Date,filesep,EmbryoName,'\*.tif']);
    %     ZoomImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,DZoom(end).name],2);


        %Check if we need to account for the image shift due to stitching
        %when figuring out the imaging region.

        if yShift<0


            ImageCenter=[Rows1x/2,Columns1x/2];

            %This is for the full image
            TopLeft=[ImageCenter(1)-Rows/Zoom/2-ShiftRow,...
                Columns1x*3/2-xShift-Columns1x/Zoom/2-ShiftColumn];

            BottomRight=[ImageCenter(1)+Rows/Zoom/2-ShiftRow,...
                Columns1x*3/2-xShift+Columns1x/Zoom/2-ShiftColumn];

            %This is for the half image
            TopLeftHalf=[ImageCenter(1)-Rows1x/Zoom/2+ShiftRow,...
                ImageCenter(2)-Columns1x/Zoom/2+ShiftColumn];
            BottomRightHalf=[ImageCenter(1)+Rows1x/Zoom/2+ShiftRow,...
                ImageCenter(2)+Columns1x/Zoom/2+ShiftColumn];

            coordAHalf=coordA+[-Columns1x+xShift,0];
            coordPHalf=coordP+[-Columns1x+xShift,0];

        else


            ImageCenter=[Rows1x/2,Columns1x/2];



            %This is for the full image
            TopLeft=[ImageCenter(1)-Rows/Zoom/2+ShiftRow,...
                ImageCenter(2)+Columns/Zoom/2+xShift/Zoom+ShiftColumn];
            BottomRight=[ImageCenter(1)+Rows/Zoom/2+ShiftRow,...
                ImageCenter(2)+Columns/Zoom*1.5+xShift/Zoom+ShiftColumn];



            %This is for the half image
            TopLeftHalf=[ImageCenter(1)-Rows/Zoom/2+ShiftRow,...
                ImageCenter(2)-Columns/Zoom/2+ShiftColumn];
            BottomRightHalf=[ImageCenter(1)+Rows/Zoom/2+ShiftRow,...
                ImageCenter(2)+Columns/Zoom/2+ShiftColumn];



            coordAHalf=coordA+[-Columns/2+xShift/2/Zoom,0];
            coordPHalf=coordP+[-Columns/2+xShift/2/Zoom,0];

       end
    elseif sum(~cellfun('isempty',strfind(lower({D.name}),'left'))&~cellfun('isempty',strfind(lower({D.name}),'surf')))


        HalfName=D(find(sum(~cellfun('isempty',strfind(lower({D.name}),'left'))&~cellfun('isempty',strfind(lower({D.name}),'surf'))))).name;
        HalfImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,HalfName],ChannelToLoad);
        [Rows1x,Columns1x]=size(HalfImage);

        if yShift<0

            ImageCenter=[Rows1x/2,Columns1x/2];

            %This is for the full image
            TopLeft=[ImageCenter(1)-Rows/Zoom/2-yShift,ImageCenter(2)-Columns/Zoom/2];
            BottomRight=[ImageCenter(1)+Rows/Zoom/2-yShift,ImageCenter(2)+Columns/Zoom/2];

            %This is for the half image
            TopLeftHalf=[ImageCenter(1)-Rows/Zoom/2+ShiftRow,...
                ImageCenter(2)-Columns/Zoom/2+ShiftColumn];
            BottomRightHalf=[ImageCenter(1)+Rows/Zoom/2+ShiftRow,...
                ImageCenter(2)+Columns/Zoom/2+ShiftColumn];

            coordAHalf=coordA+[0,yShift];
            coordPHalf=coordP+[0,yShift];

        else

            ImageCenter=[Rows1x/2,Columns1x/2];
            TopLeft=[ImageCenter(1)-Rows/Zoom/2,ImageCenter(2)-Columns/Zoom/2];
            BottomRight=[ImageCenter(1)+Rows/Zoom/2,ImageCenter(2)+Columns/Zoom/2];

            TopLeftHalf=[ImageCenter(1)-Rows/Zoom/2+ShiftRow,...
                ImageCenter(2)-Columns/Zoom/2+ShiftColumn];
            BottomRightHalf=[ImageCenter(1)+Rows/Zoom/2+ShiftRow,...
                ImageCenter(2)+Columns/Zoom/2+ShiftColumn];

            coordAHalf=coordA;
            coordPHalf=coordP;
        end


    %For the case where we needed to stitch three images
    elseif sum(~cellfun('isempty',strfind(lower({D.name}),'center'))&~cellfun('isempty',strfind(lower({D.name}),'surf')))

        HalfName=D(find(sum(~cellfun('isempty',strfind(lower({D.name}),'center'))&~cellfun('isempty',strfind(lower({D.name}),'surf'))))).name;
        HalfImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo\',HalfName],ChannelToLoad);
        [Rows1x,Columns1x]=size(HalfImage);


        ImageCenter=[Rows1x/2,Columns1x/2];



        %This is for the full image
        TopLeft=[ImageCenter(1)-Rows/Zoom/2-ShiftRow,...
            Columns1x*3/2-xShift1-xShift2-Columns/Zoom/2-ShiftColumn];

        BottomRight=[ImageCenter(1)+Rows/Zoom/2-ShiftRow,...
            Columns1x*3/2--xShift1-xShift2+Columns/Zoom/2-ShiftColumn];

        %This is for the half image
        TopLeftHalf=[ImageCenter(1)-Rows/Zoom/2+ShiftRow,...
            ImageCenter(2)-Columns/Zoom/2+ShiftColumn];
        BottomRightHalf=[ImageCenter(1)+Rows/Zoom/2+ShiftRow,...
            ImageCenter(2)+Columns/Zoom/2+ShiftColumn];

        coordAHalf=coordA+[-Columns+xShift1,0];
        coordPHalf=coordP+[-Columns+xShift1,0];





    else
        error('Problem with the surface file (or its naming) in the source data folder "FullEmbryo"')
    end







    %Plot the area where we imaged on top of the embryo
    FullEmbryo=imread([DropboxFolder,filesep,Prefix,'\APDetection\FullEmbryo.tif']);


    %Check if the embryo could actually fit in one of the images. If that's the
    %case we need to shift the the AP poisitions and boxes.
    if sum(size(HalfImage)==size(FullEmbryo))==2
        ImageCenter=[Rows1x/2,Columns1x/2];

        %This is for the full image
        TopLeft=[ImageCenter(1)-Rows/Zoom/2,ImageCenter(2)-Columns/Zoom/2]...
            -[yShift,xShift];
        BottomRight=[ImageCenter(1)+Rows/Zoom/2,ImageCenter(2)+Columns/Zoom/2]...
            -[yShift,xShift];

        %This is for the acquisition image
        TopLeftHalf=[ImageCenter(1)-Rows/Zoom/2+ShiftRow,...
            ImageCenter(2)-Columns/Zoom/2+ShiftColumn];
        BottomRightHalf=[ImageCenter(1)+Rows/Zoom/2+ShiftRow,...
            ImageCenter(2)+Columns/Zoom/2+ShiftColumn];

        coordAHalf=coordA+[xShift,yShift];
        coordPHalf=coordP+[xShift,yShift];


    end






    figure(1)
    imshow(FullEmbryo,'DisplayRange',[],'InitialMagnification',100)
    hold on
    %rectangle('Position',[TopLeft([2,1]),BottomRight([2,1])-TopLeft([2,1])],'EdgeColor','r')
    plot(coordA(1),coordA(2),'.g','MarkerSize',30)
    plot(coordP(1),coordP(2),'.r','MarkerSize',30)
    plot([coordA(1),coordP(1)],[coordA(2),coordP(2)],'-b')
    hold off
    saveas(gcf, [DropboxFolder,filesep,Prefix,'\APDetection\FullEmbryoArea.tif']);





    %Now, compare it to the surface picture
    SurfName=D(find(~cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
    SurfImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo\',SurfName],ChannelToLoad);

    %Did we already have information from manual alignment?
    if ManualAlignmentDone
        TopLeftHalf=[ImageCenter(1)-Rows1x/Zoom/2+ShiftRow,...
            ImageCenter(2)-Columns1x/Zoom/2+ShiftColumn];
        BottomRightHalf=[ImageCenter(1)+Rows1x/Zoom/2+ShiftRow,...
            ImageCenter(2)+Columns1x/Zoom/2+ShiftColumn];
    end


    %Did the user want to do manual alignment?
    if ManualAlignment
       [ShiftColumn,ShiftRow,TopLeftHalf,BottomRightHalf]=ManualAPCorrection(SurfImage,TopLeftHalf,BottomRightHalf,...
           ShiftColumn,ShiftRow,coordAHalf,coordPHalf,Rows1x,Columns1x,Zoom);
    end



    figure(2)
    imshow(SurfImage,'DisplayRange',[],'InitialMagnification',100)
    hold on
    rectangle('Position',[TopLeftHalf([2,1]),BottomRightHalf([2,1])-TopLeftHalf([2,1])],'EdgeColor','r')
    plot(coordAHalf(1),coordAHalf(2),'.g','MarkerSize',30)
    plot(coordPHalf(1),coordPHalf(2),'.r','MarkerSize',30)
    plot([coordAHalf(1),coordPHalf(1)],[coordAHalf(2),coordPHalf(2)],'-b')
    plot([1],[1],'.y','MarkerSize',50)
    hold off
    saveas(gcf, [DropboxFolder,filesep,Prefix,'\APDetection\HalfEmbryoArea.tif']);


    D=dir([SourcePath,filesep,Date,filesep,EmbryoName,'\*.tif']);
    ZoomImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,D(end).name],ChannelToLoad);



    %We have the position of the anterior and posterior in the coordinates of the 
    %field of view we took. Do the mapping of the imaging region with respect to the AP axis.

    %Convert them to the corresponding zoom factor. In order to do this I
    %need to measure A and P with respect to the center of the rectangle.
    %Notice that there might be a +/-1 issue in the positioning given that
    %Matlab starts indexing from 1. However, as shown in the images below,
    %this doesn't make any real difference.   

    coordAZoom=(coordAHalf-ImageCenter)*Zoom+[Columns/2,Rows/2]-[ShiftColumn,ShiftRow]*Zoom;
    coordPZoom=(coordPHalf-ImageCenter)*Zoom+[Columns/2,Rows/2]-[ShiftColumn,ShiftRow]*Zoom;


    figure(3)
    imshow(ZoomImage,[])
    hold on
    plot([coordAZoom(1),coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'-b')
    plot([coordAZoom(1)+1,coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'--r')
    plot([coordAZoom(1)-1,coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'-.g')
    hold off
    saveas(gcf, [DropboxFolder,filesep,Prefix,'\APDetection\ZoomedEmbryoAP.tif']);


    %With AP coordinates in hand we can now determine the AP position of
    %all particles. Look I my notes in "Calculating AP positions" in Notability
    %for details of the calculation.



    %Angle between the x-axis and the AP-axis
    APAngle=atan((coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1)));
    APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);


    APPosImage=zeros(size(ZoomImage));
    [Rows,Columns]=size(ZoomImage);

    for i=1:Rows
        for j=1:Columns
            Angle=atan((i-coordAZoom(2))./(j-coordAZoom(1)));
            Distance=sqrt((coordAZoom(2)-i).^2+(coordAZoom(1)-j).^2);
            APPosition=Distance.*cos(Angle-APAngle);
            APPosImage(i,j)=APPosition/APLength;
        end
    end


    %Bin the pixels along the AP axis
    APResolution=0.025;
    APbinID=0:APResolution:1;


    APPosBinImage=zeros(size(APPosImage));
    for i=1:(length(APbinID)-1)
        FilteredMask=(APbinID(i)<=APPosImage)&(APbinID(i+1)>APPosImage);

        APPosBinImage=APPosBinImage+FilteredMask*i;
    end

    ZoomOverlay=cat(3,mat2gray(ZoomImage)/2+mat2gray(APPosBinImage)/2,...
        mat2gray(ZoomImage)/2,mat2gray(ZoomImage)/2);

    figure(4)
    imshow(ZoomOverlay)
    hold on
    plot([coordAZoom(1),coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'-b')
    plot([coordAZoom(1)+1,coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'--r')
    plot([coordAZoom(1)-1,coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'-.g')
    hold off


    if exist([DropboxFolder,filesep,Prefix,'\Particles.mat'])
        for i=1:length(Particles)
            %Angle between the x-axis and the particle using the A position as a
            %zero
            Angles=atan((Particles(i).yPos-coordAZoom(2))./(Particles(i).xPos-coordAZoom(1)));
            %Distance between the points and the A point
            Distances=sqrt((coordAZoom(2)-Particles(i).yPos).^2+(coordAZoom(1)-Particles(i).xPos).^2);
            APPositions=Distances.*cos(Angles-APAngle);
            Particles(i).APpos=APPositions/APLength;
        end
    end



    
    
    %Save AP detection information
    if exist('xShift')
        save([DropboxFolder,filesep,Prefix,'\APDetection.mat'],'coordA','coordP',...
            'xShift','yShift','coordAZoom','coordPZoom') 
    else
        save([DropboxFolder,filesep,Prefix,'\APDetection.mat'],'coordA','coordP',...
            'xShift1','yShift1','xShift2','yShift2','coordAZoom','coordPZoom') 
    end

    if ManualAlignment|ManualAlignmentDone
        save([DropboxFolder,filesep,Prefix,'\APDetection.mat'],'ManualAlignmentDone',...
            'ShiftColumn','ShiftRow','-append')
    end
    

end




%Save particle information
if exist([DropboxFolder,filesep,Prefix,'\Particles.mat'])
    if exist('Threshold1')
        save([DropboxFolder,filesep,Prefix,'\Particles.mat'],'Particles','fad','fad2',...
            'Threshold1','Threshold2');
    else
        save([DropboxFolder,filesep,Prefix,'\Particles.mat'],'Particles','fad','fad2');
    end
end
    
 


