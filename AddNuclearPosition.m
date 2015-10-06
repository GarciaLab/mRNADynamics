function AddNuclearPosition(varargin)

%V2: Changed this function to use a correlation in order to center the
%images.

%Find out which computer this is. That will determine the folder structure.
[ret, name] = system('hostname');  
if ret ~= 0,  
   if ispc  
      name = getenv('COMPUTERNAME');  
   else  
      name = getenv('HOSTNAME');  
   end  
end  
name = lower(name); 


if strcmp(name(1:end-1),'phy-tglab2')
    SourcePath='D:\Hernan\LivemRNA\Analysis\Bcd-GFP';
    FISHPath='D:\Hernan\FISHDrosophila';
    DropboxFolder='D:\Hernan\Dropbox\LivemRNAData';
elseif strcmp(name(1:end-1),'hernanx200')
    SourcePath='C:\Users\hgarcia\Documents\My Papers\LivemRNA\Analysis\Bcd-GFP';
    FISHPath='C:\Users\hgarcia\Documents\My Papers\FISHDrosophila';
    DropboxFolder='C:\Users\hgarcia\Documents\Dropbox\LivemRNAData';
elseif strcmp(name(1:end-1),'albert-pc')
    SourcePath='C:\Users\Albert\Documents\Princeton\Gregor Lab\Data Analysis\LivemRNA\Analysis\Bcd-GFP';
    FISHPath='C:\Users\Albert\Documents\Princeton\Gregor Lab\Data Analysis\FISHDrosophila';
    DropboxFolder='C:\Users\Albert\Dropbox\LivemRNAData';
else    
    error('Include the folders for this computer in the code')
end



SkipAlignment=0;
ManualAlignment=0;

if ~isempty(varargin)
    Prefix=varargin{1};
    for i=2:length(varargin)
        if strcmp(varargin{2},'SkipAlignment')
            display('Skipping alignment step')
            SkipAlignment=1;
        elseif strcmp(varargin{2},'ManualAlignment')
            ManualAlignment=1;
        else
            error('Input not recognized')
        end
    end
            
else
    FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
    Dashes=strfind(FolderTemp,'\');
    Prefix=FolderTemp((Dashes(end)+1):end);
end



load([DropboxFolder,filesep,Prefix,'\Nuclei.mat'])




if isfield(Nuclei,'APpos')
    warning('Particles.mat already has AP positions stored.')
    %pause
end



%Figure out how our field of view maps to the AP coordinates.

%First, figure out how our field of view maps to the stitched embryo
%image

%Were the images taken on the left or right half of the embryo?


%Find out the date it was taken
Dashes=findstr(Prefix,'-');
Date=Prefix(1:Dashes(3)-1);
EmbryoName=Prefix(Dashes(3)+1:end);

D=dir([SourcePath,filesep,Date,filesep,EmbryoName,'\*.tif']);

%Get the information about the zoom
ImageInfo = imfinfo([SourcePath,filesep,Date,filesep,EmbryoName,filesep,D(1).name]);



%Get the information about the AP axis as well as the image shifts
%used for the stitching of the two halves of the embryo
load([DropboxFolder,filesep,Prefix,'\APDetection.mat'])


%Figure out the zoom factor
Zoom=ExtractInformationField(ImageInfo(1),'state.acq.zoomFactor=');
Zoom=str2num(Zoom);

if Zoom==4
    Rows=256;
    Columns=512;
elseif Zoom==8
    Rows=256;
    Columns=256;
elseif Zoom==2
    Rows=512;
    Columns=512;
else
    error('Include zoom information in AddParticlePosition.m')
end




%Get the surface image in the zoomed case
D=dir([SourcePath,filesep,Date,'\BcdGFP-HisRFP\AveragedData\',Prefix,'-His_*.tif']);
ZoomImage=imread([SourcePath,filesep,Date,'\BcdGFP-HisRFP\AveragedData\',D(end).name]);



%Get the surface image in the 1x
D=dir([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo\*.tif']);
SurfName=D(find(~cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
SurfImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo\',SurfName],2);



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




figure(5)
contourf(abs(C((CRows-1)/2-RowsZoom/2:(CRows-1)/2+RowsZoom/2,...
    (CColumns-1)/2-ColumnsZoom/2:(CColumns-1)/2+ColumnsZoom/2)))


%Get the folder information again
D=dir([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo\*.tif']);

if sum(~cellfun('isempty',strfind(lower({D.name}),'right'))&~cellfun('isempty',strfind(lower({D.name}),'surf')))

    %Load the half image at the midsaggital plane
    HalfName=D(find(sum(cellfun('isempty',strfind(lower({D.name}),'right'))&cellfun('isempty',strfind(lower({D.name}),'surf'))))).name;
    %HalfImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo\',HalfName]);
    %[Rows1x,Columns1x]=size(HalfImage);
    
    %Load the half image at the surface
    HalfNameSurf=D(find(sum(cellfun('isempty',strfind(lower({D.name}),'right'))&cellfun('isempty',strfind(lower({D.name}),'surf'))))).name;
    HalfImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo\',HalfName]);
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
            Columns1x*3/2-xShift-Columns/Zoom/2-ShiftColumn];

        BottomRight=[ImageCenter(1)+Rows/Zoom/2-ShiftRow,...
            Columns1x*3/2-xShift+Columns/Zoom/2-ShiftColumn];

        %This is for the half image
        TopLeftHalf=[ImageCenter(1)-Rows/Zoom/2+ShiftRow,...
            ImageCenter(2)-Columns/Zoom/2+ShiftColumn];
        BottomRightHalf=[ImageCenter(1)+Rows/Zoom/2+ShiftRow,...
            ImageCenter(2)+Columns/Zoom/2+ShiftColumn];

        coordAHalf=coordA+[-Columns+xShift,0];
        coordPHalf=coordP+[-Columns+xShift,0];

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
    HalfImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo\',HalfName]);
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

             

        coordAHalf=coordA+[-Columns/2+xShift/2/Zoom,0];
        coordPHalf=coordP+[-Columns/2+xShift/2/Zoom,0];
        

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


    

else
    error('Problem with the surface file (or its naming) in the source data folder "FullEmbryo"')
end




%Plot the area where we imaged on top of the embryo
FullEmbryo=imread([DropboxFolder,filesep,Prefix,'\FullEmbryo.tif']);


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
saveas(gcf, [DropboxFolder,filesep,Prefix,'\FullEmbryoArea.tif']);




%Now, compare it to the surface picture
SurfName=D(find(~cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
SurfImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo\',SurfName],2);

figure(2)
imshow(SurfImage,'DisplayRange',[0,1500],'InitialMagnification',100)
hold on
rectangle('Position',[TopLeftHalf([2,1]),BottomRightHalf([2,1])-TopLeftHalf([2,1])],'EdgeColor','r')
plot(coordAHalf(1),coordAHalf(2),'.g','MarkerSize',30)
plot(coordPHalf(1),coordPHalf(2),'.r','MarkerSize',30)
plot([coordAHalf(1),coordPHalf(1)],[coordAHalf(2),coordPHalf(2)],'-b')
plot([1],[1],'.y','MarkerSize',50)
hold off
saveas(gcf, [DropboxFolder,filesep,Prefix,'\HalfEmbryoArea.tif']);






D=dir([SourcePath,filesep,Date,filesep,EmbryoName,'\*.tif']);
ZoomImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,D(end).name],2);



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
saveas(gcf, [DropboxFolder,filesep,Prefix,'\ZoomedEmbryoAP.tif']);


%With AP coordinates in hand we can now determine the AP position of
%all particles.

for i=1:length(Nuclei)
    Nuclei(i).APpos=(Nuclei(i).x-coordAZoom(1))/(coordPZoom(1)-coordAZoom(1));
end


save([DropboxFolder,filesep,Prefix,'\Nuclei.mat'],'Nuclei');
save([DropboxFolder,filesep,Prefix,'\APDetection.mat'],'coordA','coordP',...
    'xShift','yShift','coordAZoom','coordPZoom')    ;
    
    
    
    

