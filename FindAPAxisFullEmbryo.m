function [coordA,coordP,xShift,yShift]=FindAPAxisFullEmbryo(varargin)

%When we could fit the whole embryo in one image we don't do stitching.
%Instead we do a correlation between the imaging field and the full embryo
%to determine the shift and then find the AP axis.


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
    SourcePath='D:\Hernan\LivemRNA\Analysis\MS2\MCPNoNLS+MS2';
    FISHPath='D:\Hernan\FISHDrosophila';
    DropboxFolder='D:\Hernan\Dropbox\LivemRNAData';
elseif strcmp(name(1:end-1),'hernanx200')
    SourcePath='C:\Users\hgarcia\Documents\My Papers\LivemRNA\Analysis\MS2\MCPNoNLS+MS2';
    FISHPath='C:\Users\hgarcia\Documents\My Papers\FISHDrosophila';
    DropboxFolder='C:\Users\hgarcia\Documents\Dropbox\LivemRNAData';
elseif strcmp(name(1:end-1),'albert-pc')
    SourcePath='C:\Users\Albert\Documents\Princeton\Gregor Lab\Data Analysis\LivemRNA\Analysis\MS2\MCPNoNLS+MS2';
    FISHPath='C:\Users\Albert\Documents\Princeton\Gregor Lab\Data Analysis\FISHDrosophila';
    DropboxFolder='C:\Users\Albert\Dropbox\LivemRNAData';
elseif strcmp(name(1:end-1),'phy-tglab11')
    SourcePath='Z:\LivemRNA\Analysis\MS2\MCPNoNLS+MS2';
    FISHPath='Z:\FISHDrosophila';
    DropboxFolder='C:\Users\bothma\Dropbox\LivemRNAData';
elseif strcmp(name(1:end-1),'bothma-desktop')
    SourcePath='Z:\LivemRNA\Analysis\MS2\MCPNoNLS+MS2';
    FISHPath='Z:\FISHDrosophila';
    DropboxFolder='C:\Users\bothma\Dropbox\LivemRNAData';
else    
    error('Include the folders for this computer in the code')
end




for i=1:length(varargin)
    if isnumeric(varargin{i})
        if varargin{i}==1
            FipAP=1;
        end
    elseif ischar(varargin{i})
        Prefix=varargin{i};
    end
end
   
if ~exist('Prefix')
    FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
    Dashes=strfind(FolderTemp,'\');
    Prefix=FolderTemp((Dashes(end)+1):end);
end


%If the AP axis hasn't been specified check if it was specificed in the
%image names. Otherwise assume AP orientation

%Find out the date it was taken
Dashes=findstr(Prefix,'-');
Date=Prefix(1:Dashes(3)-1);
EmbryoName=Prefix(Dashes(3)+1:end);

D=dir([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo\*.tif']);
   

%Find the right and left files that do not correspond to the surface image
RightFileIndex=find(~cellfun('isempty',strfind(lower({D.name}),'right'))&...
    cellfun('isempty',strfind(lower({D.name}),'surf')));
LeftFileIndex=find(~cellfun('isempty',strfind(lower({D.name}),'left'))&...
    cellfun('isempty',strfind(lower({D.name}),'surf')));
SurfaceFileIndex=find(~cellfun('isempty',strfind(lower({D.name}),'surf')));


if (length(RightFileIndex)>1) | (length(LeftFileIndex)>1)
    error('Too many left/right files in FullEmbryo folder')
end

%Determine which image has the full embryo. This corresponds to the one
%that doesn't have surf in the name.
if ~isempty(strfind(D(SurfaceFileIndex).name,'right'))
    FullImageIndex=LeftFileIndex;
    AcqImageIndex=RightFileIndex;
elseif ~isempty(strfind(D(SurfaceFileIndex).name,'left'))
    FullImageIndex=RightFileIndex;
    AcqImageIndex=LeftFileIndex;
else
    error('Some image missing?')
end



%See if we don't want the default AP orientation
if ~exist('FlipAP')
    if strcmp(D(RightFileIndex).name,'PA')|strcmp(D(LeftFileIndex).name,'PA')
        FlipAP=1;
    else
        FlipAP=0;
    end
end

    
%Get the necessary information to load the corresponding flat field image

%Get the structure with the acquisition information
ImageInfo = imfinfo([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo',filesep,D(LeftFileIndex).name]);


%Get the flat-field information

%Figure out the zoom factor
Zoom=ExtractInformationField(ImageInfo(1),'state.acq.zoomFactor=');
%Look for the file
FFDir=dir([SourcePath,filesep,Date,'\FF',Zoom(1),'x*.*']);
%If there's more than one match then ask for help
if length(FFDir)==1
    FFFile=FFDir(1).name;
elseif isempty(FFDir)
    display('Warning, no flat field file found. Press any key to proceed without it');
    FFImage=ones(ImageInfo(1).Height,ImageInfo(1).Width);
    pause
else
    FFFile=uigetfile([Folder,'\..\FF',Zoom(1),'x*.*'],'Select flatfield file');
end


%Now do the image correlation. We'll grab a small area around the middle of
%the imaging field and correlate it to the larger, full embryo image.

AcqImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo',filesep,D(AcqImageIndex).name]);
FullImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo',filesep,D(FullImageIndex).name]);

%Crop the imaging field by 50%
[Rows,Columns]=size(AcqImage);
AcqImageCrop=AcqImage(Rows/4:Rows/4*3,Columns/4:Columns/4*3);

%Calculate the correlation marrix and find the maximum
C = normxcorr2(AcqImageCrop, FullImage);
[Max2,MaxRows]=max(C);
[Dummy,MaxColumn]=max(Max2);
MaxRow=MaxRows(MaxColumn);

[CRows,CColumns]=size(C);



ShiftRow=MaxRow-(CRows/2+1);
ShiftColumn=MaxColumn-(CColumns/2+1);

%Create an overlay to make sure things make sense

ImShifted=uint16(zeros(size(FullImage)));

RowRange=(Rows/4:Rows/4*3)+ShiftRow;
ColumnRange=(Columns/4:Columns/4*3)+ShiftColumn;

ImShifted(RowRange,ColumnRange)=AcqImageCrop;

figure(1)
ImOverlay=cat(3,mat2gray(FullImage)+mat2gray(ImShifted),mat2gray(FullImage),mat2gray(FullImage));
imshow(ImOverlay)

xShift=-ShiftColumn;
yShift=-ShiftRow;


%Save it to the Dropbox folder
imwrite(FullImage,[DropboxFolder,filesep,Prefix,'\FullEmbryo.tif'],'compression','none');


%Now, use them to find the embryo mask
embMask = getEmbryoMask(FullImage, 20);


%This code came from Michael's code
CC=bwconncomp(embMask);
if CC.NumObjects~=1
    error('Failed to calculate embryo mask. Found more than one object in mask.');
end
    


% Rotate the mask to determine the AP axis as the extremal points of the mask
Props=regionprops(CC,'Orientation');
angle=Props.Orientation; % Angle is in DEGREES!


I_mask_rot=imrotate(embMask,-angle);
rotMatrix = [cosd(angle) sind(angle)
            -sind(angle) cosd(angle)];

        
CC=bwconncomp(I_mask_rot);
Props=regionprops(CC,'Centroid','MajorAxisLength', 'MinorAxisLength','Extrema');
% After rotation, the major axis is aligned with x axis



% for future diagnostic figures
majorAxisBegin = Props.Centroid + [Props.MajorAxisLength/2,0];
majorAxisEnd = Props.Centroid - [Props.MajorAxisLength/2,0];
minorAxisBegin = Props.Centroid + [0, Props.MinorAxisLength/2];
minorAxisEnd = Props.Centroid - [0, Props.MinorAxisLength/2];

ext=Props.Extrema;
coordP_rot=(ext(3,:)+ext(4,:))/2;
coordA_rot=(ext(7,:)+ext(8,:))/2;

if FlipAP
    temp = coordA_rot;
    coordA_rot = coordP_rot;
    coordP_rot = temp;
end


% coordA and coordP are the coordinates on the rotated image
% We should rotate them back to the coordinates of the original picture
% Remember that rotation was performed about the center of the image

%coordinates of the center of the rotated image
center_rot = 1/2*[size(I_mask_rot,2) size(I_mask_rot,1)];
%coordinates of the center of the original image
center = 1/2*[size(embMask,2) size(embMask,1)];

coordA = center + (rotMatrix * (coordA_rot-center_rot)')';
coordP = center + (rotMatrix * (coordP_rot-center_rot)')';

%Save the AP and shift information
save([DropboxFolder,filesep,Prefix,'\APDetection.mat'],'coordA','coordP',...
    'xShift','yShift');



% Save diagnostic figures to check the quality of axis determination
diagFigure = figure;
imagesc(I_mask_rot);
colormap(gray);
axis image
title('Anterior (green), posterior (red); rotated')
hold on
plot(coordA_rot(1),coordA_rot(2),'g.','MarkerSize',20);
plot(coordP_rot(1),coordP_rot(2),'r.','MarkerSize',20);
plot([majorAxisBegin(1),majorAxisEnd(1)],[majorAxisBegin(2),majorAxisEnd(2)],'b-');
plot([minorAxisBegin(1),minorAxisEnd(1)],[minorAxisBegin(2),minorAxisEnd(2)],'b-');
hold off
saveas(gcf, [DropboxFolder,filesep,Prefix,'\APMask.tif']);


clf
imagesc(FullImage)
axis image
axis off
title('Anterior (green), posterior (red); original')
hold on
plot(coordA(1),coordA(2),'g.','MarkerSize',20);
plot(coordP(1),coordP(2),'r.','MarkerSize',20);
hold off
saveas(gcf, [DropboxFolder,filesep,Prefix,'\APEmbryo.tif']);
close(diagFigure);



