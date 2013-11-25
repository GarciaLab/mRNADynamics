function [coordA,coordP,xShift,yShift]=ManualFindAPAxis(varargin)

%This function finds the AP axis of the embryo using the images taken at
%1x. The inputs are optional. 1 tells the code to flip the AP axis. The
%prefix can also be input.

%a z: Move one image up and down
%, .: Move one to the right and left
%n m: Increase or decrease contrast
%f: Flipped the images. This is useful if the 1x pictures weren't taken in
%the right orientation.

%xShift and yShift are the shifts used to stitch the images.


%Load the folder information

% ES 2013-10-29: Required for multiple users to be able to analyze data on
% one computer
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders(varargin{1});


for i=1:length(varargin)
    if isnumeric(varargin{i})
        if varargin{i}==1
            FlipAP=1;
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

%Create the output folder
mkdir([DropboxFolder,filesep,Prefix,'\APDetection'])


%See if we're dealing with a Bcd case
if (~isempty(findstr(Prefix,'Bcd')))&(isempty(findstr(Prefix,'BcdE1')))
    SourcePath=[SourcePath,'\..\..\Bcd-GFP'];
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


%See if there's a center image. This is for the case where we needed three
%images to fit the whole embryo.
CenterFileIndex=find(~cellfun('isempty',strfind(lower({D.name}),'center'))&...
    cellfun('isempty',strfind(lower({D.name}),'surf')));

if (length(RightFileIndex)>1) | (length(LeftFileIndex)>1)
    error('Too many left/right files in FullEmbryo folder')
end

if ~isempty(CenterFileIndex)
    error('Write the program to deal with three images')
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

if isempty(Zoom)
    Dtemp=dir([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo\temp\*.tif']);
    LeftFileIndex=find(~cellfun('isempty',strfind(lower({Dtemp.name}),'left'))&...
        cellfun('isempty',strfind(lower({Dtemp.name}),'surf')));
    ImageInfo = imfinfo([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo\temp',filesep,Dtemp(LeftFileIndex).name]);
    Zoom=ExtractInformationField(ImageInfo(1),'state.acq.zoomFactor=');
end

%Look for the file
FFDir=dir([SourcePath,filesep,Date,'\FF',Zoom(1:end-1),'x*.*']);
%If there's more than one match then ask for help
if length(FFDir)==1
    FFFile=FFDir(1).name;
elseif isempty(FFDir)
    display('Warning, no flat field file found. Press any key to proceed without it');
    FFImage=ones(ImageInfo(1).Height,ImageInfo(1).Width);
    pause
else
    FFFile=uigetfile([Folder,'\..\FF',Zoom(1:end-1),'x*.*'],'Select flatfield file');
end




%Prepare the images for stichting (crop, etc.)
   
leftfilename=[SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo',filesep,D(LeftFileIndex).name];
rightfilename=[SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo',filesep,D(RightFileIndex).name];
FFfilename=[SourcePath,filesep,Date,filesep,FFFile];

%Parameters:
Margin=10;      %Pixels to get rid of at the left and margins
    

% Read in flat field, filter, convert from uint16 to double
FF=imread(FFfilename);
FF=imfilter(double(FF), fspecial('disk', 30), 'replicate', 'same');
FF=FF/mean(FF(:));
FF=(FF-1)*1+1;
   



%See if we have a histone channel. It assumes that it's the second one
%if so.


ImInfo=imfinfo(leftfilename);

if length(ImInfo)==2
   % Read in left and right images
   left=imread(leftfilename,2);
   right=imread(rightfilename,2);
else
   % Read in left and right images
   left=imread(leftfilename);
   right=imread(rightfilename);
end

% x and y limits
x_min=280; x_max=450; %y_min=-40; y_max=40;
y_min=-50; y_max=50;


%Crop the images
FF=FF(Margin+1:end-Margin,Margin+1:end-Margin);
left=left(Margin+1:end-Margin,Margin+1:end-Margin);
right=right(Margin+1:end-Margin,Margin+1:end-Margin);
   
DisplayRange(1)=min([min(min(left)),min(min(right))]);
DisplayRange(2)=max([max(max(left)),max(max(right))]);



xo1=250;
yo1=-50;
Flipped=0;  %Flag indicating whether we've flipped the images
EmbryoFigure=figure;
cc=1;

%Plot the image and overlay with the different particles found

while (cc~=13)
    figure(EmbryoFigure)
    imm1 = imstitch(left,right, xo1, yo1,[1 2]);
    imshow(imm1,DisplayRange)
    title([num2str(xo1),'/',num2str(yo1)])
    
    figure(EmbryoFigure)
    ct=waitforbuttonpress;
    cc=get(EmbryoFigure,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    %Left-right, fine
    if (ct~=0)&(cc=='.')
        xo1=xo1+1;
    elseif (ct~=0)&(cc==',')
        xo1=xo1-1;
    %Left-right, coarse
    elseif (ct~=0)&(cc=='>')
        xo1=xo1+10;
    elseif (ct~=0)&(cc=='<')
        xo1=xo1-10;        
        
    %Up-down, fine   
    elseif (ct~=0)&(cc=='a')
        yo1=yo1+1;
    elseif (ct~=0)&(cc=='z')
        yo1=yo1-1;
        
    %Up-down, coarse   
    elseif (ct~=0)&(cc=='A')
        yo1=yo1+10;
    elseif (ct~=0)&(cc=='Z')
        yo1=yo1-10;
        
        
    %Change the display range
    elseif (ct~=0)&(cc=='m')
        DisplayRange(2)=DisplayRange(2)/1.5;
    elseif (ct~=0)&(cc=='n')
        DisplayRange(2)=DisplayRange(2)*1.5;
    elseif (ct~=0)&(cc=='r')
        DisplayRange(1)=min([min(min(left)),min(min(right))]);
        DisplayRange(2)=max([max(max(left)),max(max(right))]);
        
        
    %Rotate images by 90 degrees and flip horizontally
    elseif (ct~=0)&(cc=='f')
        if Flipped==0
            display('Swapping and rotating images')
            
            LeftTemp=imrotate(left,-90);
            LeftTemp=flipdim(LeftTemp,2);
            
            RightTemp=imrotate(right,-90);
            RightTemp=flipdim(RightTemp,2);
            
            right=LeftTemp;
            left=RightTemp;
      
            Flipped=1;
            
        else
            display('Images have been flipped already')
        end
        
        
        
    
    
    end
    
end

[h, w] = size(left);
imm2=imm1(:,1:2*w+1-xo1);

%Rotate the image again if it was flipped
if Flipped==1
    imm2=flipdim(imm2,2);
    imm2=imrotate(imm2,90);
    xShift= w+yo1;
    yShift= -h+xo1;
else
    xShift=xo1;
    yShift=yo1;
end
    
    
APImage=imm2;


 
%Save it to the Dropbox folder
imwrite(uint16(APImage),[DropboxFolder,filesep,Prefix,'\APDetection\FullEmbryo.tif'],'compression','none');


%Now, use them to find the embryo mask
embMask = getEmbryoMask(APImage, 20);


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
if isempty(CenterFileIndex)
    save([DropboxFolder,filesep,Prefix,'\APDetection.mat'],'coordA','coordP',...
        'xShift','yShift');
else
    save([DropboxFolder,filesep,Prefix,'\APDetection.mat'],'coordA','coordP',...
        'xShift1','yShift1','xShift2','yShift2');
end



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
saveas(gcf, [DropboxFolder,filesep,Prefix,'\APDetection\APMask.tif']);

clf
imagesc(APImage)
axis image
axis off
title('Anterior (green), posterior (red); original')
hold on
plot(coordA(1),coordA(2),'g.','MarkerSize',20);
plot(coordP(1),coordP(2),'r.','MarkerSize',20);
hold off
saveas(gcf, [DropboxFolder,filesep,Prefix,'\APDetection\APEmbryo.tif']);
close(diagFigure);



