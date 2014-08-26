function [coordA,coordP,xShift,yShift]=ManualFindAPAxis(varargin)

%This function finds the AP axis of the embryo using the images taken at
%1x. The inputs are optional. 1 tells the code to flip the AP axis. The
%prefix can also be input.

%a z: Move one image up and down
%, .: Move one to the right and left
%n m: Increase or decrease contrast
%s: Switch between left/center and center/left. This only works when having
%three images

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
   


%Are we dealing with a top-down or bottom-up orientation?
if sum(~cellfun('isempty',strfind(lower({D.name}),'right')))&...
        sum(~cellfun('isempty',strfind(lower({D.name}),'left')))
    LeftRight=1;
elseif sum(~cellfun('isempty',strfind(lower({D.name}),'top')))&...
        sum(~cellfun('isempty',strfind(lower({D.name}),'bottom')))
    LeftRight=0;
end

%See if there's a center image. This is for the case where we needed three
%images to fit the whole embryo.
CenterFileIndex=find(~cellfun('isempty',strfind(lower({D.name}),'center'))&...
    cellfun('isempty',strfind(lower({D.name}),'surf')));
if ~isempty(CenterFileIndex)
    ThreeImages=1;
    centerfilename=[SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo',filesep,D(CenterFileIndex).name];
else
    ThreeImages=0;
end



if LeftRight
    %Find the right and left files that do not correspond to the surface image
    RightFileIndex=find(~cellfun('isempty',strfind(lower({D.name}),'right'))&...
        cellfun('isempty',strfind(lower({D.name}),'surf')));
    LeftFileIndex=find(~cellfun('isempty',strfind(lower({D.name}),'left'))&...
        cellfun('isempty',strfind(lower({D.name}),'surf')));
else
    RightFileIndex=find(~cellfun('isempty',strfind(lower({D.name}),'bottom'))&...
        cellfun('isempty',strfind(lower({D.name}),'surf')));
    LeftFileIndex=find(~cellfun('isempty',strfind(lower({D.name}),'top'))&...
        cellfun('isempty',strfind(lower({D.name}),'surf')));
end



if (length(RightFileIndex)>1) | (length(LeftFileIndex)>1)
    error('Too many left/right files in FullEmbryo folder')
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
    if LeftRight
        LeftFileIndex=find(~cellfun('isempty',strfind(lower({Dtemp.name}),'left'))&...
            cellfun('isempty',strfind(lower({Dtemp.name}),'surf')));
    else
        LeftFileIndex=find(~cellfun('isempty',strfind(lower({Dtemp.name}),'top'))&...
            cellfun('isempty',strfind(lower({Dtemp.name}),'surf')));
    end
    ImageInfo = imfinfo([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo\temp',filesep,Dtemp(LeftFileIndex).name]);
    Zoom=ExtractInformationField(ImageInfo(1),'state.acq.zoomFactor=');
end


%Prepare the images for stichting (crop, etc.)
   
leftfilename=[SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo',filesep,D(LeftFileIndex).name];
rightfilename=[SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo',filesep,D(RightFileIndex).name];

%Parameters:
Margin=0;      %Pixels to get rid of at the left and margins
    

%See if we have a histone channel. It assumes that it's the second one
%if so.


ImInfo=imfinfo(leftfilename);

if length(ImInfo)==2
   % Read in left and right images
   left=imread(leftfilename,2);
   right=imread(rightfilename,2);
   if ThreeImages
       center=imread(centerfilename,2);
   end
else
   % Read in left and right images
   left=imread(leftfilename);
   right=imread(rightfilename);
   if ThreeImages
       center=imread(centerfilename);
   end
end

%Rotate the images if we're in top-down mode
if ~LeftRight
    left=imrotate(left,90);
    right=imrotate(right,90);
    
    if ThreeImages
        center=imrotate(center,90);
    end
end



% x and y limits
x_min=280; x_max=450; %y_min=-40; y_max=40;
y_min=-50; y_max=50;


%Crop the images
left=left(Margin+1:end-Margin,Margin+1:end-Margin);
right=right(Margin+1:end-Margin,Margin+1:end-Margin);
   
xo1=250;
yo1=-50;
xoDisplay=xo1;
yoDisplay=yo1;

if ThreeImages
    xo2=250;
    yo2=-50;
    leftDisplay=left;
    rightDisplay=center;
else
    leftDisplay=left;
    rightDisplay=right;
end

EmbryoFigure=figure;
cc=1;

%Sitch the image and figure out the display range
imm1 = imstitch(left,right, xoDisplay, yoDisplay,[1 2]);
DisplayRange(1)=min(min(imadjust(mat2gray(imm1))));
DisplayRange(2)=max(max(imadjust(mat2gray(imm1))));

%Are we aligning left-center or center-right? This is only relevant if we
%have three images.
CurrentPair=1;


%Plot the image and overlay with the different particles found
while (cc~=13)
    
    if ThreeImages
       1+1; 
    end
    
    
    figure(EmbryoFigure)
    imm1 = imstitch(leftDisplay,rightDisplay, xoDisplay, yoDisplay,[1 2]);
    imshow(imadjust(mat2gray(imm1)),DisplayRange)
    title([num2str(xoDisplay+2*Margin),'/',num2str(yoDisplay)])
    
    figure(EmbryoFigure)
    ct=waitforbuttonpress;
    cc=get(EmbryoFigure,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    %Left-right, fine
    if (ct~=0)&(cc=='.')
        xoDisplay=xoDisplay+1;
    elseif (ct~=0)&(cc==',')
        xoDisplay=xoDisplay-1;
    %Left-right, coarse
    elseif (ct~=0)&(cc=='>')
        xoDisplay=xoDisplay+10;
    elseif (ct~=0)&(cc=='<')
        xoDisplay=xoDisplay-10;        
        
    %Up-down, fine   
    elseif (ct~=0)&(cc=='a')
        yoDisplay=yoDisplay+1;
    elseif (ct~=0)&(cc=='z')
        yoDisplay=yoDisplay-1;
        
    %Up-down, coarse   
    elseif (ct~=0)&(cc=='A')
        yoDisplay=yoDisplay+10;
    elseif (ct~=0)&(cc=='Z')
        yoDisplay=yoDisplay-10;
        
        
    %Change the display range
    elseif (ct~=0)&(cc=='m')
        DisplayRange(2)=DisplayRange(2)/1.5;
    elseif (ct~=0)&(cc=='n')
        DisplayRange(2)=DisplayRange(2)*1.5;
    elseif (ct~=0)&(cc=='r')
        DisplayRange(1)=min([min(min(left)),min(min(right))]);
        DisplayRange(2)=max([max(max(left)),max(max(right))]);
        
        
    %Switch between embryo pairs
    elseif (ct~=0)&(cc=='s')&ThreeImages
        if CurrentPair==1
            CurrentPair=2;
            
            %Save the stichting information
            xo1=xoDisplay;
            yo1=yoDisplay;
            xoDisplay=xo2;
            yoDisplay=yo2;
            
            leftDisplay=center;
            rightDisplay=right;
            
        elseif CurrentPair==2;
            CurrentPair=1;
            
            %Save the stichting information
            xo2=xoDisplay;
            yo2=yoDisplay;
            xoDisplay=xo1;
            yoDisplay=yo1;
            
            leftDisplay=left;
            rightDisplay=center;
        end
        
    %Debug mode
    elseif (ct~=0)&(cc=='9')
        keyboard;
    
    end
    
end

if ~ThreeImages
    xo1=xoDisplay;
    yo1=yoDisplay;    
    
    [h, w] = size(left);
    imm2=imm1(:,1:2*w+1-xo1);

    %Rotate the image again if it was flipped
    if ~LeftRight
        imm2=imrotate(imm2,-90);
        xShift= yo1;
        yShift= xo1;
    else
        xShift=xo1+2*Margin;
        yShift=yo1;
    end


    APImage=imm2;
else
    if CurrentPair==1
        xo1=xoDisplay;
        yo1=yoDisplay;
    elseif CurrentPair==2;
        xo2=xoDisplay;
        yo2=yoDisplay;
    end 
    
    %Stitch all images together
    
    %Stitch left and center together
    [h, w] = size(left);
    imm2 = zeros(h, 2*w);

    %imm2(:, 1:w) = LeftImage;
    imm2(:, 1+Margin:w-Margin) = left(:,1+Margin:end-Margin);

    imm2(:, w+1-xo1+Margin:2*w-xo1-Margin) = circshift(center(1:h, 1+Margin:end-Margin), [yo1, 0]);

    imm2=circshift(imm2(:,:),[-yo1, 0]);



    %Stitch left+center with right
    [h,w]=size(imm2);
    [hRight,wRight]=size(right);

    APImage = zeros(h, w+wRight);
    APImage(:, 1:w)=imm2;
    APImage(:, w+1-xo2-xo1+Margin:wRight+w-xo2-xo1-Margin) =...
       circshift(right(1:h,1+Margin:end-Margin), [yo2, 0]);

    APImage=APImage(:,1:w+wRight+1-xo2-xo1);

    if ~LeftRight
        %Rotate the image back
        APImage=imrotate(APImage,-90);

        %Swap the shifts
        yShift1=xo1;
        xShift1=yo1;
        yShift2=xo2;
        xShift2=yo2;
    else
        yShift1=yo1;
        xShift1=xo1;
        yShift2=yo2;
        xShift2=xo2;
    end
    
    
end


 
%Save it to the Dropbox folder
imwrite(uint16(APImage),[DropboxFolder,filesep,Prefix,'\APDetection\FullEmbryo.tif'],'compression','none');


%Now, use them to find the embryo mask
embMask = getEmbryoMask(APImage, 20);


%This code came from Michael's code
% ES 2014-03-01: I changed this to handle cases where the threshold is too
% high but the embryo is still stitchable.
CC=bwconncomp(embMask);
if CC.NumObjects~=1
    disp('Failed to calculate embryo mask. Found more than one object in mask. Re-running with a lower threshold.');
    embMask = GetEmbryoMaskWithThreshold(APImage, 20, 5);
    CC = bwconncomp(embMask);
    if CC.NumObjects~=1
        error('Failed to calculate embryo mask. Found more than one object in mask. Script aborting.');
    end
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
        'xShift','yShift','Margin');
else
    save([DropboxFolder,filesep,Prefix,'\APDetection.mat'],'coordA','coordP',...
        'xShift1','yShift1','xShift2','yShift2','Margin');
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



