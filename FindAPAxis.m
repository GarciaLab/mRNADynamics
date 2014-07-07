function [coordA,coordP,xShift,yShift]=FindAPAxis(varargin)

%This function finds the AP axis of the embryo using the images taken at
%1x. The inputs are optional. 1 tells the code to flip the AP axis. The
%prefix can also be input.

%xShift and yShift are the shifts used to stitch the images.
%xShift is the displacement of the right image with respect to the left
%image. Positive xShift moves the right image towards the left.
%yShift is the displacement of the left image with respect to the right
%image. Positive yShift moves the left image up.


%Parameters
Margin=0;

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
if (~isempty(findstr(Prefix,'Bcd')))&...
        (isempty(findstr(Prefix,'BcdE1'))&(isempty(findstr(Prefix,'NoBcd'))))
    SourcePath=[SourcePath,'\..\..\Bcd-GFP'];
end


%If the AP axis hasn't been specified check if it was specificed in the
%image names. Otherwise assume AP orientation

%Find out the date it was taken
Dashes=findstr(Prefix,'-');
Date=Prefix(1:Dashes(3)-1);
EmbryoName=Prefix(Dashes(3)+1:end);

D=dir([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo\*.tif']);
   
%Determine wether we're dealing with a left-right embryo or a top-bottom
%one
LeftRight=1;

if ~isempty(strfind(lower(D(1).name),'bottom'))|~isempty(strfind(lower(D(1).name),'top'))
    LeftRight=0;
    error('The part of the code to do the automatic stitching on top/bottom embryos needs to be implemented')
end



%Find the right and left files that do not correspond to the surface image
if LeftRight
    RightFileIndex=find(~cellfun('isempty',strfind(lower({D.name}),'right'))&...
        cellfun('isempty',strfind(lower({D.name}),'surf')));
    LeftFileIndex=find(~cellfun('isempty',strfind(lower({D.name}),'left'))&...
        cellfun('isempty',strfind(lower({D.name}),'surf')));
else
    LeftFileIndex=find(~cellfun('isempty',strfind(lower({D.name}),'top'))&...
        cellfun('isempty',strfind(lower({D.name}),'surf')));
    RightFileIndex=find(~cellfun('isempty',strfind(lower({D.name}),'bottom'))&...
        cellfun('isempty',strfind(lower({D.name}),'surf')));
end
    
    
    
%See if there's a center image. This is for the case where we needed three
%images to fit the whole embryo.
CenterFileIndex=find(~cellfun('isempty',strfind(lower({D.name}),'center'))&...
    cellfun('isempty',strfind(lower({D.name}),'surf')));

if (length(RightFileIndex)>1) | (length(LeftFileIndex)>1)
    error('Too many left/right files in FullEmbryo folder')
end


%See if we don't want the default AP orientation
if ~exist('FlipAP', 'var')
    if ~isempty(strfind(D(RightFileIndex).name,'PA')) || ~isempty(strfind(D(LeftFileIndex).name,'PA'))
        % ES 2014-01-11: originally this was "strcmp" and not "strfind",
        % but I don't think "strcmp" could work
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



%Sticht the images
if isempty(CenterFileIndex)
    [APImage,xShift,yShift]=EmbryoStitchNoMargin([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo',filesep,D(LeftFileIndex).name],...
       [SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo',filesep,D(RightFileIndex).name],...
       [SourcePath,filesep,Date,filesep,FFFile],[],[],Margin);
else
    [APImage1,xShift1,yShift1]=EmbryoStitchNoMargin([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo',filesep,D(LeftFileIndex).name],...
       [SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo',filesep,D(CenterFileIndex).name],...
       [SourcePath,filesep,Date,filesep,FFFile],[],[],Margin);

    [APImage2,xShift2,yShift2]=EmbryoStitchNoMargin([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo',filesep,D(CenterFileIndex).name],...
       [SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo',filesep,D(RightFileIndex).name],...
       [SourcePath,filesep,Date,filesep,FFFile],[],[],Margin);

  
   
   %Now try to put everything together
   LeftImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo',filesep,D(LeftFileIndex).name]);
   CenterImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo',filesep,D(CenterFileIndex).name]);
   RightImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,'\FullEmbryo',filesep,D(RightFileIndex).name]);
   
   %Get the FF information
   

    %If it's there, load the image and smooth it with a Gaussian.
    filtStd=30;         %This came from the FISH code.
    FFImage=double(imread([SourcePath,filesep,Date,filesep,FFFile],1));
    FFImage=imfilter(FFImage,fspecial('gaussian',2*filtStd,filtStd),'symmetric');
    FFImage=imdivide(FFImage,double(max(FFImage(:))));

    LeftImage=imdivide(double(LeftImage),FFImage);
    CenterImage=imdivide(double(CenterImage),FFImage);
    RightImage=imdivide(double(RightImage),FFImage);

   
    %Stitch left and center together
    [h, w] = size(LeftImage);
    imm2 = zeros(h, 2*w);

    %imm2(:, 1:w) = LeftImage;
    imm2(:, 1+Margin:w-Margin) = LeftImage(:,1+Margin:end-Margin);

    %imm2(:, w+1-xShift1:2*w-xShift1) = circshift(CenterImage(1:h, :), [yShift1, 0]);

    imm2(:, w+1-xShift1+Margin:2*w-xShift1-Margin) = circshift(CenterImage(1:h, 1+Margin:end-Margin), [yShift1, 0]);

    imm2=circshift(imm2(:,:),[-yShift1, 0]);



    %Stitch left+center with right
    [h,w]=size(imm2);
    [hRight,wRight]=size(RightImage);

    APImage = zeros(h, w+wRight);
    APImage(:, 1:w)=imm2;
    APImage(:, w+1-xShift2-xShift1+Margin:wRight+w-xShift2-xShift1-Margin) =...
       circshift(RightImage(1:h,1+Margin:end-Margin), [yShift2, 0]);

    APImage=APImage(:,1:w+wRight+1-xShift2-xShift1);

  
end    
    
   
%Save it to the Dropbox folder
imwrite(uint16(APImage),[DropboxFolder,filesep,Prefix,'\APDetection\FullEmbryo.tif'],'compression','none');


%Now, use them to find the embryo mask
embMask = getEmbryoMask(APImage, 20);


%This code came from Michael's code
% ES 2014-01-07: I changed this to handle cases where the threshold is too
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



