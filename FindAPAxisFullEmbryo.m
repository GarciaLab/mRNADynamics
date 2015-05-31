function [coordA,coordP,xShift,yShift]=FindAPAxisFullEmbryo(varargin)

%When we could fit the whole embryo in one image we don't do stitching.
%Instead we do a correlation between the imaging field and the full embryo
%to determine the shift and then find the AP axis.


%Load the folder information
[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
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
    Dashes=strfind(FolderTemp,filesep);
    Prefix=FolderTemp((Dashes(end)+1):end);
end


%If the AP axis hasn't been specified check if it was specificed in the
%image names. Otherwise assume AP orientation

%Find out the date it was taken
Dashes=findstr(Prefix,'-');
Date=Prefix(1:Dashes(3)-1);
EmbryoName=Prefix(Dashes(3)+1:end);

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



%Determine whether we're dealing with 2-photon data from Princeton or LSM
%data. 2-photon data uses TIF files. In LSM mode multiple files will be
%combined into one.
DTIF=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*.tif']);
DLSM=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*.lsm']);
DLIF=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*.lif']);

if (length(DTIF)>0)
    display('2-photon @ Princeton data mode')
    D=DTIF;
    FileMode='TIF';
elseif (length(DLIF)>0)
    display('LIF export mode')
    D=DLIF;
    FileMode='LIFExport';
elseif (length(DLSM)>0)
    display('LSM mode')
    error('Add support for LSM files')
    D=DLSM;
    FileMode='LSM';
else
    error('File type not recognized')
end


%Find the right and left files that do not correspond to the surface image
MidFileIndex=find(~cellfun('isempty',strfind(lower({D.name}),'mid')));

if (length(MidFileIndex)>1) %| (length(SurfaceFileIndex)>1)
    error('Too many left/right files in FullEmbryo folder')
end


%See if we don't want the default AP orientation
if ~exist('FlipAP')
    if strcmp(D(MidFileIndex).name,'PA')%|strcmp(D(SurfaceFileIndex).name,'PA')
        FlipAP=1;
    else
        FlipAP=0;
    end
end

    

if strcmp(FileMode,'TIF')
    MidImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,D(MidFileIndex).name],2);
elseif strcmp(FileMode,'LIFExport')
    
    %Figure out which channel to use
    HisChannel=find(~cellfun(@isempty,strfind(lower({Channel1{1},Channel2{1}}),'mcherry'))|...
        ~cellfun(@isempty,strfind(lower({Channel1{1},Channel2{1}}),'his')));
    
    if isempty(HisChannel)
        error('LIF Mode error: Channel name not recognized. Check MovieDatabase.XLSX')
    end
    
    %Rotates the full embryo image to match the rotation of the zoomed
    %time series
    zoom_angle = 0;
    full_embryo_angle = 0;
    
    LIFMid=bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,D(MidFileIndex).name]);
    %MidImage=LIFMid{1}{HisChannel,1};
    %By looking at the last image we make sure we're avoiding the
    %individual tiles if we're dealing with tile scan
    MidImage=LIFMid{end,1}{HisChannel,1};
    
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
            filesep, 'MetaData', filesep,'*Mid*.xml']);
        xml_file2 = xml_file_path2(1).name;
        xDoc2 = searchXML([SourcePath, filesep, Date, filesep, EmbryoName, filesep,'FullEmbryo', filesep,...
                'MetaData', filesep, xml_file2]);
        full_embryo_angle = str2double(evalin('base','rot'));
    else 
        warning('No full embryo metadata found.')
    end
    
    evalin('base','clear rot')
    MidImage = imrotate(MidImage, -zoom_angle + full_embryo_angle);
end

%Get the necessary information to load the corresponding flat field image
% 
% %Get the structure with the acquisition information
% ImageInfo = imfinfo([SourcePath,filesep,Date,filesep,EmbryoName,filesep,...
%     'FullEmbryo',filesep,D(SurfaceFileIndex).name]);
% 
% 
% %Get the flat-field information
% 
% %Figure out the zoom factor
% Zoom=ExtractInformationField(ImageInfo(1),'state.acq.zoomFactor=');
% if ~isempty(Zoom)
%     %Look for the file
%     FFDir=dir([SourcePath,filesep,Date,filesep,'FF',Zoom(1),'x*.*']);
%     %If there's more than one match then ask for help
%     if length(FFDir)==1
%         FFFile=FFDir(1).name;
%     elseif isempty(FFDir)
%         display('Warning, no flat field file found. Press any key to proceed without it');
%         FFImage=ones(ImageInfo(1).Height,ImageInfo(1).Width);
%         pause
%     else
%         FFFile=uigetfile([Folder,filesep,'..',filesep,'FF',Zoom(1),'x*.*'],'Select flatfield file');
%     end
% else
%     warning('No zoom information found. Cannot process flat field.')
% end
% 
% %Now do the image correlation. We'll grab a small area around the middle of
% %the imaging field and correlate it to the larger, full embryo image.
% 
% % AcqImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,D(AcqImageIndex).name]);
% % FullImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,D(FullImageIndex).name]);

% 
% %Crop the imaging field by 50%
% [Rows,Columns]=size(AcqImage);
% AcqImageCrop=AcqImage(Rows/4:Rows/4*3,Columns/4:Columns/4*3);
% 
% %Calculate the correlation marrix and find the maximum
% C = normxcorr2(AcqImageCrop, FullImage);
% [Max2,MaxRows]=max(C);
% [Dummy,MaxColumn]=max(Max2);
% MaxRow=MaxRows(MaxColumn);
% 
% [CRows,CColumns]=size(C);
% 
% 
% 
% ShiftRow=MaxRow-(CRows/2+1);
% ShiftColumn=MaxColumn-(CColumns/2+1);
% 
% %Create an overlay to make sure things make sense
% 
% ImShifted=uint16(zeros(size(FullImage)));
% 
% RowRange=(Rows/4:Rows/4*3)+ShiftRow;
% ColumnRange=(Columns/4:Columns/4*3)+ShiftColumn;
% 
% ImShifted(RowRange,ColumnRange)=AcqImageCrop;
% 
% figure(1)
% ImOverlay=cat(3,mat2gray(FullImage)+mat2gray(ImShifted),mat2gray(FullImage),mat2gray(FullImage));
% imshow(ImOverlay)
% 
% xShift=-ShiftColumn;
% yShift=-ShiftRow;
% 
% 
% %Save it to the Dropbox folder
% imwrite(FullImage,[DropboxFolder,filesep,Prefix,filesep,'FullEmbryo.tif'],'compression','none');
% 


%Save it to the Dropbox folder
mkdir([DropboxFolder,filesep,Prefix,filesep,'APDetection'])
imwrite(uint16(MidImage),[DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo.tif'],'compression','none');

%Now, use them to find the embryo mask
embMask = getEmbryoMaskLive(MidImage, 50);


%This code came from Michael's code
diagFigure = figure;
CC=bwconncomp(embMask);
if CC.NumObjects~=1
    warning('Failed to calculate embryo mask. Found more than one object in mask. Assigning arbitrary A and P positions.');
    coordA=[1,1];
    coordP=[1,1];
else
    
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
    
    % Save diagnostic figures to check the quality of axis determination
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
    saveas(gcf, [DropboxFolder,filesep,Prefix,filesep,'APMask.tif']);

end
    



%Save the AP and shift information
save([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'],'coordA','coordP');




clf
imagesc(MidImage)
axis image
axis off
title('Anterior (green), posterior (red); original')
hold on
plot(coordA(1),coordA(2),'g.','MarkerSize',20);
plot(coordP(1),coordP(2),'r.','MarkerSize',20);
hold off
saveas(gcf, [DropboxFolder,filesep,Prefix,filesep,'APEmbryo.tif']);
close(diagFigure);



