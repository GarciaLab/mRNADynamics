function [coordA,coordP,xShift,yShift]=FindAPAxisFullEmbryo(varargin)

%When we could fit the whole embryo in one image we don't do stitching.
%Instead we do a correlation between the imaging field and the full embryo
%to determine the shift and then find the AP axis.


%Parameters:
%First, the prefix.
%There after:
%FlipAP- Switches anterior and posterior poles
%CorrectAxis- Runs a correction script after automatic detection

CorrectAxis = 0;

%Load the folder information
[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
    DetermineLocalFolders(varargin{1});

Prefix=varargin{1};

for i=2:length(varargin)
    if isnumeric(varargin{i})
        if varargin{i}==1
            FlipAP=1;
        end
    elseif strcmp(varargin{i},'CorrectAxis')
        CorrectAxis = 1;
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



%Determine whether we're dealing with 2-photon data from Princeton, LSM, or
%LIF data. 2-photon data uses TIF files. In LSM mode multiple files will be
%combined into one.
DTIF=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*.tif']);
DLSM=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*.lsm']);
DCZI=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*.czi']);
DLIF=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*.lif']);

if (length(DTIF)>0)&(length(DLIF)==0)
    display('2-photon @ Princeton data mode')
    D=DTIF;
    FileMode='TIF';
elseif (length(DLIF)>0)
    display('LIF export mode')
    D=DLIF;
    FileMode='LIFExport';
elseif (length(DLSM)>0)|(length(DCZI)>0)
    display('LSM mode')
    D=[DLSM,DCZI];
    FileMode='LSM';
else
    error('File type not recognized')
end


% Identify the midsagittal image
MidFileIndex=find(~cellfun('isempty',strfind(lower({D.name}),'mid')));
SurfFileIndex=find(~cellfun('isempty',strfind(lower({D.name}),'surf')));

if (length(MidFileIndex)>1)
    error('Too many midsagittal files in FullEmbryo folder')
end


%See if we don't want the default AP orientation
if ~exist('FlipAP')
    if strcmp(D(MidFileIndex).name,'PA')
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
    LIFSurf=bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,D(SurfFileIndex).name]);

    %By looking at the last image we make sure we're avoiding the
    %individual tiles if we're dealing with tile scan
    MidImage=LIFMid{end,1}{HisChannel,1};
    SurfImage=LIFSurf{end,1}{HisChannel,1};
    if size(MidImage) ~= size(SurfImage)
            MidImage = imresize(MidImage,length(SurfImage)/length(MidImage));
    end
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
        evalin('base','clear rot')
        xDoc2 = searchXML([SourcePath, filesep, Date, filesep, EmbryoName, filesep,'FullEmbryo', filesep,...
                'MetaData', filesep, xml_file2]);
         full_embryo_angle = str2double(evalin('base','rot'));
    else 
        warning('No full embryo metadata found.')
    end
    
    evalin('base','clear rot')
    MidImage = imrotate(MidImage, -zoom_angle + full_embryo_angle);
elseif strcmp(FileMode,'LSM')
    
    %Figure out which channel to use
    HisChannel=find(~cellfun(@isempty,strfind(lower({Channel1{1},Channel2{1}}),'mcherry'))|...
        ~cellfun(@isempty,strfind(lower({Channel1{1},Channel2{1}}),'his')));
    
    if isempty(HisChannel)
        error('LSM Mode error: Channel name not recognized. Check MovieDatabase.XLSX')
    end
    
    %Rotates the full embryo image to match the rotation of the zoomed
    %time series
    zoom_angle = 0;
    full_embryo_angle = 0;
    
    LSMMid=bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,D(MidFileIndex).name]);
    LSMMeta=LSMMid{:,4};
    LSMMeta2=LSMMid{:,2};
    
    %The first image in a tile scane on a size scope seems to be the one
    %that preserves pixel size.
    %individual tiles if we're dealing with tile scan. Also, in CZI files,
    %this seems to ensure a high-contrast image as well.
    MidImage=LSMMid{end,1}{HisChannel,:};
    
    %Figure out the rotation of the full embryo image
    %This works for LSM files
    full_embryo_angle = LSMMeta2.get('Recording Rotation #1');
    %If full_embryo_angle is empty, chances are we have a CZI file
    if isempty(full_embryo_angle)
        full_embryo_angle=str2num(LSMMeta2.get('Global HardwareSetting|ParameterCollection|RoiRotation #1'));
    else
        error('Could not extract rotation of FullEmbryo images')
    end
    
    
    %Figure out the rotation of the zoomed-in image. We need to check for
    %both LSM and CZI files.
    DLSMZoom=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'*.lsm']);
    DLSMZoom=[DLSMZoom,...
        dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'*.czi'])];
    LSMZoom=bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,...
        DLSMZoom(1).name]);
    LSMMetaZoom2=LSMZoom{:,2};
    
    %This works for LSM files
    zoom_angle=LSMMetaZoom2.get('Recording Rotation #1');
    %If full_embryo_angle is empty, chances are we have a CZI file
    if isempty(zoom_angle)
        zoom_angle=str2num(LSMMetaZoom2.get('Global HardwareSetting|ParameterCollection|RoiRotation #1'));
    else
        error('Could not extract rotation of FullEmbryo images')
    end

    MidImage = imrotate(MidImage, -zoom_angle + full_embryo_angle);
end


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

if CorrectAxis
    CorrectAPAxis(Prefix);
end
