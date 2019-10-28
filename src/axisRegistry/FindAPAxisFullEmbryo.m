function [coordA,coordP,xShift,yShift]=FindAPAxisFullEmbryo(Prefix, varargin)
%%
%When we could fit the whole embryo in one image we don't do stitching.
%Instead we do a correlation between the imaging field and the full embryo
%to determine the shift and then find the AP axis.


%Parameters:
%First, the prefix.
%There after:
%FlipAP- Switches anterior and posterior poles
%CorrectAxis- Runs a correction script after automatic detection

CorrectAxis = 1;
optionalResults = '';

for i=1:length(varargin)
    if isnumeric(varargin{i})
        if varargin{i}==1
            FlipAP=1;
        end
    elseif strcmp(varargin{i},'CorrectAxis')
        CorrectAxis = 1;
    elseif strcmp(varargin{i},'optionalResults')
        optionalResults = varargin{i+1};
    end
end
   
[SourcePath, ~, DefaultDropboxFolder, DropboxFolder, ~, ~,...
~, ~] = DetermineAllLocalFolders(Prefix, optionalResults);

%If the AP axis hasn't been specified check if it was specificed in the
%image names. Otherwise assume AP orientation

%Find out the date it was taken
Dashes=findstr(Prefix,'-');
Date=Prefix(1:Dashes(3)-1);
EmbryoName=Prefix(Dashes(3)+1:end);

%Figure out what type of experiment we have. Note: we name the var "DateFromDateColumn" to avoid shadowing previously defined "Date" var.
[DateFromDateColumn, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF, Channel3,~,~, ~, ~]...
    = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);

%Determine whether we're dealing with 2-photon data from Princeton, LSM, or
%LIF data. 2-photon data uses TIF files. In LSM mode multiple files will be
%combined into one.
DTIF=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*.tif']);
DLSM=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*.lsm']);
DCZI=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*.czi']);
DLIF=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*.lif']);
DSPIN=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*.nd']);     %Nikon spinning disk . CS20170911

if (~isempty(DTIF))&(isempty(DLIF))&(isempty(DSPIN))        %CS20170911
    disp('2-photon @ Princeton data mode')
    D=DTIF;
    FileMode='TIF';
elseif (~isempty(DLIF))
    disp('LIF export mode')
    D=DLIF;
    FileMode='LIFExport';
elseif (~isempty(DLSM))|(~isempty(DCZI))
    disp('LSM mode')
    D=[DLSM,DCZI];
    FileMode='LSM';
elseif (~isempty(DSPIN))        %CS20170911
    disp('Nikon spinning disk mode with .nd files')
    D=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo', filesep,'*.nd']);     %spinning disk with .nd output files
    FileMode='DSPIN';
else
    error('File type not recognized')
end


% Identify the midsagittal image
MidFileIndex=find(~cellfun('isempty',strfind(lower({D.name}),'mid')));
SurfFileIndex=find(~cellfun('isempty',strfind(lower({D.name}),'surf')));

if (length(MidFileIndex)>1)& strcmpi(FileMode, ~'DSPIN')
    error('Too many midsagittal files in FullEmbryo folder')
end



%See if we don't want the default AP orientation
if ~exist('FlipAP', 'var')
    if strcmp(FileMode, 'TIF')||strcmp(FileMode, 'LIFExport')||strcmp(FileMode, 'LSM')
        if strcmp(D(MidFileIndex).name,'PA')
            FlipAP=1;
        else
            FlipAP=0;
        end
    end
    if strcmp(FileMode,'DSPIN')     %CS20170912
        FlipAP=0;
    end
end
     
    ChannelToLoadTemp =contains([Channel1,Channel2,Channel3],'nuclear','IgnoreCase',true);

      if sum(ChannelToLoadTemp) && sum(ChannelToLoadTemp)==1
            HisChannel=find(ChannelToLoadTemp);
        elseif sum(ChannelToLoadTemp) && length(ChannelToLoadTemp)>=2
            ChannelToLoad=find(ChannelToLoadTemp);
            HisChannel = ChannelToLoad(1);
      end
    if isempty(HisChannel)
        error('LIF Mode error: Channel name not recognized. Check MovieDatabase.XLSX')
    end


if strcmp(FileMode,'TIF')
    MidImage=imread([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,D(MidFileIndex).name],2);
elseif strcmp(FileMode,'LIFExport')
       
   
    
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
    if isfolder([SourcePath, filesep, Date, filesep, EmbryoName, filesep, 'MetaData'])
        xml_file_path = dir([SourcePath, filesep, Date, filesep, EmbryoName, filesep, 'MetaData', filesep, '*.xml']);
        xml_file = xml_file_path(1).name;
        xDoc = searchXML([SourcePath, filesep, Date, filesep, EmbryoName, filesep, 'MetaData', filesep, xml_file]);
        zoom_angle = str2double(evalin('base','rot'));
    else 
        warning('No time series metadata found.')
    end
    if isfolder([SourcePath, filesep, Date, filesep, EmbryoName, filesep, 'FullEmbryo', filesep...
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
    
   
    
    %Rotates the full embryo image to match the rotation of the zoomed
    %time series
    zoom_angle = NaN;        %Set defaults to NaN so we can tell if the 
    full_embryo_angle = NaN;     %angles were successfully extracted
    
    LSMMid=bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,D(MidFileIndex).name]);
    LSMSurf = bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,D(SurfFileIndex).name]);
    LSMMeta=LSMMid{:,4};
    LSMMeta2=LSMMid{:,2};
    
    %Look for the image with the largest size. In this way, we avoid
    %loading individual tiles in the case of a tile scan.
    for i=1:size(LSMMid,1)
        SizesImages(i)=size(LSMMid{i,1}{1,1},1);
    end
    [~,ImageCellToUse]=max(SizesImages);
    
    %individual tiles if we're dealing with tile scan. Also, in CZI files,
    %this seems to ensure a high-contrast image as well.
    MidImage=LSMMid{ImageCellToUse,1}{HisChannel,:};
    SurfImage = LSMSurf{ImageCellToUse,1}{HisChannel,:};
    
    %Figure out the rotation of the full embryo image
    %This works for LSM files
    full_embryo_angle = LSMMeta2.get('Recording Rotation #1');
    %If full_embryo_angle is empty, chances are we have a CZI file
    if isempty(full_embryo_angle)
        full_embryo_angle=str2num(LSMMeta2.get('Global HardwareSetting|ParameterCollection|RoiRotation #1'));
    %Check if angle was successfully extracted
    elseif isnan(full_embryo_angle)
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
    %Check if angle was successfully extracted
    elseif isnan(zoom_angle)
        error('Could not extract rotation of FullEmbryo images')
    end

    MidImage = imrotate(MidImage, -zoom_angle + full_embryo_angle);
    SurfImage = imrotate(SurfImage, -zoom_angle + full_embryo_angle);

elseif strcmp(FileMode, 'DSPIN')        %CS20170911 This is really long-winded atm! Need to simplify.
%%  
    %Find and open the mid-saggistal .nd files in the FullEmbryo folder
    if find(~cellfun('isempty',strfind({D.name},'ANT_mid')))>0
        SurfFileIndexAnt = find(~cellfun('isempty',strfind({D.name},'ANT_mid')));
        SurfFileIndexPost = find(~cellfun('isempty',strfind({D.name},'POST_mid')));
        
        SurfImageAnt = bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,D(SurfFileIndexAnt).name]);
        SurfAntMeta = SurfImageAnt{4};
        
        SurfImagePost = bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,D(SurfFileIndexPost).name]);
        SurfPostMeta = SurfImagePost{4};
        
    else
        error('FullEmbryo images mislabeled? Should be of form ANT_mid, POST_mid, ANT_surf, POST_surf')
    end
   
    %Get some useful parameters
    try
        NSlices = SurfAntMeta.getPixelsSizeZ(0).getValue();
        NChannels = SurfAntMeta.getPixelsSizeC(0).getValue();
        PizelSize = SurfAntMeta.getPixelsPhysicalSizeX(0).getValue();
    catch
        NSlices = str2double(SurfAntMeta.getPixelsSizeZ(0));
        NChannels = str2double(SurfAntMeta.getPixelsSizeC(0));
        PizelSize = str2double(SurfAntMeta.getPixelsPhysicalSizeX(0).value);%Size in microns of one pixel. NB for some of the OME metadata you need to add .value on end to get it out
    end
    
    %Figure out which channel to use (His channel)
    HisChannel=find(~cellfun(@isempty,strfind(lower({Channel1{1},Channel2{1}}),'mcherry'))|...
        ~cellfun(@isempty,strfind(lower({Channel1{1},Channel2{1}}),'his')));
    
    %Get the stack of images in the histone channel for ANT and POST images and do a maximum projection
    temp = SurfImageAnt{1}(:,2);
    IndexAntPost = find(~cellfun('isempty',strfind({temp{1:end}},['C=', num2str(HisChannel), '/', num2str(NChannels)])));
    AntImStack = SurfImageAnt{1}(IndexAntPost, 1);
    PostImStack = SurfImagePost{1}(IndexAntPost, 1);
    
    %Max Project His channel for ANT and POST images
    for i = 1:NSlices
        AntImStack2(:,:,i) = AntImStack{i};
        PostImStack2(:,:,i) = PostImStack{i};
    end
    AntMaxProj = max(AntImStack2, [], 3); 
    PostMaxProj = max(PostImStack2, [], 3);
    
    %Write these to files because EmbryoStitchNoMargin (which will stitch them together) wants it that way
    imwrite(AntMaxProj, [SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'AntMaxProj.tif']); 
    imwrite(PostMaxProj, [SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'PostMaxProj.tif']); 
    
    %While at it, EmbryoStitchNoMargin also needs a flatfield. I don't take
    %one at this magnification, so make it up
    FFImage = ones(512,672); %make this automatically find the size later
    imwrite(FFImage, [SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'FF20x.tif']);
    Margin = 0;

    %Define inputs to EmbryoStitchNoMargin, which will stitch the ANT and POST files together and record the xShift and yShift needed to do so
    FFFile=[SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'FF20x.tif'];  
    LeftImageFile = [SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'AntMaxProj.tif'];
    RightImageFile = [SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'PostMaxProj.tif']; 
    
    if FlipAP==0    %If FlipAP not specified in input, assume Anterior image on Left and Posterior on Right. 
        [APImage,xShift,yShift]=EmbryoStitchNoMargin(LeftImageFile,RightImageFile, FFFile, [], [], Margin);
    else          %Otherwise assume the opposite. 
        [APImage,xShift,yShift]=EmbryoStitchNoMargin(RightImageFile, LeftImageFile, FFFile, [], [], Margin);
    end
    
    %Rename for use later in code
    MidImage = APImage; 
    %Display image to check it's stitched correctly
    imshow(imadjust(mat2gray(MidImage,[0 65535])),'DisplayRange',[],'InitialMagnification',100);
    
    %While we're at it, AddParticlePosition also needs the surface embryo,
    %so generate a stiched together max projection of the ANT and POST
    %surface images as well. 
    if find(~cellfun('isempty',strfind({D.name},'ANT_surf')))>0
        SurfFileIndexAnt = find(~cellfun('isempty',strfind({D.name},'ANT_surf.nd')));
        SurfFileIndexPost = find(~cellfun('isempty',strfind({D.name},'POST_surf.nd')));
        SurfImageAnt = bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,D(SurfFileIndexAnt).name]);
        SurfAntMeta = SurfImageAnt{4};
        SurfImagePost = bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,D(SurfFileIndexPost).name]);
        SurfPostMeta = SurfImagePost{4};   
    else
        error('FullEmbryo images mislabeled? Should be of form ANT_mid, POST_mid, ANT_surf, POST_surf')
    end
   
    %Get some useful parameters
    try
        NSlices = SurfAntMeta.getPixelsSizeZ(0).getValue();
        NChannels = SurfAntMeta.getPixelsSizeC(0).getValue();
        PizelSize = SurfAntMeta.getPixelsPhysicalSizeX(0).getValue();
    catch
        NSlices = str2double(SurfAntMeta.getPixelsSizeZ(0));
        NChannels = str2double(SurfAntMeta.getPixelsSizeC(0));
        PizelSize = str2double(SurfAntMeta.getPixelsPhysicalSizeX(0).value);%Size in microns of one pixel. NB for some of the OME metadata you need to add .value on end to get it out
    end
    
    %Get the stack of images in the histone channel for ANT and POST images and do a maximum projection
    temp = SurfImageAnt{1}(:,2);
    IndexAntPost = find(~cellfun('isempty',strfind({temp{1:end}},['C=', num2str(HisChannel), '/', num2str(NChannels)])));
    AntImStack = SurfImageAnt{1}(IndexAntPost, 1);
    PostImStack = SurfImagePost{1}(IndexAntPost, 1);
    
    %Max Project His channel for ANT and POST images
    for i = 1:NSlices
        AntImStack2(:,:,i) = AntImStack{i};
        PostImStack2(:,:,i) = PostImStack{i};
    end
    AntMaxProj = max(AntImStack2, [], 3); 
    PostMaxProj = max(PostImStack2, [], 3);
    
    %Write these to files because EmbryoStitchNoMargin (which will stitch them together) wants it that way
    %The FlatField has already been written, above
    imwrite(AntMaxProj, [SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'AntMaxProj.tif']); 
    imwrite(PostMaxProj, [SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'PostMaxProj.tif']); 

    %Define inputs to EmbryoStitchNoMargin, which will stitch the ANT and POST files together and record the xShift and yShift needed to do so
    LeftImageFile = [SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'AntMaxProj.tif'];
    RightImageFile = [SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'PostMaxProj.tif']; 
    
    if FlipAP==0    %If FlipAP not specified in input, assume Anterior image on Left and Posterior on Right. 
        [APImage,xShift,yShift]=EmbryoStitchNoMargin(LeftImageFile,RightImageFile, FFFile, [], [], Margin);
    else          %Otherwise assume the opposite. 
        [APImage,xShift,yShift]=EmbryoStitchNoMargin(RightImageFile, LeftImageFile, FFFile, [], [], Margin);
    end
    
    %Rename for use later in code
    SurfImage = APImage; 
    
end

%%
%Save it to the Dropbox folder
mkdir([DropboxFolder,filesep,Prefix,filesep,'APDetection'])
imwrite(uint16(MidImage),[DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo.tif'],'compression','none');
imwrite(uint16(SurfImage),[DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryoSurf.tif'],'compression','none');

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

if CC.NumObjects~=1 || CorrectAxis
    CorrectAPAxis(Prefix);
end
