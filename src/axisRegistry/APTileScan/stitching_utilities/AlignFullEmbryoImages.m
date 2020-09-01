function AlignFullEmbryoImages(Prefix)
% AlignFullEmbryoImages.m
% author: Gabriella Martini
% date created: 8/22/20
% date last modified: 8/23/20


%% 

[SourcePath, ~, DefaultDropboxFolder, DropboxFolder, ~, ~,...
    ~, ~] = DetermineAllLocalFolders(Prefix);
stitchingDataFolder = [DropboxFolder,filesep,Prefix,filesep,'FullEmbryoStitching'];

if exist([stitchingDataFolder, filesep,'SurfTileArray.mat'], 'file')
     SurfData = load([stitchingDataFolder, filesep, 'SurfTileArray.mat']);
     SurfTileArray = SurfData.tile_array;
else
     error('No surface TileArray data stored. Seed a new TileArray using "NewTileArrayFromMetadata".')  
end



if exist([stitchingDataFolder, filesep,'MidTileArray.mat'], 'file')
     MidData = load([stitchingDataFolder, filesep, 'MidTileArray.mat']);
     MidTileArray = MidData.tile_array;
else
     error('No midsaggital plane TileArray data stored. Seed a new TileArray using "NewTileArrayFromMetadata".')  
end


%% 



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


[FileMode, EmbryoName, projectDate] = getMicroscope(Prefix);
rawPrefixPath = [SourcePath,filesep,projectDate,filesep,EmbryoName,filesep];
fullEmbryoPath = [rawPrefixPath, 'FullEmbryo', filesep];


%Datatype is hardcoded in, unlike in FindAPAxisFullEmbryo
DLIF=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*.lif']);
FileMode='LIFExport';

DSTITCH=dir([stitchingDataFolder, filesep, '*.tif']);

%% 


% Identify the midsagittal image
MidLifFileIndex=find(~cellfun('isempty',strfind(lower({DLIF.name}),'midtile')));
SurfLifFileIndex=find(~cellfun('isempty',strfind(lower({DLIF.name}),'surftile')));
MidFileIndex=find(~cellfun('isempty',strfind(lower({DSTITCH.name}),'midtilestitch_max.tif')));
SurfFileIndex=find(~cellfun('isempty',strfind(lower({DSTITCH.name}),'surftilestitch_max.tif')));



%% 

MatchLifFileIndex=find(~cellfun('isempty',strfind(lower({DLIF.name}),'fullmatchstack.lif')));
if ~isempty(MatchLifFileIndex)
    MatchFilePresent = true;
else
    MatchFilePresent = false;
end

MatchPostLifFileIndex=find(~cellfun('isempty',strfind(lower({DLIF.name}),'fullmatchpost.lif')));
if ~isempty(MatchPostLifFileIndex)
    PostFilePresent = true;
else
    PostFilePresent = false;
end






%% %Rotates the full embryo image and match stack image to match the rotation of zoomed data
mid_embryo_angle = 0;
surf_embryo_angle = 0;




LIFMid=bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,DLIF(MidLifFileIndex).name]);
LIFSurf=bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,DLIF(SurfLifFileIndex).name]);

MIDMeta = LIFMid{:, 4};
SURFMeta = LIFSurf{:,4};

% Gets rotation information for full embryo images. 
if isfolder([fullEmbryoPath,'MetaData'])
    xml_file_path = dir([fullEmbryoPath,'MetaData', filesep,'*Mid*.xml']);
    xml_file = xml_file_path(1).name;
    xDoc = xmlread([fullEmbryoPath,'MetaData', filesep, xml_file]);
    xml_struct = xml2struct([fullEmbryoPath,'MetaData', filesep, xml_file]);
    attach_idx = size(xml_struct.Data.Image.Attachment, 2);
    mid_embryo_angle = str2double(xml_struct.Data.Image.Attachment{1,attach_idx}.ATLConfocalSettingDefinition.Attributes.RotatorAngle);
else
    warning('No full embryo metadata found.')
end

evalin('base','clear rot')
if isfolder([fullEmbryoPath,'MetaData'])
    xml_file_path2 = dir([fullEmbryoPath,'MetaData', filesep,'*Surf*.xml']);
    xml_file2 = xml_file_path2(1).name;
    xDoc2 = xmlread([fullEmbryoPath,'MetaData', filesep, xml_file2]);
    xml_struct2 = xml2struct([fullEmbryoPath,'MetaData', filesep, xml_file2]);
    attach_idx = size(xml_struct2.Data.Image.Attachment, 2);
    surf_embryo_angle = str2double(xml_struct2.Data.Image.Attachment{1,attach_idx}.ATLConfocalSettingDefinition.Attributes.RotatorAngle);
else
    warning('No full embryo metadata found.')
end

MidImageNoPadding=imread([stitchingDataFolder,filesep,DSTITCH(MidFileIndex).name]);
SurfImageNoPadding=imread([stitchingDataFolder,filesep,DSTITCH(SurfFileIndex).name]);
MidFilteredImageNoPadding = imread([stitchingDataFolder,filesep,'MidTileStitch_MaxFiltered.tif']);
SurfFilteredImageNoPadding = imread([stitchingDataFolder,filesep,'SurfTileStitch_MaxFiltered.tif']);

if mid_embryo_angle ~= 0
    MidImageNoPadding = imrotate(MidImageNoPadding, mid_embryo_angle);
    MidFilteredImageNoPadding = imrotate(MidFilteredImageNoPadding, mid_embryo_angle);
end
if surf_embryo_angle ~= 0
    SurfImageNoPadding = imrotate(SurfImageNoPadding, surf_embryo_angle);
    SurfFilteredImageNoPadding = imrotate(SurfFilteredImageNoPadding, surf_embryo_angle);
end

%% 

%% 




[SurfNTiles, SurfNFrames, SurfNSlices, SurfNPlanes, SurfNChannels,...
    SurfFrame_Times] = getFrames(SURFMeta);

[MidNTiles, MidNFrames, MidNSlices, MidNPlanes, MidNChannels,...
    MidFrame_Times] = getFrames(MIDMeta);

PixelSize = double(MIDMeta.getPixelsPhysicalSizeX(0).value);% units: microns
PixelSize_m = double(PixelSize)*10^(-6);
%% 

if MatchFilePresent
    match_stack_angle = 0;
    LIFMatch=bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,DLIF(MatchLifFileIndex).name]);
    MatchMeta = LIFMatch{:, 4};
    [MatchNTiles, MatchNFrames, MatchNSlices, MatchNPlanes, MatchNChannels,...
        MatchFrame_Times] = getFrames(MatchMeta);
    InitialStackTime = [];
    zPosition = [];
    MatchFrameInfo = recordFrameInfo(MatchNFrames, MatchNSlices, InitialStackTime, MatchMeta, zPosition);
    [coatChannel, histoneChannel, fiducialChannel, inputProteinChannel, FrameInfo] =...
    LIFExportMode_interpretChannels(ExperimentType, Channel1, Channel2, Channel3, MatchFrameInfo);
    framesIndex  = 1;
    matchStack = zeros(size(LIFMatch{1}{1,1}, 1), size(LIFMatch{1}{1,1}, 2), MatchNSlices, 'uint16');
    z = 1;
    firstImage = fiducialChannel;
    lastImage = framesIndex * MatchNSlices * MatchNChannels;
    for imagesIndex = firstImage:MatchNChannels:lastImage
        matchStack(:, :, z) = LIFMatch{1}{imagesIndex, 1};
        z = z + 1;
    end
    FullStackMaxProj = max(matchStack, [], 3);
    MidMatchMaxProj = max(matchStack(:,:,(MatchNSlices-4):MatchNSlices), [], 3);
    SurfMatchMaxProj = max(matchStack(:,:,1:5), [], 3);
    try 
        PixelSize_matchstack = double(MatchMeta.getPixelsPhysicalSizeX(0).value);% units: microns
    catch 
        PixelSize_matchstack = double(MatchMeta.getPixelsPhysicalSizeX(0));% units: microns
    end
    if isfolder([fullEmbryoPath,'MetaData'])
        xml_file_path3 = dir([fullEmbryoPath,'MetaData', filesep,'*FullMatchStack*.xml']);
        xml_file3 = xml_file_path3(1).name;
        xDoc3 = xmlread([fullEmbryoPath,'MetaData', filesep, xml_file3]);
        xml_struct3 = xml2struct([fullEmbryoPath,'MetaData', filesep, xml_file3]);
        attach_idx = size(xml_struct3.Data.Image.Attachment, 2);
        match_stack_angle = str2double(xml_struct3.Data.Image.Attachment{1,attach_idx}.ATLConfocalSettingDefinition.Attributes.RotatorAngle);
    else
        warning('No full embryo metadata found.')
    end
    

    
    MidMatchRot = imrotate(MidMatchMaxProj, match_stack_angle);
    SurfMatchRot = imrotate(SurfMatchMaxProj, match_stack_angle);
    
    ZoomRatio = PixelSize(1)/PixelSize_matchstack(1);
    
    if ZoomRatio <= 1
        im1_mid = imresize(MidMatchRot, 1/ZoomRatio);
        im1_surf = imresize(SurfMatchRot, 1/ZoomRatio);
        im2_mid = MidImageNoPadding;
        im2_surf = SurfImageNoPadding;
    else
        %Enlarge the zoomed out image so we can do the cross-correlation
        im1_mid = MidMatchRot;
        im1_surf = SurfMatchRot;
        im2_mid = imresize(MidImageNoPadding, ZoomRatio);
        im2_surf = imresize(SurfImageNoPadding, ZoomRatio);
    end

    C_mid = gather(normxcorr2(gpuArray(im1_mid), gpuArray(im2_mid)));
    [MaxRow_mid, MaxColumn_mid] = find(ismember(C_mid, max(C_mid(:)))); 
    C_surf= normxcorr2(im1_surf, im2_surf);
    C_surf = gather(normxcorr2(gpuArray(im1_surf), gpuArray(im2_surf)));
    [MaxRow_surf, MaxColumn_surf] = find(ismember(C_surf, max(C_surf(:)))); 

 
    %% 
    
    if ZoomRatio <= 1
        DeltaR = MaxRow_surf-MaxRow_mid;
        DeltaC = MaxColumn_surf-MaxColumn_mid;
    else
       DeltaR = round((MaxRow_surf-MaxRow_mid)/ZoomRatio, 0);
       DeltaC = round((MaxColumn_surf-MaxColumn_mid)/ZoomRatio, 0);
    end

    if DeltaR > 0
        if DeltaC > 0
            MidImageTempPadding = zeros(size(MidImageNoPadding, 1) + DeltaR,...
                                         size(MidImageNoPadding, 2) + DeltaC,...   
                                         'uint16');
            MidImageTempPadding((DeltaR+1):size(MidImageTempPadding, 1),...
                (DeltaC+1):size(MidImageTempPadding, 2)) = MidImageNoPadding;
            SurfImageTempPadding = SurfImageNoPadding;
            
            MidImageFilteredTempPadding = zeros(size(MidImageNoPadding, 1) + DeltaR,...
                                         size(MidImageNoPadding, 2) + DeltaC,...   
                                         'uint16');
            MidImageFilteredTempPadding((DeltaR+1):size(MidImageTempPadding, 1),...
                (DeltaC+1):size(MidImageTempPadding, 2)) = MidFilteredImageNoPadding;
            SurfImageFilteredTempPadding = SurfFilteredImageNoPadding;
        else
            MidImageTempPadding = zeros(size(MidImageNoPadding, 1) + DeltaR,...
                                         size(MidImageNoPadding, 2),...   
                                         'uint16');
            MidImageTempPadding((DeltaR+1):size(MidImageTempPadding, 1),...
                1:size(MidImageTempPadding, 2)) = MidImageNoPadding;
            SurfImageTempPadding = zeros(size(SurfImageNoPadding, 1),...
                                         size(SurfImageNoPadding, 2)-DeltaC,...   
                                         'uint16');
            SurfImageTempPadding(1:size(SurfImageTempPadding, 1),...
                (1-DeltaC):size(SurfImageTempPadding, 2)) = SurfImageNoPadding;
            
            
            MidImageFilteredTempPadding = zeros(size(MidImageNoPadding, 1) + DeltaR,...
                                         size(MidImageNoPadding, 2),...   
                                         'uint16');
            MidImageFilteredTempPadding((DeltaR+1):size(MidImageTempPadding, 1),...
                1:size(MidImageTempPadding, 2)) = MidFilteredImageNoPadding;
            SurfImageFilteredTempPadding = zeros(size(SurfImageNoPadding, 1),...
                                         size(SurfImageNoPadding, 2)-DeltaC,...   
                                         'uint16');
            SurfImageFilteredTempPadding(1:size(SurfImageTempPadding, 1),...
                (1-DeltaC):size(SurfImageTempPadding, 2)) = SurfFilteredImageNoPadding;
        end
    else
        if DeltaC > 0
            MidImageTempPadding = zeros(size(MidImageNoPadding, 1),...
                                         size(MidImageNoPadding, 2)+DeltaC,...   
                                         'uint16');
            MidImageTempPadding(1:size(MidImageTempPadding, 1),...
                (DeltaC+1):size(MidImageTempPadding, 2)) = MidImageNoPadding;
            SurfImageTempPadding = zeros(size(SurfImageNoPadding, 1)-DeltaR,...
                                         size(SurfImageNoPadding, 2),...   
                                         'uint16');
            SurfImageTempPadding((1-DeltaR):size(SurfImageTempPadding, 1),...
                1:size(SurfImageTempPadding, 2)) = SurfImageNoPadding;
            
            MidImageFilteredTempPadding = zeros(size(MidImageNoPadding, 1),...
                                         size(MidImageNoPadding, 2)+DeltaC,...   
                                         'uint16');
            MidImageFilteredTempPadding(1:size(MidImageTempPadding, 1),...
                (DeltaC+1):size(MidImageTempPadding, 2)) = MidFilteredImageNoPadding;
            SurfImageFilteredTempPadding = zeros(size(SurfImageNoPadding, 1)-DeltaR,...
                                         size(SurfImageNoPadding, 2),...   
                                         'uint16');
            SurfImageFilteredTempPadding((1-DeltaR):size(SurfImageTempPadding, 1),...
                1:size(SurfImageTempPadding, 2)) = SurfFilteredImageNoPadding;
        else
            MidImageTempPadding = MidImageNoPadding;
            SurfImageTempPadding = zeros(size(SurfImageNoPadding, 1) - DeltaR,...
                                         size(SurfImageNoPadding, 2) - DeltaC,...   
                                         'uint16');
            SurfImageTempPadding((-DeltaR+1):size(SurfImageTempPadding, 1),...
                (-DeltaC+1):size(SurfImageTempPadding, 2)) = SurfImageNoPadding;
            
            MidImageFilteredTempPadding = MidFilteredImageNoPadding;
            SurfImageFilteredTempPadding = zeros(size(SurfImageNoPadding, 1) - DeltaR,...
                                         size(SurfImageNoPadding, 2) - DeltaC,...   
                                         'uint16');
            SurfImageFilteredTempPadding((-DeltaR+1):size(SurfImageTempPadding, 1),...
                (-DeltaC+1):size(SurfImageTempPadding, 2)) = SurfFilteredImageNoPadding;
        end 
    end
    
    full_embryo_nrows = max([size(MidImageTempPadding, 1), size(SurfImageTempPadding, 1)]);
    full_embryo_ncols = max([size(MidImageTempPadding, 2), size(SurfImageTempPadding, 2)]);
    MidMaxImage = zeros(full_embryo_nrows, full_embryo_ncols, 'uint16');
    SurfMaxImage =  zeros(full_embryo_nrows, full_embryo_ncols, 'uint16');
    MidFilteredImage = zeros(full_embryo_nrows, full_embryo_ncols, 'uint16');
    SurfFilteredImage = zeros(full_embryo_nrows, full_embryo_ncols, 'uint16');
    MidMaxImage(1:size(MidImageTempPadding, 1), 1:size(MidImageTempPadding, 2)) = MidImageTempPadding;
    SurfMaxImage(1:size(SurfImageTempPadding, 1), 1:size(SurfImageTempPadding, 2)) = SurfImageTempPadding;
    MidFilteredImage(1:size(MidImageFilteredTempPadding, 1), 1:size(MidImageFilteredTempPadding, 2)) = MidImageFilteredTempPadding;
    SurfFilteredImage(1:size(SurfImageFilteredTempPadding, 1), 1:size(SurfImageFilteredTempPadding, 2)) = SurfImageFilteredTempPadding;
    
    if PostFilePresent 
        post_match_angle = 0;
        LIFPostMatch=bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,DLIF(MatchPostLifFileIndex).name]);
        PostMatchMeta = LIFPostMatch{:, 4};
        [PostMatchNTiles, PostMatchNFrames, PostMatchNSlices, PostMatchNPlanes, PostMatchNChannels,...
            PostMatchFrame_Times] = getFrames(PostMatchMeta);
        PostInitialStackTime = [];
        PostzPosition = [];
        PostMatchFrameInfo = recordFrameInfo(PostMatchNFrames, PostMatchNSlices, PostInitialStackTime, PostMatchMeta, PostzPosition);
        [coatChannel, histoneChannel, fiducialChannel, inputProteinChannel, FrameInfo] =...
        LIFExportMode_interpretChannels(ExperimentType, Channel1, Channel2, Channel3, PostMatchFrameInfo);
        framesIndex  = 1;
        PostMatchStack = zeros(size(LIFPostMatch{1}{1,1}, 1), size(LIFPostMatch{1}{1,1}, 2), PostMatchNSlices, 'uint16');
        z = 1;
        firstImage = fiducialChannel;
        lastImage = framesIndex * PostMatchNSlices * PostMatchNChannels;
        for imagesIndex = firstImage:PostMatchNChannels:lastImage
            PostMatchStack(:, :, z) = LIFPostMatch{1}{imagesIndex, 1};
            z = z + 1;
        end
        FullPostStackMaxProj = max(PostMatchStack, [], 3);
        MidPostMatchMaxProj = max(PostMatchStack(:,:,(PostMatchNSlices-4):PostMatchNSlices), [], 3);
        SurfPostMatchMaxProj = max(PostMatchStack(:,:,1:4), [], 3);
        try 
            PixelSize_postmatchstack = double(PostMatchMeta.getPixelsPhysicalSizeX(0).value);% units: microns
        catch 
            PixelSize_postmatchstack = double(PostMatchMeta.getPixelsPhysicalSizeX(0));% units: microns
        end
        if isfolder([fullEmbryoPath,'MetaData'])
            xml_file_path4 = dir([fullEmbryoPath,'MetaData', filesep,'*FullMatchPost*.xml']);
            xml_file4 = xml_file_path4(1).name;
            xDoc4 = xmlread([fullEmbryoPath,'MetaData', filesep, xml_file4]);
            xml_struct4 = xml2struct([fullEmbryoPath,'MetaData', filesep, xml_file4]);
            attach_idx = size(xml_struct4.Data.Image.Attachment, 2);
            post_match_angle = str2double(xml_struct4.Data.Image.Attachment{1,attach_idx}.ATLConfocalSettingDefinition.Attributes.RotatorAngle);
        else
            warning('No full embryo metadata found.')
        end



        MidPostMatchRot = imrotate(MidPostMatchMaxProj, post_match_angle);
        SurfPostMatchRot = imrotate(SurfPostMatchMaxProj, post_match_angle);

        PostZoomRatio = PixelSize(1)/PixelSize_matchstack(1);

        if PostZoomRatio <= 1
            im1_mid_post = imresize(MidPostMatchRot, 1/ZoomRatio);
            im1_surf_post = imresize(SurfPostMatchRot, 1/ZoomRatio);
        else
            %Enlarge the zoomed out image so we can do the cross-correlation
            im1_mid_post = MidPostMatchRot;
            im1_surf_post = SurfPostMatchRot;
            
        end

        C_surf_post= normxcorr2(im1_surf_post, im2_surf);
        C_surf_post = gather(normxcorr2(gpuArray(im1_surf_post), gpuArray(im2_surf)));
        [MaxRow_surf_post, MaxColumn_surf_post] = find(ismember(C_surf_post, max(C_surf_post(:)))); 
        RowMinsTemp = [1, MaxRow_surf-(size(im1_surf, 1)-1),  MaxRow_surf_post-(size(im1_surf_post, 1)-1)];
        RowMaxsTemp = [size(im2_surf, 1), MaxRow_surf,  MaxRow_surf_post];
        ColumnMinsTemp = [1, MaxColumn_surf-(size(im1_surf, 2)-1),  MaxColumn_surf_post-(size(im1_surf_post, 2)-1)];
        ColumnMaxsTemp = [size(im2_surf, 2), MaxColumn_surf,  MaxColumn_surf_post];
        if (all(RowMinsTemp > 0) & all(ColumnMinsTemp > 0) & (RowMaxsTemp(1) == max(RowMaxsTemp)) & (ColumnMaxsTemp(1) == max(ColumnMaxsTemp)))
            MidEmbryoImageV2 = im2_surf;
            MidEmbryoImageV2(RowMinsTemp(2):RowMaxsTemp(2), ColumnMinsTemp(2):ColumnMaxsTemp(2)) = ...
                MidEmbryoImageV2(RowMinsTemp(2):RowMaxsTemp(2), ColumnMinsTemp(2):ColumnMaxsTemp(2)) + im1_mid;
            MidEmbryoImageV2(RowMinsTemp(3):RowMaxsTemp(3), ColumnMinsTemp(3):ColumnMaxsTemp(3)) = ...
                MidEmbryoImageV2(RowMinsTemp(3):RowMaxsTemp(3), ColumnMinsTemp(3):ColumnMaxsTemp(3)) + im1_mid_post;
        else
            error('Stack match image boundaries not supported.');
        end
        
        if ZoomRatio > 1
            MidEmbryoImageV2 = imresize(MidEmbryoImageV2, 1/ZoomRatio);
        end


        if DeltaR > 0
            if DeltaC > 0
                MidEmbryoImageV2TempPadding = MidEmbryoImageV2;
            else

                MidEmbryoImageV2TempPadding = zeros(size(MidEmbryoImageV2, 1),...
                                             size(MidEmbryoImageV2, 2)-DeltaC,...   
                                             'uint16');
                MidEmbryoImageV2TempPadding(1:size(MidEmbryoImageV2TempPadding, 1),...
                    (1-DeltaC):size(MidEmbryoImageV2TempPadding, 2)) = MidEmbryoImageV2;

            end
        else
            if DeltaC > 0
                MidEmbryoImageV2TempPadding = zeros(size(MidEmbryoImageV2, 1)-DeltaR,...
                                             size(MidEmbryoImageV2, 2),...   
                                             'uint16');
                MidEmbryoImageV2TempPadding((1-DeltaR):size(MidEmbryoImageV2TempPadding, 1),...
                    1:size(MidEmbryoImageV2TempPadding, 2)) = MidEmbryoImageV2;

            else
                MidEmbryoImageV2TempPadding = zeros(size(MidEmbryoImageV2, 1) - DeltaR,...
                                             size(MidEmbryoImageV2, 2) - DeltaC,...   
                                             'uint16');
                MidEmbryoImageV2TempPadding((-DeltaR+1):size(MidEmbryoImageV2TempPadding, 1),...
                    (-DeltaC+1):size(MidEmbryoImageV2TempPadding, 2)) = MidEmbryoImageV2;
            end 
        end


        MidEmbryoMaxImageV2 =  zeros(full_embryo_nrows, full_embryo_ncols, 'uint16');
        MidEmbryoMaxImageV2(1:size(MidEmbryoImageV2TempPadding, 1), 1:size(MidEmbryoImageV2TempPadding, 2)) = MidEmbryoImageV2TempPadding;
       
        
    end
    
else
    surf_ypos = [];
    surf_xpos = [];
    for i=0:(SurfNTiles-1)
        surf_xpos(length(surf_xpos)+1) = -1*double(SURFMeta.getPlanePositionX(i,0).value);
        surf_ypos(length(surf_ypos)+1) = double(SURFMeta.getPlanePositionY(i,0).value);
    end
    surfxbound = find(surf_xpos == min(surf_xpos));
    surfybound = find(surf_ypos == min(surf_ypos));
    surfreftile = intersect(surfxbound, surfybound);

    mid_ypos = [];
    mid_xpos = [];
    for i=0:(MidNTiles-1)
        mid_xpos(length(mid_xpos)+1) = -1*double(MIDMeta.getPlanePositionX(i,0).value);
        mid_ypos(length(mid_ypos)+1) = double(MIDMeta.getPlanePositionY(i,0).value);
    end
    midxbound = find(mid_xpos == min(mid_xpos));
    midybound = find(mid_ypos == min(mid_ypos));
    midreftile = intersect(midxbound, midybound);
    DeltaR = round((surf_xpos(surfreftile)-mid_xpos(midreftile))/PixelSize_m, 0);
    DeltaC = round((surf_ypos(surfreftile)-mid_ypos(midreftile))/PixelSize_m, 0);
    
    if ~all(size(MidImageNoPadding) == size(SurfImageNoPadding))
        nrows = max(size(MidImageNoPadding, 1), size(SurfImageNoPadding, 1));
        ncols = max(size(MidImageNoPadding, 2), size(SurfImageNoPadding, 2));
        MidMaxImage = zeros(nrows, ncols, 'uint16');
        SurfMaxImage = zeros(nrows, ncols, 'uint16');
        MidFilteredImage = zeros(nrows, ncols, 'uint16');
        SurfFilteredImage = zeros(nrows, ncols, 'uint16');
        [hm, wm] = size(MidImageNoPadding, 1:2);
        [hs, ws] = size(SurfImageNoPadding, 1:2);
        if DeltaR > 0
            if DeltaC > 0
                MidMaxImage = MidImageNoPadding;
                SurfMaxImage(DeltaR + 1:DeltaR + hs, DeltaC + 1:DeltaC + ws) = SurfImageNoPadding;
                MidFilteredImage = MidImageFilteredNoPadding;
                SurfFilteredImage(DeltaR + 1:DeltaR + hs, DeltaC + 1:DeltaC + ws) = SurfImageFilteredNoPadding;
            else
                MidMaxImage(1:hm, DeltaC + 1:DeltaC+wm) = MidImageNoPadding;
                SurfMaxImage(DeltaR + 1:DeltaR + hs, 1:ws) = SurfImageNoPadding;
                MidFilteredImage(1:hm, DeltaC + 1:DeltaC+wm) = MidImageFilteredNoPadding;
                SurfFilteredImage(DeltaR + 1:DeltaR + hs, 1:ws) = SurfImageFilteredNoPadding;

            end
        else
            if DeltaC > 0
                MidMaxImage(DeltaR+1:DeltaR+hm, 1:wm) = MidImageNoPadding;
                SurfMaxImage(1:hs, DeltaC+1:DeltaC+ws) = SurfImageNoPadding;
                MidFilteredImage(DeltaR+1:DeltaR+hm, 1: wm) = MidImageFilteredNoPadding;
                SurfFilteredImage(1:hs, DeltaC+1:DeltaC+ws) = SurfImageFilteredNoPadding;
            else
                MidMaxImage(DeltaR + 1:DeltaR + hm, DeltaC + 1:DeltaC + wm) = MidImageNoPadding;
                SurfMaxImage(1:hs, 1:ws) = SurfImageNoPadding;
                MidFilteredImage(DeltaR + 1:DeltaR + hm, DeltaC + 1:DeltaC + wm) = MidFilteredImageNoPadding;
                SurfFilteredImage(1:hs, 1:ws) = SurfFilteredImageNoPadding;
            end
        end



        %=MidImage = imresize(MidImage,length(SurfImage)/length(MidImage));
    else
        MidMaxImage = MidImageNoPadding;
        SurfMaxImage = SurfImageNoPadding;
        MidFilteredImage = MidImageFilteredNoPadding;
        SurfFilteredImage = SurfImageFilteredNoPadding;
    end
    % 
    


end
%% 

%% Save stitched data in tif files

outputFolder = [DropboxFolder,filesep,Prefix,filesep,'FullEmbryoStitching'];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
APDetectionFolder = [DropboxFolder,filesep,Prefix,filesep,'APDetection'];
if ~exist(APDetectionFolder, 'dir')
    mkdir([DropboxFolder,filesep,Prefix,filesep,'APDetection']);
end

MidMaxProjFileName = 'MidTileStitch_MaxPadded.tif';
MidGaussFileName = 'MidTileStitch_MaxFilteredPadded.tif';
% MidZStacksFileName = 'MidTileStitchPadded.tif';
SurfMaxProjFileName = 'SurfTileStitch_MaxPadded.tif';
SurfGaussFileName = 'SurfTileStitch_MaxFilteredPadded.tif';
% SurfZStacksFileName = 'SurfTileStitchPadded.tif';
scale_factor = 6*10^4/max(max(MidMaxImage));
t1 =  cat(3, MidMaxImage, zeros(size(MidMaxImage), 'uint16'), MidMaxImage);
t1 = t1*scale_factor;
scale_factor_2 = 6*10^4/max(max(SurfMaxImage));
t2 =  cat(3, SurfMaxImage, SurfMaxImage, zeros(size(SurfMaxImage), 'uint16'));
t2 = t2*scale_factor_2;
testoverlay = zeros([size(SurfMaxImage) 3], 'uint16');
testoverlay = testoverlay + t1;
testoverlay = testoverlay + t2;
figure(1)
image(testoverlay)

    

imwrite(MidMaxImage,[outputFolder, filesep,MidMaxProjFileName],'compression','none');
imwrite(MidFilteredImage,[outputFolder, filesep,MidGaussFileName],'compression','none');
% imwrite(MidImage,[outputFolder, filesep,MidZStacksFileName],'compression','none');
imwrite(SurfMaxImage,[outputFolder, filesep,SurfMaxProjFileName],'compression','none'); 
imwrite(SurfFilteredImage,[outputFolder, filesep,SurfGaussFileName],'compression','none');
% imwrite(SurfImage,[outputFolder, filesep,SurfZStacksFileName],'compression','none');

if PostFilePresent
    figure(2)
    imagesc(MidEmbryoMaxImageV2)
    MidMaxProjFileNameV2 = 'MidTileStitch_MaxPaddedV2.tif';
    imwrite(MidEmbryoMaxImageV2,[outputFolder, filesep,MidMaxProjFileNameV2],'compression','none');
end
