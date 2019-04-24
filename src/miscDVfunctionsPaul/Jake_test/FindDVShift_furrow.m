function  [DV_shift] = FindDVShift_furrow(varargin)
%% Initialization


%% Part 1: Read image data

Prefix=varargin{1};
[SourcePath, ~, DefaultDropboxFolder, DropboxFolder, ~, ~,...
~, ~] = DetermineAllLocalFolders(Prefix);

%Figure out what type of experiment we have. Note: we name the var "DateFromDateColumn" to avoid shadowing previously defined "Date" var.
[DateFromDateColumn, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder)

%Find out the date it was taken
Dashes=findstr(Prefix,'-');
Date=Prefix(1:Dashes(3)-1);
EmbryoName=Prefix(Dashes(3)+1:end);

%Find Zoom Factor
    FileMode='LIFExport';
    %Find the zoomed movie pixel size
    D=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'*.',FileMode(1:3)]);

    %Load only the metadata from the zoomed images
    MetaReader=bfGetReader([SourcePath,filesep,Date,filesep,EmbryoName,filesep,D(end).name]);
    MetaZoom=MetaReader.getMetadataStore();
    try
        PixelSizeZoom=str2double(MetaZoom.getPixelsPhysicalSizeX(0).value);
    catch
        PixelSizeZoom=str2double(MetaZoom.getPixelsPhysicalSizeX(0));
    end

    %Find the full embryo pixel size and load the image
    D=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*surf*.',FileMode(1:3)]);

    ImageTemp=bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,D(end).name]);
    MetaFullEmbryo= ImageTemp{:, 4};
    PixelSizeFullEmbryo=str2double(MetaFullEmbryo.getPixelsPhysicalSizeX(0) );
    try
        PixelSizeFullEmbryo=str2double(MetaFullEmbryo.getPixelsPhysicalSizeX(0).value);
    catch
        PixelSizeFullEmbryo=str2double(MetaFullEmbryo.getPixelsPhysicalSizeX(0));
    end

    %Zoom factor
    MovieZoom=PixelSizeFullEmbryo(1)/PixelSizeZoom(1);
    SurfZoom=1;     %We'll call the zoom of the full embryo image 1


% Read full embryo image (surf, mid)
FullEmbryo=imread([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo.tif']);
%FullEmbryoSurf=imread([DropboxFolder,filesep,Prefix,filesep,'DV',filesep,'surf_max.tif']);
FullEmbryoSurf = imread('embryo2_surf_max.lif_furrow_ch01.tif');

%% Part 4: Calculate DV position

ImMid = FullEmbryo;
ImSurf = FullEmbryoSurf;

%{
%Let's define vector for our clicks. The first will be our x-values.
v_APx=[v_AP(1,1)];
%The second vector will be for our y-values
v_APy=[v_AP(1,2)];
%}

imshow(ImSurf,[]);
hold on
furrow=ginput(1); %Click once on the tip of the ventral furrow

%Get the information about the AP axis as well as the image shifts
%used for the stitching of the two halves of the embryo
load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])
APx = [coordA(1),coordP(1)];
APy = [coordA(2),coordP(2)];

%Find the equation of the AP axis
mAP=(APy(1)-APy(2))/(APx(1)-APx(2));
bAP=APy(1)-mAP*APx(1);

%Find the equation of ventral furrow
mAP_v = mAP;
bAP_v = furrow(2)-mAP_v*furrow(1);

APx_v = APx;
APy_v = [mAP_v*APx_v(1)+bAP_v,mAP_v*APx_v(2)+bAP_v];

plot(APx,APy,'o-r','LineWidth',4);
plot(APx_v,APy_v,'o-r','LineWidth',4);

%DV_shift = (bAP_v-bAP);

end

