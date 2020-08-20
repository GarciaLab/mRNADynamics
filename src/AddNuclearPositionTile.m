function AddNuclearPositionTile(Prefix)
% author: Gabriella Martini
% date created: 1/29/19
% date last modified: 12/31/19
%Default set of variables to save

SelectChannel=0;
SkipAlignment=false;
ManualAlignment=false;


saveVars={'coordA','coordP','coordAZoom','coordPZoom'};
%Get the relevant folders for this data set
[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath,...
configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);
%Get the relevant folders for this data set
[RawDynamicsPath, ~, DefaultDropboxFolder, DropboxFolder, ~, PreProcPath,...
    configValues,...
    movieDatabasePath, movieDatabase] = DetermineAllLocalFolders(Prefix);
%Find out the date it was taken
stitchingDataFolder = [DropboxFolder,filesep,Prefix,filesep,'FullEmbryoStitching'];

[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoopEnd, APResolution,...
    Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
    nc9, nc10, nc11, nc12, nc13, nc14, CF,Channel3,prophase,metaphase, anaphase] = getExperimentDataFromMovieDatabase(Prefix, movieDatabase);

Dashes=findstr(Prefix,'-');
Date=Prefix(1:Dashes(3)-1);
EmbryoName=Prefix(Dashes(3)+1:end);


%% 


%Load the AP detection information generated by AddParticlePosition. This
%means that all alignments and manual corrections have to be done at that
%level.
load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])

%Get the surface image in the zoomed case by looking at the last
%frame of our movie
%See if we had any lineage/nuclear information
hisDir=dir([PreProcPath,filesep,Prefix,filesep,'*-His_*']);
if ~isempty(hisDir)
    histoneChannelPresent = true;
else
    histoneChannelPresent = false;
end
ZoomImage=imread([PreProcPath,filesep,Prefix,filesep,hisDir(end-1).name]);
[FileMode, EmbryoName, projectDate] = getMicroscope(Prefix);
rawPrefixPath = [RawDynamicsPath,filesep,projectDate,filesep,EmbryoName,filesep];
fullEmbryoPath = [rawPrefixPath, 'FullEmbryo', filesep];
%% 
ChannelToLoadTemp=contains([Channel1,Channel2,Channel3],'nuclear','IgnoreCase',true);
if sum(ChannelToLoadTemp) && sum(ChannelToLoadTemp)==1
    ChannelToLoad=find(ChannelToLoadTemp);
elseif sum(ChannelToLoadTemp) && length(ChannelToLoadTemp)>=2
    ChannelToLoad=find(ChannelToLoadTemp);
    ChannelToLoad = ChannelToLoad(1);
else
    error('No histone channel found. Was it defined in MovieDatabase as :Nuclear or :InvertedNuclear?')
end
%% %Datatype is hardcoded in, unlike in FindAPAxisFullEmbryo
DLIF=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*.lif']);
FileMode='LIFExport';

DSTITCH=dir([DropboxFolder,filesep,Prefix,filesep,'StitchedEmbryoImages\*.tif']);

MidLifFileIndex=find(~cellfun('isempty',strfind(lower({DLIF.name}),'midtile')));
SurfLifFileIndex=find(~cellfun('isempty',strfind(lower({DLIF.name}),'surftile')));


LIFMid=bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,DLIF(MidLifFileIndex).name]);
LIFSurf=bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,DLIF(SurfLifFileIndex).name]);

MIDMeta = LIFMid{:, 4};
SURFMeta = LIFSurf{:,4};

%% 

MidImage=imread([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo.tif']);
SurfImage=imread([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryoSurf.tif']);

MidImageStack=imread([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryoStack.tif']);
SurfImageStack=imread([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryoSurfStack.tif']);
%% 

%Rotates the full embryo image to match the rotation of the zoomed
%time series
zoom_angle = 0;
full_embryo_angle = 0;

[SurfNSeries, SurfNFrames, SurfNSlices, SurfNPlanes, SurfNChannels, SurfFrame_Times] = getFrames(SURFMeta);
SurfNTiles = SurfNSeries;
[MidNSeries, MidNFrames, MidNSlices, MidNPlanes, MidNChannels, MidFrame_Times] = getFrames(MIDMeta);
MidNTiles = MidNSeries;


PixelSize = double(MIDMeta.getPixelsPhysicalSizeX(1).value);% units: microns


%% 
SurfName=[];

%Figure out the different channels
%NuclearChannel=contains([Channel1,Channel2,Channel3],'nuclear','IgnoreCase',true);


if any(contains([Channel1,Channel2,Channel3],'nuclear','IgnoreCase',true))
    if SelectChannel
        HisChannel = ChannelToLoad;
        InvertHis=0;
    elseif ~any(contains([Channel1,Channel2,Channel3],'inverted','IgnoreCase',true))
        HisChannel = ChannelToLoad;
        InvertHis=0;
    elseif any(contains([Channel1,Channel2,Channel3],'inverted','IgnoreCase',true))
        HisChannel= ChannelToLoad;
        InvertHis=1;
    else
        error('You should check the MovieDatabase.csv to see whether ":Nuclear" or "invertedNuclear" is plugged into your channels.')
    end
end

%% 

%Find the zoomed movie pixel size
DZOOM=dir([rawPrefixPath,'*.',FileMode(1:3)]);

%Load only the metadata from the zoomed images
MetaReader=bfGetReader([rawPrefixPath,DZOOM(end).name]);
MetaZoom=MetaReader.getMetadataStore();

PixelSizeZoom=str2double(MetaZoom.getPixelsPhysicalSizeX(0).value);
PixelSizeFullEmbryoSurf = PixelSize;

NChannelsMeta=SURFMeta.getChannelCount(0);
NSlices=str2double(SURFMeta.getPixelsSizeZ(0));

% if InvertHis
%     SurfImage=SurfImage(:,:,HisChannel+round(NSlices/2)*NChannelsMeta);
% else
%     SurfImage=max(SurfImage,[],3);
% end
%% 

%SurfImage = imrotate(SurfImage, -zoom_angle + full_embryo_angle);
% 
% mkdir([DropboxFolder,filesep,Prefix, filesep, 'DV']);
% maxSurfSavePath = [DropboxFolder,filesep,Prefix, filesep, 'DV', filesep, 'surf_max.tif'];
% imwrite(SurfImage,maxSurfSavePath);


%clear ImageTemp

%Zoom factor
MovieZoom=PixelSizeFullEmbryoSurf(1)/PixelSizeZoom(1);
SurfZoom=1;     %We'll call the zoom of the full embryo image 1

%Get the size of the zoom image
Columns = str2double(MetaZoom.getPixelsSizeX(0));
Rows = str2double(MetaZoom.getPixelsSizeY(0));
surf_size = size(SurfImage);
SurfColumns= surf_size(2);
SurfRows=surf_size(1);
ZoomRatio = MovieZoom / SurfZoom;

%% 
%Load the schnitzcells
load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])

if isfield(schnitzcells,'APpos')
    warning([Prefix,'_lin.mat already has AP positions stored.'])
end
if isfield(schnitzcells,'DVpos')
    warning([Prefix,'_lin.mat already has DV positions stored.'])
end

DHis=dir([PreProcPath,filesep,Prefix,filesep,Prefix,'-His*.tif']);
if ~isempty(DHis)
    ZoomImage=imread([PreProcPath,filesep,Prefix,filesep,DHis(end-1).name]);
else
    disp('Did you run ExportDataForLivemRNA again, after editing the MovieDatabase.csv with ":Nuclear" ("invertedNuclear")?')
    % This might be the case, for instance, if you're just trying
    % to find AP information about an image without using FISH
    % code. In that case, just extract the nuclei from the last
    % raw image.
    DGFP = dir([rawPrefixPath, '*.tif']);
    ImageInfo = imfinfo([rawPrefixPath, DGFP(end).name]);
    NumFramesAndSlices = length(ImageInfo)/2;
    RawImage3M = NaN(Rows, Columns, NumFramesAndSlices);
    for lImageIndex = 1:NumFramesAndSlices
        RawImage3M(:, :, lImageIndex) = imread([rawPrefixPath,DGFP(end).name],'Index', 2*(lImageIndex-1) + ChannelToLoad);
    end
    ZoomImage = max(RawImage3M, [], 3) / 256;
end

%% 
%Do a correlation between the zoomed in and zoomed out surface images
%to figure out the shift.


if ~SkipAlignment && histoneChannelPresent

    ShiftColumn=0;
    ShiftRow=0;
    NucMaskZoomIn = false(size(ZoomImage));
    NucMaskZoomOut = false(size(SurfImage));

    if ZoomRatio < 24  %ZoomRatio > 1 && ZoomRatio < 24. AR 12/4/17- where did this number come from

        %Enlarge the zoomed out image so we can do the cross-correlation
        SurfImageResized=imresize(SurfImage, ZoomRatio);

        %Calculate the correlation matrix and find the maximum
        im1 = ZoomImage;
        im2 = SurfImageResized;

%             if InvertHis
%                 im1 = imcomplement(im1);
%                 im2 = imcomplement(im2);
%             end

        %             try
        %                 C = gather(normxcorr2(gpuArray(im1), gpuArray(im2)));
        %             catch
        C = normxcorr2(im1, im2);
        %             end

        [Max2,MaxRows]=max(C);
        [~,MaxColumn]=max(Max2);
        MaxRow=MaxRows(MaxColumn);
        [CRows,CColumns]=size(C);


        %This shift is now converted to the zoom out distances. If we
        %want to translate to the zoomed in coordinates we need to
        %multiply again by ZoomRatio.
        ShiftRow=round((MaxRow-(CRows/2+1))/ZoomRatio);
        ShiftColumn=round((MaxColumn-(CColumns/2+1))/ZoomRatio);

        %How well did we do with the alignment?
        [RowsResized,ColumnsResized]=size(SurfImageResized);
        [RowsZoom,ColumnsZoom]=size(ZoomImage);

        %If manual alignment was done before then load the results
        if exist('ManualAlignmentDone', 'var')
            if ManualAlignmentDone
                disp('Manual alignment results saved.')
                if ~yToManualAlignmentPrompt
                    Answer=input('Would you like to use them (y/n)? ','s');
                else
                    Answer = 'y';
                end
                if strcmpi(Answer,'y')
                    load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'],'ShiftRow','ShiftColumn')
                elseif strcmpi(Answer,'n')
                    disp('Deleting manual alignment results')
                    clear ManualAlignmentDone
                else
                    msgbox('Error: Answer not recognized');
                    error('Answer not recognized')
                end
            end
        end

        %If we can, we'll crop the surface image to center the overlay
        %on the area of acquisition. If not, we'll just show both images
        %without cropping.

        if ((round(RowsResized/2-RowsZoom/2+ShiftRow*ZoomRatio))>0)&...
                (round(RowsResized/2+RowsZoom/2-1+ShiftRow*ZoomRatio)<=RowsResized)&...
                (round(ColumnsResized/2-ColumnsZoom/2+ShiftColumn*ZoomRatio)>0)&...
                (round(ColumnsResized/2+ColumnsZoom/2-1+ShiftColumn*ZoomRatio)<=ColumnsResized)
            RowsResizedRange=...
                round(RowsResized/2-RowsZoom/2+ShiftRow*ZoomRatio):...
                round(RowsResized/2+RowsZoom/2-1+ShiftRow*ZoomRatio);
            ColumnsResizedRange=...
                round(ColumnsResized/2-ColumnsZoom/2+ShiftColumn*ZoomRatio):...
                round(ColumnsResized/2+ColumnsZoom/2-1+ShiftColumn*ZoomRatio);
        else
            RowsResizedRange=1:RowsResized;
            ColumnsResizedRange=1:ColumnsResized;
        end


        %Make an overlay of the zoomed in and zoomed out real
        %images as well as of a quickly segmented nuclear mask. If this
        %fails, we'll swtich to ManualAlignment.

        try
            %Real image overlay
            %Crop the zoomed out image to match the zoomed in one
            SurfImageResizeZoom=...
                SurfImageResized(RowsResizedRange,ColumnsResizedRange);
            ImOverlay=cat(3,mat2gray(SurfImageResizeZoom),...
                +mat2gray(ZoomImage),zeros(size(SurfImageResizeZoom)));

            %Nuclear mask overlay
            NucMaskZoomOut= GetNuclearMask(SurfImage,2.5,0);
            NucMaskZoomOutResized=imresize(NucMaskZoomOut, ZoomRatio);
            NucMaskZoomOutResizedCropped=...
                NucMaskZoomOutResized(RowsResizedRange,ColumnsResizedRange);
            NucMaskZoomIn=GetNuclearMask(ZoomImage,8,2);
            ImOverlayMask=cat(3,mat2gray(NucMaskZoomOutResizedCropped),...
                +mat2gray(NucMaskZoomIn),zeros(size(NucMaskZoomOutResizedCropped)));

            alOvFig = figure(1);
            imOv = subplot(2,1,1);
            imshow(ImOverlay, 'Parent', imOv)
            imOvMask = subplot(2,1,2);
            imshow(ImOverlayMask, 'Parent', imOvMask)
            saveas(alOvFig, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'AlignmentOverlay.tif']);




            %Show the correlation image, but crop it a little bit
            contFig = figure(2);
            contAxes = axes(contFig);
            %HG: Note that I changed the ranges here by two pixels at
            %least.
            lX = uint16((CRows+1)/2-RowsZoom+1);
            uX = uint16((CRows-1)/2+RowsZoom);
            lY = uint16((CColumns+1)/2-ColumnsZoom+1);
            uY = uint16((CColumns-1)/2+ColumnsZoom);

            contourf(contAxes,imresize(abs(C(lX:uX,lY:uY)),.5));
            saveas(contFig, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'AlignmentCorrelation.tif']);
        catch
            warning('Could not generate correlation image. Switching to manual alignment')
            ManualAlignment=true;
        end

    else
        warning('Not able to do the cross correlation. Switching to manual alignment mode.')     
        ManualAlignment=true;
    end


    %See if we need the manual alignment
    if ManualAlignment
        %See if we need to load the previous manual alignment results
        [ShiftColumn,ShiftRow]=ManualAPCorrection(SurfImage,ZoomImage,C,ZoomRatio,ShiftRow,ShiftColumn,...
            FullEmbryo, ZoomRatio, SurfRows,Rows, Columns, coordA, coordP, SurfColumns);
        ManualAlignmentDone = true;
    end

else
    warning('Not able to do the cross correlation. Assuming no shift between surface-level and movie-level images.')
end
%% 

ImageCenter=[SurfRows/2,SurfColumns/2];
    
%This is for the acquisition/zoom image
TopLeft=[ImageCenter(1)-Rows/ZoomRatio/2+ShiftRow,...
    ImageCenter(2)-Columns/ZoomRatio/2+ShiftColumn];
BottomRight=[ImageCenter(1)+Rows/ZoomRatio/2+ShiftRow,...
    ImageCenter(2)+Columns/ZoomRatio/2+ShiftColumn];
coordAZoom=(coordA-[TopLeft(2),TopLeft(1)])*ZoomRatio;
coordPZoom=(coordP-[TopLeft(2),TopLeft(1)])*ZoomRatio;
if exist('coordD', 'var')
     coordDZoom=(coordD-[TopLeft(2),TopLeft(1)])*ZoomRatio;
    coordVZoom=(coordV-[TopLeft(2),TopLeft(1)])*ZoomRatio;
end

fullFigure = figure(7);
fullAxes = axes(fullFigure);
imshow(imadjust(mat2gray(SurfImage)),'DisplayRange',[],'InitialMagnification',100,'Parent', fullAxes)
%imshow(imadjust(mat2gray(SurfImageFlat)),'DisplayRange',[],'InitialMagnification',100,'Parent', fullAxes)
hold(fullAxes, 'on')
rectangle('Position',[TopLeft([2,1]),BottomRight([2,1])-TopLeft([2,1])],'EdgeColor','r')
plot(coordA(1),coordA(2),'.g','MarkerSize',30)
plot(coordP(1),coordP(2),'.r','MarkerSize',30)
plot([coordA(1),coordP(1)],[coordA(2),coordP(2)],'-b')
hold(fullAxes, 'off')
saveas(fullFigure, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryoArea.tif']);

zoomImageFigure = figure(9);
zoomImageAxes = axes(zoomImageFigure);
imshow(imadjust(ZoomImage),[], 'Parent', zoomImageAxes)
%imshow(NucMaskZoomIn)
%imshow(NucMaskZoomOutResizedCropped)
hold(zoomImageAxes,'on')
plot(zoomImageAxes,[coordAZoom(1),coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'-b')
plot(zoomImageAxes,[coordAZoom(1)+1,coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'--r')
plot(zoomImageAxes,[coordAZoom(1)-1,coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'-.g')
hold(zoomImageAxes,'off')
saveas(zoomImageFigure, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'ZoomedEmbryoAP.tif']);


%With AP coordinates in hand we can now determine the AP position of
%all particles. Look I my notes in "Calculating AP positions" in Notability
%for details of the calculation.

%Angle between the x-axis and the AP-axis
APAngle=atan2((coordPZoom(2)-coordAZoom(2)),(coordPZoom(1)-coordAZoom(1)));
APLength = norm(coordAZoom-coordPZoom);
    if exist('coordD', 'var')
        DVLength = norm(coordDZoom-coordVZoom);
    else
        DVLength = APLength/2; %rough estimate but surprisingly accurate
    end
saveVars = [saveVars, 'APLength', 'DVLength'];

APPosImage=zeros(size(ZoomImage));
[Rows,Columns]=size(ZoomImage);

for y=1:Rows
    for x=1:Columns
        Angle=atan2((y-coordAZoom(2)),(x-coordAZoom(1)));

        Distance=sqrt((coordAZoom(2)-y).^2+(coordAZoom(1)-x).^2);
        APPosition=Distance.*cos(Angle-APAngle);
        APPosImage(y,x)=APPosition/APLength;
    end
end


APbinID = 0:APResolution:1;


APPosBinImage=zeros(size(APPosImage));
for i=1:(length(APbinID)-1)
    FilteredMask=(APbinID(i)<=APPosImage)&(APbinID(i+1)>APPosImage);

    APPosBinImage=APPosBinImage+FilteredMask*i;
end


zoomOverlayFigure = figure(10);
zoomOverlayAxes = axes(zoomOverlayFigure);
ZoomOverlay=cat(3,mat2gray(ZoomImage)/2+mat2gray(APPosBinImage)/2,...
    mat2gray(ZoomImage)/2,mat2gray(ZoomImage)/2);
imshow(ZoomOverlay, 'Parent', zoomOverlayAxes)
hold(zoomOverlayAxes, 'on')
plot(zoomOverlayAxes,[coordAZoom(1),coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'-b')
plot(zoomOverlayAxes,[coordAZoom(1)+1,coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'--r')
plot(zoomOverlayAxes,[coordAZoom(1)-1,coordPZoom(1)],[coordAZoom(2),coordPZoom(2)],'-.g')
hold(zoomOverlayAxes,'off')

% DV_correction = 0;
% for i=1:length(Particles{ChN})
%     %Angle between the x-axis and the particle using the A position as a
%     %zero
%     Angles=atan2((Particles{ChN}(i).yPos-coordAZoom(2)),...
%         (Particles{ChN}(i).xPos-coordAZoom(1)));
% 
%     %Distance between the points and the A point
%     Distances=sqrt((coordAZoom(2)-Particles{ChN}(i).yPos).^2+(coordAZoom(1)-Particles{ChN}(i).xPos).^2);
%     APPositions=Distances.*cos(Angles-APAngle);
%     Particles{ChN}(i).APpos=APPositions/APLength;
% end
%% 

%% 


for i=1:length(schnitzcells)
    %Angle between the x-axis and the particle using the A position as a
    %zero
    Angles=atan2((double(schnitzcells(i).ceny)-coordAZoom(2)),...
        (double(schnitzcells(i).cenx)-coordAZoom(1)));

    %Distance between the points and the A point
    Distances=sqrt((coordAZoom(2)-double(schnitzcells(i).ceny)).^2+(coordAZoom(1)-double(schnitzcells(i).cenx)).^2);
    APPositions=Distances.*cos(Angles-APAngle);
    schnitzcells(i).APpos=APPositions/APLength;

    %Determine the distance perpendicular to the AP axis. This is a
    %proxy for a DV axis.

    schnitzcells(i).DVpos=Distances.*sin(Angles-APAngle);

end





save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells');

end 
    
    
    
