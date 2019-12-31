function FindAPAxisTile(Prefix, varargin)
% author: Gabriella Martini
% date created: 12/28/19
% date last modified: 12/30/19


%% Parse Inputs
if ~exist('Prefix')
    FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
    Dashes=strfind(FolderTemp,filesep);
    Prefix=FolderTemp((Dashes(end)+1):end);
end 


% Parse inputs
x = 1;
while x <= length(varargin)
    switch varargin{x}
        case{'NIterations'}
            NIterations = varargin{x+1};
            x = x + 1;
            fprintf('Number of Iterations: %d\n', NIterations);
        case {'FullyAutomate'}
            FullyAutomate = true;
            fprintf('Stitching fully automated.\n')
        case {'StitchManually'}
            StitchManually = true;
            fprintf('Stitching to be performed manually.\n')
        case {'MaxStep'}
            MaxStep = varargin{x+1};
            x = x+1;
            fprintf('Max Step Size to be used in stitching loop: %d\n', MaxStep)
        case{'MaxOverlap'}
            MaxOverlap = varargin{x+1};
            x = x+1;
            fprintf('Max overlap between adjacent tiles to be used in stitching loop: %d\n', RowMaxOverlap)
    end
    x = x +1;
end





%% 


[SourcePath, ~, DefaultDropboxFolder, DropboxFolder, ~, ~,...
~, ~] = DetermineAllLocalFolders(Prefix);

stitchingDataFolder = [DropboxFolder,filesep,Prefix,filesep,'FullEmbryoStitching'];
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

%Datatype is hardcoded in, unlike in FindAPAxisFullEmbryo
DLIF=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*.lif']);
FileMode='LIFExport';

DSTITCH=dir([DropboxFolder,filesep,Prefix,filesep,'StitchedEmbryoImages\*.tif']);

%% 

% Identify the midsagittal image
MidLifFileIndex=find(~cellfun('isempty',strfind(lower({DLIF.name}),'midtile')));
SurfLifFileIndex=find(~cellfun('isempty',strfind(lower({DLIF.name}),'surftile')));

%% 



%Rotates the full embryo image to match the rotation of the zoomed
%time series
zoom_angle = 0;
full_embryo_angle = 0;

LIFMid=bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,DLIF(MidLifFileIndex).name]);
LIFSurf=bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,DLIF(SurfLifFileIndex).name]);

MIDMeta = LIFMid{:, 4};
SURFMeta = LIFSurf{:,4};
%% 
[SurfNSeries, SurfNFrames, SurfNSlices, SurfNPlanes, SurfNChannels, SurfFrame_Times] = getFrames(SURFMeta);
SurfNTiles = SurfNSeries;
[MidNSeries, MidNFrames, MidNSlices, MidNPlanes, MidNChannels, MidFrame_Times] = getFrames(MIDMeta);
MidNTiles = MidNSeries;


PixelSize = double(MIDMeta.getPixelsPhysicalSizeX(1).value);% units: microns
PixelSize_m = double(PixelSize)*10^(-6);
%% 


MidFilteredImage = imread([stitchingDataFolder, filesep, 'MidTileStitch_MaxFilteredPadded.tif']);
MidImageStack = imread([stitchingDataFolder, filesep, 'MidTileStitchPadded.tif']);
MidImage =imread([stitchingDataFolder, filesep, 'MidTileStitch_MaxPadded.tif']);
SurfFilteredImage = imread([stitchingDataFolder, filesep, 'SurfTileStitch_MaxFilteredPadded.tif']);
SurfImageStack = imread([stitchingDataFolder, filesep, 'SurfTileStitchPadded.tif']);
SurfImage = imread([stitchingDataFolder, filesep, 'SurfTileStitch_MaxPadded.tif']);
%% 
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
        filesep, 'MetaData', filesep,'*SurfTile*.xml']);
    xml_file2 = xml_file_path2(1).name;
    evalin('base','clear rot')
    xDoc2 = searchXML([SourcePath, filesep, Date, filesep, EmbryoName, filesep,'FullEmbryo', filesep,...
            'MetaData', filesep, xml_file2]);
     full_embryo_angle = str2double(evalin('base','rot'));
else 
    warning('No full embryo metadata found.')
end

MidImage = imrotate(MidImage, -zoom_angle + full_embryo_angle);
SurfImage = imrotate(SurfImage, -zoom_angle + full_embryo_angle);
MidImageStack = imrotate(MidImageStack, -zoom_angle + full_embryo_angle);
SurfImageStack = imrotate(SurfImageStack, -zoom_angle + full_embryo_angle);
MidFilteredImage = imrotate(MidFilteredImage, -zoom_angle + full_embryo_angle);
SurfFilteredImage = imrotate(SurfFilteredImage, -zoom_angle + full_embryo_angle);

%Save it to the Dropbox folder
if ~exist([DropboxFolder,filesep,Prefix,filesep,'APDetection'], 'dir')
    mkdir([DropboxFolder,filesep,Prefix,filesep,'APDetection'])
end
imwrite(uint16(MidImage),[DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo.tif'],'compression','none');
imwrite(uint16(SurfImage),[DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryoSurf.tif'],'compression','none');

imwrite(uint16(MidImageStack),[DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryoStack.tif'],'compression','none');
imwrite(uint16(SurfImageStack),[DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryoSurfStack.tif'],'compression','none')
    


%%


APImage = MidFilteredImage;

APImageFig=figure;
apAx = axes(APImageFig);
%Range for image display
DisplayRange=[min(min(APImage)),max(max(APImage))];

coordA=[1,1];
coordP=size(MidImage);
%Now, do the correction
cc=1;

while (cc~='x')
    
    
    imshow(imadjust(APImage), 'Parent', apAx)
    %imshow(APImage,DisplayRange)
    
    axis image
    axis off
    title('Anterior (green), posterior (red); original')
    hold on
    
    try
        plot(coordA(1),coordA(2),'g.','MarkerSize',20);
    catch
        %not sure what happened here.
    end
    
    try
        plot(coordP(1),coordP(2),'r.','MarkerSize',20);
    catch
        %not sure what happened here.
    end
    
    if exist('coordV', 'var')
        plot(coordV(1),coordV(2),'m.','MarkerSize',20);
    end
    if exist('coordD', 'var')
        plot(coordD(1),coordD(2),'y.','MarkerSize',20);  
    end
    
    hold off
    
    figure(APImageFig)
    ct=waitforbuttonpress;
    cc=get(APImageFig,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    
    if (ct~=0)&(cc=='c')        %Clear all AP information
        coordA=[];
        coordP=[];
    elseif (ct~=0)&(cc=='a')	%Select anterior end
        [coordAx,CoordAy]=ginputc(1,'Color',[1,1,1]);
        coordA = [coordAx,CoordAy];
    elseif (ct~=0)&(cc=='p')    %Select posterior end
        [coordPx,CoordPy]=ginputc(1,'Color',[1,1,1]);
        coordP = [coordPx,CoordPy];
    elseif (ct~=0)&(cc=='.')    %Increase contrast
        DisplayRange(2)=DisplayRange(2)/2;
        
    elseif (ct~=0)&(cc==',')    %Decrease contrast
        DisplayRange(2)=DisplayRange(2)*2;
        
    elseif (ct~=0)&(cc=='r')    %Reset the contrast
        DisplayRange=[min(min(APImage)),max(max(APImage))];
        
    elseif (ct==0)&(strcmp(get(APImageFig,'SelectionType'),'alt')) %Delete the point that was clicked on
        cc=1;
        
        [~,MinIndex]=min((cm(1,1)-[coordA(1),coordP(1)]).^2+(cm(1,2)-[coordA(2),coordP(2)]).^2);
        
        if MinIndex==1
            coordA=[];
        elseif MinIndex==2
            coordP=[];
        end
    elseif (ct~=0)&(cc=='m')        %Manual stitching mode
        %ManualStitch
    elseif (ct~=0)&(cc=='s')
        coordPTemp=coordA;
        coordA=coordP;
        coordP=coordPTemp;
    elseif (ct~=0)&(cc=='v')
        [coordVx,CoordVy]=ginputc(1,'Color','m');
        coordV = [coordVx,CoordVy];
        saveVars = [saveVars, 'coordV'];
        dv = true;
    elseif (ct~=0)&(cc=='d')
        [coordDx,CoordDy]=ginputc(1,'Color','y');
        coordD = [coordDx,CoordDy];
        saveVars = [saveVars, 'coordD'];
        dv = true;
        
        
    end
end

%Save the information
saveVars = {};
saveVars = [saveVars, 'coordA', 'coordP'];
if exist('xShift', 'var')
    saveVars = [saveVars, 'xShift'];
end
if exist('xShift1', 'var')
    saveVars = [saveVars, 'xShift1'];
end

save([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'],saveVars{:});
%Redo the diagnostic plots

diagFigure = figure;
imagesc(APImage)
axis image
axis off
title('Anterior (green), posterior (red); corrected')
hold on
plot(coordA(1),coordA(2),'g.','MarkerSize',20);
plot(coordP(1),coordP(2),'r.','MarkerSize',20);
% if dv
%     plot(coordV(1),coordV(2),'m.','MarkerSize',20);
%     plot(coordD(1),coordD(2),'y.','MarkerSize',20);
% end
hold off
saveas(gcf, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'APEmbryo-Manual.tif']);

close all;

end
