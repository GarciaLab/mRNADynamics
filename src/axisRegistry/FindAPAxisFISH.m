function FindAPAxisFISH(Prefix, varargin)

%This function allows to check and manually determine the AP axis


%Manual:

%c: clear all AP information
%a: Select anterior end
%p: Select posterior end
%.: Increase contrast
%,: Decrease contrast
%r: Reset the contrast
%s: Swap AP positions
%right click: Delete the point that was clicked on
%m: Manual stitching mode
%x: Exit and save


close all

%Load the folder information

[RawDataPath,ProcPath,DropboxFolder,MS2CodePath, PreProcPath]=...
    DetermineLocalFolders(Prefix);




saveVars = {};
dv = false;
% 
% %Load the current AP detection data
% load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])
% 
% %Load the full embryo image
% APImage=readTiffStack([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo.tif']);

optionalResults = '';

[rawDataPath,ProcPath,DropboxFolder,MS2CodePath, PreProcPath,...
    rawDataFolder, Prefix, ExperimentType,Channel1,Channel2,OutputFolder,...
    Channel3, spotChannels, MovieDataBaseFolder, movieDatabase]...
    = readMovieDatabase(Prefix, optionalResults);

DataFolder = [DropboxFolder, filesep, Prefix];
FilePrefix = [Prefix, '_'];

[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution, ...
    Channel1, Channel2, Objective, Power, DataFolderColumnValue, ~, Comments, ...
    nc9, nc10, nc11, nc12, nc13, nc14, CF, Channel3, prophase, metaphase] =...
    getExperimentDataFromMovieDatabase(Prefix, movieDatabase);

load([DataFolder, filesep, 'FrameInfo.mat'], 'FrameInfo')

[xSize, ySize, pixelSize, zStep, snippet_size, LinesPerFrame, PixelsPerLine,...
    nEmbryos] = getFrameInfoParams(FrameInfo);

coordA = {};
coordP = {};
coordD = {};
coordV = {};


Channels={Channel1{1},Channel2{1}, Channel3{1}};

InputChannelTemp1 = strfind({lower(Channel1{1}),lower(Channel2{1}), lower(Channel3{1})},'input');
InputChannelTemp2=~cellfun(@isempty,InputChannelTemp1);
InputChannel = find(InputChannelTemp2);
   nameSuffix=['_ch',iIndex(InputChannel,2)];
           NumberSlices=FrameInfo(1).NumberSlices;
        NumberSlices2 = NumberSlices + 2;
inputImage=zeros(LinesPerFrame,PixelsPerLine,NumberSlices2);

APImageStack = zeros(LinesPerFrame,PixelsPerLine,nEmbryos);
%Load the z-stack for this frame

for e = 1:nEmbryos
    parfor z=1:NumberSlices2   
        inputImage(:,:,z)=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(e,3),'_z',iIndex(z,2),nameSuffix,'.tif']);
    end
    APImageStack(:,:,e) = max(inputImage, [], 3);
end

    
% for f = 1:nEmbryos
%     APImageStack(:,:,f)=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'-','His_',iIndex(f,3),'.tif']);
% end
% APImageStack = inputStack;

currentEmbryo = 1;

APImage = APImageStack(:,:,currentEmbryo);

APImageFig=figure;
apAx = axes(APImageFig);
%Range for image display
DisplayRange=[min(min(APImage)),max(max(APImage))/2];


%Now, do the correction
cc=1;

while (cc~='x')
    
    APImage = APImageStack(:,:,currentEmbryo);
    APImage = mat2gray(APImage,double(DisplayRange));
    imshow(APImage, 'Parent', apAx)
    %imshow(APImage,DisplayRange)
    
    axis image
    axis off
    title(['Anterior (green), posterior (red); original. Embryo ', num2str(currentEmbryo)])
    hold on
    
    try
        plot(coordA{currentEmbryo}(1),coordA{currentEmbryo}(2),'g.','MarkerSize',20);
    end
    try
        plot(coordP{currentEmbryo}(1),coordP{currentEmbryo}(2),'r.','MarkerSize',20);
    end
    
    try
        plot(coordV{currentEmbryo}(1),coordV{currentEmbryo}(2),'m.','MarkerSize',20);
    end
    try
        plot(coordD{currentEmbryo}(1),coordD{currentEmbryo}(2),'y.','MarkerSize',20);
    end
    
    hold off
    
    figure(APImageFig)
    
    ct=waitforbuttonpress;
    cc=get(APImageFig,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    if (ct~=0)&(cc=='c')        %Clear all AP information
        coordA{currentEmbryo}=[];
        coordP{currentEmbryo}=[];
    elseif (ct~=0)&(cc=='.')&currentEmbryo <nEmbryos	%move to next embryo
        currentEmbryo = currentEmbryo+1;
     elseif (ct~=0)&(cc==',')&currentEmbryo >1%%move to previous embryo
         currentEmbryo = currentEmbryo - 1;
     elseif (ct~=0)&(cc=='k') %%move to previous embryo
           currentEmbryo = inputdlg('Embryo to jump to:', ...
                'Move to embryo');
            currentEmbryo = str2double(currentEmbryo{1});
    elseif (ct~=0)&(cc=='a')	%Select anterior end
        [coordAx,coordAy]=ginputc(1,'Color',[1,1,1]);
        coordA{currentEmbryo} = [coordAx,coordAy];
    elseif (ct~=0)&(cc=='p')    %Select posterior end
        [coordPx,coordPy]=ginputc(1,'Color',[1,1,1]);
        coordP{currentEmbryo} = [coordPx,coordPy];
    elseif (ct~=0)&(cc=='v')
        [coordVx,coordVy]=ginputc(1,'Color','m');
        coordV{currentEmbryo} = [coordVx,coordVy];
        saveVars = [saveVars, 'coordV'];
        dv = true;
    elseif (ct~=0)&(cc=='d')
        [coordDx,coordDy]=ginputc(1,'Color','y');
        coordD{currentEmbryo} = [coordDx,coordDy];
        saveVars = [saveVars, 'coordD'];
        dv = true;
    elseif (ct~=0)&(cc=='g')    %Increase contrast
        DisplayRange(2)=DisplayRange(2)/2;
        
    elseif (ct~=0)&(cc=='b')    %Decrease contrast
        DisplayRange(2)=DisplayRange(2)*2;
        
    elseif (ct~=0)&(cc=='r')    %Reset the contrast
        DisplayRange=[min(min(APImage)),max(max(APImage))];
   elseif (ct~=0)&(cc=='0')    
        keyboard;
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
% 
% diagFigure = figure;
% imagesc(APImage)
% axis image
% axis off
% title('Anterior (green), posterior (red); corrected')
% hold on
% plot(coordA(1),coordA(2),'g.','MarkerSize',20);
% plot(coordP(1),coordP(2),'r.','MarkerSize',20);
% if dv
%     plot(coordV(1),coordV(2),'m.','MarkerSize',20);
%     plot(coordD(1),coordD(2),'y.','MarkerSize',20);
% end
% hold off
% saveas(gcf, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'APEmbryo-Manual.tif']);
% 
% close(diagFigure);


end