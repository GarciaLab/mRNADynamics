function CorrectAPAxis(varargin)

%This function allows you to check and manually determine the AP axis


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
    DetermineLocalFolders(varargin{1});


if ~isempty(varargin)
    Prefix=varargin{1};
else
    FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
    Dashes=strfind(FolderTemp,filesep);
    Prefix=FolderTemp((Dashes(end)+1):end);
end

saveVars = {};
dv = false;

%Load the current AP detection data
load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])

%Load the full embryo image
APImage=imread([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo.tif']);



APImageFig=figure;
apAx = axes(APImageFig);
%Range for image display
DisplayRange=[min(min(APImage)),max(max(APImage))];


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
if dv
    plot(coordV(1),coordV(2),'m.','MarkerSize',20);
    plot(coordD(1),coordD(2),'y.','MarkerSize',20);
end
hold off
saveas(gcf, [DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'APEmbryo-Manual.tif']);

close(diagFigure);

end

