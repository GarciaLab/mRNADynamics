function CorrectAPAxis(varargin)

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


close all

%Load the folder information

[Dummy,XLS]=xlsread('ComputerFolders.xlsx');

%Find out which computer this is. That will determine the folder structure.
[ret, name] = system('hostname');  
if ret ~= 0,  
   if ispc  
      name = getenv('COMPUTERNAME');  
   else  
      name = getenv('HOSTNAME');  
   end  
end  
name = lower(name); 


%Find which computer we are dealing with:
ComputerColumn=find(strcmp(XLS(1,:),name(1:end-1)));

%Now load the corresponding folders
SourceRow=find(strcmp(XLS(:,1),'SourcePath'));
FISHRow=find(strcmp(XLS(:,1),'FISHPath'));
DropboxRow=find(strcmp(XLS(:,1),'DropboxFolder'));
SchnitzRow=find(strcmp(XLS(:,1),'SchnitzcellsFolder'));



%Assign the folders
SourcePath=XLS{SourceRow,ComputerColumn};
FISHPath=XLS{FISHRow,ComputerColumn};
DropboxFolder=XLS{DropboxRow,ComputerColumn};
SchnitzcellsFolder=XLS{SchnitzRow,ComputerColumn};



if ~isempty(varargin)
    Prefix=varargin{1};
else
    FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
    Dashes=strfind(FolderTemp,'\');
    Prefix=FolderTemp((Dashes(end)+1):end);
end



%Load the current AP detection data
load([DropboxFolder,filesep,Prefix,'\APDetection.mat'])

%Load the full embryo image
APImage=imread([DropboxFolder,filesep,Prefix,'\APDetection\FullEmbryo.tif']);



APImageFig=figure;
%Range for image display
DisplayRange=[min(min(APImage)),max(max(APImage))];


%Now, do the correction
cc=1;
while (cc~=13)

    
    
    imshow(APImage,DisplayRange)
    axis image
    axis off
    title('Anterior (green), posterior (red); original')
    hold on
    try
        plot(coordA(1),coordA(2),'g.','MarkerSize',20);        
    end
    try
        plot(coordP(1),coordP(2),'r.','MarkerSize',20);
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
        coordA=ginput(1);
        
    elseif (ct~=0)&(cc=='p')    %Select posterior end
        coordP=ginput(1);
    elseif (ct~=0)&(cc=='.')    %Increase contrast
        DisplayRange(2)=DisplayRange(2)/2;
        
    elseif (ct~=0)&(cc==',')    %Decrease contrast
        DisplayRange(2)=DisplayRange(2)*2;
        
    elseif (ct~=0)&(cc=='r')    %Reset the contrast
        DisplayRange=[min(min(APImage)),max(max(APImage))];

    elseif (ct==0)&(strcmp(get(APImageFig,'SelectionType'),'alt')) %Delete the point that was clicked on
        cc=1;
    
        [Dummy,MinIndex]=min((cm(1,1)-[coordA(1),coordP(1)]).^2+(cm(1,2)-[coordA(2),coordP(2)]).^2);
            
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
    end
end
            
%Save the information
if exist('xShift')
    save([DropboxFolder,filesep,Prefix,'\APDetection.mat'],'coordA','coordP',...
        'xShift','yShift');
else
   save([DropboxFolder,filesep,Prefix,'\APDetection.mat'],'coordA','coordP',...
        'xShift1','yShift1','xShift2','yShift2');
end
    
%Redo the diagnostic plots

diagFigure = figure;
imagesc(APImage)
axis image
axis off
title('Anterior (green), posterior (red); corrected')
hold on
plot(coordA(1),coordA(2),'g.','MarkerSize',20);
plot(coordP(1),coordP(2),'r.','MarkerSize',20);
hold off
saveas(gcf, [DropboxFolder,filesep,Prefix,'\APDetection\APEmbryo-Manual.tif']);
close(diagFigure);

