function CheckNucleiSegmentation(varargin)

%This code allows you to check the nuclear segmentation performed by
%Timon's code implemented in SegmentNucleiTimon.m

%To do:
%1) Allow to edit the size and angle of an ellipse

%Usage:

%.  - Move a frame forward and keep in the new frame only areas that
%     overlap with the previous ones.
%,  - Move a frame backwards
%j  - Jump to a frame
%s  - Save current analysis
%m  - Increase contrast
%n  - Decrease contrast
%r  - Reset contrast setting


%right click  - delete region
%left click - add region with default nc radius





close all


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


%Load the folder information
[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders;


if ~isempty(varargin)
    Prefix=varargin{1};
else
    FolderTemp=uigetdir(DefaultDropboxFolder,'Choose folder with files to analyze');
    Dashes=strfind(FolderTemp,'\');
    Prefix=FolderTemp((Dashes(end)+1):end);
end


[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders(Prefix);


%Set the source folders
Folder=[FISHPath,'\Analysis\',Prefix,'_\preanalysis\'];
FileName=['CompactResults_',Prefix,'_.mat'];

%Set the destination folders
OutputFolder=[DropboxFolder,filesep,Prefix];
FilePrefix=FileName(16:end-4);
DataFolder=[Folder,'..\..\..\Data\',FilePrefix(1:end-1)];


%Find out how many frames we have
D=dir([FISHPath,'\Data\',Prefix,'\*-His*.tif']);
if length(D)==0
    warning('The name format is a mess. I had to do this for KITP')
    D=dir([FISHPath,'\Data\',Prefix,'\*_His*.tif']);
end
TotalFrames=length(D);



%Load the information about the nc from the XLS file
[Num,Txt]=xlsread([DefaultDropboxFolder,'\HGMovieDatabaseV2.xlsx']);
XLSHeaders=Txt(1,:);
Txt=Txt(2:end,:);

%Find the different columns.
DataFolderColumn=find(strcmp(XLSHeaders,'DataFolder'));
nc9Column=find(strcmp(XLSHeaders,'nc9'));
nc10Column=find(strcmp(XLSHeaders,'nc10'));
nc11Column=find(strcmp(XLSHeaders,'nc11'));
nc12Column=find(strcmp(XLSHeaders,'nc12'));
nc13Column=find(strcmp(XLSHeaders,'nc13'));
nc14Column=find(strcmp(XLSHeaders,'nc14'));
CFColumn=find(strcmp(XLSHeaders,'CF'));
Channel2Column=find(strcmp(XLSHeaders,'Channel2'));

%Convert the prefix into the string used in the XLS file
Dashes=findstr(Prefix,'-');

%Find the corresponding entry in the XLS file
if (~isempty(findstr(Prefix,'Bcd')))&(isempty(findstr(Prefix,'BcdE1')))&...
        (isempty(findstr(Prefix,'NoBcd')))
    XLSEntry=find(strcmp(Txt(:,DataFolderColumn),...
        [Date,'\BcdGFP-HisRFP']));
else
    XLSEntry=find(strcmp(Txt(:,DataFolderColumn),...
        [Prefix(1:Dashes(3)-1),filesep,Prefix(Dashes(3)+1:end)]));
end


if strcmp(Txt(XLSEntry,Channel2Column),'His-RFP')
    nc9=Num(XLSEntry,nc9Column-6);
    nc10=Num(XLSEntry,nc10Column-6);
    nc11=Num(XLSEntry,nc11Column-6);
    nc12=Num(XLSEntry,nc12Column-6);
    nc13=Num(XLSEntry,nc13Column-6);
    nc14=Num(XLSEntry,nc14Column-6);
    %This is in case the last column for CF is all nan and is not part of
    %the Num matrix
    if size(Num,2)==CFColumn-6    
        CF=Num(XLSEntry,CFColumn-6);
    else
        CF=nan;
    end
else
    error('nc information not define in HGMovieDatabase.xlsx')
end



%Get the nuclei segmentation data
load([DropboxFolder,filesep,Prefix,'\Ellipses.mat']);


% 
% 
% 
% %Compare the results
% Image=imread([FISHPath,'\Data\',Prefix,filesep,Prefix,'-His_',iIndex(nc14+15,3),'.tif']);
% figure(1)
% imshow(Image,[])
% hold on
% PlotHandle=[];
% for i=1:length(Ellipses)
%     PlotHandle=[PlotHandle,ellipse(Ellipses(i,3),Ellipses(i,4),...
%         Ellipses(i,5),Ellipses(i,1)+1,Ellipses(i,2)+1)];
% end
% hold off
% set(PlotHandle,'Color','r')
        



%Get the information about the Histone channel images
D=dir([FISHPath,'\Data\',Prefix,'\*-His*.tif']);
if length(D)==0
    warning('The name format is a mess. I had to do this for KITP')
    D=dir([FISHPath,'\Data\',Prefix,'\*_His*.tif']);
end
TotalFrames=length(D);


%Get information about the image size
HisImage=imread([FISHPath,'\Data\',Prefix,filesep,D(1).name]);
[Rows,Cols]=size(HisImage);
DisplayRange=[min(min(HisImage)),max(max(HisImage))];



%Make a vector containing the nc corresponding to each frame
for i=1:length(D)
    if i<nc9
        nc(i)=8;
    elseif (i>nc9)&(i<nc10)
        nc(i)=9;
    elseif (i>nc10)&(i<nc11)
        nc(i)=10;
    elseif (i>nc11)&(i<nc12)
        nc(i)=11;
    elseif (i>nc12)&(i<nc13)
        nc(i)=12;
    elseif (i>nc13)&(i<nc14)
        nc(i)=13;
    elseif i>nc14
        nc(i)=14;
    end
end


%Set the figure sizes for Albert or the other computers
if strcmp(name(1:end-1),'albert-pc')
    Overlay=figure;
    set(gcf,'Position',[10         113        1376         872])

    OriginalImage=figure;
    set(gcf,'Position',[1398         659         512         326])
else
    Overlay=figure;
    set(gcf,'Position',[6   603   676   342])

    OriginalImage=figure;
    set(gcf,'Position',[ 704   382   512   256])
end


CurrentFrame=1;
cc=1;

while (cc~=13)
    
    %Load the image
    HisImage=imread([FISHPath,'\Data\',Prefix,filesep,D(CurrentFrame).name]);
    
    
    %Get the information about the centroids
    [NCentroids,Dummy]=size(Ellipses{CurrentFrame});

    
    figure(Overlay)
    imshow(HisImage,DisplayRange,'Border','Loose')
    hold on
    PlotHandle=[];
    for i=1:NCentroids
        PlotHandle=[PlotHandle,ellipse(Ellipses{CurrentFrame}(i,3),...
            Ellipses{CurrentFrame}(i,4),...
            Ellipses{CurrentFrame}(i,5),Ellipses{CurrentFrame}(i,1)+1,...
            Ellipses{CurrentFrame}(i,2)+1)];
    end
    hold off
    set(PlotHandle,'Color','r')
    
     

    FigureTitle=['Frame: ',num2str(CurrentFrame),'/',num2str(TotalFrames),...
        ', nc: ',num2str(nc(CurrentFrame))];
    title(FigureTitle)
  
    
    figure(OriginalImage)
    imshow(HisImage,DisplayRange,'Border','Tight')
    
    figure(Overlay)
    ct=waitforbuttonpress;
    cc=get(Overlay,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    


    if (ct~=0)&(cc=='.')&(CurrentFrame<TotalFrames)
        CurrentFrame=CurrentFrame+1;
        %DisplayRange=[min(min(HisImage)),max(max(HisImage))];
    elseif (ct~=0)&(cc==',')&(CurrentFrame>1)
        CurrentFrame=CurrentFrame-1;
        %DisplayRange=[min(min(HisImage)),max(max(HisImage))];
    elseif (ct~=0)&(cc=='s')
        save([DropboxFolder,filesep,Prefix,'\Ellipses.mat'],'Ellipses')
        display('Ellipses saved.')
    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'normal'))
        cc=1;
        if (cm(1,2)>0)&(cm(1,1)>0)&(cm(1,2)<=Rows)&(cm(1,1)<=Cols)
            
            %Add a circle to this location with the mean radius of the
            %ellipses found in this frame

            %(x, y, a, b, theta, maxcontourvalue, time,
            %particle_id)
            MeanRadius=mean((Ellipses{CurrentFrame}(:,3)+Ellipses{CurrentFrame}(:,4))/2);
            Ellipses{CurrentFrame}(end+1,:)=...
                [cm(1,1),cm(1,2),MeanRadius,MeanRadius,0,0,0,0];
        end
            
            
            
        
    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'alt'))
        cc=1;
        if (cm(1,2)>0)&(cm(1,1)>0)&(cm(1,2)<=Rows)&(cm(1,1)<=Cols)
            %Find out which ellipses we clicked on so we can delete it
            
            %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
            Distances=sqrt((Ellipses{CurrentFrame}(:,1)-cm(1,1)).^2+...
                (Ellipses{CurrentFrame}(:,2)-cm(1,2)).^2);
            [MinValue,MinIndex]=min(Distances);

            Ellipses{CurrentFrame}=[Ellipses{CurrentFrame}(1:MinIndex-1,:);...
                Ellipses{CurrentFrame}(MinIndex+1:end,:)];
        end
    
    elseif (ct~=0)&(cc=='j')
        iJump=input('Frame to jump to: ');
        if (floor(iJump)>0)&(iJump<=TotalFrames)
            CurrentFrame=iJump;
        end
        
    elseif (ct~=0)&(cc=='m')    %Increase contrast
        DisplayRange(2)=DisplayRange(2)/1.5;
        
    elseif (ct~=0)&(cc=='n')    %Decrease contrast
        DisplayRange(2)=DisplayRange(2)*1.5;
        
    elseif (ct~=0)&(cc=='r')    %Reset the contrast
        DisplayRange=[min(min(HisImage)),max(max(HisImage))];
    end
end



save([DropboxFolder,filesep,Prefix,'\Ellipses.mat'],'Ellipses')
display('Ellipses saved.')






