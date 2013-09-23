function CheckDivisionTimes(varargin)

%The idea is to determine division times with a higher spatial resolution
%by doing it per AP bin.

close all

%Find out which computer this is. That will determine the folder structure.
%Information about about folders

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


%Find out the date it was taken
Dashes=findstr(Prefix,'-');
Date=Prefix(1:Dashes(3)-1);
EmbryoName=Prefix(Dashes(3)+1:end);


%See if we're dealing with a Bcd case
if (~isempty(findstr(Prefix,'Bcd')))&(isempty(findstr(Prefix,'BcdE1')))&...
        (isempty(findstr(Prefix,'NoBcd')))
    SourcePath=[SourcePath,'\..\..\Bcd-GFP'];
end



D=dir([SourcePath,filesep,Date,filesep,EmbryoName,'\*.tif']);

%Get the information about the zoom
ImageInfo = imfinfo([SourcePath,filesep,Date,filesep,EmbryoName,filesep,D(1).name]);

%Figure out the zoom factor
Zoom=ExtractInformationField(ImageInfo(1),'state.acq.zoomFactor=');
Zoom=str2num(Zoom);

if Zoom==4
    Rows=256;
    Columns=512;
elseif Zoom==8
    Rows=256;
    Columns=256;
elseif Zoom==16
    Rows=128;
    Columns=128;
elseif Zoom==2
    Rows=256;
    Columns=256;
else
    error('Include zoom information in CheckDivisionTimes.m')
end



%Get the surface image in the zoomed case
if (~isempty(findstr(Prefix,'Bcd')))&(isempty(findstr(Prefix,'BcdE1')))&...
        (isempty(findstr(Prefix,'NoBcd')))
    D=dir([SourcePath,filesep,Date,'\BcdGFP-HisRFP\AveragedData\*His_*.tif']);
    ZoomImage=imread([SourcePath,filesep,Date,'\BcdGFP-HisRFP\AveragedData\',D(end-10).name]);
else
    D=dir([FISHPath,'\Data\',Prefix,filesep,Prefix,'-His*.tif']);
    ZoomImage=imread([FISHPath,'\Data\',Prefix,filesep,D(end-10).name]);
end



%Get the information about the AP axis as well as the image shifts
%used for the stitching of the two halves of the embryo
load([DropboxFolder,filesep,Prefix,'\APDetection.mat'])

%Angle between the x-axis and the AP-axis
APAngle=atan((coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1)));
APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);


APPosImage=zeros(size(ZoomImage));
[Rows,Columns]=size(ZoomImage);

for i=1:Rows
    for j=1:Columns
        Angle=atan((i-coordAZoom(2))./(j-coordAZoom(1)));
        Distance=sqrt((coordAZoom(2)-i).^2+(coordAZoom(1)-j).^2);
        APPosition=Distance.*cos(Angle-APAngle);
        APPosImage(i,j)=APPosition/APLength;
    end
end


%Bin the pixels along the AP axis
APResolution=0.025;
APbinID=0:APResolution:1;


APPosBinImage=zeros(size(APPosImage));
for i=1:(length(APbinID)-1)
    FilteredMask=(APbinID(i)<=APPosImage)&(APbinID(i+1)>APPosImage);
    
    APPosBinImage=APPosBinImage+FilteredMask*i;
end


%Load the information about the nc from the XLS file
[Num,Txt]=xlsread([DropboxFolder,'\HGMovieDatabaseV2.xlsx']);
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



%Load the division information if it's already there
if exist([DropboxFolder,filesep,Prefix,'\APDivision.mat'])
    load([DropboxFolder,filesep,Prefix,'\APDivision.mat'])
else
    %Matrix where we'll store the information about the divisions
    APDivision=zeros(14,length(APbinID));
end


%Show the frames on top of the AP bins
Overlay=figure;



CurrentFrame=1;
CurrentNC=13;

cc=1;



while (cc~=13)
    
    figure(Overlay)
    
    %Load the image
    %Get the surface image in the zoomed case
    if (~isempty(findstr(Prefix,'Bcd')))&(isempty(findstr(Prefix,'BcdE1')))&...
            (isempty(findstr(Prefix,'NoBcd')))
        HisImage=imread([SourcePath,filesep,Date,'\BcdGFP-HisRFP\AveragedData\',D(CurrentFrame).name]);

    else
        HisImage=imread([FISHPath,'\Data\',Prefix,filesep,D(CurrentFrame).name]);
    end
    
    
    
    
    %Generate the image with the information about divisions
    %Green: This AP bin divided in the current frame
    %Blue: This AP bin has been determined, but divided in another frame
    %than the current one
    %Red: The division of this AP bin has not been determined.
    
    BlueImage=zeros(size(APPosBinImage));
    GreenImage=zeros(size(APPosBinImage));
    RedImage=zeros(size(APPosBinImage));
    
    for i=1:length(APDivision(CurrentNC,:))
        if APDivision(CurrentNC,i)
            if APDivision(CurrentNC,i)==CurrentFrame
                GreenImage(APPosBinImage==i)=i;
            else
                BlueImage(APPosBinImage==i)=i;
            end
        else
            RedImage(APPosBinImage==i)=i;
        end
    end
           
                
                
    
    
    %Combine the images
    HisOverlay=cat(3,mat2gray(HisImage)+mat2gray(RedImage),...
        mat2gray(HisImage)+mat2gray(GreenImage)/2,...
        mat2gray(HisImage)+mat2gray(BlueImage));
    
    imshow(HisOverlay)
    
    set(gcf,'Name',(['Frame: ',num2str(CurrentFrame),'/',num2str(length(D)),...
        '. Current nc:',num2str(CurrentNC)]));
    
    
    figure(Overlay)
    ct=waitforbuttonpress;
    cc=get(Overlay,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    %Move frames
    if (ct~=0)&(cc=='.')&(CurrentFrame<length(D))
        CurrentFrame=CurrentFrame+1;
    elseif (ct~=0)&(cc==',')&(CurrentFrame>1)
        CurrentFrame=CurrentFrame-1;
        
    %Move nc
    elseif (ct~=0)&(cc=='m')&(CurrentNC<14)
        CurrentNC=CurrentNC+1;
        eval(['CurrentFrame=nc',num2str(CurrentNC)]);
    elseif (ct~=0)&(cc=='n')&(CurrentNC>8)&eval(['nc',num2str(CurrentNC-1),'~=0'])
        CurrentNC=CurrentNC-1;    
        eval(['CurrentFrame=nc',num2str(CurrentNC)]);
        
        
        
    %Save
    elseif (ct~=0)&(cc=='s')
        save([DropboxFolder,filesep,Prefix,'\APDivision.mat'],'APDivision')
        display('Data saved')
    %Select a time for division
    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'normal'))
        cc=1;
        if (cm(1,2)>0)&(cm(1,1)>0)&(cm(1,2)<=Rows)&(cm(1,1)<=Columns)
            APDivision(CurrentNC,APPosBinImage(round(cm(1,2)),round(cm(1,1))))=CurrentFrame;
        end
    end
end


save([DropboxFolder,filesep,Prefix,'\APDivision.mat'],'APDivision')
display('Data saved')

