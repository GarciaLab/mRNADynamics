function CheckDivisionTimesDV(varargin)

%The idea is to determine division times with a higher spatial resolution
%by doing it per DV bin.

%m n: Move between nuclear cycles
%, .: Move between frames
%Click: Division of clicked DV bin in current frame
%r  : Reset the information for the current nuclear cycle
%s  : Save the information
%x  : Save and quit


close all

%Find out which computer this is. That will determine the folder structure.
%Information about about folders

[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath,...
    configValues, movieDatabasePath] = DetermineAllLocalFolders(varargin{1});


if ~isempty(varargin)
    Prefix=varargin{1};
    
else
    FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
    Dashes=strfind(FolderTemp,filesep);
    Prefix=FolderTemp((Dashes(end)+1):end);
end


%Find out the date it was taken
Dashes=findstr(Prefix,'-');
Date=Prefix(1:Dashes(3)-1);
EmbryoName=Prefix(Dashes(3)+1:end);


D=dir([PreProcPath,filesep,Prefix,filesep,Prefix,'-His*.tif']);
ZoomImage=imread([PreProcPath,filesep,Prefix,filesep,D(end).name]);



%Get the information about the DV axis as well as the image shifts
%used for the stitching of the two halves of the embryo
load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])

if ~exist('coordPZoom', 'var')
    warning('AddParticlePosition should have been run first. Running it now.')
    AddParticlePosition(Prefix, 'ManualAlignment')
    load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])
end

%Angle between the x-axis and the DV-axis
APAngle=atan((coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1)));
%Correction for if DVAngle is in quadrants II or III
if coordPZoom(1)-coordAZoom(1) < 0
    APAngle = APAngle + pi;
end
APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);

%need to double-check this.
DVAngle = APAngle + (pi/2);

DVPosImage=zeros(size(ZoomImage));
[Rows,Columns]=size(ZoomImage);

for i=1:Rows
    for j=1:Columns
        Angle=atan((i-coordAZoom(2))./(j-coordAZoom(1)));
        % Correction for if Angle is in quadrant II
        if (j-coordAZoom(1) < 0)
            Angle = Angle + pi;
        end
        
        
        Distance=sqrt((coordAZoom(2)-i).^2+(coordAZoom(1)-j).^2);
        DVPosition=Distance.*cos(Angle-DVAngle);
        DVPosImage(i,j)=DVPosition/DVLength;
    end
end

[DateFromDateColumn, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, DVResolution,...
    Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
    nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);

ncs=[nc9,nc10,nc11,nc12,nc13,nc14];

DVbinID=0:DVResolution:1;


DVPosBinImage=zeros(size(DVPosImage));
for i=1:(length(DVbinID)-1)
    FilteredMask=(DVbinID(i)<=DVPosImage)&(DVbinID(i+1)>DVPosImage);
    
    DVPosBinImage=DVPosBinImage+FilteredMask*i;
end


%Load the division information if it's already there
if exist([DropboxFolder,filesep,Prefix,filesep,'DVDivision.mat'], 'file')
    load([DropboxFolder,filesep,Prefix,filesep,'DVDivision.mat'], 'DVDivision')
    %Check if we changed the number of DV bins
    if size(DVDivision,2)~=length(DVbinID)
        DVDivision=zeros(14,length(DVbinID));
    end
else
    %Matrix where we'll store the information about the divisions
    DVDivision=zeros(14,length(DVbinID));
end


%Show the frames on top of the DV bins
figureOverlay=figure;
axOverlay = axes(figureOverlay);


try
    CurrentFrame = min(ncs(ncs~=0));
catch
    CurrentFrame=1;
end

CurrentNC=find(ncs,1)+8;

cc=1;



while (cc~='x')
    
    
    if ~isnan(CurrentFrame)
        HisImage=imread([PreProcPath,filesep,Prefix,filesep,D(CurrentFrame).name]);
    end
    
    
    
    
    %Generate the image with the information about divisions
    %Green: This DV bin divided in the current frame
    %Blue: This DV bin has been determined, but divided in another frame
    %than the current one
    %Red: The division of this DV bin has not been determined.
    
    BlueImage=zeros(size(DVPosBinImage));
    GreenImage=zeros(size(DVPosBinImage));
    RedImage=zeros(size(DVPosBinImage));
    
    for i=1:length(DVDivision(CurrentNC,:))
        if DVDivision(CurrentNC,i)
            if DVDivision(CurrentNC,i)==CurrentFrame
                GreenImage(DVPosBinImage==i)=i;
            else
                BlueImage(DVPosBinImage==i)=i;
            end
        else
            RedImage(DVPosBinImage==i)=i;
        end
    end
    
    
    
    
    
    %Combine the images
    HisOverlay=cat(3,mat2gray(HisImage)+mat2gray(RedImage),...
        mat2gray(HisImage)+mat2gray(GreenImage)/2,...
        mat2gray(HisImage)+mat2gray(BlueImage));
    
    imshow(HisOverlay, 'Parent', axOverlay)
    
    set(figureOverlay,'Name',(['Frame: ',num2str(CurrentFrame),'/',num2str(length(D)),...
        '. Current nc:',num2str(CurrentNC)]));
    
    
    %     figure(Overlay)
    ct=waitforbuttonpress;
    cc=get(figureOverlay,'currentcharacter');
    cm=get(axOverlay,'CurrentPoint');
    
    %Move frames
    if (ct~=0)&(cc=='.')&(CurrentFrame<length(D))
        CurrentFrame=CurrentFrame+1;
    elseif (ct~=0)&(cc==',')&(CurrentFrame>1)
        CurrentFrame=CurrentFrame-1;
    elseif (ct~=0)&(cc=='>')& (CurrentFrame+10)< length(D)
        CurrentFrame=CurrentFrame+10;
    elseif (ct~=0)&(cc=='<')&( CurrentFrame-10) > 1
        CurrentFrame=CurrentFrame-10;
        %Move nc
    elseif (ct~=0)&(cc=='m')&(CurrentNC<14) & ~isnan(ncs(CurrentNC+1-8))
        CurrentNC=CurrentNC+1;
        eval(['CurrentFrame=nc',num2str(CurrentNC)]);
    elseif (ct~=0)&(cc=='n')&(CurrentNC>8)&eval(['nc',num2str(CurrentNC-1),'~=0'])
        CurrentNC=CurrentNC-1;
        eval(['CurrentFrame=nc',num2str(CurrentNC)]);
        
        %Reset the information
    elseif (ct~=0)&(cc=='r')
        DVDivision(CurrentNC,:)=0;
        
        %Save
    elseif (ct~=0)&(cc=='s')
        save([DropboxFolder,filesep,Prefix,filesep,'DVDivision.mat'],'DVDivision')
        disp('DVDivision.mat saved.');
        %Select a time for division
    elseif (ct==0)&(strcmp(get(figureOverlay,'SelectionType'),'normal'))
        cc=1;
        if (cm(1,2)>0)&(cm(1,1)>0)&(cm(1,2)<=Rows)&(cm(1,1)<=Columns)
            DVDivision(CurrentNC,DVPosBinImage(round(cm(1,2)),round(cm(1,1))))=CurrentFrame;
        end
        
        %Debug mode
    elseif (ct~=0)&(cc=='9')
        keyboard;
    end
end


save([DropboxFolder,filesep,Prefix,filesep,'DVDivision.mat'],'DVDivision');
disp('DVDivision.mat saved.');
close(figureOverlay);