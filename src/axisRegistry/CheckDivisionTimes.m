function CheckDivisionTimes(Prefix, varargin)

%The idea is to determine division times with a higher spatial resolution
%by doing it per AP bin.

%m n: Move between nuclear cycles
%, .: Move between frames
%Click: Division of clicked AP bin in current frame
%r  : Reset the information for the current nuclear cycle
%s  : Save the information
%x  : Save and quit


%Find out which computer this is. That will determine the folder structure.
%Information about about folders

liveExperiment = LiveExperiment(Prefix);

[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath,...
    configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);


lazy = false;
for arg = 1:length(varargin)
    if strcmpi(varargin{arg}, 'lazy')
        lazy = true;
    end
end


%Find out the date it was taken
Dashes=findstr(Prefix,'-');
Date=Prefix(1:Dashes(3)-1);
EmbryoName=Prefix(Dashes(3)+1:end);


%See if we're dealing with a Bcd case
if (~isempty(findstr(Prefix,'Bcd')))&(isempty(findstr(Prefix,'BcdE1')))&...
        (isempty(findstr(Prefix,'NoBcd')))&(isempty(findstr(Prefix,'Bcd1x')))
    SourcePath=[SourcePath,filesep,'..',filesep,'..',filesep,'Bcd-GFP'];
end



try 
    hisMat = getHisMat(liveExperiment); 
    ZoomImage = hisMat(end);
    nFrames = size(hisMat, 3);
catch
    
D=dir([PreProcPath,filesep,Prefix,filesep,Prefix,'-His*.tif']);
ZoomImage=imread([PreProcPath,filesep,Prefix,filesep,D(end).name]);
nFrames = length(D);
end




%Get the information about the AP axis as well as the image shifts
%used for the stitching of the two halves of the embryo
load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])

if ~lazy
    if ~exist('coordPZoom', 'var')
        warning('AddParticlePosition should have been run first. Running it now.')
        AddParticlePosition(Prefix, 'ManualAlignment')
        load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])
    end
end

%Angle between the x-axis and the AP-axis
APAngle=atan((coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1)));
%Correction for if APAngle is in quadrants II or III
if coordPZoom(1)-coordAZoom(1) < 0
    APAngle = APAngle + pi;
end
APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);


APPosImage=zeros(size(ZoomImage));
[Rows,Columns]=size(ZoomImage);

for i=1:Rows
    for j=1:Columns
        Angle=atan((i-coordAZoom(2))./(j-coordAZoom(1)));
        % Correction for if Angle is in quadrant II
        if (j-coordAZoom(1) < 0)
            Angle = Angle + pi;
        end
        
        
        Distance=sqrt((coordAZoom(2)-i).^2+(coordAZoom(1)-j).^2);
        APPosition=Distance.*cos(Angle-APAngle);
        APPosImage(i,j)=APPosition/APLength;
    end
end

[DateFromDateColumn, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
    Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
    nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);

ncs=[nc9,nc10,nc11,nc12,nc13,nc14];

APbinID=0:APResolution:1;


APPosBinImage=zeros(size(APPosImage));
for i=1:(length(APbinID)-1)
    FilteredMask=(APbinID(i)<=APPosImage)&(APbinID(i+1)>APPosImage);
    
    APPosBinImage=APPosBinImage+FilteredMask*i;
end


%Load the division information if it's already there
if exist([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat'], 'file')
    load([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat'], 'APDivision')
    %Check if we changed the number of AP bins
    if size(APDivision,2)~=length(APbinID)
        APDivision=zeros(14,length(APbinID));
    end
else
    %Matrix where we'll store the information about the divisions
    APDivision=zeros(14,length(APbinID));
end


%lazy mode
if lazy
    apsPresent = unique(APPosBinImage(:));
    for nc = 1:length(ncs)
        for ap = apsPresent
            if ncs(nc)~=0
                APDivision(nc+8, ap) = ncs(nc);
            end
        end
    end
end

%Show the frames on top of the AP bins
figureOverlay=figure;
axOverlay = axes(figureOverlay);


try
    CurrentFrame = min(ncs(ncs~=0));
catch
    CurrentFrame=1;
end

CurrentNC=find(ncs,1)+8;

cc=1;


if ~lazy
    while (cc~='x')
    
    
    if ~isnan(CurrentFrame)
        try
            HisImage = hisMat(:, :, CurrentFrame);
        catch
            HisImage=imread([PreProcPath,filesep,Prefix,filesep,D(CurrentFrame).name]);
        end
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
    
    imshow(HisOverlay, 'Parent', axOverlay)
    
    set(figureOverlay,'Name',(['Frame: ',num2str(CurrentFrame),'/',num2str(nFrames),...
        '. Current nc:',num2str(CurrentNC)]));
    
    
    %     figure(Overlay)
    ct=waitforbuttonpress;
    cc=get(figureOverlay,'currentcharacter');
    cm=get(axOverlay,'CurrentPoint');
    
    %Move frames
    if (ct~=0)&(cc=='.')&(CurrentFrame<nFrames)
        CurrentFrame=CurrentFrame+1;
    elseif (ct~=0)&(cc==',')&(CurrentFrame>1)
        CurrentFrame=CurrentFrame-1;
    elseif (ct~=0)&(cc=='>')& (CurrentFrame+10)< nFrames
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
        APDivision(CurrentNC,:)=0;
        
        %Save
    elseif (ct~=0)&(cc=='s')
        save([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat'],'APDivision')
        disp('APDivision.mat saved.');
        %Select a time for division
    elseif (ct==0)&(strcmp(get(figureOverlay,'SelectionType'),'normal'))
        cc=1;
        if (cm(1,2)>0)&(cm(1,1)>0)&(cm(1,2)<=Rows)&(cm(1,1)<=Columns)
            APDivision(CurrentNC,APPosBinImage(round(cm(1,2)),round(cm(1,1))))=CurrentFrame;
        end
        
        %Debug mode
    elseif (ct~=0)&(cc=='9')
        keyboard;
    end
end
end


save([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat'],'APDivision');
disp('APDivision.mat saved.');
close(figureOverlay);