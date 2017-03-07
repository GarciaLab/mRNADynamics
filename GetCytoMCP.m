function [Mean,SD,Median,Max,...
    MeanAPProfile,SDAPProfile,SEAPProfile]=GetCytoMCP(Prefix)




%Are the observed differences in offset related to the total amount of
%fluorescent protein that the mother deposited in the embryo?

%close all

%Get the folders
[SourcePath,ProcPath,DefaultDropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;
[SourcePath,ProcPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);

%Load FrameInfo
load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])

%Find out the image size
ImageSize=[FrameInfo(1).LinesPerFrame,FrameInfo(1).PixelsPerLine];
%How many slices do we have?
ZSlices=FrameInfo(1).NumberSlices;

%Folder for report figures
mkdir([DropboxFolder,filesep,Prefix,filesep,'CytoFluo'])


if ~exist([ProcPath,filesep,Prefix,filesep,'CytoImages.mat'])

    %Get the flat field and smooth it with a Gaussian.
    
    FFDir=dir([PreProcPath,filesep,Prefix,filesep,'*FF.tif']);
    if ~isempty(FFDir)
        FFImage=double(imread([PreProcPath,filesep,Prefix,filesep,FFDir(1).name],1));
        filtStd=30;         %This came from the FISH code, in pixels.
        FFImage=imfilter(FFImage,fspecial('gaussian',2*filtStd,filtStd),'symmetric');
        FFImage=imdivide(FFImage,double(max(FFImage(:))));
    else
        warning('No flat field found')
        FFImage=ones(ImageSize);
    end
    

    %Get the ellipses
    load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'])
 
    %This is the structuring element we will use to dilate the mask. 1.5um
    %might be a good choice.
    SESize=round(1.5/FrameInfo(1).PixelSize);
    StrucElement=strel('disk',SESize);

    h=waitbar(0,'Calculating the maximum projection');
    for i=1:length(Ellipses)
        waitbar(i/length(Ellipses),h)

        %Create a mask of the ellipses

        %Initialize the image. In order to avoid problems with the drawing of
        %ellipses we'll make the image larger and then reduce it.
        IncreasePixels=200;
        Mask=zeros(ImageSize+IncreasePixels*2);

        %How many ellipses do we have?
        [NEllipses,Dummy]=size(Ellipses{i});

        for j=1:NEllipses
            %Put ellipses at their corresponding positions. THe value of each
            %ellipse corresponds to its index number
            Mask=ellipseMatrix(round(Ellipses{i}(j,2))+IncreasePixels,...
                round(Ellipses{i}(j,1))+IncreasePixels,...
                Ellipses{i}(j,3), Ellipses{i}(j,4), Ellipses{i}(j,5), Mask, 1);
        end

        %Go back to the original image size
        Mask=Mask(IncreasePixels+1:end-IncreasePixels,IncreasePixels+1:end-IncreasePixels);
        Mask=imdilate(Mask,StrucElement);
        Mask=~Mask;


        %Now, get the information for each Z slice. Note that I changed
        %this to account for the fact that we now have blank images at the
        %beginning and end of the stack
        for j=2:ZSlices-1
            MCPImage(:,:,j-1)=double(imread([PreProcPath,filesep,Prefix,filesep,filesep,Prefix,'_',iIndex(i,3),'_z',...
                iIndex(j,2),'.tif']));
            MCPImage(:,:,j-1)=imdivide(MCPImage(:,:,j-1),FFImage);
        end

        %Find the maximum and the mean
        MaxImageTemp=max(MCPImage,[],3);
        MeanImageTemp=mean(MCPImage,3);

        MaxImage(:,:,i)=immultiply(MaxImageTemp,Mask);
        MeanImage(:,:,i)=immultiply(MeanImageTemp,Mask);

    end
    close(h)

    %Save to the FISH path so that we don't overwhelm the Dropbox folder!
    save([ProcPath,filesep,Prefix,'_',filesep,'CytoImages.mat'],'MaxImage','MeanImage')
else
    display('Using saved CytoImages.mat located in the PreProcessed folder.')
    
    load([ProcPath,filesep,Prefix,filesep,'CytoImages.mat']);
end

%Get the information
for i=1:size(MeanImage,3)
    MeanImageTemp=MeanImage(:,:,i);
    MeanImageTemp=MeanImageTemp(:);
    MeanImageTemp=MeanImageTemp(MeanImageTemp~=0);
    
    Mean(i)=mean(MeanImageTemp);
    SD(i)=std(MeanImageTemp);
    Median(i)=median(MeanImageTemp);
    
    MaxImageTemp=MaxImage(:,:,i);
    MaxImageTemp=MaxImageTemp(:);
    MaxImageTemp=MaxImageTemp(MaxImageTemp~=0);
    Max(i)=mean(MaxImageTemp);
end


%Generate reporting plots


%Get profile information as a function of AP.

%Get the AP information
load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])
    
%Now, assign an AP value to each pixel

%Angle between the x-axis and the AP-axis
APAngle=atan2((coordPZoom(2)-coordAZoom(2)),(coordPZoom(1)-coordAZoom(1)));
APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);

APPosImage=zeros(ImageSize);
Rows=ImageSize(1);
Cols=ImageSize(2);

for i=1:Rows
    for j=1:Cols
        Angle=atan2((i-coordAZoom(2)),(j-coordAZoom(1)));
        Distance=sqrt((coordAZoom(2)-i).^2+(coordAZoom(1)-j).^2);
        APPosition=Distance.*cos(Angle-APAngle);
        APPosImage(i,j)=APPosition/APLength;
    end
end



%Bin the pixels along the AP axis
[XLSNum,XLSTxt,XLSRaw]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
Dashes=findstr(Prefix,'-');
PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));
    if isempty(PrefixRow)
        PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'/',Prefix(Dashes(3)+1:end)]));
        if isempty(PrefixRow)
            error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
        end
    end
APResolutionColumn = find(strcmp(XLSRaw(1,:),'APResolution'));
APResolution = XLSRaw{PrefixRow,APResolutionColumn};
APbinID=0:APResolution:1;

MeanAPProfile=nan(length(APbinID),size(MeanImage,3));
SDAPProfile=nan(length(APbinID),size(MeanImage,3));
SEAPProfile=nan(length(APbinID),size(MeanImage,3));



for j=1:size(MeanImage,3)
    MeanImageTemp=MeanImage(:,:,j);
    
    for i=1:(length(APbinID)-1)
        FilteredMask=(APbinID(i)<=APPosImage)&(APbinID(i+1)>APPosImage);

        %Check that this filter mask region applies to the current data set
        if sum(sum(FilteredMask))
            FilteredMeanImage=immultiply(MeanImageTemp,FilteredMask);

            %Check that we will have enough statistics
            if sum(sum(FilteredMeanImage>0))>10
            
                MeanAPProfile(i,j)=mean(FilteredMeanImage(FilteredMeanImage>0));
                SDAPProfile(i,j)=std(FilteredMeanImage(FilteredMeanImage>0));
                SEAPProfile(i,j)=std(FilteredMeanImage(FilteredMeanImage>0))/...
                    sqrt(sum(sum(FilteredMeanImage>0)));
            end
        end
    end
end




colormap(jet(128));
cmap=colormap;

for j=1:size(MeanImage,3)
    ColorTime(j,:)= cmap(round((j-1)/...
        (size(MeanImage,3)-1)*127+1),:);
end

figure(1)
hold on
for i=1:size(MeanImage,3)
    plot(APbinID,MeanAPProfile(:,i),'-','color',ColorTime(i,:))
end
hold off
box on
h = colorbar;
xlabel('AP position (x/L)')
ylabel('Cyto fluorescence (au)')
ylabel(h,'Time (frames)')
%Figure out the x axis range
APbins=find(sum(~isnan(MeanAPProfile)'));
if length(APbins)>1
    xlim([APbinID(APbins(1)),APbinID(APbins(end))])
else
    xlim([APbinID(APbins)*0.8,APbinID(APbins)*1.2])
end
saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'CytoFluo',filesep,'CytoFluoAPTime.tif'])


figure(2)
hold on
for i=1:20:size(MeanImage,3)
    errorbar(APbinID,MeanAPProfile(:,i),SDAPProfile(:,i),'.-','color',ColorTime(i,:))
end
hold off
box on
h = colorbar;
xlabel('AP position (x/L)')
ylabel('Cyto fluorescence (au)')
ylabel(h,'Time (frames)')
if length(APbins)>1
    xlim([APbinID(APbins(1)),APbinID(APbins(end))])
else
    xlim([APbinID(APbins)*0.8,APbinID(APbins)*1.2])
end
saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'CytoFluo',filesep,'CytoFluoAPTime-ErrorBars.tif'])




%Create an image of fluorescence as a function of time by averaging over
%the whole field of view


figure(3)
errorbar(1:length(Mean),Mean,SD,'.-k')
xlabel('Frame')
ylabel('Mean cyto fluorescence (\pm SD)')
saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'CytoFluo',filesep,'CytoFluoTime.tif'])
