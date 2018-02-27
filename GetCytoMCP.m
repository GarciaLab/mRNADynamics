function [Mean,SD,Median,Max,...
    MeanAPProfile,SDAPProfile,SEAPProfile]=GetCytoMCP(Prefix)




%Are the observed differences in offset related to the total amount of
%fluorescent protein that the mother deposited in the embryo?

%Get the folders
[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, SchnitzcellsFolder,...
configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);

%Load FrameInfo
load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])

%Do we have multiple channels? In that case, we will extract the
%cytoplasmic fluorescence from each channel.
NChannels=...
    length(dir([PreProcPath,filesep,Prefix,filesep,Prefix,'_001_z01*.tif']));

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
        for ChN=1:NChannels
            try
                FFImage{ChN}=double(imread([PreProcPath,filesep,Prefix,filesep,FFDir(1).name],ChN));
            catch
                FFImage{ChN}=double(imread([PreProcPath,filesep,Prefix,filesep,FFDir(1).name],1));
            end
            filtStd=30;         %This came from the FISH code, in pixels.
            FFImage{ChN}=imfilter(FFImage{ChN},fspecial('gaussian',2*filtStd,filtStd),'symmetric');
            FFImage{ChN}=imdivide(FFImage{ChN},double(max(FFImage{ChN}(:))));
        end
    else
        warning('No flat field found')
        for ChN=1:NChannels
            FFImage{ChN}=ones(ImageSize);
        end
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


        for ChN=1:NChannels
            %Now, get the information for each Z slice. Note that I changed
            %this to account for the fact that we now have blank images at the
            %beginning and end of the stack
            
            for j=2:3%ZSlices-1
                %I need to do this because the naming convention can be
                %different when I have only one channel.
%                 if NChannels>1
%                     FileName=[Prefix,'_',iIndex(i,3),'_z',iIndex(j,2),...
%                         '_ch',iIndex(ChN,2),'.tif'];
%                 else
%                     FileName=[Prefix,'_',iIndex(i,3),'_z',iIndex(j,2),'.tif'];
%                 end
                % YJK : Now, since I changed the ExportDataForFISH so that it will export
                % the data as ch01 even in case there is only one channel (protein) ( this
                % change was for the TrackNuclei), we don't need this separation.
                FileName=[Prefix,'_',iIndex(i,3),'_z',iIndex(j,2),...
                         '_ch',iIndex(ChN,2),'.tif'];
               
                MCPImage(:,:,j-1)=double(imread([PreProcPath,filesep,Prefix,filesep,filesep,FileName]));
                MCPImage(:,:,j-1)=imdivide(MCPImage(:,:,j-1),FFImage{ChN});
            end

            %Find the maximum and the mean
            MaxImageTemp=max(MCPImage,[],3);
            MeanImageTemp=mean(MCPImage,3);
            MedianImageTemp=median(MCPImage,3);

            MaxImage{ChN}(:,:,i)=immultiply(MaxImageTemp,Mask);
            MeanImage{ChN}(:,:,i)=immultiply(MeanImageTemp,Mask);
            MedianImage{ChN}(:,:,i)=immultiply(MedianImageTemp,Mask);
        end
    end
    close(h)

    %Save to the FISH path so that we don't overwhelm the Dropbox folder!
    save([ProcPath,filesep,Prefix,'_',filesep,'CytoImages.mat'],'MaxImage','MeanImage','MedianImage')
else
    display('Using saved CytoImages.mat located in the PreProcessed folder.')

    load([ProcPath,filesep,Prefix,filesep,'CytoImages.mat']);
end
    
    
%Get the information
for ChN=1:NChannels
    for i=1:size(MeanImage,3)
        MeanImageTemp=MeanImage{ChN}(:,:,i);
        MeanImageTemp=MeanImageTemp(:);
        MeanImageTemp=MeanImageTemp(MeanImageTemp~=0);

        Mean{ChN}(i)=mean(MeanImageTemp);
        SD{ChN}(i)=std(MeanImageTemp);
        Median{ChN}(i)=median(MeanImageTemp);

        MaxImageTemp=MaxImage{ChN}(:,:,i);
        MaxImageTemp=MaxImageTemp(:);
        MaxImageTemp=MaxImageTemp(MaxImageTemp~=0);
        Max{ChN}(i)=mean(MaxImageTemp);
    end
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
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder)

APbinID=0:APResolution:1;


for ChN=1:NChannels
    MeanAPProfile{ChN}=nan(length(APbinID),size(MeanImage{ChN},3));
    SDAPProfile{ChN}=nan(length(APbinID),size(MeanImage{ChN},3));
    SEAPProfile{ChN}=nan(length(APbinID),size(MeanImage{ChN},3));
end


for ChN=1:NChannels
    for j=1:size(MeanImage{ChN},3)
        MeanImageTemp=MeanImage{ChN}(:,:,j);

        for i=1:(length(APbinID)-1)
            FilteredMask=(APbinID(i)<=APPosImage)&(APbinID(i+1)>APPosImage);

            %Check that this filter mask region applies to the current data set
            if sum(sum(FilteredMask))
                FilteredMeanImage=immultiply(MeanImageTemp,FilteredMask);

                %Check that we will have enough statistics
                if sum(sum(FilteredMeanImage>0))>10
                    MeanAPProfile{ChN}(i,j)=mean(FilteredMeanImage(FilteredMeanImage>0));
                    SDAPProfile{ChN}(i,j)=std(FilteredMeanImage(FilteredMeanImage>0));
                    SEAPProfile{ChN}(i,j)=std(FilteredMeanImage(FilteredMeanImage>0))/...
                        sqrt(sum(sum(FilteredMeanImage>0)));
                end
            end
        end
    end
end



colormap(jet(128));
cmap=colormap;

for ChN=1:NChannels

    for j=1:size(MeanImage{ChN},3)
        ColorTime(j,:)= cmap(round((j-1)/...
            (size(MeanImage{ChN},3)-1)*127+1),:);
    end

    figure((i-1)*3+1)
    hold on
    for i=1:size(MeanImage{ChN},3)
        plot(APbinID,MeanAPProfile{ChN}(:,i),'-','color',ColorTime(i,:))
    end
    hold off
    box on
    h = colorbar;
    xlabel('AP position (x/L)')
    ylabel('Cyto fluorescence (au)')
    ylabel(h,'Time (frames)')
    %Figure out the x axis range
    APbins=find(sum(~isnan(MeanAPProfile{ChN})'));
    if length(APbins)>1
        xlim([APbinID(APbins(1)),APbinID(APbins(end))])
    else
        xlim([APbinID(APbins)*0.8,APbinID(APbins)*1.2])
    end
    title(['Channel ',num2str(ChN)])
    saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'CytoFluo',filesep,...
        'CytoFluoAPTime_ch',iIndex(ChN,2),'.tif'])


    figure((i-1)*3+2)
    hold on
    for i=1:20:size(MeanImage{ChN},3)
        errorbar(APbinID,MeanAPProfile{ChN}(:,i),SDAPProfile{ChN}(:,i),'.-','color',ColorTime(i,:))
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
    title(['Channel ',num2str(ChN)])
    saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'CytoFluo',filesep,...
        'CytoFluoAPTime-ErrorBars_ch',iIndex(ChN,2),'.tif'])



    %Create an image of fluorescence as a function of time by averaging over
    %the whole field of view
    figure((i-1)*3+3)
    errorbar(1:length(Mean{ChN}),Mean{ChN},SD{ChN},'.-k')
    xlabel('Frame')
    ylabel('Mean cyto fluorescence (\pm SD)')
        title(['Channel ',num2str(ChN)])
    saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'CytoFluo',filesep,...
        'CytoFluoTime_ch',iIndex(ChN,2),'.tif'])
end