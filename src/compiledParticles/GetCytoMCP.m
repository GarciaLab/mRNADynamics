function [Mean,SD,Median,Max,...
    MeanAPProfile,SDAPProfile,SEAPProfile]=GetCytoMCP(Prefix)

%Are the observed differences in offset related to the total amount of
%fluorescent protein that the mother deposited in the embryo?


liveExperiment = LiveExperiment(Prefix);

%Get the folders
[SourcePath, ProcPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath,...
    configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);

[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoopEnd, APResolution,...
    Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
    nc9, nc10, nc11, nc12, nc13, nc14, CF,Channel3,prophase,metaphase, anaphase] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);


%Load FrameInfo
load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo')

%Do we have multiple channels? In that case, we will extract the
%cytoplasmic fluorescence from each channel.

NChannels = length(liveExperiment.spotChannels);

%Find out the image size
ImageSize=[FrameInfo(1).LinesPerFrame,FrameInfo(1).PixelsPerLine];
%How many slices do we have?
ZSlices=FrameInfo(1).NumberSlices;

%Folder for depositing figures
mkdir([DropboxFolder,filesep,Prefix,filesep,'CytoFluo'])

if ~exist([ProcPath,filesep,Prefix,filesep,'CytoImages.mat'], 'file')
    %
    %     %Get the flat field and smooth it with a Gaussian.
    %     FFDir=dir([PreProcPath,filesep,Prefix,filesep,'*FF.tif']);
    %     if ~isempty(FFDir)
    %         for ChN=1:NChannels
    %             try
    %                 FFImage{ChN}=double(imread([PreProcPath,filesep,Prefix,filesep,FFDir(1).name],ChN));
    %             catch
    %                 FFImage{ChN}=double(imread([PreProcPath,filesep,Prefix,filesep,FFDir(1).name],1));
    %             end
    %             filtStd=30;         %This came from the FISH code, in pixels.
    %             FFImage{ChN}=imfilter(FFImage{ChN},fspecial('gaussian',2*filtStd,filtStd),'symmetric');
    %             FFImage{ChN}=imdivide(FFImage{ChN},double(max(FFImage{ChN}(:))));
    %         end
    %     else
    %         warning('No flat field found')
    %         for ChN=1:NChannels
    %             FFImage{ChN}=ones(ImageSize);
    %         end
    %     end
    for ChN=1:NChannels
        FFImage{ChN}=ones(ImageSize);
    end
    
    %Get the ellipses
    load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'], 'Ellipses')
    
    %This is the structuring element we will use to dilate the mask. 1.5um
    %might be a good choice.
    SESize=round(1.5/FrameInfo(1).PixelSize);
    StrucElement=strel('disk',SESize);
    
    h=waitbar(0,'Calculating the maximum projection');
    
    %     argin = tryGPU;
    
    for frame=1:length(Ellipses)
        waitbar(frame/length(Ellipses),h)
        
        %Create a mask of the ellipses
        
        %Initialize the image. In order to avoid problems with the drawing of
        %ellipses we'll make the image larger and then reduce it.
        IncreasePixels=200;
        %         Mask=zeros(ImageSize+IncreasePixels*2, argin{:});
        Mask=zeros(ImageSize+IncreasePixels*2);
        
        %How many ellipses do we have?
        [NEllipses,~]=size(Ellipses{frame});
        
        for j=1:NEllipses
            %Put ellipses at their corresponding positions. THe value of each
            %ellipse corresponds to its index number
            Mask=ellipseMatrix(round(Ellipses{frame}(j,2))+IncreasePixels,...
                round(Ellipses{frame}(j,1))+IncreasePixels,...
                Ellipses{frame}(j,3), Ellipses{frame}(j,4), Ellipses{frame}(j,5), Mask, 1);
        end
        
        %Go back to the original image size
        Mask=Mask(IncreasePixels+1:end-IncreasePixels,IncreasePixels+1:end-IncreasePixels);
        % Mask's underlying class needs to be logical or uint8 or else
        % morphopInputParser will throw an error.
        Mask=logical(Mask);
        Mask=imdilate(Mask,StrucElement);
        Mask=~Mask;
        
        
        for ChN=1:NChannels
            %Now, get the information for each Z slice. Note that I changed
            %this to account for the fact that we now have blank images at the
            %beginning and end of the stack
            
            %             for j=2:ZSlices-1
            %
            %                 FileName=[Prefix,'_',iIndex(frame,3),'_z',iIndex(j,2),...
            %                          '_ch',iIndex(ChN,2),'.tif'];
            %
            %
            % %                 if ~isempty(argin)
            % %                     tempMCPImage = gpuArray(imread([PreProcPath,filesep,Prefix,filesep,FileName]));
            % %                     tempMCPImage = double(gather(tempMCPImage));
            % %                 else
            % %                     tempMCPImage = double(imread([PreProcPath,filesep,Prefix,filesep,FileName]));
            % %                 end
            %                                 tempMCPImage = double(imread([PreProcPath,filesep,Prefix,filesep,FileName]));
            %
            %                 tempMCPImage = double(imread([PreProcPath,filesep,Prefix,filesep,FileName]));
            %                 MCPImage(:,:,j-1)=imdivide(tempMCPImage,FFImage{ChN});
            %
            %             end
            %
            MCPImage = getMovieFrame(liveExperiment, frame, liveExperiment.spotChannels(ChN));
            
            %Find the maximum and the mean
            MaxImageTemp=max(MCPImage,[],3);
            MeanImageTemp=mean(MCPImage,3);
            MedianImageTemp=median(MCPImage,3);
            
            MaxImage{ChN}(:,:,frame)=immultiply(MaxImageTemp,Mask);
            MeanImage{ChN}(:,:,frame)=immultiply(MeanImageTemp,Mask);
            MedianImage{ChN}(:,:,frame)=immultiply(MedianImageTemp,Mask);
        end
    end
    close(h)
    
    %Save to the processed path so that we don't overwhelm the Dropbox folder!
    save([ProcPath,filesep,Prefix,'_',filesep,'CytoImages.mat'],'MaxImage','MeanImage','MedianImage')
else
    disp('Using saved CytoImages.mat located in the PreProcessed folder.')
    
    load([ProcPath,filesep,Prefix,filesep,'CytoImages.mat']);
end


%Get the information
for ChN = 1:NChannels
    for frame=1:size(MeanImage,3)
        MeanImageTemp=MeanImage{ChN}(:,:,frame);
        MeanImageTemp=MeanImageTemp(:);
        MeanImageTemp=MeanImageTemp(MeanImageTemp~=0);
        
        Mean{ChN}(frame)=mean(MeanImageTemp);
        SD{ChN}(frame)=std(MeanImageTemp);
        Median{ChN}(frame)=median(MeanImageTemp);
        
        MaxImageTemp=MaxImage{ChN}(:,:,frame);
        MaxImageTemp=MaxImageTemp(:);
        MaxImageTemp=MaxImageTemp(MaxImageTemp~=0);
        Max{ChN}(frame)=mean(MaxImageTemp);
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

for frame=1:Rows
    for j=1:Cols
        Angle=atan2((frame-coordAZoom(2)),(j-coordAZoom(1)));
        Distance=sqrt((coordAZoom(2)-frame).^2+(coordAZoom(1)-j).^2);
        APPosition=Distance.*cos(Angle-APAngle);
        APPosImage(frame,j)=APPosition/APLength;
    end
end

APbinID=0:APResolution:1;


for ChN=1:NChannels
    MeanAPProfile{ChN}=nan(length(APbinID),size(MeanImage{ChN},3));
    SDAPProfile{ChN}=nan(length(APbinID),size(MeanImage{ChN},3));
    SEAPProfile{ChN}=nan(length(APbinID),size(MeanImage{ChN},3));
end


for ChN=1:NChannels
    for j=1:size(MeanImage{ChN},3)
        MeanImageTemp=MeanImage{ChN}(:,:,j);
        
        for frame=1:(length(APbinID)-1)
            FilteredMask=(APbinID(frame)<=APPosImage)&(APbinID(frame+1)>APPosImage);
            
            %Check that this filter mask region applies to the current data set
            if sum(sum(FilteredMask))
                FilteredMeanImage=immultiply(MeanImageTemp,FilteredMask);
                
                %Check that we will have enough statistics
                if sum(sum(FilteredMeanImage>0))>10
                    MeanAPProfile{ChN}(frame,j)=mean(FilteredMeanImage(FilteredMeanImage>0));
                    SDAPProfile{ChN}(frame,j)=std(FilteredMeanImage(FilteredMeanImage>0));
                    SEAPProfile{ChN}(frame,j)=std(FilteredMeanImage(FilteredMeanImage>0))/...
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
    
    figure((frame-1)*3+1)
    hold on
    for frame=1:size(MeanImage{ChN},3)
        plot(APbinID,MeanAPProfile{ChN}(:,frame),'-','color',ColorTime(frame,:))
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
    
    
    figure((frame-1)*3+2)
    hold on
    for frame=1:20:size(MeanImage{ChN},3)
        errorbar(APbinID,MeanAPProfile{ChN}(:,frame),SDAPProfile{ChN}(:,frame),'.-','color',ColorTime(frame,:))
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
    figure((frame-1)*3+3)
    errorbar(1:length(Mean{ChN}),Mean{ChN},SD{ChN},'.-k')
    xlabel('Frame')
    ylabel('Mean cyto fluorescence (\pm SD)')
    title(['Channel ',num2str(liveExperiment.spotChannels(ChN))])
    saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'CytoFluo',filesep,...
        'CytoFluoTime_ch',iIndex(liveExperiment.spotChannels(ChN),2),'.tif'])
    
end