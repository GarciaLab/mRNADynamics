function RotateRawDataImages(Prefix)
% Midpoints are in columns first (x) then rows (y) - y is smaller when
% higher up
ShowPlots = false;
liveExperiment = LiveExperiment(Prefix);
PixelSize_um = liveExperiment.pixelSize_um;
xSize = liveExperiment.xDim;
ySize = liveExperiment.yDim;
FrameInfo = getFrameInfo(liveExperiment);
nSlices=FrameInfo(1).NumberSlices;
load([liveExperiment.resultsFolder, filesep, 'MarkAndFindInfo.Mat'], 'MarkAndFindInfo');
CompiledEmbryoPath = [liveExperiment.resultsFolder, filesep, 'CompiledEmbryos.Mat'];
load(CompiledEmbryoPath);

FrameInfo = getFrameInfo(liveExperiment);
NEmbryos = MarkAndFindInfo.NSeries;


ProcPath = liveExperiment.userProcFolder;
DropboxFolder = liveExperiment.userResultsFolder;

Channel1 = liveExperiment.Channel1;
Channel2 = liveExperiment.Channel2;
Channel3 = liveExperiment.Channel3;
Channel4 = liveExperiment.Channel4;
Channel5 = liveExperiment.Channel5;


xSize = liveExperiment.xDim;
ySize = liveExperiment.yDim;

Channels = {Channel1, Channel2, Channel3, Channel4, Channel5};






close all
nEmbryos = length(FrameInfo);
%Get information about the image size
% HisImage=imread([PreProcPath,filesep,Prefix,filesep,D(1).name]);


%%


xSize_Temp = 2*xSize;
ySize_Temp = 2*ySize;




for embryoIndex = 1:nEmbryos
    if CompiledEmbryos.Approved(embryoIndex) & CompiledEmbryos.CoordAs(embryoIndex,1) > 0 & CompiledEmbryos.CoordPs(embryoIndex,1) > 0 & CompiledEmbryos.CoordDs(embryoIndex,1) > 0
        for channelIndex = 1:length(Channels)
            EmbryoChStackFile = [liveExperiment.preFolder, filesep, Prefix,...
                '_Position',iIndex(embryoIndex,3),'_001_ch', iIndex(channelIndex,2),'.tif'];
            
            
            EmbryoImage = imreadStack2(EmbryoChStackFile, liveExperiment.yDim, liveExperiment.xDim,...
                liveExperiment.nFrames);
            
            RotatedEmbryoMat = zeros(ySize/2, xSize,nSlices, 'double');
            
            
            DisplayRange=[min(EmbryoImage(:)),max(EmbryoImage(:))];
            TempImage = zeros(ySize_Temp,xSize_Temp,nSlices, 'double');
            
            TempImage(uint16(ySize/2-CompiledEmbryos.yShift(embryoIndex)+1):uint16(ySize/2-CompiledEmbryos.yShift(embryoIndex)+ySize),...
                uint16(xSize/2+CompiledEmbryos.xShift(embryoIndex)+1):uint16(xSize/2+CompiledEmbryos.xShift(embryoIndex)+xSize),:)= EmbryoImage;
            
            if ShowPlots
                TempFig = figure(2);
                set(TempFig,'units', 'normalized', 'position',[0.01, .1, .7, .7]);
                
                TempAxes = axes(TempFig,'Units', 'normalized', 'Position', [0 0 1 1]);
                imTempImage = imshow(TempImage(:,:,3),DisplayRange,'Parent',TempAxes);
                axis off
            end
            
            TempImage2 = imrotate(TempImage, -CompiledEmbryos.APRotationAngles(embryoIndex),'bilinear', 'crop');
            
            theta = -CompiledEmbryos.APRotationAngles(embryoIndex)*pi/180;
            RotMat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
            
            
            if ShowPlots
                TempFig2 = figure(3);
                set(TempFig2,'units', 'normalized', 'position',[0.01, .1, .7, .7]);
                
                
                TempAxes2 = axes(TempFig2,'Units', 'normalized', 'Position', [0 0 1 1]);
                imTempImage2 = imshow(TempImage2(:,:,3),DisplayRange,'Parent',TempAxes2);
                axis off
                hold on
                
                hold off
            end
            %%
            TempImage3 =  TempImage2(round(3*xSize/4)+1:round(5*xSize/4),xSize/2+1:3*xSize/2,:);
            
            
            
            if CompiledEmbryos.FlippedOrientation(embryoIndex)
                TempImage3 = flipud(TempImage3);
                
            end
            
            if ShowPlots
                TempFig3 = figure(4);
                set(TempFig3,'units', 'normalized', 'position',[0.01, .1, .7, .7]);
                
                TempAxes3 = axes(TempFig3,'Units', 'normalized', 'Position', [0 0 1 1]);
                imTempImage3 = imshow(TempImage3(:,:,3),DisplayRange,'Parent',TempAxes3);
                axis off
                
            end
            RotatedEmbryoChStackFile = [liveExperiment.preFolder, filesep, Prefix,...
                '_Rotated_Position',iIndex(embryoIndex,3),'_001_ch', iIndex(channelIndex,2),'.tif'];
            TempImage3 = uint16(TempImage3);
            for zindex = 1:nSlices
                if zindex == 1
                    imwrite(TempImage3(:,:,zindex), RotatedEmbryoChStackFile);
                else
                    imwrite(TempImage3(:,:,zindex), RotatedEmbryoChStackFile, 'WriteMode', 'append');
                end
            end
            
        end
    else
        BlankImage = zeros(ySize/2, xSize, nSlices, 'uint16');
        for channelIndex = 1:length(Channels)
            RotatedEmbryoChStackFile = [liveExperiment.preFolder, filesep, Prefix,...
                '_Rotated_Position',iIndex(embryoIndex,3),'_001_ch', iIndex(channelIndex,2),'.tif'];
            for zindex = 1:nSlices
                if zindex == 1
                    imwrite(BlankImage(:,:,zindex), RotatedEmbryoChStackFile);
                else
                    imwrite(BlankImage(:,:,zindex), RotatedEmbryoChStackFile, 'WriteMode', 'append');
                end
            end
        end
    end
end



