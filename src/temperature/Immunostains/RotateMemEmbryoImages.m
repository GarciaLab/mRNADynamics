function RotateMemEmbryoImages(Prefix)
% Midpoints are in columns first (x) then rows (y) - y is smaller when
% higher up
liveExperiment = LiveExperiment(Prefix);
FrameInfo = getFrameInfo(liveExperiment);
load([liveExperiment.resultsFolder, filesep, 'MarkAndFindInfo.Mat'], 'MarkAndFindInfo');
PreProcFolder = liveExperiment.preFolder;


NEmbryos = MarkAndFindInfo.NSeries;
MembraneZoomResultsFolder = [liveExperiment.resultsFolder, filesep, 'ZoomMembraneInfo'];


MembraneZoomPixelPath = [MembraneZoomResultsFolder, filesep, 'MembraneZoomPixelSize.mat'];
load(MembraneZoomPixelPath,'PixelSize_um');

%%

MembraneZoomImageSizePath = [MembraneZoomResultsFolder, filesep, 'MembraneZoomImageSize.mat'];
load(MembraneZoomImageSizePath);
NEmbryos = MarkAndFindInfo.NSeries;

%%


RotatedMembraneMat = zeros(ySize/2, xSize, zSize, NEmbryos,'uint16');





CompiledEmbryoPath = [liveExperiment.resultsFolder, filesep, 'CompiledEmbryos.Mat'];
load(CompiledEmbryoPath);

%%


xSize_Temp = 2*xSize;
ySize_Temp = 2*ySize;



BlankImage  = zeros(ySize/2,xSize,'uint16');
for embryoIndex = 1:NEmbryos
    if CompiledEmbryos.Approved(embryoIndex) & CompiledEmbryos.MemCoordAs(embryoIndex,1) > 0 & CompiledEmbryos.MemCoordPs(embryoIndex,1) > 0 & CompiledEmbryos.MemCoordDs(embryoIndex,1) > 0
        
        NewName = [Prefix, '_Position',iIndex(embryoIndex,3),...
            '_ZoomMembrane', '.tif'];
        EmbryoImage = zeros(ySize, xSize, zSize, 'uint16');
        for k = 1:zSize
            EmbryoImage(:,:,k) = imread([PreProcFolder, filesep, NewName], k);
        end
        DisplayRange=[min(EmbryoImage(:)),max(EmbryoImage(:))];
        for zIndex = 1:zSize
            TempImage = zeros(ySize_Temp,xSize_Temp,'uint16');
            
            TempImage(uint16(ySize/2-CompiledEmbryos.MemyShift(embryoIndex)+1):uint16(ySize/2-CompiledEmbryos.MemyShift(embryoIndex)+ySize),...
                uint16(xSize/2+CompiledEmbryos.MemxShift(embryoIndex)+1):uint16(xSize/2+CompiledEmbryos.MemxShift(embryoIndex)+xSize))= EmbryoImage(:,:,zIndex);
            
            %         if ShowPlots
            %             TempFig = figure(2);
            %             set(TempFig,'units', 'normalized', 'position',[0.01, .1, .7, .7]);
            %
            %             TempAxes = axes(TempFig,'Units', 'normalized', 'Position', [0 0 1 1]);
            %             imTempImage = imshow(TempImage,DisplayRange,'Parent',TempAxes);
            %             axis off
            %         end
            
            TempImage2 = imrotate(TempImage, -CompiledEmbryos.MemAPRotationAngles(embryoIndex),'bilinear', 'crop');
            
            theta = -CompiledEmbryos.MemAPRotationAngles(embryoIndex)*pi/180;
            RotMat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
            
            
            %             if ShowPlots
            %                 TempFig2 = figure(3);
            %                 set(TempFig2,'units', 'normalized', 'position',[0.01, .1, .7, .7]);
            %
            %
            %                 TempAxes2 = axes(TempFig2,'Units', 'normalized', 'Position', [0 0 1 1]);
            %                 imTempImage2 = imshow(TempImage2,DisplayRange,'Parent',TempAxes2);
            %                 axis off
            %                 hold on
            %
            %                 hold off
            %             end
            %%
            TempImage3 =  TempImage2(round(3*xSize/4)+1:round(5*xSize/4),xSize/2+1:3*xSize/2);
            
            
            
            if CompiledEmbryos.FlippedOrientation(embryoIndex)
                TempImage3 = flipud(TempImage3);
                
            end
            
            %             if ShowPlots
            %                 TempFig3 = figure(4);
            %                 set(TempFig3,'units', 'normalized', 'position',[0.01, .1, .7, .7]);
            %
            %                 TempAxes3 = axes(TempFig3,'Units', 'normalized', 'Position', [0 0 1 1]);
            %                 imTempImage3 = imshow(TempImage3,DisplayRange,'Parent',TempAxes3);
            %                 axis off
            %
            %             end
            
            
            RotatedMembraneMat(:,:,zIndex, embryoIndex) = TempImage3;
            RotatedNewName = [Prefix, '_Position',iIndex(embryoIndex,3),...
                '_RotatedZoomMembrane', '.tif'];
            if zIndex == 1
                imwrite(squeeze(RotatedMembraneMat(:,:,zIndex, embryoIndex)), [PreProcFolder, filesep, RotatedNewName]);
            else
                imwrite(squeeze(RotatedMembraneMat(:,:,zIndex, embryoIndex)), [PreProcFolder, filesep, RotatedNewName],...
                    'WriteMode', 'append');
            end
        end
        
        
    else
        for zIndex = 1:zSize
            RotatedNewName = [Prefix, '_Position',iIndex(embryoIndex,3),...
                '_RotatedZoomMembrane', '.tif'];
            if zIndex == 1
                imwrite(BlankImage, [PreProcFolder, filesep, RotatedNewName]);
            else
                imwrite(BlankImage, [PreProcFolder, filesep, RotatedNewName],...
                    'WriteMode', 'append');
            end
        end
    end
    
end








