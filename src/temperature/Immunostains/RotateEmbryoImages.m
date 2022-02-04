function RotateEmbryoImages(Prefix)
% Midpoints are in columns first (x) then rows (y) - y is smaller when
% higher up
ShowPlots = false;
liveExperiment = LiveExperiment(Prefix);
PixelSize_um = liveExperiment.pixelSize_um;
xSize = liveExperiment.xDim;
ySize = liveExperiment.yDim;
FrameInfo = getFrameInfo(liveExperiment);
load([liveExperiment.resultsFolder, filesep, 'MarkAndFindInfo.Mat'], 'MarkAndFindInfo');

NEmbryos = MarkAndFindInfo.NSeries;

%%
MembraneMat = getMembraneMat(liveExperiment);
HisMat = getHisMat(liveExperiment);


RotatedMembraneMat = zeros(ySize/2, xSize, NEmbryos,'double');
RotatedHisMat = zeros(ySize/2, xSize, NEmbryos,'double');

if isfile([liveExperiment.preFolder, filesep, Prefix, '-CustomHis.tif'])
    CustomHisMat = imreadStack2([liveExperiment.preFolder, filesep,...
        liveExperiment.Prefix, '-CustomHis.tif'], liveExperiment.yDim, liveExperiment.xDim, liveExperiment.nFrames);
    
    RotatedCustomHisMat = zeros(ySize/2, xSize, NEmbryos,'double');
end

if isfile([liveExperiment.preFolder, filesep, Prefix, '-CustomMembrane.tif'])
    CustomMembraneMat = imreadStack2([liveExperiment.preFolder, filesep,...
        liveExperiment.Prefix, '-CustomMembrane.tif'], liveExperiment.yDim, liveExperiment.xDim, liveExperiment.nFrames);
    
    RotatedCustomMembraneMat = zeros(ySize/2, xSize, NEmbryos,'double');
end


CompiledEmbryoPath = [liveExperiment.resultsFolder, filesep, 'CompiledEmbryos.Mat'];
load(CompiledEmbryoPath);

%%
xSize = size(MembraneMat,2);
ySize = size(MembraneMat,1);
NEmbryos = size(MembraneMat,3);

xSize_Temp = 2*xSize;
ySize_Temp = 2*ySize;


for embryoIndex = 1:NEmbryos
    if CompiledEmbryos.Approved(embryoIndex) & CompiledEmbryos.CoordAs(embryoIndex,1) > 0 & CompiledEmbryos.CoordPs(embryoIndex,1) > 0 & CompiledEmbryos.CoordDs(embryoIndex,1) > 0
        EmbryoImage = MembraneMat(:,:,embryoIndex);
        DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
        TempImage = zeros(ySize_Temp,xSize_Temp,'double');
        
        TempImage(uint16(ySize/2-CompiledEmbryos.yShift(embryoIndex)+1):uint16(ySize/2-CompiledEmbryos.yShift(embryoIndex)+ySize),...
            uint16(xSize/2+CompiledEmbryos.xShift(embryoIndex)+1):uint16(xSize/2+CompiledEmbryos.xShift(embryoIndex)+xSize))= EmbryoImage;
        
        if ShowPlots
            TempFig = figure(2);
            set(TempFig,'units', 'normalized', 'position',[0.01, .1, .7, .7]);
            
            TempAxes = axes(TempFig,'Units', 'normalized', 'Position', [0 0 1 1]);
            imTempImage = imshow(TempImage,DisplayRange,'Parent',TempAxes);
            axis off
        end
        
        TempImage2 = imrotate(TempImage, -CompiledEmbryos.APRotationAngles(embryoIndex),'bilinear', 'crop');
        
        theta = -CompiledEmbryos.APRotationAngles(embryoIndex)*pi/180;
        RotMat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
        
        
        if ShowPlots
            TempFig2 = figure(3);
            set(TempFig2,'units', 'normalized', 'position',[0.01, .1, .7, .7]);
            
            
            TempAxes2 = axes(TempFig2,'Units', 'normalized', 'Position', [0 0 1 1]);
            imTempImage2 = imshow(TempImage2,DisplayRange,'Parent',TempAxes2);
            axis off
            hold on
            
            hold off
        end
        %%
        TempImage3 =  TempImage2(round(3*xSize/4)+1:round(5*xSize/4),xSize/2+1:3*xSize/2);
        
        
        
        if CompiledEmbryos.FlippedOrientation(embryoIndex)
            TempImage3 = flipud(TempImage3);
            
        end
        
        if ShowPlots
            TempFig3 = figure(4);
            set(TempFig3,'units', 'normalized', 'position',[0.01, .1, .7, .7]);
            
            TempAxes3 = axes(TempFig3,'Units', 'normalized', 'Position', [0 0 1 1]);
            imTempImage3 = imshow(TempImage3,DisplayRange,'Parent',TempAxes3);
            axis off
            
        end
        
        
        RotatedMembraneMat(:,:,embryoIndex) = TempImage3;
        
        
        EmbryoImage = HisMat(:,:,embryoIndex);
        DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
        TempImage = zeros(ySize_Temp,xSize_Temp,'double');
        
        TempImage(uint16(ySize/2-CompiledEmbryos.yShift(embryoIndex)+1):uint16(ySize/2-CompiledEmbryos.yShift(embryoIndex)+ySize),...
            uint16(xSize/2+CompiledEmbryos.xShift(embryoIndex)+1):uint16(xSize/2+CompiledEmbryos.xShift(embryoIndex)+xSize))= EmbryoImage;
        
        
        TempImage2 = imrotate(TempImage, -CompiledEmbryos.APRotationAngles(embryoIndex),'bilinear', 'crop');
        
        theta = -CompiledEmbryos.APRotationAngles(embryoIndex)*pi/180;
        RotMat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
        
        
        %%
        TempImage3 =  TempImage2(round(3*xSize/4)+1:round(5*xSize/4),xSize/2+1:3*xSize/2);
        
        
        if CompiledEmbryos.FlippedOrientation(embryoIndex)
            TempImage3 = flipud(TempImage3);
            
        end
        
        
        RotatedHisMat(:,:,embryoIndex) = TempImage3;
        
        if exist('RotatedCustomHisMat', 'var')
            EmbryoImage = CustomHisMat(:,:,embryoIndex);
            DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
            TempImage = zeros(ySize_Temp,xSize_Temp,'double');
            
            TempImage(uint16(ySize/2-CompiledEmbryos.yShift(embryoIndex)+1):uint16(ySize/2-CompiledEmbryos.yShift(embryoIndex)+ySize),...
                uint16(xSize/2+CompiledEmbryos.xShift(embryoIndex)+1):uint16(xSize/2+CompiledEmbryos.xShift(embryoIndex)+xSize))= EmbryoImage;
            
            
            TempImage2 = imrotate(TempImage, -CompiledEmbryos.APRotationAngles(embryoIndex),'bilinear', 'crop');
            
            theta = -CompiledEmbryos.APRotationAngles(embryoIndex)*pi/180;
            RotMat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
            
            
            %%
            TempImage3 =  TempImage2(round(3*xSize/4)+1:round(5*xSize/4),xSize/2+1:3*xSize/2);
            
            
            if CompiledEmbryos.FlippedOrientation(embryoIndex)
                TempImage3 = flipud(TempImage3);
                
            end
            
            
            RotatedCustomHisMat(:,:,embryoIndex) = TempImage3;
        end
        
        if exist('RotatedCustomMembraneMat', 'var')
            EmbryoImage = CustomMembraneMat(:,:,embryoIndex);
            DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
            TempImage = zeros(ySize_Temp,xSize_Temp,'double');
            
            TempImage(uint16(ySize/2-CompiledEmbryos.yShift(embryoIndex)+1):uint16(ySize/2-CompiledEmbryos.yShift(embryoIndex)+ySize),...
                uint16(xSize/2+CompiledEmbryos.xShift(embryoIndex)+1):uint16(xSize/2+CompiledEmbryos.xShift(embryoIndex)+xSize))= EmbryoImage;
            
            
            TempImage2 = imrotate(TempImage, -CompiledEmbryos.APRotationAngles(embryoIndex),'bilinear', 'crop');
            
            theta = -CompiledEmbryos.APRotationAngles(embryoIndex)*pi/180;
            RotMat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
            
            
            %%
            TempImage3 =  TempImage2(round(3*xSize/4)+1:round(5*xSize/4),xSize/2+1:3*xSize/2);
            
            
            if CompiledEmbryos.FlippedOrientation(embryoIndex)
                TempImage3 = flipud(TempImage3);
                
            end
            
            
            RotatedCustomMembraneMat(:,:,embryoIndex) = TempImage3;
        end
    end
end




RotatedHisFile = [liveExperiment.preFolder, filesep, Prefix, '-His_Rotated.tif'];
RotatedMemFile = [liveExperiment.preFolder, filesep, Prefix, '-Membrane_Rotated.tif'];



saveNuclearProjection(RotatedHisMat,RotatedHisFile);
saveNuclearProjection(RotatedMembraneMat,RotatedMemFile);

if exist('RotatedCustomHisMat', 'var')
    RotatedCustomHisFile = [liveExperiment.preFolder, filesep,Prefix, '-CustomHis_Rotated.tif'];
    saveNuclearProjection(RotatedCustomHisMat,RotatedCustomHisFile);
end

if exist('RotatedCustomMembraneMat', 'var')
    RotatedCustomMembraneFile = [liveExperiment.preFolder, filesep,Prefix, '-CustomMembrane_Rotated.tif'];
    saveNuclearProjection(RotatedCustomMembraneMat,RotatedCustomMembraneFile);
end
