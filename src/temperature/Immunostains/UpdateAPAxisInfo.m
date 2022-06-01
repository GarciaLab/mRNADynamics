function CompiledEmbryos = UpdateAPAxisInfo(Prefix, CompiledEmbryos)

liveExperiment = LiveExperiment(Prefix);
if isfield(CompiledEmbryos, 'MemCoordAs')
MembraneZoomResultsFolder = [liveExperiment.resultsFolder, filesep, 'ZoomMembraneInfo'];


MembraneZoomPixelPath = [MembraneZoomResultsFolder, filesep, 'MembraneZoomPixelSize.mat'];

load(MembraneZoomPixelPath,'PixelSize_um');
 MembraneZoomImageSizePath = [MembraneZoomResultsFolder, filesep, 'MembraneZoomImageSize.mat'];
 load(MembraneZoomImageSizePath);
 
MemPixelSize_um = PixelSize_um;
MemxSize = xSize;
MemySize = ySize;
MemzSize = zSize;
end

PixelSize_um = liveExperiment.pixelSize_um;

MembraneMat = getMembraneMat(liveExperiment);
xSize = size(MembraneMat, 2);
ySize = size(MembraneMat, 1);
NEmbryos = size(MembraneMat, 3);

    
try
    FrameInfo = getFrameInfo(liveExperiment);
catch
    FrameInfo = {};
end
if isfile([liveExperiment.resultsFolder, filesep, 'MarkAndFindInfo.Mat'])
    load([liveExperiment.resultsFolder, filesep, 'MarkAndFindInfo.Mat'], 'MarkAndFindInfo');
    NEmbryos = MarkAndFindInfo.NSeries;
else
    MarkAndFindInfo = {};
    NEmbryos = size(CompiledEmbryos.CoordAs, 1);
end



%%
CompiledEmbryos.APLengths = NaN(1, NEmbryos);
CompiledEmbryos.DVLengths = NaN(1, NEmbryos);
CompiledEmbryos.APRotationAngles = NaN(1, NEmbryos);
CompiledEmbryos.FlippedOrientation = zeros(1, NEmbryos,'logical');
CompiledEmbryos.APSlopes = NaN(1, NEmbryos);

CompiledEmbryos.APIntercepts = NaN(1, NEmbryos);
CompiledEmbryos.Midpoints = NaN(NEmbryos,2);
CompiledEmbryos.DVSlopes = NaN(1, NEmbryos);
CompiledEmbryos.DVIntercepts = NaN(1, NEmbryos);
CompiledEmbryos.Checked = zeros(1, NEmbryos,'logical');

CompiledEmbryos.xShift = NaN(1, NEmbryos);
CompiledEmbryos.yShift = NaN(1, NEmbryos);

for embryoIndex =1:NEmbryos
    if ~CompiledEmbryos.Approved(embryoIndex) | ((CompiledEmbryos.CoordAs(embryoIndex,1) > 0) & (CompiledEmbryos.CoordPs(embryoIndex,1) > 0)  & (CompiledEmbryos.CoordDs(embryoIndex,1) > 0))
        CompiledEmbryos.Checked(embryoIndex) = true;
    end
    if (CompiledEmbryos.CoordAs(embryoIndex,1) > 0) & (CompiledEmbryos.CoordPs(embryoIndex,1) > 0) & (CompiledEmbryos.CoordDs(embryoIndex,1) > 0)
        CompiledEmbryos.APLengths(embryoIndex) = sqrt((CompiledEmbryos.CoordAs(embryoIndex,1)-CompiledEmbryos.CoordPs(embryoIndex,1))^2+...
            (CompiledEmbryos.CoordAs(embryoIndex,2)-CompiledEmbryos.CoordPs(embryoIndex,2))^2)*PixelSize_um;
        

        CompiledEmbryos.APSlopes(embryoIndex) = (CompiledEmbryos.CoordAs(embryoIndex,2)-CompiledEmbryos.CoordPs(embryoIndex,2))/(CompiledEmbryos.CoordAs(embryoIndex,1)-CompiledEmbryos.CoordPs(embryoIndex,1));
        CompiledEmbryos.APIntercepts(embryoIndex) = CompiledEmbryos.CoordPs(embryoIndex,2)-CompiledEmbryos.APSlopes(embryoIndex) *CompiledEmbryos.CoordPs(embryoIndex,1);
        theta = 180/pi*atan(-(CompiledEmbryos.CoordAs(embryoIndex,2)-CompiledEmbryos.CoordPs(embryoIndex,2)) / (CompiledEmbryos.CoordAs(embryoIndex,1)-CompiledEmbryos.CoordPs(embryoIndex,1)));
        if CompiledEmbryos.CoordAs(embryoIndex,1) > CompiledEmbryos.CoordPs(embryoIndex,1)
            CompiledEmbryos.APRotationAngles(embryoIndex) = 180+theta;
            
            if CompiledEmbryos.CoordDs(embryoIndex,2) < CompiledEmbryos.APSlopes(embryoIndex)* CompiledEmbryos.CoordDs(embryoIndex,1)+CompiledEmbryos.APIntercepts(embryoIndex)
                CompiledEmbryos.FlippedOrientation(embryoIndex) = true;
            end
            
        elseif CompiledEmbryos.CoordAs(embryoIndex,1) < CompiledEmbryos.CoordPs(embryoIndex,1)
            CompiledEmbryos.APRotationAngles(embryoIndex) = mod(theta,360);
            
            if CompiledEmbryos.CoordDs(embryoIndex,2) > CompiledEmbryos.APSlopes(embryoIndex)* CompiledEmbryos.CoordDs(embryoIndex,1)+CompiledEmbryos.APIntercepts(embryoIndex)
                CompiledEmbryos.FlippedOrientation(embryoIndex) = true;
            end
            
        else
            if CompiledEmbryos.CoordAs(embryoIndex,2) > CompiledEmbryos.CoordPs(embryoIndex,2)
                CompiledEmbryos.APRotationAngles(embryoIndex) = 270;
                if CompiledEmbryos.CoordDs(embryoIndex,1) <  CompiledEmbryos.CoordAs(embryoIndex,1)
                    CompiledEmbryos.FlippedOrientation(embryoIndex) = true;
                end
            else
                CompiledEmbryos.APRotationAngles(embryoIndex) = 90;
                if CompiledEmbryos.CoordDs(embryoIndex,1) >  CompiledEmbryos.CoordAs(embryoIndex,1)
                    CompiledEmbryos.FlippedOrientation(embryoIndex) = true;
                end
            end
        end
        
        CompiledEmbryos.Midpoints(embryoIndex,:) = [(CompiledEmbryos.CoordAs(embryoIndex,1)+CompiledEmbryos.CoordPs(embryoIndex,1))/2, (CompiledEmbryos.CoordAs(embryoIndex,2)+CompiledEmbryos.CoordPs(embryoIndex,2))/2];
        CompiledEmbryos.DVSlopes(embryoIndex) = -1/CompiledEmbryos.APSlopes(embryoIndex) ;
        CompiledEmbryos.DVIntercepts(embryoIndex) = CompiledEmbryos.Midpoints(embryoIndex,2)-CompiledEmbryos.DVSlopes(embryoIndex)*CompiledEmbryos.Midpoints(embryoIndex,1);
        DVSlopeTemp = (CompiledEmbryos.CoordDs(embryoIndex,2)-CompiledEmbryos.CoordVs(embryoIndex,2))/(CompiledEmbryos.CoordDs(embryoIndex,1)-CompiledEmbryos.CoordVs(embryoIndex,1));
        DVInterceptTemp = CompiledEmbryos.CoordDs(embryoIndex,2)-DVSlopeTemp*CompiledEmbryos.CoordDs(embryoIndex,1);
        DVtheta = atan(DVSlopeTemp);
        APtheta = atan(CompiledEmbryos.APSlopes(embryoIndex));
        DeltaTheta = APtheta-DVtheta;

        CompiledEmbryos.DVLengths(embryoIndex) = sqrt((CompiledEmbryos.CoordDs(embryoIndex,1)-CompiledEmbryos.CoordVs(embryoIndex,1))^2+(CompiledEmbryos.CoordDs(embryoIndex,2)-CompiledEmbryos.CoordVs(embryoIndex,2))^2)*PixelSize_um*abs(sin(DeltaTheta));


    end
end

%% Add xShift & yShift Info
xSize_Temp = 2*xSize;
ySize_Temp = 2*ySize;


for embryoIndex = 1:NEmbryos
    if CompiledEmbryos.Checked(embryoIndex)
        EmbryoMidpoint = CompiledEmbryos.Midpoints(embryoIndex,:);
        EmbryoMidpoint_Temp = zeros(1,2,'double');
        EmbryoMidpoint_Temp(1) = uint16(round(EmbryoMidpoint(1)-xSize/2+xSize_Temp/2));
        EmbryoMidpoint_Temp(2) = uint16(round(EmbryoMidpoint(2)-ySize/2+ySize_Temp/2));
        CompiledEmbryos.xShift(embryoIndex) = round(xSize-EmbryoMidpoint_Temp(1));
        CompiledEmbryos.yShift(embryoIndex) = round(EmbryoMidpoint_Temp(2)-ySize);
    end
end

CompiledEmbryos.RotatedCoordAs = TransformToRotatedImageCoordinates(CompiledEmbryos.CoordAs, Prefix, CompiledEmbryos);
CompiledEmbryos.RotatedCoordPs = TransformToRotatedImageCoordinates(CompiledEmbryos.CoordPs, Prefix, CompiledEmbryos);
CompiledEmbryos.RotatedCoordDs = TransformToRotatedImageCoordinates(CompiledEmbryos.CoordDs, Prefix, CompiledEmbryos);
CompiledEmbryos.RotatedCoordVs = TransformToRotatedImageCoordinates(CompiledEmbryos.CoordVs, Prefix, CompiledEmbryos);

%%
if isfield(CompiledEmbryos, 'MemCoordAs')
CompiledEmbryos.MemAPLengths = NaN(1, NEmbryos);
CompiledEmbryos.MemDVLengths = NaN(1, NEmbryos);
CompiledEmbryos.MemAPRotationAngles = NaN(1, NEmbryos);
CompiledEmbryos.MemAPSlopes = NaN(1, NEmbryos);

CompiledEmbryos.MemAPIntercepts = NaN(1, NEmbryos);
CompiledEmbryos.MemMidpoints = NaN(NEmbryos,2);
CompiledEmbryos.MemDVSlopes = NaN(1, NEmbryos);
CompiledEmbryos.MemDVIntercepts = NaN(1, NEmbryos);


CompiledEmbryos.MemxShift = NaN(1, NEmbryos);
CompiledEmbryos.MemyShift = NaN(1, NEmbryos);

for embryoIndex =1:NEmbryos

    if (CompiledEmbryos.MemCoordAs(embryoIndex,1) > 0) & (CompiledEmbryos.MemCoordPs(embryoIndex,1) > 0) & (CompiledEmbryos.MemCoordDs(embryoIndex,1) > 0)
        
        CompiledEmbryos.MemAPLengths(embryoIndex) = sqrt((CompiledEmbryos.MemCoordAs(embryoIndex,1)-CompiledEmbryos.MemCoordPs(embryoIndex,1))^2+...
            (CompiledEmbryos.MemCoordAs(embryoIndex,2)-CompiledEmbryos.MemCoordPs(embryoIndex,2))^2)*MemPixelSize_um;
        

        CompiledEmbryos.MemAPSlopes(embryoIndex) = (CompiledEmbryos.MemCoordAs(embryoIndex,2)-CompiledEmbryos.MemCoordPs(embryoIndex,2))/(CompiledEmbryos.MemCoordAs(embryoIndex,1)-CompiledEmbryos.MemCoordPs(embryoIndex,1));
        CompiledEmbryos.MemAPIntercepts(embryoIndex) = CompiledEmbryos.MemCoordPs(embryoIndex,2)-CompiledEmbryos.MemAPSlopes(embryoIndex) *CompiledEmbryos.MemCoordPs(embryoIndex,1);
        theta = 180/pi*atan(-(CompiledEmbryos.MemCoordAs(embryoIndex,2)-CompiledEmbryos.MemCoordPs(embryoIndex,2)) / (CompiledEmbryos.MemCoordAs(embryoIndex,1)-CompiledEmbryos.MemCoordPs(embryoIndex,1)));
        MemFlippedO = false;
        if CompiledEmbryos.MemCoordAs(embryoIndex,1) > CompiledEmbryos.MemCoordPs(embryoIndex,1)
            if CompiledEmbryos.MemCoordDs(embryoIndex,2) < CompiledEmbryos.MemAPSlopes(embryoIndex)* CompiledEmbryos.MemCoordDs(embryoIndex,1)+CompiledEmbryos.MemAPIntercepts(embryoIndex)
                MemFlippedO = true;
            end
            
        elseif CompiledEmbryos.MemCoordAs(embryoIndex,1) < CompiledEmbryos.MemCoordPs(embryoIndex,1)
            if CompiledEmbryos.MemCoordDs(embryoIndex,2) > CompiledEmbryos.MemAPSlopes(embryoIndex)* CompiledEmbryos.MemCoordDs(embryoIndex,1)+CompiledEmbryos.MemAPIntercepts(embryoIndex)
                 MemFlippedO = true;
            end
            
        else
            if CompiledEmbryos.MemCoordAs(embryoIndex,2) > CompiledEmbryos.MemCoordPs(embryoIndex,2)
                CompiledEmbryos.MemAPRotationAngles(embryoIndex) = 270;
                if CompiledEmbryos.MemCoordDs(embryoIndex,1) <  CompiledEmbryos.MemCoordAs(embryoIndex,1)
                       MemFlippedO = true;
                end
            else
                CompiledEmbryos.MemAPRotationAngles(embryoIndex) = 90;
                if CompiledEmbryos.MemCoordDs(embryoIndex,1) >  CompiledEmbryos.MemCoordAs(embryoIndex,1)
                     MemFlippedO = true;
                end
            end
        end
        
        if (~CompiledEmbryos.FlippedOrientation(embryoIndex) & MemFlippedO) | (CompiledEmbryos.FlippedOrientation(embryoIndex) & ~MemFlippedO) 
            CoordDHold = CompiledEmbryos.MemCoordDs(embryoIndex,:);
            CoordVHold = CompiledEmbryos.MemCoordVs(embryoIndex,:);
            CompiledEmbryos.MemCoordDs(embryoIndex,:) = CoordVHold;
            CompiledEmbryos.MemCoordVs(embryoIndex,:) = CoordDHold;
            CompiledEmbryos.MemAPLengths(embryoIndex) = sqrt((CompiledEmbryos.MemCoordAs(embryoIndex,1)-CompiledEmbryos.MemCoordPs(embryoIndex,1))^2+...
                (CompiledEmbryos.MemCoordAs(embryoIndex,2)-CompiledEmbryos.MemCoordPs(embryoIndex,2))^2)*MemPixelSize_um;
            
            
            CompiledEmbryos.MemAPSlopes(embryoIndex) = (CompiledEmbryos.MemCoordAs(embryoIndex,2)-CompiledEmbryos.MemCoordPs(embryoIndex,2))/(CompiledEmbryos.MemCoordAs(embryoIndex,1)-CompiledEmbryos.MemCoordPs(embryoIndex,1));
            CompiledEmbryos.MemAPIntercepts(embryoIndex) = CompiledEmbryos.MemCoordPs(embryoIndex,2)-CompiledEmbryos.MemAPSlopes(embryoIndex) *CompiledEmbryos.MemCoordPs(embryoIndex,1);
            theta = 180/pi*atan(-(CompiledEmbryos.MemCoordAs(embryoIndex,2)-CompiledEmbryos.MemCoordPs(embryoIndex,2)) / (CompiledEmbryos.MemCoordAs(embryoIndex,1)-CompiledEmbryos.MemCoordPs(embryoIndex,1)));
        end
        
        
        
        if CompiledEmbryos.MemCoordAs(embryoIndex,1) > CompiledEmbryos.MemCoordPs(embryoIndex,1)
            CompiledEmbryos.MemAPRotationAngles(embryoIndex) = 180+theta;
            
        elseif CompiledEmbryos.MemCoordAs(embryoIndex,1) < CompiledEmbryos.MemCoordPs(embryoIndex,1)
            CompiledEmbryos.MemAPRotationAngles(embryoIndex) = mod(theta,360);
            
        else
            if CompiledEmbryos.MemCoordAs(embryoIndex,2) > CompiledEmbryos.MemCoordPs(embryoIndex,2)
                CompiledEmbryos.MemAPRotationAngles(embryoIndex) = 270;
            else
                CompiledEmbryos.MemAPRotationAngles(embryoIndex) = 90;
            end
        end
        
        CompiledEmbryos.MemMidpoints(embryoIndex,:) = [(CompiledEmbryos.MemCoordAs(embryoIndex,1)+CompiledEmbryos.MemCoordPs(embryoIndex,1))/2, (CompiledEmbryos.MemCoordAs(embryoIndex,2)+CompiledEmbryos.MemCoordPs(embryoIndex,2))/2];
        CompiledEmbryos.MemDVSlopes(embryoIndex) = -1/CompiledEmbryos.MemAPSlopes(embryoIndex) ;
        CompiledEmbryos.MemDVIntercepts(embryoIndex) = CompiledEmbryos.MemMidpoints(embryoIndex,2)-CompiledEmbryos.MemDVSlopes(embryoIndex)*CompiledEmbryos.MemMidpoints(embryoIndex,1);
        MemDVSlopeTemp = (CompiledEmbryos.MemCoordDs(embryoIndex,2)-CompiledEmbryos.MemCoordVs(embryoIndex,2))/(CompiledEmbryos.MemCoordDs(embryoIndex,1)-CompiledEmbryos.MemCoordVs(embryoIndex,1));
        MemDVInterceptTemp = CompiledEmbryos.MemCoordDs(embryoIndex,2)-MemDVSlopeTemp*CompiledEmbryos.MemCoordDs(embryoIndex,1);
        MemDVtheta = atan(MemDVSlopeTemp);
        MemAPtheta = atan(CompiledEmbryos.MemAPSlopes(embryoIndex));
        MemDeltaTheta = MemAPtheta-MemDVtheta;

        CompiledEmbryos.MemDVLengths(embryoIndex) = sqrt((CompiledEmbryos.MemCoordDs(embryoIndex,1)-CompiledEmbryos.MemCoordVs(embryoIndex,1))^2+(CompiledEmbryos.MemCoordDs(embryoIndex,2)-CompiledEmbryos.MemCoordVs(embryoIndex,2))^2)*MemPixelSize_um*abs(sin(MemDeltaTheta));


    end
end

%% Add xShift & yShift Info
MemxSize_Temp = 2*MemxSize;
MemySize_Temp = 2*MemySize;


for embryoIndex = 1:NEmbryos
    if CompiledEmbryos.Checked(embryoIndex)
        EmbryoMidpoint = CompiledEmbryos.MemMidpoints(embryoIndex,:);
        EmbryoMidpoint_Temp = zeros(1,2,'double');
        EmbryoMidpoint_Temp(1) = uint16(round(EmbryoMidpoint(1)-MemxSize/2+MemxSize_Temp/2));
        EmbryoMidpoint_Temp(2) = uint16(round(EmbryoMidpoint(2)-MemySize/2+MemySize_Temp/2));
        CompiledEmbryos.MemxShift(embryoIndex) = round(MemxSize-EmbryoMidpoint_Temp(1));
        CompiledEmbryos.MemyShift(embryoIndex) = round(EmbryoMidpoint_Temp(2)-MemySize);
    end
end

CompiledEmbryos.MemRotatedCoordAs = TransformToMemRotatedImageCoordinates(CompiledEmbryos.MemCoordAs, Prefix, CompiledEmbryos, MembraneMat);
CompiledEmbryos.MemRotatedCoordPs = TransformToMemRotatedImageCoordinates(CompiledEmbryos.MemCoordPs, Prefix, CompiledEmbryos, MembraneMat);
CompiledEmbryos.MemRotatedCoordDs = TransformToMemRotatedImageCoordinates(CompiledEmbryos.MemCoordDs, Prefix, CompiledEmbryos, MembraneMat);
CompiledEmbryos.MemRotatedCoordVs = TransformToMemRotatedImageCoordinates(CompiledEmbryos.MemCoordVs, Prefix, CompiledEmbryos, MembraneMat);

end
