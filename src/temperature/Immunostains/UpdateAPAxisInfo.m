function CompiledEmbryos = UpdateAPAxisInfo(Prefix, CompiledEmbryos)

liveExperiment = LiveExperiment(Prefix);
PixelSize_um = liveExperiment.pixelSize_um;
xSize = liveExperiment.xDim;
ySize = liveExperiment.yDim;
FrameInfo = getFrameInfo(liveExperiment);
load([liveExperiment.resultsFolder, filesep, 'MarkAndFindInfo.Mat'], 'MarkAndFindInfo');

NEmbryos = MarkAndFindInfo.NSeries;

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

