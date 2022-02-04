function Coords = TransformToOriginalImageCoordinates(RotatedCoords, Prefix,CompiledEmbryos, liveExperiment, EmbryoLabel)

if ~exist('liveExperiment', 'var')
    liveExperiment = LiveExperiment(Prefix);
end
xSize = liveExperiment.xDim;
ySize = liveExperiment.yDim;

N = size(RotatedCoords, 1);

Coords = zeros(size(RotatedCoords),'double');


xSize_Temp = 2*xSize;
ySize_Temp = 2*ySize;
%%

if ~exist('EmbryoLabel', 'var')
    for embryoIndex = 1:N
        if RotatedCoords(embryoIndex,1) > 0
            
            Coord_Temp = RotatedCoords(embryoIndex,:);
            
            if CompiledEmbryos.FlippedOrientation(embryoIndex)
                Coord_Temp(2) = ySize/2-Coord_Temp(2)+1;
            end
            
            Coord_Temp2 = [Coord_Temp(1)+xSize/2 Coord_Temp(2)+3*ySize/4  ];
            
            Coord_Temp3 = [Coord_Temp2(1)-xSize_Temp/2-0.5 Coord_Temp2(2)-ySize_Temp/2-0.5 ];
            
            theta = CompiledEmbryos.APRotationAngles(embryoIndex)*pi/180;
            RotMat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
            Coord_Temp4 = (RotMat*Coord_Temp3.').';
            
            Coord_Temp5 = [Coord_Temp4(1)+xSize_Temp/2+0.5 Coord_Temp4(2)+ySize_Temp/2+0.5 ];
            
            
            
            Coords(embryoIndex,:)  = [Coord_Temp5(1)-xSize/2-CompiledEmbryos.xShift(embryoIndex) Coord_Temp5(2)-ySize/2+CompiledEmbryos.yShift(embryoIndex)];
            
            
            
        end
    end
    
else
    for coordIndex = 1:N
            Coord_Temp = RotatedCoords(coordIndex,:);
            
            if CompiledEmbryos.FlippedOrientation(EmbryoLabel)
                Coord_Temp(2) = ySize/2-Coord_Temp(2)+1;
            end
            
            Coord_Temp2 = [Coord_Temp(1)+xSize/2 Coord_Temp(2)+3*ySize/4  ];
            
            Coord_Temp3 = [Coord_Temp2(1)-xSize_Temp/2-0.5 Coord_Temp2(2)-ySize_Temp/2-0.5 ];
            
            theta = CompiledEmbryos.APRotationAngles(EmbryoLabel)*pi/180;
            RotMat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
            Coord_Temp4 = (RotMat*Coord_Temp3.').';
            
            Coord_Temp5 = [Coord_Temp4(1)+xSize_Temp/2+0.5 Coord_Temp4(2)+ySize_Temp/2+0.5 ];
            
            
            
            Coords(coordIndex,:)  = [Coord_Temp5(1)-xSize/2-CompiledEmbryos.xShift(EmbryoLabel) Coord_Temp5(2)-ySize/2+CompiledEmbryos.yShift(EmbryoLabel)];
            
            
            
        end
    end
end
