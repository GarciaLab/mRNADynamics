function RotatedCoords = TransformToRotatedImageCoordinates(Coords, Prefix,CompiledEmbryos, MembraneMat)


liveExperiment = LiveExperiment(Prefix);
if exist('MembraneMat', 'var')
    xSize = size(MembraneMat, 2);
    ySize = size(MembraneMat, 1);
    NEmbryos = size(MembraneMat, 3);
else
    MembraneMat = getMembraneMat(liveExperiment);
    xSize = size(MembraneMat, 2);
    ySize = size(MembraneMat, 1);
    NEmbryos = size(MembraneMat, 3);
end
try 
    FrameInfo = getFrameInfo(liveExperiment);
catch
    FrameInfo = [];
end

N = size(Coords, 1);

RotatedCoords = zeros(size(Coords),'double');


xSize_Temp = 2*xSize;
ySize_Temp = 2*ySize;
%%

for embryoIndex = 1:N
    if Coords(embryoIndex,1) > 0

        Coord_Temp = [Coords(embryoIndex,1)+xSize/2+CompiledEmbryos.xShift(embryoIndex) Coords(embryoIndex,2)+ySize/2-CompiledEmbryos.yShift(embryoIndex)];
     

        Coord_Temp2 = [Coord_Temp(1)-xSize_Temp/2-0.5 Coord_Temp(2)-ySize_Temp/2-0.5 ];
        
        theta = -CompiledEmbryos.APRotationAngles(embryoIndex)*pi/180;
        RotMat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
        Coord_Temp3 = (RotMat*Coord_Temp2.').';

        Coord_Temp4 = [Coord_Temp3(1)+xSize_Temp/2+0.5 Coord_Temp3(2)+ySize_Temp/2+0.5 ];

        %%
        if ~isempty(FrameInfo)
            Coord_Temp5 = [Coord_Temp4(1)-xSize/2 Coord_Temp4(2)-3*ySize/4  ];
        else
            Coord_Temp5 = [Coord_Temp4(1)-xSize/4 Coord_Temp4(2)-3*ySize/4  ];
        end
        
        if CompiledEmbryos.FlippedOrientation(embryoIndex)
            Coord_Temp5(2) = ySize/2-Coord_Temp5(2)+1;
        end

        
        RotatedCoords(embryoIndex,:) = Coord_Temp5;

  
    end
end