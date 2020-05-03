function [Ellipses, schnitzcells] = makeEllipsesFromSchnitzcells(schnitzcells, nFrames)

Ellipses = cell(nFrames, 1);
schnitzcellsFieldsToTranscribe = ["cenx"
    "ceny"
    "smaj"
    "smin"
    "orientationAngle"];


for frame = 1:nFrames
   
    ellipseFrame = [];
    
    for schnitzIndex = 1:length(schnitzcells)
        
        frameIndex = find(schnitzcells(schnitzIndex).frames == frame);
        
        if ~isempty(frameIndex)
            
            ellipseRow = [];
            
            %take care of Ellipses columns 1 to 5
            for k = 1:length(schnitzcellsFieldsToTranscribe)
                schnitzField = schnitzcells(schnitzIndex).(schnitzcellsFieldsToTranscribe{k});
                ellipseRow = [ellipseRow, schnitzField(frameIndex)];
            end
            
            %expand Ellipses to have columns 6, 7, 8, 9
            ellipseRow =[ellipseRow,...
                0, 0, 0, uint16(schnitzIndex)];
            
            %add this Ellipse to the current frame of Ellipses
            if isempty(ellipseFrame)
                ellipseFrame = ellipseRow;
            else
                ellipseFrame = [ellipseFrame;...
                    ellipseRow];
            end
            
            %add cellno field to make sure schnitzcells and Ellipses
            %correspond correctly
            schnitzcells(schnitzIndex).cellno(frameIndex) =...
                uint16(size(ellipseFrame, 1));
        end
    end
    
    Ellipses{frame} = ellipseFrame;
end

end