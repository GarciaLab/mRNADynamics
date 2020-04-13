function [Ellipses, schnitzcells] = makeEllipsesFromSchnitzcells(schnitzcells, nFrames)

Ellipses = cell(nFrames, 1);
names = string(fieldnames(schnitzcells));
names = names(names~='frames'...
& names~='Fluo'...
& names~='cellno');


for f = 1:nFrames
   
    ellipseFrame = [];
    
    for s = 1:length(schnitzcells)
        
        frameIndex = find(schnitzcells(s).frames == f);
        
        if ~isempty(frameIndex)
            
            ellipseRow = [];
            
            for k = 1:length(names)
                schnitzField = schnitzcells(s).(names{k});
                ellipseRow = [ellipseRow, schnitzField(frameIndex)];
            end
            
            ellipseRow =[ellipseRow,...
                0, 0, 0, s];
            
            if isempty(ellipseFrame)
                ellipseFrame = ellipseRow;
            else
                ellipseFrame = [ellipseFrame;...
                    ellipseRow];
            end
            
            schnitzcells(s).cellno(frameIndex) =...
                size(ellipseFrame, 1);
        end
    end
    
    Ellipses{f} = ellipseFrame;
end

end