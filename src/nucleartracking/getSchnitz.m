function out = getSchnitz(ellipse, schnitzcells, frame, varargin)


tolerance = 1; %2 pixel tolerance for distance
out = [];

if ~isempty(varargin)
    cellno = varargin{1};
end

ellipseCentroid = double([ellipse(1), ellipse(2)]);
foundIt = false;
nSchnitz = length(schnitzcells);
schnitzIndex = 0;

while ~foundIt && schnitzIndex < nSchnitz
    
    schnitzIndex = schnitzIndex+1;
    schnitzFrameIndex = find(schnitzcells(schnitzIndex).frames == frame);
    
    %check for repeated frames.
    assert(length(schnitzFrameIndex) == 1 || isempty(schnitzFrameIndex));
    
    if ~isempty(schnitzFrameIndex)
        
        foundIt = abs(double(schnitzcells(schnitzIndex).ceny(schnitzFrameIndex)) - ellipseCentroid(2)) < tolerance...
            &...
            abs(double(schnitzcells(schnitzIndex).cenx(schnitzFrameIndex)) - ellipseCentroid(1)) < tolerance;
        
%         foundIt = double(schnitzcells(s).cellno) == double(cellno);
        
    end
    
    if foundIt
        out = uint16(schnitzIndex);
    end

end

end