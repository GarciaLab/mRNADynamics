function val = getSchnitz(ellipse, schnitzcells, frame, varargin)


val = [];

if ~isempty(varargin)
    cellno = varargin{1};
end
center = round([ellipse(1), ellipse(2)]);
foundIt = false;
len = length(schnitzcells);
s = 0;

while ~foundIt && s < len
    
    s = s+1;
    schnitzFrameIndex = find(schnitzcells(s).frames == frame);
    
    %check for repeated frames.
    assert(length(schnitzFrameIndex) == 1 || isempty(schnitzFrameIndex));
    
    if ~isempty(schnitzFrameIndex)
        foundIt = round(schnitzcells(s).ceny(schnitzFrameIndex)) == center(2) &...
            round(schnitzcells(s).cenx(schnitzFrameIndex)) == center(1);
%         foundIt = double(schnitzcells(s).cellno) == double(cellno);
        
    end
    
    if foundIt
        val = uint16(s);
    end

end

end