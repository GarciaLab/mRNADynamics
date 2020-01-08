function val = getSchnitz(ellipse, schnitzcells, frame)


val = [];

center = [ellipse(1), ellipse(2)];
foundIt = false;
len = length(schnitzcells);
s = 0;

while ~foundIt && s < len
    
    s = s+1;
    schnitzFrameIndex = find(schnitzcells(s).frames == frame);
    if ~isempty(schnitzFrameIndex)
        foundIt = schnitzcells(s).ceny(schnitzFrameIndex) == center(2) & schnitzcells(s).cenx(schnitzFrameIndex) == center(1);
    end
    
    if foundIt
        val = uint16(s);
    end

end

end