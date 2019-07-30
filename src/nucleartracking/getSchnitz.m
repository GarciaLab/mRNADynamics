function val = getSchnitz(ellipse, schnitzcells, frame)


val = [];

center = [ellipse(1), ellipse(2)];
foundIt = false;
s = 0;
while ~foundIt & s < length(schnitzcells)
    s = s+1;
    schnitz = schnitzcells(s);
    schnitzFrameIndex = find(schnitz.frames == frame);
    if ~isempty(schnitzFrameIndex)
        foundIt = schnitz.ceny(schnitzFrameIndex) == center(2) & schnitz.cenx(schnitzFrameIndex) == center(1);
    end
    
    if foundIt
        val = uint16(s);
%         'hooray'
%     else
%         'uh oh empty'
    end

end

end