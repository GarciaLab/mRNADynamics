function val = getSchnitz(ellipse, schnitzcells, frame)


val = [];

center = [ellipse(1), ellipse(2)];
foundIt = false;
s = 0;
while ~foundIt & s <= length(schnitzcells)
    s = s+1;
    schnitz = schnitzcells(s);
    schnitzInd = find(schnitz.frames == frame);
    if ~isempty(schnitzInd)
        foundIt = schnitz.ceny(schnitzInd) == center(2) & schnitz.cenx(schnitzInd) == center(1);
    end
    
    if foundIt
%         Ellipses{frame}(ellipseIndex, 9) = s;
        val = s;
    end

end

end