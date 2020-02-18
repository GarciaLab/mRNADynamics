function Ellipses = addSchnitzIndexToEllipses(Ellipses, schnitzcells)

for frame = 1:length(Ellipses)
    
    ellipseFrame = cell2mat(Ellipses(frame));
    for ellipseIndex = 1:size(ellipseFrame, 1)
        ellipse = ellipseFrame(ellipseIndex, :);
        schnitz = getSchnitz(ellipse, schnitzcells, frame, ellipseIndex);
        if ~isempty(schnitz)
            Ellipses{frame}(ellipseIndex, 9) = schnitz;
        end
    end
    
end


end