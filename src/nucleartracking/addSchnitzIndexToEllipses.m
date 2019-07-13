function Ellipses = addSchnitzIndexToEllipses(Ellipses, schnitzcells)

for frame = 1:length(Ellipses)
    
    ellipseFrame = cell2mat(Ellipses(frame));
    for ellipseIndex = 1:size(ellipseFrame, 1)
        ellipse = ellipseFrame(ellipseIndex, :);
        Ellipses{frame}(ellipseIndex, 9) = getSchnitz(ellipse, schnitzcells, frame);
    end
    
end


end