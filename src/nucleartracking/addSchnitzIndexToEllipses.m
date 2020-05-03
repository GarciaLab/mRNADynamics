function [Ellipses, schnitzcells] = addSchnitzIndexToEllipses(Ellipses, schnitzcells)

ellipsesOld = Ellipses;
schnitzcellsOld = schnitzcells;

for frame = 1:length(Ellipses)
    
    ellipseFrame = cell2mat(Ellipses(frame));
    if ~isempty(ellipseFrame)
        for ellipseIndex = 1:size(ellipseFrame, 1)
            ellipse = ellipseFrame(ellipseIndex, :);
            schnitzIndex = getSchnitz(ellipse, schnitzcells, frame, ellipseIndex);
            if ~isempty(schnitzIndex)
                Ellipses{frame}(ellipseIndex, 9) = uint16(schnitzIndex);
                schnitzcells(schnitzIndex).cellno(schnitzcells(...
                    schnitzIndex).frames == frame) = uint16(ellipseIndex);
            end
        end
    end
    
end

ellipsesSizeUnchanged(ellipsesOld, Ellipses);
schnitzcellsSizeUnchanged(schnitzcellsOld, schnitzcells); 

end