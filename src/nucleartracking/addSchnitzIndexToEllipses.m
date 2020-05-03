function [Ellipses, schnitzcells] = addSchnitzIndexToEllipses(Ellipses, schnitzcells)


disp('Ensuring consistency between Ellipses and schnitzcells...');

ellipsesOld = Ellipses;
schnitzcellsOld = schnitzcells;

for frame = 1:length(Ellipses)
    
    ellipseFrame = cell2mat(Ellipses(frame));
    if ~isempty(ellipseFrame)
        
        %reset the schnitzIndex column of Ellipses to be safe. 
        Ellipses{frame}(:, 9) = [];
        
        for ellipseIndex = 1:size(ellipseFrame, 1)
            
            ellipse = ellipseFrame(ellipseIndex, :);
            schnitzIndex = getSchnitz(ellipse, schnitzcells, frame, ellipseIndex);
            
            if ~isempty(schnitzIndex)
                
                assert(schnitzIndex <= length(schnitzcells));
                
                Ellipses{frame}(ellipseIndex, 9) = uint16(schnitzIndex);
                schnitzcells(schnitzIndex).cellno(schnitzcells(...
                    schnitzIndex).frames == frame) = uint16(ellipseIndex);
                
                
            end
        end
    end
    
end

ellipsesSizeUnchanged(ellipsesOld, Ellipses);
schnitzcellsSizeUnchanged(schnitzcellsOld, schnitzcells); 

disp('Schnitzcells and ellipses successfully inter-referenced.');


end