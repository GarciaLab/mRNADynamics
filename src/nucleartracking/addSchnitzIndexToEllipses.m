function [Ellipses, schnitzcells] = addSchnitzIndexToEllipses(Ellipses, schnitzcells)


disp('Ensuring consistency between Ellipses and schnitzcells...');

ellipsesOld = Ellipses;
schnitzcellsOld = schnitzcells;

for frame = 1:length(Ellipses)
    
    ellipseFrame = cell2mat(Ellipses(frame));
    if ~isempty(ellipseFrame)
        
        %reset the schnitzIndex column of Ellipses to be safe. 
        if size(Ellipses, 2) == 9
            Ellipses{frame}(:, 9) = [];
        end
        
        for ellipseIndex = 1:size(ellipseFrame, 1)
            
            ellipse = ellipseFrame(ellipseIndex, :);
            schnitzIndex = getSchnitz(ellipse, schnitzcells, frame, ellipseIndex);
            
            if ~isempty(schnitzIndex)
                
                assert(schnitzIndex <= length(schnitzcells));               
                Ellipses{frame}(ellipseIndex, 9) = uint16(schnitzIndex);
                schnitzcells(schnitzIndex).cellno(schnitzcells(...
                    schnitzIndex).frames == frame) = uint16(ellipseIndex);
            else
                Ellipses{frame}(ellipseIndex, 9) = 0;
                
            end
        end
        
        %make sure no ellipses are assigned to the same schnitz in each
        %frame. multiple zeros are okay. 
        schnitzReferences = Ellipses{frame}(:, 9);
        nonZeroSchnitzReferences = schnitzReferences(schnitzReferences ~= 0); 
        assert( length(unique(nonZeroSchnitzReferences) ) == length(nonZeroSchnitzReferences) );
        
    end
    
end

ellipsesSizeUnchanged(ellipsesOld, Ellipses);
schnitzcellsSizeUnchanged(schnitzcellsOld, schnitzcells); 

disp('Schnitzcells and ellipses successfully inter-referenced.');


end