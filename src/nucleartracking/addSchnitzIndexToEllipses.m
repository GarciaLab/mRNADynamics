function [Ellipses, schnitzcells] = addSchnitzIndexToEllipses(Ellipses, schnitzcells)


disp('Ensuring consistency between Ellipses and schnitzcells...');

ellipsesOld = Ellipses;
schnitzcellsOld = schnitzcells;

%reset the cellno field of schnitzcells

% schnitzcells = rmfield(schnitzcells, 'cellno');


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
%                 schnitzcells(schnitzIndex).cellno(schnitzcells(...
%                     schnitzIndex).frames == frame) = uint16(ellipseIndex);
            else
                Ellipses{frame}(ellipseIndex, 9) = 0;
                
            end
        end
        
        %make sure no ellipses are assigned to the same schnitz in each
        %frame. multiple zeros are okay. 
        Ellipses{frame} = sortrows(Ellipses{frame}, 9);
        schnitzReferences = Ellipses{frame}(:, 9);
        nonZeroSchnitzReferences = schnitzReferences(schnitzReferences ~= 0); 
        uniqueReferences =  length(unique(nonZeroSchnitzReferences) ) == length(nonZeroSchnitzReferences);
        if ~uniqueReferences
%             Ellipses{frame}(hist(schnitzReferences,unique(schnitzReferences))>1, :) = []; %#ok<HIST>
%             warning('Overlapping ellipses found. Removing extra ellipses');
            error(['Overlapping ellipses found in frame ',num2str(frame),...
                'please correct and try again.'])
            %it's worth noting that there's a good chance that after all of
            %these procedures are done, schnitzcells.cellno references are
            %mixed up and not useful. that said, i don't think they were
            %useful to begin with. if you want useful cellnos in
            %schnitzcells, my advice is to append an extra loop after all
            %this work on Ellipses is done to fix the cellnos in
            %Schnitzcells (ie loop over schnitzcells and compare to
            %references in the ninth column of Ellipses). -AR 9/4/20
        end
                
    end
    
end


% %validation. this fails often, and that's upsetting. But it doesn't 
%lead to obvious problems down the road, so it's easier to ignore it for
%now.
% for s = 1:length(schnitzcells)
%     assert(length(schnitzcells(s).cellno) == length(schnitzcells(s).frames));
% end

% 
% ellipsesSizeUnchanged(ellipsesOld, Ellipses);
% schnitzcellsSizeUnchanged(schnitzcellsOld, schnitzcells); 

% disp('Schnitzcells and ellipses successfully inter-referenced.');


end