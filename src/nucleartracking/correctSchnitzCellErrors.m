function [valid_schnitzes, Ellipses] = correctSchnitzCellErrors(schnitzcells, Ellipses)

num_frames = zeros(1, length(schnitzcells));

% secondtrigger = false;

for i=1:length(schnitzcells)
    
    num_frames(i) = length(schnitzcells(i).frames);
    
    for j=1:num_frames(i)
        
        cellno_found = false;
        f = schnitzcells(i).frames(j);
        aCell = schnitzcells(i).cellno(j);
        
        if aCell > 0
            
            x0 = round(Ellipses{f}(aCell, 1));
            y0 = round(Ellipses{f}(aCell, 2));
            
            if (x0 ~= schnitzcells(i).cenx(j)) || (y0 ~= schnitzcells(i).ceny(j))
                
                Distances=sqrt((Ellipses{f}(:,1) - double(schnitzcells(i).cenx(j))).^2 +...
                    (Ellipses{f}(:,2) - double(schnitzcells(i).ceny(j)).^2));
                schnitzcells(i).cellno(j) = find(Distances == min(Distances), 1);
                cellno_found = true;
                
            else
                cellno_found = true;
            end
            
        else
            
            Distances=sqrt((Ellipses{f}(:,1)- double(schnitzcells(i).cenx(j))).^2 +...
                (Ellipses{f}(:,2)- double(schnitzcells(i).ceny(j))).^2);
            schnitzcells(i).cellno(j) = find(Distances == min(Distances), 1);
            
            cellno_found = true;
            
        end
        
        
        if ~cellno_found
            schnitzcells(i).cellno(j) = NaN;
        else
            Ellipses{f}(schnitzcells(i).cellno(j), 9) = i;
        end
        
        
        %         if (i == 10) & (j == 10)
        %             disp('pause')
        %             secondtrigger = true;
        %         end
        %         if secondtrigger & (Ellipses{23}(schnitzcells(10).cellno(10), 9) ~= 10)
        %             disp('second pause')
        %         end
    end
    
end

valid_schnitzes = schnitzcells(num_frames > 0);