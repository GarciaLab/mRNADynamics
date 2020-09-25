function [broken] = checkSchnitzCellErrors(schnitzcells, Ellipses)
broken = zeros(1, length(schnitzcells));
for i=1:length(schnitzcells)
    for j=1:length(schnitzcells(i).frames)
        f = schnitzcells(i).frames(j);
        aCell = schnitzcells(i).cellno(j);
        if aCell > 0
            x0 = round(Ellipses{f}(aCell, 1));
            y0 = round(Ellipses{f}(aCell, 2));
            if (x0 ~= schnitzcells(i).cenx(j)) | (y0 ~= schnitzcells(i).ceny(j))
                Distances =  sqrt((Ellipses{f}(:,1)-schnitzcells(i).cenx(j)).^2+(Ellipses{f}(:,2)-schnitzcells(i).ceny(j)).^2); 
                if find(Distances == min(Distances), 1) ~= schnitzcells(i).cellno(j)
                    broken(i) = j;
                end
                %schnitzcells(i).Flag = 10;
            end
            if Ellipses{f}(aCell, 9) ~= i
                broken(i) = -1;
            end
        else
            broken(i) = j;
        end
    end

end