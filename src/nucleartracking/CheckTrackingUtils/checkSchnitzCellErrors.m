function [broken] = checkSchnitzCellErrors(schnitzcells, Ellipses)
broken = zeros(1, length(schnitzcells));
for i=1:length(schnitzcells)
    for j=1:length(schnitzcells(i).frames)
        f = schnitzcells(i).frames(j);
        aCell = schnitzcells(i).cellno(j);
        x0 = round(Ellipses{f}(aCell, 1));
        y0 = round(Ellipses{f}(aCell, 2));
        if (x0 ~= schnitzcells(i).cenx(j)) | (y0 ~= schnitzcells(i).ceny(j))
            broken(i) = j;
            %schnitzcells(i).Flag = 10;
        end
    end

end