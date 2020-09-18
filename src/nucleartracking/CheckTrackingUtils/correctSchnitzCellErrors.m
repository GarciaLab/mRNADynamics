function valid_schnitzes = correctSchnitzCellErrors(schnitzcells, Ellipses)
num_frames = zeros(1, length(schnitzcells));
for i=1:length(schnitzcells)
    num_frames(i) = length(schnitzcells(i).frames);
    for j=1:num_frames(i)
        cellno_found = false;
        f = schnitzcells(i).frames(j);
        aCell = schnitzcells(i).cellno(j);
        x0 = round(Ellipses{f}(aCell, 1));
        y0 = round(Ellipses{f}(aCell, 2));
        if (x0 ~= schnitzcells(i).cenx(j)) | (y0 ~= schnitzcells(i).ceny(j))
            for k = 1:size(Ellipses{f}, 1)
                x02 = round(Ellipses{f}(k, 1));
                y02 = round(Ellipses{f}(k, 2));
                if (x02 == schnitzcells(i).cenx(j)) & (y02 == schnitzcells(i).ceny(j))
                    schnitzcells(i).cellno(j) = k;
                    cellno_found = true;
                    break
                end
                
            end
        else
            cellno_found = true;
        end
        if ~cellno_found
            schnitzcells(i).cellno(j) = NaN;
        end
    end

end

valid_schnitzes = schnitzcells(num_frames > 0);