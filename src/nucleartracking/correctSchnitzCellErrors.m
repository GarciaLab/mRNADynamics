function [valid_schnitzes, Ellipses] = correctSchnitzCellErrors(schnitzcells, Ellipses, Prefix)
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
            
            if (x0 ~= schnitzcells(i).cenx(j)) | (y0 ~= schnitzcells(i).ceny(j))
                Distances=sqrt((Ellipses{f}(:,1)-double(schnitzcells(i).cenx(j))).^2+(Ellipses{f}(:,2)-double(schnitzcells(i).ceny(j))).^2); 
                schnitzcells(i).cellno(j) = find(Distances == min(Distances), 1);
                cellno_found = true;

            else
                cellno_found = true;
            end
        else
            Distances=sqrt((Ellipses{f}(:,1)-schnitzcells(i).cenx(j)).^2+(Ellipses{f}(:,2)-schnitzcells(i).ceny(j)).^2); 
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
%% 

if exist('Prefix', 'var')
    lE= LiveExperiment(Prefix);
    FrameInfo = getFrameInfo(lE);
    ncs = [lE.nc9, lE.nc10, lE.nc11, lE.nc12, lE.nc13, lE.nc14, length(FrameInfo)];
    NCsForFrames = zeros(1, length(FrameInfo));
    for nc=9:14
        if ncs(nc-7) > 0
            if nc == 14
                NCsForFrames(ncs(nc-8):end) = nc;
            elseif ncs(nc-8) > 0
                NCsForFrames(ncs(nc-8):(ncs(nc-7)-1)) = nc;
            else
                NCsForFrames(1:(ncs(nc-7)-1)) = nc;
            end
        end
    end
    for i=1:length(valid_schnitzes)
        valid_schnitzes(i).cycle = floor(mode(NCsForFrames(valid_schnitzes(i).frames)));
        
    end
end








