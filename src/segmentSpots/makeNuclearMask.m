function nuclearMask = makeNuclearMask(ellipsesFrame, dim)

%make a mask from the ellipses structure. this is currently used for
%segmenting loci in segmentSpots.


radScale = 1.3;

nuclearMask = false(dim(1), dim (2));

for e = 1:size(ellipsesFrame, 1)
    
    ceny = ellipsesFrame(e, 1);
    cenx = ellipsesFrame(e, 2);
    rad = ellipsesFrame(e,3)*radScale;
    xrange = max(ceil(cenx-rad), 1) : min(ceil(cenx+rad),dim(1));
    yrange = max(ceil(ceny-rad), 1) : min(ceil(ceny+rad),dim(2));
    for x = xrange
        for y = yrange
            nuclearMask(x, y) = norm([x, y] - [cenx, ceny]) < rad;
        end
    end
end




end