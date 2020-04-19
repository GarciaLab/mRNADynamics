function circleMask = filledVisCircles(im, centers, radii)

dim=[size(im, 1), size(im,2)];
circleMask=zeros(dim(1), dim(2));
for c = 1:length(radii)
    rad = radii(c); cenx = centers(c, 1); ceny = centers(c, 2);
    xrange = max(ceil(cenx-rad), 1) : min(ceil(cenx+rad),dim(1));
    yrange = max(ceil(ceny-rad), 1) : min(ceil(ceny+rad),dim(2));
    for x = xrange
        for y = yrange
            circleMask(y, x) = norm([x, y] - [cenx, ceny]) < rad;
        end
    end
    
end

% imshow(circleMask,[])

end