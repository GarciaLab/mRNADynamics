function mask = makeMask(centers,radii,imageSize)


mask = zeros(imageSize);
for i = 1:length(radii)
    circle = [centers(i,:),radii(i)];
    [xx,yy] = ndgrid((1:imageSize(1))-circle(1),(1:imageSize(2))-circle(2));
    mask = mask + (xx.^2 + yy.^2)<circle(3)^2;
end
mask = uint16(mask);

end
