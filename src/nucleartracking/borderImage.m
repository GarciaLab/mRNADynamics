function border = borderImage(bw)

xDim = size(bw, 2);
yDim = size(bw, 1);

border= false(yDim, xDim);

for x = 1:xDim
    for y = 1:yDim
        if x==1 || y==1 || x==xDim || y==yDim
            border(x,y) = true;
        end
    end
end

end