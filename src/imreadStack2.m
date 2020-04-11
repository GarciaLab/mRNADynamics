function stack = imreadStack2(file, yDim, xDim, nPages)

    if contains(file, '.tif')
    t = Tiff(file,'r');

    stack = zeros(yDim, xDim, nPages);


    for z = 1:nPages-1

        stack(:, :, z) = read(t); 
        nextDirectory(t);

    end
    close(t);
    else
        stack = load(file); 
        varName = fieldnames(stack);
        stack = stack.(varName{1});
    end

end