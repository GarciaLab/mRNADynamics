function stack = imreadStack2(file, yDim, xDim, nPages)

    if contains(file, '.tif')
        if false
            t = Tiff(file,'r');

            stack = zeros(yDim, xDim, nPages);


            for z = 0:nPages

                stack(:, :, z+1) = read(t);
                nextDirectory(t);

            end
            close(t);
        end
        
        stack = imreadStack(file);
        
        
    else
        stack = load(file);
        varName = fieldnames(stack);
        stack = stack.(varName{1});
    end
    

end