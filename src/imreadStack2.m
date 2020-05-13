function stack = imreadStack2(file, yDim, xDim, nPages)

if contains(file, '.tif')
    if false
        t = Tiff(file,'r');
        
        stack = zeros(yDim, xDim, nPages);
        
        for z = 1:nPages - 1
            
            stack(:, :, z) = read(t);
            nextDirectory(t);
        end
        stack(:, :, nPages) = read(t);
        
        close(t);
        
    else
        
        stack = imreadStack(file);
        
    end
    
elseif contains(file, '.mat')
    
    stack = load(file);
    varName = fieldnames(stack);
    stack = stack.(varName{1});
    
end


end