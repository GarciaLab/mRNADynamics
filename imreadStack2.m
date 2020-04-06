function stack = imreadStack2(file, yDim, xDim, nPages)

% tic
t = Tiff(file,'r');

stack = zeros(yDim, xDim, nPages);


parfor z = 1:nPages-1
    
    stack(:, :, z) = read(t); 
    nextDirectory(t);

end
close(t);
% toc

end