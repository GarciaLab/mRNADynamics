function livemRNAImageMatSaver(file, movieMat, fileTypeIn)

%by default, save the files as .mat. %save as
%tif plane or stack as preferred
if nargin==2, fileType = '.mat';
else  fileType = fileTypeIn; end

[precision, movieMat] = getImagePrecision(movieMat);

varStr = var2str(movieMat);


%% Saving 

switch fileType
    
    case '.mat'
        
        if whos(varStr).bytes < 2E9
            
            save(file, varStr, '-v6');
            
        else
            %in case the image is too big for v6, we'll use the clunky
            %v7.3 filetype. %it's inefficient out of the box, 
            %but we can configure it to chunk entire images into
            %contiguous blocks of memory. note that after
            %chunking, loading and saving are
            %still quiteinefficient relative
            %to v6 (but more compact). 
            imageDims = size(movieMat);
            chunkSize = [imageDims(1), imageDims(2),...
                ones(1, numel(imageDims)-2)];
            
            varMatic = newmatic(file ,true,...
                newmatic_variable(varStr, precision,...
                imageDims, chunkSize));
            
            varMatic.(varStr) = movieMat;
            
        end
        
    case '.tif'
                 
      error('Saving as .tif not supported yet.');
      
    otherwise
        
        error('Unrecognized file type.');
        
end


end