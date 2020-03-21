function livemRNAImageMatSaver(file, movieMat, fileTypeIn)

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