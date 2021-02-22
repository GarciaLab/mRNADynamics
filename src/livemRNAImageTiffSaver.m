function livemRNAImageTiffSaver(Prefix, image, channelIn)

[precision, image] = getImagePrecision(image);

saveAsStack = true;

if strcmpi(precision, 'double')
    error('Error trying to save double precision file to Tif.')
end

liveExperiment = LiveExperiment(Prefix);

preFolder = liveExperiment.preFolder;

numDims = numel(size(image));
nCh = 1;

if numDims == 3
    nFrames = size(image, 3);
    channel = channelIn;
    error('Saving z-projections isn''t supported yet.')
elseif numDims == 4
    zSize = size(image, 3);
    nFrames = size(image, 4);
    channel = channelIn;
elseif numDims == 5
    zSize = size(image, 3);
    nFrames = size(image, 4);
    nCh = size(image, 5);
end

varStr = var2str(image);

%% Saving
for ch = 1:nCh
    
    for f = 1:nFrames
        
        if numDims == 4
            imageStack = image(:, :, :, f);
        elseif numDims ==5
            imageStack = image(:, :, :, f, channel(ch));
        end
        
        if ~saveAsStack
            
            for zIndex = 1:zSize
                singlePlaneFile= [preFolder, filesep,...
                    Prefix, '_', iIndex(f, 3), '_z', iIndex(zIndex, 2), ...
                    '_ch', iIndex(channel(ch), 2), '.tif'];
                imwrite(imageStack(:, :, zIndex), singlePlaneFile);
            end
            
        else
            
            stackFile = [preFolder, filesep, Prefix, '_',...
                iIndex(f, 3), '_ch', iIndex(channel(ch), 2), '.tif'];
            
            imwrite(imageStack(:, :, 1), stackFile);
            
            for z = 2:size(imageStack, 3)
                imwrite(imageStack(:, :, z),...
                    stackFile, 'WriteMode', 'append');
            end
        end
        
    end %frame loop
    
end %channel loop


end
