% Generates difference of Gaussian images
function [sigmas] = generateDifferenceOfGaussianImages(FISHPath, customFilter, nCh, ExperimentType, FrameInfo,  coatChannel, numFrames, displayFigures, zSize, PreProcPath, Prefix, filterType, highPrecision, sigmas, nWorkers, app, kernelSize)

DogOutputFolder = [FISHPath, filesep, Prefix, '_', filesep, 'dogs'];
mkdir(DogOutputFolder)

if contains(filterType, '_3D','IgnoreCase',true)
    filter3D = true;
    filterType = filterType(1:end-3);
else
    filter3D = false;
end

ps = parallel.Settings;
ps.Pool.AutoCreate = false;

if nWorkers > 1 && ~displayFigures
    maxWorkers = nWorkers;
    ps.Pool.AutoCreate = true;
    
    try
        parpool(maxWorkers);
    catch
        
        try
            parpool; % in case there aren't enough cores on the computer
        catch
            % parpool throws an error if there's a pool already running.
        end
        
    end
    
else
    
    try %#ok<TRYNC>
        poolobj = gcp('nocreate');
        delete(poolobj);
    end
    
end

clear rawdir;

pixelSize = FrameInfo(1).PixelSize * 1000; %nm
zStep = FrameInfo(1).ZStep*1000; %nm
close all;

if isempty(kernelSize)
    filterSize = round(2000 / pixelSize); %2000nm seems to be a good size empirically -AR
else
    filterSize = kernelSize;
end

if ~customFilter
    sigmas = {1, round(42000 / pixelSize)}; %42000nm seems to be a good size empirically -AR
end

if filter3D
    imStack = zeros(FrameInfo(1).LinesPerFrame,FrameInfo(1).PixelsPerLine, zSize-2);
else
    imStack = [];
end

for channelIndex = 1:nCh
    
    h = waitbar(0, ['Filtering images: Channel ', num2str(channelIndex)]);
    
    % (MT, 2018-02-11) Added support for lattice imaging, maybe
    % temporary - FIX LATER
    
    if strcmpi(ExperimentType, 'inputoutput') || strcmpi(ExperimentType, 'lattice')
        nameSuffix = ['_ch', iIndex(coatChannel, 2)];
    else
        nameSuffix = ['_ch', iIndex(channelIndex, 2)];
    end
    
    if displayFigures && isempty(app)
        filterFig = figure();
        filterAxes = axes(filterFig);
    end
    if ~filter3D
        parfor current_frame = 1:numFrames
            
            %         waitbar(current_frame / numFrames, h);
            
            for zIndex = 1:zSize
                generateDoGs(DogOutputFolder, PreProcPath, Prefix, current_frame, nameSuffix, filterType, sigmas, filterSize, ...
                    highPrecision, zIndex, displayFigures, app, numFrames);
            end
        end
        
    else
        for currentFrame = 1:numFrames
            for zIndex = 2:zSize-1
                imPath = [PreProcPath, filesep, Prefix, filesep, Prefix, '_', iIndex(currentFrame, 3), '_z', ...
                    iIndex(zIndex, 2), nameSuffix, '.tif'];
                imStack(:,:,zIndex-1) = double(imread(imPath));
            end
            generateDoGs(DogOutputFolder, PreProcPath, Prefix, currentFrame, nameSuffix, filterType, sigmas, filterSize, ...
                highPrecision, zIndex, displayFigures, app, numFrames, imStack, zSize, zStep);
        end
    end
    
    
    close(h);
    
end

end

function generateDoGs(DogOutputFolder, PreProcPath, Prefix, current_frame, nameSuffix, filterType, sigmas, filterSize, highPrecision, zIndex, displayFigures, app, numFrames, varargin)

if ~isempty(varargin)
    im = varargin{1};
    zSize = varargin{2};
    zStep = varargin{3};
    dim = 3;
else
    fileName = [PreProcPath, filesep, Prefix, filesep, Prefix, '_', iIndex(current_frame, 3), '_z', ...
        iIndex(zIndex, 2), nameSuffix, '.tif'];
    im = double(imread(fileName));
    dim = 2;
end

if sum(im(:)) ~= 0
    
    if strcmpi(filterType, 'Difference_of_Gaussian')
        dog = filterImage(im, filterType, sigmas, 'filterSize',filterSize, 'zStep', zStep);
        if highPrecision
            dog = (dog + 100) * 10;
        end
%         if dim == 3
%             kernelSize = 3*sigmas{2}+1;
% %             dog = padarray(dog(kernelSize:end - kernelSize - 1, kernelSize:end - kernelSize - 1, :), [kernelSize, kernelSize], 0, 'both');
%         end
    else
        dog = (filterImage(im, filterType, sigmas, 'filterSize',filterSize) + 100) * 10;
    end
    
else
    dog = im;
end

if dim == 2
    dog = padarray(dog(filterSize:end - filterSize - 1, filterSize:end - filterSize - 1), [filterSize, filterSize], 0,'both');
    dog_name = ['DOG_', Prefix, '_', iIndex(current_frame, 3), '_z', iIndex(zIndex, 2), nameSuffix, '.tif'];
    dog_full_path = [DogOutputFolder, filesep, dog_name];
    imwrite(uint16(dog), dog_full_path)
elseif dim == 3
    
    dog = cat(3, zeros(size(dog, 1), size(dog, 2)), dog);
    dog(:, :, zSize) = zeros(size(dog, 1), size(dog, 2));
    for z = 1:zSize
        dog_name = ['DOG_', Prefix, '_', iIndex(current_frame, 3), '_z', iIndex(z, 2), nameSuffix, '.tif'];
        dog_full_path = [DogOutputFolder, filesep, dog_name];
        imwrite(uint16(dog(:,:, z)), dog_full_path);
    end
end

if displayFigures && dim == 2
    
    if ~ isempty(app)
        ax = app{1};
    else
        ax = gca;
    end
    if dim == 2
        imshow(dog, [median(dog(:)), max(dog(:))], 'Parent', ax, 'InitialMagnification', 'fit');
    else
        imshow(dog(:, :, zIndex), [median(dog(:)), max(dog(:))], 'Parent', ax, 'InitialMagnification', 'fit');
    end
    title(ax, [nameSuffix(2:end), ' frame: ', num2str(current_frame), '/', num2str(numFrames), ' z: ', num2str(zIndex)], 'Interpreter', 'none')
    pause(.05)
end
end
