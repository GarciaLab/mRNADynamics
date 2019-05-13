% Generates difference of Gaussian images
function dogs = generateDifferenceOfGaussianImages(ProcPath, customFilter, nCh, ExperimentType, FrameInfo,  coatChannel, numFrames, displayFigures, zSize, PreProcPath, Prefix, filterType, highPrecision, sigmas, nWorkers, app, kernelSize, noSave)

dogs = [];

DogOutputFolder = [ProcPath, filesep, Prefix, '_', filesep, 'dogs'];
mkdir(DogOutputFolder)

if contains(filterType, '_3D','IgnoreCase',true)
    filter3D = true;
    filterType = filterType(1:end-3);
else
    filter3D = false;
end

pixelSize = FrameInfo(1).PixelSize * 1000; %nm
zStep = FrameInfo(1).ZStep*1000; %nm

if isempty(kernelSize)
    filterSize = round(2000 / pixelSize); %2000nm seems to be a good size empirically -AR
else
    filterSize = kernelSize;
end

if filter3D
    imStack = zeros(FrameInfo(1).LinesPerFrame,FrameInfo(1).PixelsPerLine, zSize-2);
else
    imStack = [];
end

for channelIndex = 1:nCh
    
    if ~filter3D
        h = waitbar(0, ['Filtering images: Channel ', num2str(channelIndex)]);
    else
        h = [];
    end
    
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
       for current_frame = 1:numFrames
            
            %             waitbar(current_frame / numFrames, h);
            
            for zIndex = 1:zSize
                generateDoGs(DogOutputFolder, PreProcPath, Prefix, current_frame, nameSuffix, filterType, sigmas, filterSize, ...
                    highPrecision, zIndex, displayFigures, app, numFrames);
            end
        end
        
    else
        
        %         for currentFrame = 1:numFrames
        %              waitbar(currentFrame / numFrames, h);
        %
        %             for zIndex = 2:zSize-1
        %                 imPath = [PreProcPath, filesep, Prefix, filesep, Prefix, '_', iIndex(currentFrame, 3), '_z', ...
        %                     iIndex(zIndex, 2), nameSuffix, '.tif'];
        %                 imStack(:,:,zIndex-1) = single(imread(imPath));
        %             end
        %             generateDoGs(DogOutputFolder, PreProcPath, Prefix, currentFrame, nameSuffix, filterType, sigmas, filterSize, ...
        %                 highPrecision, zIndex, displayFigures, app, numFrames, imStack, zSize, zStep);
        %         end
        format = [FrameInfo(1).LinesPerFrame, FrameInfo(1).PixelsPerLine, zSize];
        
        if noSave
            dogs = zeros(format(1), format(2), format(3)-2, numFrames);
        end
        sigmas = {round(210/pixelSize), floor(800/pixelSize)};
        padSize = 2*sigmas{2};
        pixVol = format(1)*format(2)*format(3);
        maxGPUMem = 1.8E9;
%         maxGPUMem = evalin('base', 'maxGPUMem'); %for testing
        maxPixVol = maxGPUMem / 4; %bytes in a single
        chunkSize = floor(maxPixVol/pixVol);
        chunks = [1:chunkSize:numFrames, numFrames+1];
        for i = 1:length(chunks)-1
            g = makeGiantImage(PreProcPath, format, padSize, chunks(i), chunks(i+1)-1, Prefix, coatChannel);
            gt = permute(g, [2 1 3]);
            gdog = filterImage(gt, filterType, sigmas, 'zStep', zStep);
            gdogt = permute(gdog, [2 1 3]);
            dogs(:,:,:,chunks(i):chunks(i+1)-1) = extractFromGiant(gdogt, format, padSize, chunks(i), chunks(i+1)-1, Prefix, coatChannel, ProcPath, noSave);
        end
    end
    
%     imshow(dogs(:,:, 5, 5),[]);
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
    zStep = NaN;
end

if sum(im(:)) ~= 0
    
    if strcmpi(filterType, 'Difference_of_Gaussian')
        dog = filterImage(im, filterType, sigmas, 'filterSize',filterSize, 'zStep', zStep);
        if highPrecision
            dog = (dog + 100) * 10;
        end
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
