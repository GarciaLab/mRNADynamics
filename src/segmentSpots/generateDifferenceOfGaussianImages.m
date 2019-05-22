% Generates difference of Gaussian images
function dogs = generateDifferenceOfGaussianImages(ProcPath, ExperimentType, FrameInfo, spotChannels,...
    numFrames, displayFigures, zSize, PreProcPath, Prefix, filterType, highPrecision,...
    sigmas, app, kernelSize, noSave, numType, gpu, saveAsMat, saveType)

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

extractOpts = {};
if saveAsMat | strcmpi(saveType, '.mat')
    extractOpts = [extractOpts, 'mat'];
end


nCh = length(spotChannels);
for channelIndex = 1:nCh
    
    waitbarFigure = waitbar(0, ['Filtering images: Channel ', num2str(channelIndex)]);
    
    % (MT, 2018-02-11) Added support for lattice imaging, maybe
    % temporary - FIX LATER
    
    if strcmpi(ExperimentType, 'inputoutput') || strcmpi(ExperimentType, 'lattice')
        nameSuffix = ['_ch', iIndex(spotChannels, 2)];
    else
        nameSuffix = ['_ch', iIndex(channelIndex, 2)];
    end
    
    if displayFigures && isempty(app)
        filterFig = figure();
    end
    
    q = parallel.pool.DataQueue;
    afterEach(q, @nUpdateWaitbar);
    p = 1;

    
    if ~filter3D
        parfor current_frame = 1:numFrames
                        
            
            for zIndex = 1:zSize
                generateDoGs(DogOutputFolder, PreProcPath, Prefix, current_frame, nameSuffix, filterType, sigmas, filterSize, ...
                    highPrecision, zIndex, displayFigures, app, numFrames);
            end
                    send(q, current_frame);

        end

    else
        
        format = [FrameInfo(1).LinesPerFrame, FrameInfo(1).PixelsPerLine, zSize];
        
        if noSave
            dogs = zeros(format(1), format(2), format(3)-2, numFrames);
        end
        sigmas = {round(210/pixelSize), floor(800/pixelSize)};
        padSize = 2*sigmas{2};
        pixVol = format(1)*format(2)*format(3);
%         if ~strcmpi(gpu, 'noGPU')
            maxMem = .8E9;
%         else
%             maxMem = 30E9;
%         end
        %         maxGPUMem = evalin('base', 'maxGPUMem'); %for testing
        maxPixVol = maxMem / 4; %bytes in a single
        chunkSize = floor(maxPixVol/pixVol);
        chunks = [1:chunkSize:numFrames, numFrames+1];
        
        for i = 1:length(chunks)-1
            waitbar(chunks(i) / numFrames, waitbarFigure);
            g = makeGiantImage(PreProcPath, format, padSize, chunks(i), chunks(i+1)-1, Prefix, spotChannels, numType, gpu);
            gt = permute(g, [2 1 3]);
            gdog = filterImage(gt, filterType, sigmas, 'zStep', zStep, numType);
            gdogt = permute(gdog, [2 1 3]);
            dogs(:,:,:,chunks(i):chunks(i+1)-1) = extractFromGiant(gdogt, format, padSize, chunks(i), chunks(i+1)-1, Prefix, spotChannels, ProcPath, noSave, numType, extractOpts{:});
        end
        
        %     imshow(dogs(:,:, 5, 5),[]);
    end
    
    close(waitbarFigure);
    
end

function nUpdateWaitbar(~)
    waitbar(p/numFrames, waitbarFigure);
    p = p + 1;
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

