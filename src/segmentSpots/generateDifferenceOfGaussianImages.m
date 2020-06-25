% Generates difference of Gaussian images
function generateDifferenceOfGaussianImages(ProcPath, spotChannels,...
    numFrames, displayFigures, zSize, PreProcPath, Prefix, filterType, highPrecision,...
    sigmas, app, kernelSize, noSave, numType, gpu, saveAsMat, saveType)

liveExperiment = LiveExperiment(Prefix);

FrameInfo = getFrameInfo(liveExperiment);


gpu = 'noGPU';

DogOutputFolder = strrep([ProcPath, filesep, Prefix, '_', filesep, 'dogs'], '\\', '\');
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
if saveAsMat || strcmpi(saveType, '.mat')
    extractOpts = [extractOpts, 'mat'];
end


movieMat = getMovieMat(liveExperiment);

saveAsStacks = true;

dogStr = 'dogStack_';

nCh = length(spotChannels);


format = [FrameInfo(1).LinesPerFrame, FrameInfo(1).PixelsPerLine, zSize];

%2 spot 2 color will break this
dogMat = zeros(format(1),format(2), zSize, numFrames,  'uint16');

for ch = spotChannels
    
    nameSuffix = ['_ch', iIndex(ch, 2)];

    waitbarFigure = waitbar(0, ['Filtering images: Channel ', num2str(ch)]);
    
    %     nameSuffix = ['_ch', iIndex(ch, 2)];
    
    
    if displayFigures && isempty(app)
        filterFig = figure();
    end
    
    q = parallel.pool.DataQueue;
    afterEach(q, @nUpdateWaitbar);
    p = 1;
    
    
    if strcmpi(gpu, 'noGPU')
        for currentFrame = 1:numFrames
            
            if strcmpi(filterType, 'Difference_of_Gaussian')
                dogMat(:, :,:,currentFrame) = uint16((filterImage(double(movieMat(...
                    :, :, :, currentFrame, ch)), filterType, sigmas,...
                    'filterSize',filterSize, 'zStep', zStep) + 100) * 100);
            else
                dogMat(:, :,:,currentFrame)  = uint16((filterImage(double(...
                    movieMat(:, :, :, currentFrame, ch)), filterType, sigmas,...
                    'filterSize',filterSize) + 100) * 100);
            end
            
            if saveAsStacks
                dogStack = dogMat(:, :, :, currentFrame);
                dogStackFile = [DogOutputFolder, filesep, dogStr, Prefix, '_', iIndex(currentFrame, 3),...
                    nameSuffix,'.tif'];
                for z = 1:zSize
                    slice = dogStack(:, :, z);
                 if z == 1
                    imwrite(slice, dogStackFile);
                 else
                     imwrite(slice, dogStackFile,'WriteMode', 'append');
                 end
                end
            end
            
            send(q, currentFrame);
        end
        
    else
        
        format = [FrameInfo(1).LinesPerFrame, FrameInfo(1).PixelsPerLine, zSize];
        
        dogs = zeros(format(1), format(2), format(3)-2, numFrames);
        
        sigmas = {round(210/pixelSize), floor(800/pixelSize)};
        padSize = 2*sigmas{2};
        pixVol = format(1)*format(2)*format(3);
        if ~strcmpi(gpu, 'noGPU')
            maxMem = .8E9;
        else
            maxMem = 30E9;
        end
        %         maxGPUMem = evalin('base', 'maxGPUMem'); %for testing
        maxPixVol = maxMem / 4; %bytes in a single
        chunkSize = floor(maxPixVol/pixVol);
        chunks = [1:chunkSize:numFrames, numFrames+1];
        
        for i = 1:length(chunks)-1
            waitbar(chunks(i) / numFrames, waitbarFigure);
            g = makeGiantImage(Prefix, PreProcPath, format, padSize, chunks(i), chunks(i+1)-1, spotChannels, numType, gpu);
            g = permute(g, [2 1 3]);
            g = filterImage(g, filterType, sigmas, 'zStep', zStep, numType);
            g = permute(g, [2 1 3]);
            dogs(:,:,:,chunks(i):chunks(i+1)-1) = extractFromGiant(g, format, padSize, chunks(i), chunks(i+1)-1, Prefix, spotChannels, ProcPath, noSave, numType, extractOpts{:});
        end
        
        %             imshow(dogs(:,:, 5, 5),[]); %useful line for debugging
    end
    
    try close(waitbarFigure); catch; end
    
end

    function nUpdateWaitbar(~)
        try waitbar(p/numFrames, waitbarFigure);catch; end
        p = p + 1;
    end

end


function dogFrame = generateDoGs(DogOutputFolder,...
    PreProcPath, Prefix, currentFrame, nameSuffix, filterType,...
    sigmas, filterSize, highPrecision, zIndex, ...
    displayFigures, app, numFrames, saveType, im, zSize, zStep)


dim = 3;


if strcmpi(filterType, 'Difference_of_Gaussian')
    dogFrame = filterImage(im, filterType, sigmas,...
        'filterSize',filterSize, 'zStep', zStep);
    if highPrecision
        dogFrame = (dogFrame + 100) * 100;
    end
else
    dogFrame = (filterImage(im, filterType,...
        sigmas, 'filterSize',filterSize) + 100) * 100;
end


if displayFigures && dim == 2
    
    if ~ isempty(app)
        ax = app{1};
    else
        ax = gca;
    end
    if dim == 2
        imshow(dogFrame, [median(dogFrame(:)), max(dogFrame(:))], 'Parent', ax, 'InitialMagnification', 'fit');
    else
        imshow(dogFrame(:, :, zIndex), [median(dogFrame(:)), max(dogFrame(:))], 'Parent', ax, 'InitialMagnification', 'fit');
    end
    title(ax, [nameSuffix(2:end), ' frame: ', num2str(currentFrame), '/', num2str(numFrames), ' z: ', num2str(zIndex)], 'Interpreter', 'none')
    pause(.05)
end

end