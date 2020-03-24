% Generates difference of Gaussian images
function generateDifferenceOfGaussianImages(ProcPath, spotChannels,...
    numFrames, displayFigures, zSize, PreProcPath, Prefix, filterType, highPrecision,...
    sigmas, app, kernelSize, noSave, numType, gpu, saveAsMat, saveType)

thisExperiment = liveExperiment(Prefix);

FrameInfo = getFrameInfo(thisExperiment);


gpu = 'noGPU';

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
if saveAsMat || strcmpi(saveType, '.mat')
    extractOpts = [extractOpts, 'mat'];
end

% stacksPath = [PreProcPath, filesep, Prefix, filesep, 'stacks'];

movieMat = getMovieMat(thisExperiment);

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
            
            %             im = double(permute(squeeze(movieMat(ch,:, currentFrame,:,:)), [2,3,1]));
            
            %             dogMat(currentFrame, :,:,:) = generateDoGs(DogOutputFolder, PreProcPath, Prefix, currentFrame, nameSuffix, filterType,...
            %                 sigmas, filterSize, ...
            %                 highPrecision, -1, displayFigures, app, numFrames, saveType, im, zSize, zStep);
            %
            
            if strcmpi(filterType, 'Difference_of_Gaussian')
                dogMat(:, :,:,currentFrame) = uint16((filterImage(double(squeeze(movieMat(:, :, :, currentFrame, ch))), filterType, sigmas, 'filterSize',filterSize, 'zStep', zStep) + 100) * 100);
                %
                if highPrecision
                    %                         dogMat(currentFrame, :,:,:) = (dogMat(currentFrame, :,:,:) + 100) * 100;
                    %                     end
                end
            else
                dogMat(:, :,:,currentFrame)  = uint16((filterImage(double(squeeze(movieMat(:, :, :, currentFrame, ch))), filterType, sigmas, 'filterSize',filterSize) + 100) * 100);
            end
            
            if saveAsStacks
                dogStack = squeeze(dogMat(:, :, :, currentFrame));
                dogStackFile = [DogOutputFolder, filesep, dogStr, Prefix, '_', iIndex(currentFrame, 3),...
                    nameSuffix,'.mat'];
                save(dogStackFile, 'dogStack','-v6');
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
    
    try close(waitbarFigure); end
    
end

    function nUpdateWaitbar(~)
        try waitbar(p/numFrames, waitbarFigure); end
        p = p + 1;
    end

dogFile = [ProcPath, filesep, Prefix,'_', filesep, Prefix, '_dogMat.mat'];

% if whos(var2str(dogMat)).bytes < 2E9
%     save(dogFile, 'dogMat', '-v6');
% else
%     save(dogFile, 'dogMat', '-v7.3', '-nocompression');
% end

end


function dogFrame = generateDoGs(DogOutputFolder, PreProcPath, Prefix, currentFrame, nameSuffix, filterType,...
    sigmas, filterSize, highPrecision, zIndex, displayFigures, app, numFrames, saveType, im, zSize, zStep)


dim = 3;


if strcmpi(filterType, 'Difference_of_Gaussian')
    dogFrame = filterImage(im, filterType, sigmas, 'filterSize',filterSize, 'zStep', zStep);
    if highPrecision
        dogFrame = (dogFrame + 100) * 100;
    end
else
    dogFrame = (filterImage(im, filterType, sigmas, 'filterSize',filterSize) + 100) * 100;
end


% if dim == 2
%     dog = padarray(dog(filterSize:end - filterSize - 1, filterSize:end - filterSize - 1), [filterSize, filterSize], 0,'both');
%     dog_name = ['DOG_', Prefix, '_', iIndex(currentFrame, 3), '_z', iIndex(zIndex, 2), nameSuffix, saveType];
%     dog_full_path = [DogOutputFolder, filesep, dog_name];
%     if strcmpi(saveType, '.tif')
%         imwrite(uint16(dog), dog_full_path)
%     elseif strcmpi(saveType, '.mat')
%         save(dog_full_path, 'dog', '-v6');
%     end
% elseif dim == 3
%     
%     dog = cat(3, zeros(size(dog, 1), size(dog, 2)), dog);
%     dog(:, :, zSize) = zeros(size(dog, 1), size(dog, 2));
%     for z = 1:zSize
%         dog_name = ['DOG_', Prefix, '_', iIndex(currentFrame, 3), '_z', iIndex(z, 2), nameSuffix, saveType];
%         dog_full_path = [DogOutputFolder, filesep, dog_name];
%         if strcmpi(saveType, '.tif')
%             imwrite(uint16(dog(:,:, z)), dog_full_path);
%         elseif strcmpi(saveType, '.mat')
%             dog = dog(:,:,z);
%             save(dog_full_path, 'dog', '-v6');
%         end
%     end
% end

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