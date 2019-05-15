% Function to apply selected filters to specified frame of a time series
function featureTable = filterImageFrames(Prefix,sigmaVec,featureCell,varargin)
featureTable = [];
for i = 1:numel(varargin)
end
[~,ProcPath,DropboxFolder,~, PreProcPath,...
    ~, Prefix, ~,~,~,~,...
    ~, spotChannels] = readMovieDatabase(Prefix);

if length(spotChannels) > 1
    error('nope. talk to nick');
end

% determine size of image stack
load([DropboxFolder, filesep, Prefix, filesep, 'FrameInfo.mat']);
zDim = FrameInfo(1).NumberSlices + 2;
yDim = FrameInfo(1).LinesPerFrame;
xDim = FrameInfo(1).PixelsPerLine;
pixelSize = FrameInfo(1).PixelSize*1000;
numFrames = length(FrameInfo);
zStep = FrameInfo(1).ZStep * 1000; %nm

format = [yDim, xDim, zDim];

noSave = true;
if noSave
    dogs = zeros(format(1), format(2), format(3)-2, length(FrameInfo));
end
sigmas = {round(210/pixelSize), floor(800/pixelSize)};
padSize = 2*sigmas{2};
pixVol = format(1)*format(2)*format(3);
maxGPUMem = 1.8E9;
%         maxGPUMem = evalin('base', 'maxGPUMem'); %for testing
maxPixVol = maxGPUMem / 4; %bytes in a single
chunkSize = floor(maxPixVol/pixVol);
chunks = [1:chunkSize:numFrames, numFrames+1];

% apply filters
for k = 1:length(chunks)-1
    g = makeGiantImage(PreProcPath, format, padSize, chunks(k), chunks(k+1)-1, Prefix, spotChannels);
    gt = permute(g, [2 1 3]);
    for i = 1:numel(featureCell)
        feature = featureCell{i};
        for j = 1:numel(sigmaVec)
            if strcmpi(feature,'Difference_of_Gaussian')
                gdog = filterImage(gt, feature, {sigmaVec{j}, sigmaVec{j}*4}, 'zStep', zStep);
            else
                gdog = filterImage(gt, feature, {sigmaVec{j}}, 'zStep', zStep);
            end
            gdogt = permute(gdog, [2 1 3]);
            imshow(gdogt(:,:,5),[]);
            dogs(:,:,:,chunks(k):chunks(k+1)-1) = gather(extractFromGiant(gdogt, format, padSize, chunks(k), chunks(k+1)-1,...
                Prefix, spotChannels, ProcPath, noSave));
            
            if i == 1 && j == 1
                featureTable = dogs(:);
            else
                featureTable = [featureTable dogs(:)];
            end
        end
    end
end

end
