% Function to apply selected filters to specified frame of a time series
function pChunk = filterImageFrames(Prefix,sigmaVec,featureCell, beta, varargin)

beta = single(beta);

for i = 1:numel(varargin)
end
[~,ProcPath,DropboxFolder,~, PreProcPath,...
    ~, Prefix, ~,~,~,~,...
    ~, spotChannels] = readMovieDatabase(Prefix);

if length(spotChannels) > 1
    error('nope. talk to nick');
end

% determine size of image stack
load([DropboxFolder, filesep, Prefix, filesep, 'FrameInfo.mat'], 'FrameInfo');
zDim = FrameInfo(1).NumberSlices + 2;
yDim = FrameInfo(1).LinesPerFrame;
xDim = FrameInfo(1).PixelsPerLine;
pixelSize = FrameInfo(1).PixelSize*1000;
numFrames = length(FrameInfo);
zStep = FrameInfo(1).ZStep * 1000; %nm
% nFilters = numel(sigmaVec)*numel(featureCell);


format = [yDim, xDim, zDim];

noSave = false;
% if noSave
% %     dogs = zeros(format(1), format(2), format(3)-2, length(FrameInfo), 'single');
% %       dogs = zeros(format(1), format(2), format(3)-2, length(FrameInfo),'single', 'gpuArray');
% end
sigmas = {round(210/pixelSize), floor(800/pixelSize)};
padSize = 2*sigmas{2};
% padSize = 2*max(cell2mat(sigmaVec));
pixVol = format(1)*format(2)*format(3);
maxGPUMem = .5E9;
%         maxGPUMem = evalin('base', 'maxGPUMem'); %for testing
maxPixVol = maxGPUMem / 4; %bytes in a single
chunkSize = floor(maxPixVol/pixVol);
chunks = [1:chunkSize:numFrames, numFrames+1];

h = waitbar(0, 'Frame progress');
% dim1 = format(1);
% dim2 = (format(2)*numFrames) + (padSize*(numFrames-1));
% dim3 = format(3) - 2;
% dim4 = nFilters;
% featureTable = tall([dim1*dim2*dim3*dim4, 1]);

% pChunk = zeros(format(1), format(2), format(3)-2, length(FrameInfo),'single', 'gpuArray');
%
% cnt = 1;
% szes = 1;
mem = [];
tocs = [];
tic;
mfig = figure();
mAx = axes(mfig);
plotGPUMem;
% apply filters
for k = 1:length(chunks)-1
    h = waitbar(chunks(k) / numFrames, h);
    plotGPUMem;
    g = makeGiantImage(PreProcPath, format, padSize, chunks(k), chunks(k+1)-1, Prefix, spotChannels);
%     plotGPUMem;
    %     gt = permute(g, [2 1 3]);
    pChunk = zeros(size(g), 'like', g);
%     plotGPUMem;
    hf = waitbar(0, ['Filtering frames ', num2str(chunks(k)), ' to ', num2str(chunks(k+1)-1)]);
    m = 1;
    for i = 1:numel(featureCell)
        feature = featureCell{i};
        for j = 1:numel(sigmaVec)
            
            %             ind1 = sum(szes(cnt));
            %             cnt = cnt + 1;
            %
            waitbar((i*j)/(numel(sigmaVec)*numel(featureCell)), hf);
            if strcmpi(feature,'Difference_of_Gaussian')
                
                pChunk = pChunk + beta(m).*filterImage(g, feature, {sigmaVec{j}, sigmaVec{j}*4}, 'zStep', zStep);
                
%                 plotGPUMem;
            else
                %                 g = filterImage(g, feature, {sigmaVec{j}}, 'zStep', zStep);
%                 plotGPUMem;
                pChunk = pChunk + beta(m).*filterImage(g, feature, {sigmaVec{j}}, 'zStep', zStep);
            end
            %             gdogt = permute(gdog, [2 1 3]);
            %             imshow(gdogt(:,:,5),[]);
            %             dogs(:,:,:,chunks(k):chunks(k+1)-1) = extractFromGiant(gdogt, format, padSize, chunks(k), chunks(k+1)-1,...
            %                 Prefix, spotChannels, ProcPath, noSave);
            %             pChunk = sum(pChunk, beta(m).*extractFromGiant(gdogt, format, padSize, chunks(k), chunks(k+1)-1,...
            %                 Prefix, spotChannels, ProcPath, noSave));
            %              pChunk = pChunk + beta(m).*g;
            %             szes(cnt) = size(dogs(:));
            %             szes = [szes, format(1)*format(2)*(format(3)-2)*(chunks(k+1)-chunks(k))];
            %             ind2 = szes(cnt);
            %             featureTable(ind1:ind2) = dogs(:,:,:,chunks(k):chunks(k+1)-1);
            % %             featureTable(1:nFilters*(chunks(k+1)-1 - chunks(k))*format(1)*(format(3)-1)) = dogs(:);
            m = m + 1;
        end
    end
    
    %     szes = [szes, format(1)*format(2)*(format(3)-2)*(chunks(k+1)-chunks(k))];
    %             ind2 = szes(cnt);
    %             featureTable(ind1:ind2) = dogs(:,:,:,chunks(k):chunks(k+1)-1);
%     plotGPUMem;

    pChunk = exp(-pChunk)./ (1 +exp(-pChunk));
%     imshow(pChunk(:, 1:1000, 5), [nanmedian(pChunk(:)), max(pChunk(:))]);
%     plotGPUMem;
    %             probStack = reshape(pChunk(:,2),[yDim, xDim, zMax, chunks(k+1)-chunks(k)]);
    extractFromGiant(pChunk, format, padSize, chunks(k), chunks(k+1)-1,...
        Prefix, spotChannels, ProcPath, noSave, true);
%     plotGPUMem;
    close(hf);
    
end

    function plotGPUMem()
        memi = gpuDevice;
        mem = [mem, memi.AvailableMemory / (1E9)];
        tocs = [tocs, toc];
        plot(mAx, tocs,mem);
        ylabel(mAx, 'Available GPU memory (GB)');
        xlabel(mAx, 'time (s)');
        standardizeFigure(mAx, []);
        drawnow;
    end
    close(h);
end
