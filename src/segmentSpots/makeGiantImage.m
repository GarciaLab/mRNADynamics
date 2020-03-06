function ggiantIm = makeGiantImage(imFile, format, padSize, frameRange, Prefix, channel, varargin)

%imDir is the folder where images are located
% format is like [numrows, numcolumns, num z slices]
%padSize should be twice filterSize

numType = 'single';
gpu = true;
processor = 'gpu';
movieMatCh = [];

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

firstFrame = frameRange(1);
lastFrame  = frameRange(end);

if gpu
    argin = tryGPU;
else
    argin = {};
    processor = 'cpu';
end

numFrames = (lastFrame - firstFrame) + 1;

dim = length(format);
frameInterval = padSize+format(2);

stacks = false;
if contains(imFile, lower('stacks'))
    stacks = true;
end

if stacks
    if dim == 2
        d = dir([imFile, '\*z*.tif']);
    elseif dim == 3
        d = dir([imFile, filesep,'*.tif']);
    end
end

nPads = numFrames;


if dim == 2
    ggiantIm = zeros(format(1), (format(2)*numFrames) + (padSize*nPads), numType, argin);
elseif dim == 3
    ggiantIm = zeros(format(1), (format(2)*numFrames) + (padSize*nPads),format(3)-2, numType,argin{:});
end

fcnt = 1;
for frame = firstFrame:lastFrame
    if isempty(movieMatCh)
        if stacks
            imPath = [imFile,filesep, d(frame).name];
            if dim == 2
                im = single(imread(imPath));
            elseif dim == 3
                im = single(readTiffStack(imPath));
            end
        else
             if dim == 2
                im = single(imread(imFile));
            elseif dim == 3
                im = single(readPlanes(imFile, format, frame, Prefix, channel));
             end
        end
    else
        im = movieMatCh(:, :, :, frame);
    end
    
%     gim = gpuArray(im);
if strcmpi(processor, 'gpu')
    im = gpuArray(im);
end

    if dim == 2
        padIm = padarray(gim, [0,padSize], 0, 'post');
    else
%         padIm = padarray(gim, [0,padSize, 0], 0, 'post');
    end
    
    ind1 = frameInterval*(fcnt-1) + 1;
    ind2 = frameInterval + (frameInterval*(fcnt-1));
    if dim == 2
%         ggiantIm(:,ind1:ind2) = padarray(gpuArray(im), [0,padSize, 0], 0, 'post');
    elseif dim == 3
        ggiantIm(:,ind1:ind2, :) = padarray(im, [0,padSize, 0], 0, 'post');
    end
    
    fcnt = fcnt + 1;
end

% padIm = padarray(gim, [0,padSize], 0, 'post');
% ggiantIm = repmat(padIm, [1,numFrames]);

% giantIm = gather(ggiantIm); clear ggiantIm;


    function imStack = readPlanes(imDir, format, frame, Prefix, channel)
       
        nameSuffix = ['_ch', iIndex(channel, 2)];
        imStack = zeros(format(1), format(2), format(3)-2);
        for z = 2:format(3)-2
            imPath = [imDir,filesep, Prefix, filesep,Prefix, '_', iIndex(frame, 3), '_z', ...
                iIndex(z, 2), nameSuffix, '.tif'];
            imStack(:,:,z) = imread(imPath);
        end
    end

end
