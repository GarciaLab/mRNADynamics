function ggiantIm = makeGiantImage(imIn, format, padSize,firstFrame, lastFrame, Prefix, channel)

%imDir is the folder where images are located
% format is like [numrows, numcolumns, num z slices]
%padSize should be twice filterSize

gpuDevice(1);

numFrames = (lastFrame - firstFrame) + 1;

dim = length(format);
frameInterval = padSize+format(1);

stacks = false;
if contains(imIn, lower('stacks'))
    stacks = true;
end

if stacks
    if dim == 2
        d = dir([imIn, '\*z*.tif']);
    elseif dim == 3
        d = dir([imIn, filesep,'*.tif']);
    end
end

nPads = numFrames;

if dim == 2
    ggiantIm = zeros(format(2), (format(1)*numFrames) + (padSize*nPads), 'gpuArray');
elseif dim == 3
    ggiantIm = zeros(format(2), (format(1)*numFrames) + (padSize*nPads),format(3)-2, 'single','gpuArray');
end

fcnt = 1;
for frame = firstFrame:lastFrame
    
    if stacks
        imPath = [imIn,filesep, d(frame).name];
        if dim == 2
            im = single(imread(imPath));
        elseif dim == 3
            im = single(readStack(imPath, format));
        end
    else
         if dim == 2
            im = single(imread(imIn));
        elseif dim == 3
            im = single(readPlanes(imIn, format, frame, Prefix, channel));
         end
    end
    
    gim = gpuArray(im);
    
    if dim == 2
        padIm = padarray(gim, [0,padSize], 0, 'post'); clear gim;
    else
        padIm = padarray(gim, [0,padSize, 0], 0, 'post'); clear gim;
    end
    
    ind1 = frameInterval*(fcnt-1) + 1;
    ind2 = frameInterval + (frameInterval*(fcnt-1));
    if dim == 2
        ggiantIm(:,ind1:ind2) = padIm; clear padIm;
    elseif dim == 3
        ggiantIm(:,ind1:ind2, :) = padIm; clear padIm;
    end
    
    fcnt = fcnt + 1;
end

% padIm = padarray(gim, [0,padSize], 0, 'post');
% ggiantIm = repmat(padIm, [1,numFrames]);

% giantIm = gather(ggiantIm); clear ggiantIm;

    function imStack = readStack(imPath, format)
        info = imfinfo(imPath);
        imStack = zeros(format(1), format(2), format(3)-2);
        for k = 2:format(3)-2
            imStack(:,:,k) = imread(imPath, k, 'Info', info);
        end
    end

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
