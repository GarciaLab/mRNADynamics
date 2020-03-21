function [Spots, dogs]...
    ...
    = segmentTranscriptionalLoci(...
    ...
    nCh, coatChannel, channelIndex, initialFrame, numFrames,...
    zSize, PreProcPath, Prefix, ProcPath, displayFigures,doFF, ffim,...
    Threshold, neighborhood, snippet_size, pixelSize, microscope,...
    Weka,filterMovieFlag, resultsFolder, gpu, saveAsMat, saveType, Ellipses)


cleanupObj = onCleanup(@myCleanupFun);

dogs = [];
DogOutputFolder=[ProcPath,filesep,'dogs',filesep];

%the underscore is to remove . and .. from the output structure
dogDir = dir([DogOutputFolder, '*_*']);

loadAsStacks = ~contains(dogDir(1).name, '_z');
Weka = startsWith(dogDir(1).name, 'prob');

waitbarFigure = waitbar(0, 'Segmenting spots');
set(waitbarFigure, 'units', 'normalized', 'position', [0.4, .15, .25,.1]);

Spots = repmat(struct('Fits', []), 1, numFrames);

if displayFigures
    %left bottom width height
    dogFig = figure('units', 'normalized', 'position',[.4, .5, .4, .4]);
    dogAx = axes(dogFig, 'Visible', 'off');
    gFig = figure('units', 'normalized', 'position',[0.01, .55, .33, .33]);
    gAx = axes(gFig);
    snipFig = figure('units', 'normalized', 'position',[0.4, .2, .1, .1]);
    snipAx = axes(snipFig);
    rawFig = figure('units', 'normalized', 'position',[.01, .1, .33, .33]);
    rawAx = axes(rawFig);
    
    graphicsHandles = [dogFig, dogAx, gFig, gAx, snipFig, snipAx, rawFig, rawAx];
else
    graphicsHandles = [];
    dogAx = 0;
end

if Weka
    MLFlag = 'ML';
    dogStr = 'prob';
    if Threshold < 5000
        warning('Increasing threshold to 5000. For Weka ML, you are thresholding on probability maps so the threshold shouldn''t be set below 50% = 5000.')
        Threshold = 5000;
    end
elseif loadAsStacks
    MLFlag = '';
    dogStr = 'dogStack_';
else
    MLFlag = '';
    dogStr = 'DOG_';
end


%Check how many coat channels we have and segment the appropriate channel
%accordingly

if filterMovieFlag
    filterType = 'Difference_of_Gaussian_3D';
    sigmas = {round(200/pixelSize),round(800/pixelSize)};
    filterOpts = {'nWorkers', 1, 'highPrecision', 'customFilter', filterType,...
        sigmas, 'double', 'keepPool', gpu};
    if saveAsMat
        filterOpts = [filterOpts, 'saveAsMat'];
    else
        filterOpts = [filterOpts, 'noSave'];
    end
    [~, dogs] = filterMovie(Prefix,'optionalResults', resultsFolder, filterOpts{:});
end

if nCh > 1
    ch = channelIndex;
else
    ch = coatChannel;
end

nameSuffix = ['_ch', iIndex(ch, 2)];

movieMat = loadMovieMat(Prefix);

yDim = size(movieMat, 1);
xDim = size(movieMat, 2);
nSlices = size(movieMat, 3);
nFrames = size(movieMat, 4);


% dogMat = loadDogMat(Prefix);

if Threshold == -1 && ~Weka
    
    
    if ~filterMovieFlag
        Threshold = determineThreshold(Prefix, ch,  'numFrames', numFrames);
        display(['Threshold: ', num2str(Threshold)])
    else
        Threshold = determineThreshold(Prefix, ch, 'noSave',  'numFrames', numFrames);
    end
    
    display(['Threshold: ', num2str(Threshold)])
    
end


q = parallel.pool.DataQueue;
afterEach(q, @nUpdateWaitbar);
p = 1;



zPadded = size(movieMat, 3) ~= zSize;


parfor currentFrame = initialFrame:numFrames 
    
    %report progress every tenth frame
    if ~mod(currentFrame, 10), disp(num2str(currentFrame)); end
    
    if loadAsStacks
        
        dogStackFile = [DogOutputFolder, filesep, dogStr, Prefix, '_', iIndex(currentFrame, 3),...
            nameSuffix];
        if exist([dogStackFile, '.mat'], 'file')
            
            dogStack = load([dogStackFile,'.mat'], 'dogStack');
            dogStack = dogStack.dogStack;
            
        elseif exist([dogStackFile, '.tif'], 'file')
            
            dogStack = imreadStack([dogStackFile, '.tif']);
            
        end
        
        
    end
    
    for zIndex = 1:zSize
        
        im = double(squeeze(movieMat(:, :, zIndex, currentFrame, ch)));
        try
            imAbove = double(sliceMovieMat(movieMat, ch, zIndex+1, currentFrame));
            imBelow= double(sliceMovieMat(movieMat, ch, zIndex-1, currentFrame));
        catch
            imAbove = nan(size(im,1),size(im,2));
            imBelow = nan(size(im,1),size(im,2));
        end
        
        
        if zPadded
            dogZ = zIndex;
        else
            dogZ = zIndex - 1;
        end
        
        if loadAsStacks
            dog = dogStack(:, :, dogZ);
        end
% =======
%             if isZPadded | ( ~isZPadded & (zIndex~=1 & zIndex~=zSize) )
%                 if strcmpi(saveType, '.tif')
%                     dogFileName = [DogOutputFolder, filesep, dogStr, Prefix, '_', iIndex(current_frame, 3), '_z', iIndex(dogZ, 2),...
%                         nameSuffix,'.tif'];
%                     dog = double(imread(dogFileName));
%                 elseif strcmpi(saveType, '.mat')
%                     dogFileName = [DogOutputFolder, filesep, dogStr, Prefix, '_', iIndex(current_frame, 3), '_z', iIndex(dogZ, 2),...
%                         nameSuffix,'.mat'];
%                     plane = load(dogFileName);
%                     try dog = plane.plane; catch dog = plane.dog; end
%                 elseif strcmpi(saveType, 'none')
%                     dog = dogs(:,:, dogZ, current_frame);
%                 end
%             else
%                 dog = false(size(im, 1), size(im, 2));
%             end
        
        
        if displayFigures
            dogO = im(:);
            lLim = median(dogO);
            uLim = max(dogO);
            if lLim ~= uLim
                imagescUpdate(dogAx, im, [lLim, uLim]);
            else
                imagescUpdate(dogAx, im, []);
            end
            drawnow;
        end
        
        % Apply flatfield correction
        if doFF && sum(size(im)==size(ffim))
            im = im.*ffim;
        end
        
        im_thresh = dog >= Threshold;
        
        % apply nuclear mask if it exists
        if ~isempty(Ellipses)
            ellipsesFrame = Ellipses{currentFrame};
            nuclearMask = makeNuclearMask(ellipsesFrame, [yDim, xDim]);
            %         immask = uint16(nuclearMask).*im;
            %         imshow(immask, [])
            im_thresh = im_thresh & nuclearMask;
        end
        
        if Weka
            se = strel('square', 3);
            im_thresh = imdilate(im_thresh, se); %thresholding from this classified probability map can produce non-contiguous, spurious Spots{channelIndex}. This fixes that and hopefully does not combine real Spots{channelIndex} from different nuclei
            im_thresh = im_thresh > 0;
        end
        
        [im_label, n_spots] = bwlabel(im_thresh);
        centroids = regionprops(im_thresh, 'centroid');
        
        temp_frames = {};
        temp_particles = cell(1, n_spots);
        
        if n_spots ~= 0
            
            for spotIndex = 1:n_spots
                centroid = round(centroids(spotIndex).Centroid);
                tic
                [temp_particles(spotIndex), Fits] = identifySingleSpot(spotIndex, {im,imAbove,imBelow}, im_label, dog, ...
                    neighborhood, snippet_size, pixelSize, displayFigures, graphicsHandles, microscope, 0, centroid,MLFlag, currentFrame, spotIndex, zIndex);
                Spots(currentFrame).Fits = [Spots(currentFrame).Fits, Fits];
            end
            
        end
        
    end
    
    send(q, currentFrame);
    
end


try close(waitbarFigure); end

    function nUpdateWaitbar(~)
        try waitbar(p/numFrames, waitbarFigure); end
        p = p + 1;
    end

end