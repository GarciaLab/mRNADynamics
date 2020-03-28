function [Spots, dogs]...
    ...
    = segmentTranscriptionalLoci(...
    ...
    ~, ~, channelIndex, initialFrame, lastFrame,...
    zSize, ~, Prefix, ~, shouldDisplayFigures,doFF, ffim,...
    Threshold, neighborhood, snippet_size, ~, microscope,...
    ~,filterMovieFlag, resultsFolder, gpu, saveAsMat, ~, shouldMaskNuclei,...
    autoThresh)


cleanupObj = onCleanup(@myCleanupFun);


thisExperiment = liveExperiment(Prefix); 
pixelSize_nm = thisExperiment.pixelSize_nm;

if shouldMaskNuclei
    if thisExperiment.hasEllipsesFile, Ellipses = getEllipses(thisExperiment); end
else Ellipses = []; end

dogs = [];

DogOutputFolder = [thisExperiment.procFolder, 'dogs', filesep];

dogDir = dir([DogOutputFolder, '*_ch0', num2str(channelIndex), '.*']);

shouldLoadAsStacks = ~contains(dogDir(1).name, '_z');

isFileProbMap = startsWith(dogDir(1).name, 'prob');

waitbarFigure = waitbar(0, 'Segmenting spots');
set(waitbarFigure, 'units', 'normalized', 'position', [0.4, .15, .25,.1]);

Spots = repmat(struct('Fits', []), 1, lastFrame);

if shouldDisplayFigures
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

if isFileProbMap
    MLFlag = 'ML';
    dogStr = 'prob';
    if isnan(Threshold) || Threshold < 5000
        warning('Increasing threshold to 5000. For Weka ML, you are thresholding on probability maps so the threshold shouldn''t be set below 50% = 5000.')
        Threshold = 5000;
    end
else
        MLFlag = '';

    if shouldLoadAsStacks
        dogStr = 'dogStack_';
    else
        
        dogStr = 'DOG_';
    end
end


%Check how many coat channels we have and segment the appropriate channel
%accordingly

if filterMovieFlag
    filterType = 'Difference_of_Gaussian_3D';
    sigmas = {round(200/thisExperiment.pixelSize_nm),...
        round(800/thisExperiment.pixelSize_nm)};
    filterOpts = {'nWorkers', 1, 'highPrecision', 'customFilter', filterType,...
        sigmas, 'double', 'keepPool', gpu};
    if saveAsMat, filterOpts = [filterOpts, 'saveAsMat'];
        else, filterOpts = [filterOpts, 'noSave']; end
    [~, dogs] = filterMovie(Prefix,'optionalResults', resultsFolder, filterOpts{:});
end

nameSuffix = ['_ch', iIndex(channelIndex, 2)];

movieMat = getMovieMat(thisExperiment);

yDim = size(movieMat, 1);
xDim = size(movieMat, 2);

% dogMat = loadDogMat(Prefix);

    
    if autoThresh
        if ~filterMovieFlag
            Threshold = determineThreshold(Prefix, channelIndex,  'numFrames', lastFrame);
            display(['Threshold: ', num2str(Threshold)])
        else
            Threshold = determineThreshold(Prefix, channelIndex, 'noSave',  'numFrames', lastFrame);
        end
    end
    
    display(['Spot intensity threshold: ', num2str(Threshold)])
        
isZPadded = size(movieMat, 3) ~= zSize;

q = parallel.pool.DataQueue;
afterEach(q, @nUpdateWaitbar);
p = 1;
for currentFrame = initialFrame:lastFrame 
    
    %report progress every tenth frame
    if ~mod(currentFrame, 10), disp(['Segmenting frame ', num2str(currentFrame), '...']); end
    
    if shouldLoadAsStacks
        
        dogStackFile = [DogOutputFolder, dogStr, Prefix, '_', iIndex(currentFrame, 3),...
            nameSuffix];
        
        if exist([dogStackFile, '.mat'], 'file')
            dogStack = load([dogStackFile,'.mat'], 'dogStack');
            dogStack = dogStack.dogStack;
        elseif exist([dogStackFile, '.tif'], 'file')   
            dogStack = imreadStack([dogStackFile, '.tif']);
        end
        
        
    end
    
    for zIndex = 1:zSize
        
        im = double(squeeze(movieMat(:, :, zIndex, currentFrame, channelIndex)));
        try
            imAbove = double(sliceMovieMat(movieMat, channelIndex, zIndex+1, currentFrame));
            imBelow= double(sliceMovieMat(movieMat, channelIndex, zIndex-1, currentFrame));
        catch
            imAbove = nan(size(im,1),size(im,2));
            imBelow = nan(size(im,1),size(im,2));
        end
        
        
        if isZPadded, dogZ = zIndex;
        else, dogZ = zIndex - 1; end
        
        if shouldLoadAsStacks, dog = dogStack(:, :, dogZ); end
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
        
        
        if shouldDisplayFigures
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
        if shouldMaskNuclei && ~isempty(Ellipses)
            nuclearMask = makeNuclearMask(Ellipses{currentFrame}, [yDim, xDim]);
            im_thresh = im_thresh & nuclearMask;
        end
        
        %probability map regions usually look different from dog regions and
        %require some mophological finesse
        if isFileProbMap
            se = strel('square', 3);
            im_thresh = imdilate(im_thresh, se); %thresholding from this classified probability map can produce non-contiguous, spurious Spots{channelIndex}. This fixes that and hopefully does not combine real Spots{channelIndex} from different nuclei
            im_thresh = im_thresh > 0;
        end
        
        [im_label, n_spots] = bwlabel(im_thresh);
        centroids = regionprops(im_thresh, 'centroid');
        
        temp_particles = cell(1, n_spots);
        
        if n_spots ~= 0
            
            for spotIndex = 1:n_spots
                centroid = round(centroids(spotIndex).Centroid);
                
                [temp_particles(spotIndex), Fits] = identifySingleSpot(spotIndex,...
                    {im,imAbove,imBelow}, im_label, dog, ...
                    neighborhood, snippet_size, pixelSize_nm, shouldDisplayFigures,...
                    graphicsHandles, microscope, 0, centroid,MLFlag, currentFrame, spotIndex, zIndex);
                
                Spots(currentFrame).Fits = [Spots(currentFrame).Fits, Fits];
                
            end
            
        end
        
    end
    
    send(q, currentFrame);
    
end


try close(waitbarFigure); end

    function nUpdateWaitbar(~)
        try waitbar(p/lastFrame, waitbarFigure); end
        p = p + 1;
    end

end