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

radiusScale = 1.3; %dilate the nuclear mask by this factor

liveExperiment = LiveExperiment(Prefix);
pixelSize_nm = liveExperiment.pixelSize_nm;


dogs = [];

DogOutputFolder = [liveExperiment.procFolder, 'dogs', filesep];

dogDir = dir([DogOutputFolder, '*_ch0', num2str(channelIndex), '.*']);

haveProbs = any(cellfun(@(x) contains(x, 'prob'), {dogDir.name}));
%stacks won't be indexed by _z
haveStacks = any(cellfun(@(x) ~contains(x, '_z0'), {dogDir.name}));
% 
% waitbarFigure = waitbar(0, 'Segmenting spots');
% set(waitbarFigure, 'units', 'normalized', 'position', [0.4, .15, .25,.1]);

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
    maskFig = figure;
    graphicsHandles = [dogFig, dogAx, gFig, gAx, snipFig, snipAx, rawFig, rawAx];
else
    graphicsHandles = [];
    dogAx = 0;
end

if haveProbs
    MLFlag = 'ML';
    dogStr = 'prob';
    if isnan(Threshold) || isempty(Threshold)
        warning('Increasing threshold to 5000. For Weka ML, you are thresholding on probability maps so the threshold shouldn''t be set below 50% = 5000.')
        Threshold = 5000;
    end
else
    MLFlag = '';
    
    if haveStacks
        dogStr = 'dogStack_';
    else      
        dogStr = 'DOG_';
    end
end


%Check how many coat channels we have and segment the appropriate channel
%accordingly

if filterMovieFlag
    filterType = 'Difference_of_Gaussian_3D';
    sigmas = {round(200/liveExperiment.pixelSize_nm),...
        round(800/liveExperiment.pixelSize_nm)};
    filterOpts = {'nWorkers', 1, 'highPrecision',...
         'noSave', 'customFilter', filterType,...
        sigmas, 'double', 'keepPool', gpu};
  
    [~, dogs] = filterMovie(Prefix,'optionalResults', resultsFolder, filterOpts{:});
end

nameSuffix = ['_ch', iIndex(channelIndex, 2)];

movieMat = getMovieMat(liveExperiment);

movieMatCh = double(movieMat(:, :, :, :, channelIndex));

yDim = liveExperiment.yDim;
xDim = liveExperiment.xDim;
zDim = liveExperiment.zDim;

if shouldMaskNuclei
    if liveExperiment.hasEllipsesFile
        Ellipses = getEllipses(liveExperiment);
        if isempty(Ellipses)
            disp('Don''t try to do anything to Ellipses')
        else
            Ellipses = filterEllipses(Ellipses, [yDim, xDim]);
        end
    else
        shouldMaskNuclei = false;
        Ellipses = {}; % for parfor compatibility, see below.
    end
else
    Ellipses = {}; % for parfor compatibility, see below.
end

   
    if autoThresh
        if ~filterMovieFlag
            Threshold = determineThreshold(Prefix,...
                channelIndex,  'numFrames', lastFrame, 'firstFrame', initialFrame);
            display(['Threshold: ', num2str(Threshold)])
        else
            Threshold = determineThreshold(Prefix,...
                channelIndex, 'noSave', 'numFrames', lastFrame);
        end
    end
    
    display(['Spot intensity threshold: ', num2str(Threshold)])
        
isZPadded = size(movieMat, 3) ~= zSize;

% q = parallel.pool.DataQueue;
% afterEach(q, @nUpdateWaitbar);
% p = 1;
parfor currentFrame = initialFrame:lastFrame
    % neeed for parfor compatibility
    size(Ellipses);
    % this is just to send a hint to parfor to classify Ellipses variable correctly.
    % otherwise, it fails saying it doesnt exist or it doesnt have the
    % required size, even if not used by the worker.
    % See: https://la.mathworks.com/matlabcentral/answers/570619-matlab-parfor-index-exceeds-the-number-of-array-elements
    
    imStack = movieMatCh(:, :, :, currentFrame);
%     if shouldMaskNuclei
%         ellipseFrame = Ellipses{currentFrame};
%     end
    
    %report progress every tenth frame
    if ~mod(currentFrame, 10), disp(['Segmenting frame ',...
            num2str(currentFrame), '...']); end
    
    if haveStacks
        dogStackFile = [DogOutputFolder, dogStr, Prefix, '_',...
            iIndex(currentFrame, 3),...
            nameSuffix];
        
        if exist([dogStackFile, '.mat'], 'file')
            dogStack = load([dogStackFile,'.mat'], 'dogStack');
            dogStack = dogStack.dogStack;
        elseif exist([dogStackFile, '.tif'], 'file')
            dogStack = imreadStack2([dogStackFile, '.tif'], yDim, xDim, zDim);
        else
            error('Cannot find any dogs in the ProcessedData folder. Are they missing or incorrectly named?')
        end
        
        
    end
    
    for zIndex = 1:zSize
        
        im = imStack(:, :, zIndex);
        
        try
            imAbove = imStack(:, :, zIndex+1);
            imBelow= imStack(:, :, zIndex-1);
        catch
            imAbove = nan(size(im,1),size(im,2));
            imBelow = nan(size(im,1),size(im,2));
        end
        
        
        if isZPadded, dogZ = zIndex;
        else, dogZ = zIndex - 1; end
        
        if haveStacks, dog = dogStack(:, :, dogZ); end
      
        
        if shouldDisplayFigures
            dogO = im(:);
            lLim = median(dogO);
            uLim = max(dogO);
            if lLim ~= uLim
                imagescUpdate(dogAx, im, [lLim, uLim], 'cmap', 'gray');
            else
                imagescUpdate(dogAx, im, [], 'cmap', 'gray');
            end
            drawnow;
        end
        
        % Apply flatfield correction
        if doFF && sum(size(im)==size(ffim))
            im = double(im).*ffim; % GM 9/3/20
            %im = im.*ffim;
        end
        
        maxMaxDog = max(max(dog));
        maxDogOrderOfMagnitude = floor(log10(maxMaxDog));
        thresholdOrderOfMagnitude = floor(log10(Threshold));
        
        if thresholdOrderOfMagnitude < maxDogOrderOfMagnitude
            fprintf('WARNING: Threshold probably too low, resulting in out of memory errors.\n');
        elseif thresholdOrderOfMagnitude > maxDogOrderOfMagnitude
            fprintf('WARNING: Threshold probably too high, nothing will be detected.\n');
        end
        
        dog_thresh = dog >= Threshold;
        
         % apply nuclear mask if it exists
         if shouldMaskNuclei

            if liveExperiment.hasCustomMasks
                nuclearMask = liveExperiment.getNuclearMask(currentFrame, zIndex);
            else    
                nuclearMask = makeNuclearMask(Ellipses{currentFrame}, [yDim xDim], radiusScale);
            end

            dog_thresh = dog_thresh & nuclearMask;
             
 %             if shouldDisplayFigures
 %                 figure(maskFig);
 %                 imshowpair(nuclearMask, dog, 'montage'); 
 %             end
             
         end
        
        %probability map regions usually look different from dog regions and
        %require some mophological finesse
        if haveProbs
            se = strel('square', 3);
            dog_thresh = imdilate(dog_thresh, se); %thresholding from this classified probability map can produce non-contiguous, spurious Spots{channelIndex}. This fixes that and hopefully does not combine real Spots{channelIndex} from different nuclei
            dog_thresh = dog_thresh > 0;
        end
        
        [dog_label, n_spots] = bwlabel(dog_thresh);
        centroids = regionprops(dog_thresh, 'centroid');
        
        if n_spots ~= 0
            
            for spotIndex = 1:n_spots
                centroid = round(centroids(spotIndex).Centroid);
                
                 Fits = identifySingleSpot(spotIndex,...
                    {im,imAbove,imBelow}, dog_label, dog, ...
                    neighborhood, snippet_size, pixelSize_nm, shouldDisplayFigures,...
                    graphicsHandles, microscope, 0,...
                    centroid,MLFlag, currentFrame, spotIndex, zIndex);
                
                Spots(currentFrame).Fits = [Spots(currentFrame).Fits, Fits];
                
            end
        end
    end
%     send(q, currentFrame);
end

% 
% try close(waitbarFigure); catch; end
% 
%     function nUpdateWaitbar(~)
%         try waitbar(p/lastFrame, waitbarFigure); catch; end
%         p = p + 1;
%     end

end