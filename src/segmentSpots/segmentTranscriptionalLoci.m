function Spots...
    ...
    = segmentTranscriptionalLoci(...
    ...
    ch_quantify, initialFrame, lastFrame,...
    zSize, ~, Prefix, ~, shouldDisplayFigures,doFF, ffim,...
    Threshold, neighborhood, snippet_size, ~, microscope,...
    ~,filterMovieFlag, resultsFolder, gpu, saveAsMat, ~, shouldMaskNuclei,...
    autoThresh, ch_segment)


cleanupObj = onCleanup(@myCleanupFun);

radiusScale = 1.3; %dilate the nuclear mask by this factor

liveExperiment = LiveExperiment(Prefix);
pixelSize_nm = liveExperiment.pixelSize_nm;


dogs = [];

DogOutputFolder = [liveExperiment.procFolder, 'dogs', filesep];

dogDir = dir([DogOutputFolder, '*_ch0', num2str(ch_segment), '.*']);

haveProbs = any(cellfun(@(x) contains(x, 'prob'), {dogDir.name}));
%stacks won't be indexed by _z
haveStacks = any(cellfun(@(x) ~contains(x, '_z'), {dogDir.name}));

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

movieMat = getMovieMat(liveExperiment);
if ~isempty(movieMat)
    movieMat_channel = movieMat(:, :, :, :, ch_quantify);
else
    movieMat_channel = [];
end

yDim = liveExperiment.yDim;
xDim = liveExperiment.xDim;
zDim = liveExperiment.zDim;

% if shouldMaskNuclei
%     if liveExperiment.hasEllipsesFile
%         Ellipses = getEllipses(liveExperiment);
%         if isempty(Ellipses)
%             disp('Don''t try to do anything to Ellipses')
%         else
%             Ellipses = filterEllipses(Ellipses, [yDim, xDim]);
%         end
%     else
%         shouldMaskNuclei = false;
%     end
% end

   
    if autoThresh
        if ~filterMovieFlag
            Threshold = determineThreshold(Prefix,...
                ch_segment,  'numFrames', lastFrame, 'firstFrame', initialFrame);
            display(['Threshold: ', num2str(Threshold)])
        else
            Threshold = determineThreshold(Prefix,...
                ch_segment, 'noSave', 'numFrames', lastFrame);
        end
    end
    
    display(['Spot intensity threshold: ', num2str(Threshold)])
        
isZPadded = size(movieMat, 3) ~= zSize;

q = parallel.pool.DataQueue;
afterEach(q, @nUpdateWaitbar);
p = 1;
parfor currentFrame = initialFrame:lastFrame
    if ~isempty(movieMat_channel)
        imStack = movieMat_channel(:, :, :, currentFrame);
    else
        imStack = getMovieFrame(liveExperiment, currentFrame, ch_quantify);
    end
    
   
    %report progress every tenth frame
    if ~mod(currentFrame, 10), disp(['Segmenting frame ',...
            num2str(currentFrame), '...']); end
    
    if haveStacks
        
        dogStackFile = [DogOutputFolder, dogStr, Prefix, '_',...
            iIndex(currentFrame, 3),...
            '_ch', iIndex(ch_segment, 2)];
        
        if exist([dogStackFile, '.mat'], 'file')
            dogStack = load([dogStackFile,'.mat'], 'dogStack');
            dogStack = dogStack.dogStack;
        elseif exist([dogStackFile, '.tif'], 'file')
            dogStack = imreadStack2([dogStackFile, '.tif'], yDim,...
                xDim, zDim);
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
        
        maxMaxDog = max(max(dog))
        maxDogOrderOfMagnitude = floor(log10(maxMaxDog))
        thresholdOrderOfMagnitude = floor(log10(Threshold))
        
        if thresholdOrderOfMagnitude < maxDogOrderOfMagnitude
            disp('WARNING: Threshold probably too low, resulting in out of memory errors.');
        elseif thresholdOrderOfMagnitude > maxDogOrderOfMagnitude
            disp('WARNING: Threshold probably too high, nothing will be detected.');
        end
        
        
        
        im_thresh = dog >= Threshold;
        
        % apply nuclear mask if it exists
%         if shouldMaskNuclei
%     
%             nuclearMask = makeNuclearMask(Ellipses{currentFrame}, [yDim xDim], radiusScale);
%             im_thresh = im_thresh & nuclearMask;
%             
% %             if shouldDisplayFigures
% %                 figure(maskFig);
% %                 imshowpair(nuclearMask, dog, 'montage'); 
% %             end
%             
%         end
        
        %probability map regions usually look different from dog regions and
        %require some mophological finesse
        if haveProbs
            se = strel('square', 3);
            im_thresh = imdilate(im_thresh, se); %thresholding from this classified probability map can produce non-contiguous, spurious Spots{channelIndex}. This fixes that and hopefully does not combine real Spots{channelIndex} from different nuclei
            im_thresh = im_thresh > 0;
        end
        
        [im_label, n_spots] = bwlabel(im_thresh);
        centroids = regionprops(im_thresh, 'centroid');
        
        
        if n_spots ~= 0
            
            for spotIndex = 1:n_spots
                centroid = round(centroids(spotIndex).Centroid);
                
                 Fits = identifySingleSpot(spotIndex,...
                    {im,imAbove,imBelow}, im_label, dog, ...
                    neighborhood, snippet_size, pixelSize_nm, shouldDisplayFigures,...
                    graphicsHandles, microscope, 0,...
                    centroid,MLFlag, currentFrame, spotIndex, zIndex);
                
                Spots(currentFrame).Fits = [Spots(currentFrame).Fits, Fits];
                
            end
            
        end
        
    end
    
    send(q, currentFrame);
    
end


try close(waitbarFigure); catch; end

    function nUpdateWaitbar(~)
        try waitbar(p/lastFrame, waitbarFigure); catch; end
        p = p + 1;
    end

end