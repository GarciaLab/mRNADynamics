

function Spots...
    ...
    = segment2DLoci(...
    ...
    ch_quantify, initialFrame, lastFrame,...
    zSize, ~, Prefix, ~, shouldDisplayFigures,doFF, ffim,...
    Threshold, neighborhood, snippet_size, ~, microscope,...
    ~, resultsFolder, gpu, saveAsMat, ~, shouldMaskNuclei,...
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
haveStacks = any(cellfun(@(x) isempty(regexp(x, '_z[0-9]')), {dogDir.name}));

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


movieMat = getMovieMat(liveExperiment);
if ~isempty(movieMat)
    movieMat_channel = movieMat(:, :, :, :, ch_quantify);
else
    movieMat_channel = [];
end

yDim = liveExperiment.yDim;
xDim = liveExperiment.xDim;
zDim = liveExperiment.zDim;



if autoThresh
    
    Threshold = determineThreshold(Prefix,...
        ch_segment,  'numFrames', lastFrame, 'firstFrame', initialFrame);
    display(['Threshold: ', num2str(Threshold)])
    
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
            dogStack = imreadStack([dogStackFile, '.tif']);
        else
            error('Cannot find any dogs in the ProcessedData folder. Are they missing or incorrectly named?')
        end
        
        
    end
    
    TempZSize = size(dogStack, 3);
    for zIndex = 1:TempZSize
        
        
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
            %             se = strel('square', 3);
            %             im_thresh = imdilate(im_thresh, se); %thresholding from this classified probability map can produce non-contiguous, spurious Spots{channelIndex}. This fixes that and hopefully does not combine real Spots{channelIndex} from different nuclei
            im_thresh = im_thresh > 0;
        end
        
        [im_label, n_spots] = bwlabel(im_thresh);
        
        if n_spots ~= 0
            spots_regionprops = regionprops('table', im_label, 'centroid', 'EquivDiameter',...
                'Area', 'Circularity', 'Eccentricity', 'MajorAxisLength', 'MinorAxisLength', 'PixelList');
            
            for spotIndex = 1:n_spots
                spot_props = [table2array(spots_regionprops(spotIndex, 'Centroid')),...
                    table2array(spots_regionprops(spotIndex, 'EquivDiameter')),...
                    table2array(spots_regionprops(spotIndex,'Area')),...
                    table2array(spots_regionprops(spotIndex, 'Circularity')),...
                    table2array(spots_regionprops(spotIndex, 'Eccentricity')),...
                    table2array(spots_regionprops(spotIndex, 'MajorAxisLength')),...
                    table2array(spots_regionprops(spotIndex, 'MinorAxisLength')),...
                    table2array(spots_regionprops(spotIndex, 'PixelList'))];
                Fits = identify2DSpot(spotIndex,...
                    im, im_label, dog, ...
                    neighborhood, snippet_size, pixelSize_nm, shouldDisplayFigures,...
                    graphicsHandles, 0,...
                    spot_props,MLFlag, currentFrame, spotIndex, zIndex);
                
                Spots(currentFrame).Fits = [Spots(currentFrame).Fits, Fits];
                
            end
            
        end
        
    end
    
    send(q, currentFrame);
    
end

% % %
try close(waitbarFigure); catch; end

    function nUpdateWaitbar(~)
        try waitbar(p/lastFrame, waitbarFigure); catch; end
        p = p + 1;
    end

end