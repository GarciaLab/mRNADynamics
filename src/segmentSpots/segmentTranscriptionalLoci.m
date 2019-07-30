function [Spots, dogs]...
    ...
    = segmentTranscriptionalLoci(...
    ...
    nCh, coatChannel, channelIndex, initialFrame, numFrames,...
    zSize, PreProcPath, Prefix, DogOutputFolder, displayFigures,doFF, ffim,...
    Threshold, neighborhood, snippet_size, pixelSize, microscope, intScale,...
    Weka, use_integral_center, filterMovieFlag, resultsFolder, gpu, saveAsMat, saveType, Ellipses)


dogs = [];

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
    chh = channelIndex;
else
    chh = coatChannel;
end

nameSuffix = ['_ch', iIndex(chh, 2)];

if Threshold == -1 && ~Weka
    
    if ~filterMovieFlag
        Threshold = determineThreshold(Prefix, chh);
        display(['Threshold: ', num2str(Threshold)])
    else
        Threshold = determineThreshold(Prefix, chh, 'noSave', dogs);
    end
    
    display(['Threshold: ', num2str(Threshold)])
    
end


q = parallel.pool.DataQueue;
afterEach(q, @nUpdateWaitbar);
p = 1;

if filterMovieFlag
    saveType = 'none';
end

isZPadded = false;

firstdogpath = [DogOutputFolder, filesep, dogStr, Prefix, '_', iIndex(1, 3), '_z', iIndex(1, 2),...
    nameSuffix];

matsPresent = exist([firstdogpath, '.mat'], 'file');
tifsPresent = exist([firstdogpath, '.tif'], 'file');

if isempty(saveType)
    if tifsPresent & ~matsPresent
        saveType = '.tif';
    elseif matsPresent & ~tifsPresent
        saveType = '.mat';
    elseif matsPresent & tifsPresent
        error('not sure which files to pick. check your processed folder.');
    end
end

firstdogpath = [firstdogpath, saveType];
if strcmpi(saveType, '.tif')
    firstDoG = imread(firstdogpath);
elseif strcmpi(saveType, '.mat')
    load(firstdogpath, 'plane');
    firstDoG = plane;
elseif strcmpi(saveType, 'none')
    firstDoG = dogs(:, :, 1, 1);
end

if sum(firstDoG(:)) == 0
    isZPadded = true;
end

parfor current_frame = initialFrame:numFrames
    
    for zIndex = 1:zSize
        
        imFileName = [PreProcPath, filesep, Prefix, filesep, Prefix, '_', iIndex(current_frame, 3), '_z', iIndex(zIndex, 2),...
            nameSuffix, '.tif'];
        im = double(imread(imFileName));
        try
            imAbove = double(imread([PreProcPath, filesep, Prefix, filesep, Prefix, '_', iIndex(current_frame, 3), '_z', iIndex(zIndex-1, 2),...
                nameSuffix, '.tif']));
            imBelow = double(imread([PreProcPath, filesep, Prefix, filesep, Prefix, '_', iIndex(current_frame, 3), '_z', iIndex(zIndex+1, 2),...
                nameSuffix, '.tif']));
        catch
            imAbove = nan(size(im,1),size(im,2));
            imBelow = nan(size(im,1),size(im,2));
        end
        
        
        if isZPadded
            dogZ = zIndex;
        else
            dogZ = zIndex - 1;
        end
        
%         try
            if isZPadded | ( ~isZPadded & (zIndex~=1 & zIndex~=zSize) )
                if strcmpi(saveType, '.tif')
                    dogFileName = [DogOutputFolder, filesep, dogStr, Prefix, '_', iIndex(current_frame, 3), '_z', iIndex(dogZ, 2),...
                        nameSuffix,'.tif'];
                    dog = double(imread(dogFileName));
                elseif strcmpi(saveType, '.mat')
                    dogFileName = [DogOutputFolder, filesep, dogStr, Prefix, '_', iIndex(current_frame, 3), '_z', iIndex(dogZ, 2),...
                        nameSuffix,'.mat'];
                    plane = load(dogFileName, 'plane');
                    dog = plane.plane;
                elseif strcmpi(saveType, 'none')
                    dog = dogs(:,:, dogZ, current_frame);
                end
            else
                dog = false(size(im, 1), size(im, 2));
            end
%         catch
%             error('Please run filterMovie to create DoG files.');
%         end
        
        
        
        
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
            ellipsesFrame = Ellipses{current_frame};
            nuclearMask = makeNuclearMask(ellipsesFrame, [size(im,1), size(im,2)]);
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
                    neighborhood, snippet_size, pixelSize, displayFigures, graphicsHandles, microscope, 0, centroid,MLFlag, intScale, current_frame, spotIndex, zIndex, use_integral_center);
                Spots(current_frame).Fits = [Spots(current_frame).Fits, Fits];
            end
            
        end
        
    end
    
    send(q, current_frame);
    
end


close(waitbarFigure);

    function nUpdateWaitbar(~)
        waitbar(p/numFrames, waitbarFigure);
        p = p + 1;
    end

end