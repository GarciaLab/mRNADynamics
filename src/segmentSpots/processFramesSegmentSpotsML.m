function allFrames = processFramesSegmentSpotsML(Prefix, Threshold, PreProcPath, dogsFolder, channelIndex, nameSuffix, initialFrame, numFrames, zSize, doFF, ffim, displayFigures, nWorkers, segmentSpotsWaitbar, allFrames, snippet_size, neighborhood, microscope, intScale, pixelSize)

  for currentFrame = initialFrame:numFrames
    currentFrameWaitbar = waitbar(currentFrame / (numFrames - initialFrame), segmentSpotsWaitbar);
    set(currentFrameWaitbar, 'units', 'normalized', 'position', [0.4, .15, .25, .1]);

    tifFileNamePrefix = [PreProcPath, filesep, Prefix, filesep, Prefix, '_', iIndex(currentFrame, 3), '_z'];
    tifFileNameSuffix = [nameSuffix, '.tif'];

    for zIndex = 1:zSize
      im = double(imread([tifFileNamePrefix, iIndex(zIndex, 2), tifFileNameSuffix]));

      try
        imAbove = double(imread([tifFileNamePrefix, iIndex(zIndex - 1, 2), tifFileNameSuffix]));
        imBelow = double(imread([tifFileNamePrefix, iIndex(zIndex + 1, 2), tifFileNameSuffix]));
      catch
        imAbove = nan(size(im, 1), size(im, 2));
        imBelow = nan(size(im, 1), size(im, 2));
      end

      fileName = [dogsFolder, filesep, 'prob', Prefix, '_', iIndex(currentFrame, 3), '_z', iIndex(zIndex, 2), tifFileNameSuffix];
      pMap = double(imread(fileName));

      if displayFigures
        fig = figure(1);
        imshow(im, []);
      else
        fig = [];
      end

      %apply flatfield correction
      if doFF && sum(size(im) == size(ffim))
        im = im .* ffim;
      end

      imThresh = pMap >= Threshold(channelIndex);
      se = strel('square', 3);
      imThresh = imdilate(imThresh, se); %thresholding from this classified probability map can produce non-contiguous, spurious Spots{channelIndex}. This fixes that and hopefully does not combine real Spots{channelIndex} from different nuclei
      imThresh = imThresh > 0;
      [imLabel, nSpots] = bwlabel(imThresh);
      centroids = regionprops(imThresh, 'centroid');

      tempFrames = {};
      tempParticles = cell(1, nSpots);

      if nSpots ~= 0

        if ~ displayFigures && nWorkers > 0

          parfor spotIndex = 1:nSpots

            try
              tempParticles(spotIndex) = processSpot(spotIndex, centroids, im, imAbove, imBelow, imLabel, pMap,...
                neighborhood, snippet_size, pixelSize, displayFigures, fig, microscope, intScale);
            catch
            end

          end

        else

          for spotIndex = 1:nSpots
            tempParticles(spotIndex) = processSpot(spotIndex, centroids, im, imAbove, imBelow, imLabel, pMap,...
              neighborhood, snippet_size, pixelSize, displayFigures, fig, microscope, intScale);
          end

        end

        for spotIndex = 1:nSpots

          if ~ isempty(tempParticles{spotIndex})
            tempFrames = [tempFrames, tempParticles(spotIndex)];
          end

        end

        allFrames{currentFrame, zIndex} = tempFrames;
      end

    end

  end

end

function tempParticlesForSpot = processSpot(spotIndex, centroids, im, imAbove, imBelow, imLabel, pMap, neighborhood, snippet_size, pixelSize, displayFigures, fig, microscope, intScale)
  centroid = round(centroids(spotIndex).Centroid);
  tempParticlesForSpot = identifySingleSpot(spotIndex, {im, imAbove, imBelow}, imLabel, pMap, ...
  neighborhood, snippet_size, pixelSize, displayFigures, fig, microscope, 0, centroid, 'ML', intScale);
end
