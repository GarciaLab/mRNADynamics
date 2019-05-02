function [all_frames, Spots] = segmentTranscriptionalLoci(nCh, coatChannel, channelIndex, all_frames, initialFrame, numFrames, zSize, PreProcPath, Prefix, DogOutputFolder, displayFigures,doFF, ffim, Threshold, neighborhood, snippet_size, pixelSize, microscope, intScale, Weka, use_integral_center)
  
  waitbarFigure = waitbar(0, 'Segmenting spots');

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
      Threshold = 5000;
  else
      MLFlag = '';
      dogStr = 'DOG_';
  end
  
  %Check how many coat channels we have and segment the appropriate channel
  %accordingly
  if nCh > 1
    nameSuffix = ['_ch', iIndex(channelIndex, 2)];
    if Threshold == -1 && ~Weka
        Threshold = determineThreshold(Prefix, channelIndex);
        display(['Threshold: ', num2str(Threshold)])
    end
  else
    nameSuffix= ['_ch', iIndex(coatChannel, 2)];
    if Threshold == -1 && ~Weka
        Threshold = determineThreshold(Prefix, coatChannel);
        display(['Threshold: ', num2str(Threshold)])
    end
  end        
  
  parfor current_frame = initialFrame:numFrames
    
    %w = waitbar(current_frame / numFrames, waitbarFigure);
    %set(w, 'units', 'normalized', 'position', [0.4, .15, .25,.1]);

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
      
      try
        dogFileName = [DogOutputFolder, filesep, dogStr, Prefix, '_', iIndex(current_frame, 3), '_z', iIndex(zIndex, 2),...
          nameSuffix,'.tif'];
        dog = double(imread(dogFileName));
      catch
        error('Please run filterMovie to create DoG files')
      end
      
      if displayFigures
          dogO = dog(:);
          lLim = median(dogO);
          uLim = max(dogO);
        if lLim ~= uLim
         imagesc(dogAx,dog, [lLim, uLim]);
         colormap(dogAx, 'gray')

        else
            imagesc(dogAx,dog);
            colormap(dogAx, 'gray')
        end
      end
            
      % Apply flatfield correction
      if doFF && sum(size(im)==size(ffim))
        im = im.*ffim;
      end
            
      im_thresh = dog >= Threshold;
      
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
        
        for spotIndex = 1:n_spots
          if ~isempty(temp_particles{spotIndex})
            temp_frames = [temp_frames, temp_particles(spotIndex)];
          end
        end

        all_frames{current_frame, zIndex} = temp_frames;
      end
      
    end

  end

  close(waitbarFigure);

end
