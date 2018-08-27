function all_frames = segmentTranscriptionalLoci(ExperimentType, coatChannel, channelIndex, all_frames, numFrames, zSize, PreProcPath, Prefix, DogOutputFolder, displayFigures, pool, doFF, ffim, Threshold, neighborhood, snippet_size, pixelSize, microscope, intScale)
  waitbarFigure = waitbar(0, 'Segmenting spots');

  % (MT, 2018-02-11) Added support for lattice imaging, maybe 
  % temporary - FIX LATER
  if strcmpi(ExperimentType, 'inputoutput') ||  strcmpi(ExperimentType, 'lattice')
    nameSuffix= ['_ch', iIndex(coatChannel, 2)];
    if Threshold == -1
        Threshold = determineThreshhold(Prefix, coatChannel);
    end
  else
    nameSuffix = ['_ch', iIndex(channelIndex, 2)];
    if Threshold == -1
        Threshold = determineThreshhold(Prefix, channelIndex);
    end
  end        
  
  for current_frame = 1:numFrames
    
    w = waitbar(current_frame / numFrames, waitbarFigure);
    set(w, 'units', 'normalized', 'position', [0.4, .15, .25,.1]);
        
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
        dogFileName = [DogOutputFolder, filesep, 'DOG_', Prefix, '_', iIndex(current_frame, 3), '_z', iIndex(zIndex, 2),...
          nameSuffix,'.tif'];
        dog = double(imread(dogFileName));
      catch
        error('Please re-run with threshold ''[]'' to create DoG files')
      end
      
      if displayFigures
        fig = figure(1);
        imshow(dog, [median(dog(:)), max(dog(:))]);
      else
        fig=[];
      end
            
      % Apply flatfield correction
      if doFF && sum(size(im)==size(ffim))
        im = im.*ffim;
      end
      
      im_thresh = dog >= Threshold;
      [im_label, n_spots] = bwlabel(im_thresh); 
      centroids = regionprops(im_thresh, 'centroid');

      temp_frames = {};
      temp_particles = cell(1, n_spots);
      
      if n_spots ~= 0
        if ~displayFigures && pool
          parfor spotIndex = 1:n_spots
            centroid = round(centroids(spotIndex).Centroid);
            temp_particles(spotIndex) = identifySingleSpot(spotIndex, {im,imAbove,imBelow}, im_label, dog, ...
              neighborhood, snippet_size, pixelSize, displayFigures, fig, microscope, 0, centroid, '', intScale);
          end
        else
          for spotIndex = 1:n_spots
            centroid = round(centroids(spotIndex).Centroid);
            temp_particles(spotIndex) = identifySingleSpot(spotIndex, {im,imAbove,imBelow}, im_label, dog, ...
              neighborhood, snippet_size, pixelSize, displayFigures, fig, microscope, 0, centroid, '', intScale);
          end
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
