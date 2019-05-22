% Generates the TIF stacks necessary for doing Weka classification.
% Recommended to run this before making a new classifier.
function generateTifsForWeka(Prefix, ExperimentType, PreProcPath, numFrames,...
    nCh,coatChannel, zSize, initialFrame)

  % Create stacks subfolder
  stacksPath = [PreProcPath, filesep, Prefix, filesep, 'stacks'];
  mkdir(stacksPath);

  % Sets waitbar for the whole tif generation process
  tifStacksWaitbar = waitbar(0, 'Making .tif stacks for Weka classification');

  % Make requisite TIF stacks for classification
  for channelIndex = 1:nCh
    if strcmpi(ExperimentType, 'inputoutput')
      nameSuffix = ['_ch', iIndex(coatChannel, 2)];
    else
      nameSuffix = ['_ch', iIndex(channelIndex, 2)];
    end

    currentFrameWaitbar = waitbar(0, ['Making ch0', num2str(channelIndex), ' .tif stacks for Weka classification']);
  
    for currentFrame = initialFrame:numFrames
      currentFrameWaitbar = waitbar(currentFrame / numFrames, currentFrameWaitbar);
      set(currentFrameWaitbar, 'units', 'normalized', 'position', [0.4, .15, .25, .1]);
      rawStackName = [stacksPath, filesep, iIndex(currentFrame, 3), nameSuffix, '.tif'];
    
      %Don't write new stacks if they're already made.
      %2018-08-22 MT: Now takes into account 1 vs 2 spot channels
      %when determining if you've already made the stacks
      if length(dir([stacksPath, filesep, '*.tif'])) ~= numFrames * nCh
        rawStackArray = [];
    
        for i = 1:zSize
          fileName = [PreProcPath, filesep, Prefix, filesep, Prefix, '_', iIndex(currentFrame, 3), '_z', iIndex(i, 2), ...
            nameSuffix, '.tif'];
          rawStackArray(:, :, i) = imread(fileName);
        end
    
        imwrite(uint16(rawStackArray(:, :, 1)), rawStackName);
    
        for k = 2:size(rawStackArray, 3)
          imwrite(uint16(rawStackArray(:, :, k)), rawStackName, 'WriteMode', 'append');
        end
    
        clear rawStackArray;
      end
    
    end
    
    close(currentFrameWaitbar);
    waitbar(channelIndex / nCh, tifStacksWaitbar);      
  end

  close(tifStacksWaitbar);
end
