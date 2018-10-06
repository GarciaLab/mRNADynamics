% Generates differences of gaussians files that will be later used to run segment spots ML, requires a classifier.
function generateDogsWeka(Prefix, FISHPath, MS2CodePath, PreProcPath, ExperimentType, coatChannel, zSize, numFrames, nCh, initialFrame, ignoreMemoryCheck, classifierPathCh1, classifierFolder)

  OutputFolder1 = [FISHPath, filesep, Prefix, '_', filesep, 'dogs'];
  mkdir(OutputFolder1)

  % Create stacks subfolder
  stacksPath = [PreProcPath, filesep, Prefix, filesep, 'stacks'];

  if isempty(classifierPathCh1)
    [classifierPathCh1, classifierFolder] = uigetfile([MS2CodePath, filesep, 'classifiers', filesep, '*.model']);
  end

  if nCh == 2
    [classifierPathCh2, ~] = uigetfile([MS2CodePath, filesep, 'classifiers', filesep, '*.model']);
  end

  evalin('base', 'clear probmaps');

  heapSize = java.lang.Runtime.getRuntime.maxMemory;
  if heapSize < 1E10 && ~ ignoreMemoryCheck
    error('Please increase your Java heap memory allocation to at least 10GB (Home -> Preferences -> General -> Java Heap Memory.');
  end

  zSize2 = zSize * 2;
  filterMovieWaitbarFigure = waitbar(0, 'Running Weka Classifier');

  %Make requisite TIF stacks for classification
  for channelIndex = 1:nCh

    try
      %this is just some function that can only be called if IJM is set up
      IJM.getIdentifier()
    catch
      addpath([MS2CodePath, filesep, 'Fiji.app', filesep, 'scripts'])
      ImageJ % Initialize IJM and MIJ
    end

    ijm = evalin('base', 'IJM');
    mij = evalin('base', 'MIJ');

    if strcmpi(ExperimentType, 'inputoutput')
      nameSuffix = ['_ch', iIndex(coatChannel, 2)];
    else
      nameSuffix = ['_ch', iIndex(channelIndex, 2)];
    end

    currentFrameWaitbar = waitbar(0, ['Making ch0', nCh, '.tif stacks for Weka classification']);

    for currentFrame = initialFrame:numFrames
      currentFrameWaitbar = waitbar(currentFrame / numFrames, currentFrameWaitbar);
      set(currentFrameWaitbar, 'units', 'normalized', 'position', [0.4, .15, .25, .1]);
      rawStackName = [stacksPath, filesep, iIndex(currentFrame, 3), nameSuffix, '.tif'];

      %Do the classification with Weka in Fiji
      mij.run('Trainable Weka Segmentation 3D', ['open=', rawStackName]);
      pause(10);

      if channelIndex == 1 || strcmpi(ExperimentType, 'inputoutput')
        trainableSegmentation.Weka_Segmentation.loadClassifier([classifierFolder, classifierPathCh1]);
      elseif channelIndex == 2
        trainableSegmentation.Weka_Segmentation.loadClassifier([classifierFolder, classifierPathCh2]);
      else
        error(['This pipeline does not support', ...
               'more than two spot channels. If you''re actually', ...
               'trying to segment 3 or more channels, talk to Armando to', ...
               'get this implemented. Otherwise you''ve reached an error.', ...
               'Check your data. This is probably not a bug in the code.']);
      end

      trainableSegmentation.Weka_Segmentation.getProbability();
      ijm.getDatasetAs('probmaps')
      pMapTemp = evalin('base', 'probmaps');
      pMap = [];

      for m = 1:2:zSize2
        pMap(:, :, ceil(m / 2)) = pMapTemp(:, :, m); %the even images in the original array are negatives of the odds
      end

      clear pMapTemp;
      pMap = permute(pMap, [2 1 3]) * 10000; %multiplying so this can be cast to uint16

      for i = 1:size(pMap, 3)
        p_name = ['prob', Prefix, '_', iIndex(currentFrame, 3), '_z', iIndex(i, 2), nameSuffix, '.tif'];
        imwrite(uint16(pMap(:, :, i)), [OutputFolder1, filesep, p_name])
      end

      mij.run('Close All');
      clear pMap;

    end

    close(currentFrameWaitbar);
    waitbar(channelIndex / nCh, filterMovieWaitbarFigure);      
  end

  close(filterMovieWaitbarFigure);

end
