function log = logSegmentSpots(DropboxFolder, Prefix, t, initialFrame, numFrames, Spots, falsePositives, Threshold, channelIndex, numShadows, intScale, fit3D)
  
  logFile = [DropboxFolder, filesep, Prefix, filesep, 'log.mat'];

  if exist(logFile, 'file')
    load(logFile);
  else 
    log = struct();
  end 

  log(end + 1).Date = date;
  
  log(end).Script = 'segmentSpots';
  log(end).runTime = t / 60; %min

  if ~isempty(initialFrame) 
    log(end).InitialFrame = initialFrame;
    log(end).LastFrame = numFrames;
    log(end).NFrames = numFrames - initialFrame + 1;
    log(end).TimePerFrame = (t / 60) / (numFrames - initialFrame + 1);
  else
    log(end).LastFrame = numFrames;
    log(end).TimePerFrame = (t / 60) / numFrames;
  end

  detectedCircles = 0;
  detectedBalls = 0;

  if iscell(Spots)

    for i = 1:length(Spots{channelIndex})

      for j = 1:length(Spots{channelIndex}(i).Fits)
        detectedCircles = detectedCircles + length(Spots{channelIndex}(i).Fits(j).z);
        detectedBalls = detectedBalls + 1;
      end 

    end 

  else 

    for i = 1:length(Spots)

      for j = 1:length(Spots(i).Fits)
        detectedCircles = detectedCircles + length(Spots(i).Fits(j).z);
        detectedBalls = detectedBalls + 1;
      end 

    end 

  end 

  display(['Detected Spots: ', num2str(detectedCircles)])
  log(end).falsePositives = falsePositives;
  log(end).totalCircles = detectedCircles;
  log(end).totalBalls = detectedBalls;
  log(end).avgZSize = detectedCircles / detectedBalls;
  log(end).Threshold = Threshold;
  log(end).numShadows = numShadows;
  log(end).intScale = intScale;
  log(end).fit3D = fit3D;

  if isfield(log, 'Classifier')
    log(end).Classifier = 'no Weka';
  end 

  save(logFile, 'log', '-v6');
end 
