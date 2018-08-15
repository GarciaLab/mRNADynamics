function log = logSegmentSpots(DropboxFolder, Prefix, t, numFrames, justDoG, Spots, filterType, sigmas, falsePositives, Threshold)
  logFile = [DropboxFolder, filesep, Prefix, filesep, 'log.mat'];

  if exist(logFile, 'file')
    load(logFile);
  else 
    log = struct();
  end 

  log(end + 1).Date = date;
  log(end).runTime = t / 60; %min
  % log(end).InitialFrame = initial_frame;
  log(end).LastFrame = numFrames;
  % log(end).NFrames = num_frames - initial_frame + 1;
  % log(end).TimePerFrame = (t/60)/(num_frames-initial_frame + 1);
  log(end).TimePerFrame = (t / 60) / numFrames;

  if ~ justDoG
    detectedCircles = 0;
    detectedBalls = 0;

    if iscell(Spots)

      for i = 1:length(Spots{q})

        for j = 1:length(Spots{q}(i).Fits)
          detectedCircles = detectedCircles + length(Spots{q}(i).Fits(j).z);
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

    if isfield(log, 'Classifier')
      log(end).Classifier = 'no Weka';
    end 

  else 
    log(end).Filter = filterType;
    log(end).sigmas = sigmas;
  end 

  save(logFile, 'log', '-v7.3');
end 
