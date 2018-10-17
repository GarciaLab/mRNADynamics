function writeFilterMovieLog(t, justTifs, weka, DropboxFolder, Prefix, initialFrame, numFrames, filterType, sigmas, classifierPathCh1)
  if ~justTifs
    logFile = [DropboxFolder, filesep, Prefix, filesep, 'log.mat'];

    if exist(logFile, 'file')
      load(logFile);
    else
      log = struct();
    end

    log(end + 1).Date = date;
    log(end).runTime = t / 60; % min
    
    if weka
      log(end).InitialFrame = initialFrame;
    end
    
    log(end).LastFrame = numFrames;
    
    if ~weka
      log(end).TimePerFrame = (t / 60) / numFrames;
      log(end).Filter = filterType;
      log(end).sigmas = sigmas;
    else 
      log(end).NFrames = numFrames - initialFrame + 1;
      log(end).TimePerFrame = (t / 60) / (numFrames - initialFrame + 1);
      log(end).Classifier = classifierPathCh1;
    end

    save(logFile, 'log', '-v7.3');
  end

end
