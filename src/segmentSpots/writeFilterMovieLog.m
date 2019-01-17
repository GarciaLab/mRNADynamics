function log = writeFilterMovieLog(t, Weka, DropboxFolder, Prefix, initialFrame, numFrames, filterType, sigmas, classifierPathCh1)
    
    logFile = [DropboxFolder, filesep, Prefix, filesep, 'log.mat'];

    if exist(logFile, 'file')
      load(logFile, 'log');
    else
      log = struct();
    end

    log(end + 1).Date = date;
    log(end).runTime = t / 60; % min
    
    if Weka
      log(end).InitialFrame = initialFrame;
    end
    
    log(end).LastFrame = numFrames;
    
    if ~Weka
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
