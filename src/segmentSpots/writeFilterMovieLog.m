function log = writeFilterMovieLog(t, Weka, DropboxFolder, Prefix, initialFrame, numFrames, filterType, sigmas, classifierPathCh1)
    
    logFile = [DropboxFolder, filesep, Prefix, filesep, 'log.mat'];

    if exist(logFile, 'file')
      load(logFile, 'log');
    else
      log = struct();
    end

    log(end + 1).Date = date;
    
    log(end).Script = 'filterMovie';
    log(end).runTime = t / 60; % min
    log(end).InitialFrame = initialFrame;    
    log(end).LastFrame = numFrames;
    log(end).TimePerFrame = (t / 60) / numFrames;
    log(end).NFrames = numFrames - initialFrame + 1;
    log(end).TimePerFrame = (t / 60) / (numFrames - initialFrame + 1);
    
    if Weka
      log(end).Classifier = classifierPathCh1;
    else 
      log(end).Filter = filterType;
      log(end).sigmas = sigmas;
    end

    save(logFile, 'log', '-v7.3');

end
