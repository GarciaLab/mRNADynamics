function Particles = initializeParticles()

    Particles = struct(...
        'Frame',{}, ...
        'Index', {}, ...
        'xPos', {}, ...
        'yPos', {}, ...
        'zPosDetrended', {}, ...
        'zPos', {}, ...
        'Nucleus', {}, ...
        'FirstFrame', {}, ...
        'LastFrame', {}, ...
        'Approved', {}, ...
        'FrameApproved', {}, ...
        'framesFull', {}, ...
        'logL', {}, ...
        'logLMean', {}, ...
        'xPosInf', {}, ...
        'yPosInf', {}, ...
        'zPosDetrendedInf', {}, ...
        'xPosSEInf', {}, ...
        'yPosSEInf', {}, ...
        'zPosSEInf', {}, ...
        'zPosInf', {}, ...
        'smoothedPredictions', {},...
        'smoothedPredictionsSE',{},...
        'obsFrameFilter', {});