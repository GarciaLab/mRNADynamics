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
        'ManuallyReviewed', {}, ...
        'framesFull', {}, ...
        'logL', {}, ...
        'logLMean', {}, ...
        'logLArray', {}, ...
        'xPosInf', {}, ...
        'yPosInf', {}, ...
        'zPosDetrendedInf', {}, ...
        'xPosSEInf', {}, ...
        'yPosSEInf', {}, ...
        'zPosSEInf', {}, ...
        'zPosInf', {}, ...
        'smoothedPredictions', {},...
        'smoothedPredictionsSE',{},...
        'obsFrameFilter', {});%,...
%         'logLDistance', {});