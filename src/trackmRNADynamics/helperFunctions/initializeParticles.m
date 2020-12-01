function Particles = initializeParticles(has3DInfo)
    if ~has3DInfo
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
            'StableParticleID',{},...
            'obsFrameFilter', {});
    else
        Particles = struct(...
            'Frame',{}, ...
            'Index', {}, ...
            'xPos', {}, ...
            'yPos', {}, ...
            'zPosDetrended', {}, ...
            'zPos', {}, ...
            'xPos3D', {}, ...
            'yPos3D', {}, ...
            'zPosDetrended3D', {}, ...
            'zPos3D', {}, ...
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
            'xPosInf3D', {}, ...
            'yPosInf3D', {}, ...
            'zPosDetrendedInf3D', {}, ...
            'xPosSEInf3D', {}, ...
            'yPosSEInf3D', {}, ...
            'zPosSEInf3D', {}, ...
            'zPosInf3D', {}, ...
            'smoothedPredictions3D', {},...
            'smoothedPredictionsSE3D',{},...
            'StableParticleID',{},...
            'obsFrameFilter', {});
    end