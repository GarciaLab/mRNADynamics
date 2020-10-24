function particleTracks = initializeParticleTracks()
    % create an empty array of tracks
    particleTracks = struct(...
        'Index', {}, ... 
        'Frame', {}, ... 
        'MeasurementVec', {}, ...   
        'zPos', {}, ...
        'Nucleus', {}, ...
        'kalmanFilter', {}, ...
        'age', {}, ...
        'totalVisibleCount', {}, ...
        'consecutiveInvisibleCount', {},...
        'firstFrame', {}, ...
        'lastFrame', {});
end