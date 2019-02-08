% This is not so much a CheckParticleTracking test "per-se", but more of a
% test of the whole pipeline up to the CheckParticleTracking point, for a
% short test movie.
function testCase = testCheckParticleTracking(testCase)
    try
        tic;
        disp(['Running CheckParticleTracking test with prefix ', testCase.Prefix]);
        fprintf('Test run started at %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));

        ExportDataForLivemRNA(testCase.Prefix);
        filterMovie(testCase.Prefix);
        segmentSpots(testCase.Prefix, testCase.DoG);
        TrackmRNADynamics(testCase.Prefix);
        CheckParticleTracking(testCase.Prefix, testCase.CheckParticleTrackingOption)
        
        elapsedTime = toc;
        fprintf('Test run for %s ended successfully at %s\n', testCase.Prefix, datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
        fprintf('Elapsed time for test was %d minutes and %f seconds\n', floor(elapsedTime/60), rem(elapsedTime,60));
    catch Exception
        disp(['Test ended with exception: ', Exception]);
    end
end
