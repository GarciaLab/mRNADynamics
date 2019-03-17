% This is not so much a CheckParticleTracking test "per-se", but more of a
% test of the whole pipeline up to the CheckParticleTracking point, for a
% short test movie.
function testCase = testCheckParticleTracking(testCase)
    tic;
    disp(['Running Integration Test with prefix ', testCase.Prefix]);
    fprintf('Test run started at %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));

    disp('ExportDataForLivemRNA....');
    ExportDataForLivemRNA(testCase.Prefix);
    disp('ExportDataForLivemRNA, completed.');

    disp('FilterMovie....');
    filterMovie(testCase.Prefix);
    assertExpectedResults(testCase, 'filterMovie', 'FrameInfo.mat')
    disp('FilterMovie, completed.');

    disp('SegmentSpots....');
    segmentSpots(testCase.Prefix, testCase.DoG);
    assertExpectedResults(testCase, 'segmentSpots', 'FrameInfo.mat')
    assertExpectedResults(testCase, 'segmentSpots', 'Spots.mat')
    disp('SegmentSpots, completed.');

    disp('TrackmRNADynamics....');
    TrackmRNADynamics(testCase.Prefix);
    assertExpectedResults(testCase, 'trackmRNADynamics', 'FrameInfo.mat')
    assertExpectedResults(testCase, 'trackmRNADynamics', 'Spots.mat')
    assertExpectedResults(testCase, 'trackmRNADynamics', 'Particles.mat')
    disp('TrackmRNADynamics, completed.');

    disp('CheckParticleTracking....');
    CheckParticleTracking(testCase.Prefix, testCase.CheckParticleTrackingOption)
    assertExpectedResults(testCase, 'checkParticleTracking', 'FrameInfo.mat')
    assertExpectedResults(testCase, 'checkParticleTracking', 'Spots.mat')
    assertExpectedResults(testCase, 'checkParticleTracking', 'Particles.mat')
    disp('CheckParticleTracking, completed.');

    elapsedTime = toc;
    fprintf('Test run for %s ended successfully at %s\n', testCase.Prefix, datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
    fprintf('Elapsed time for test was %d minutes and %f seconds\n', floor(elapsedTime/60), rem(elapsedTime,60));
end

function assertExpectedResults(testCase, stepName, fileName)
  
  % because current directory could have change during execution of tests,
  % we position ourselves again in the root of the project (mRNADynamics folder)
  [~, ~, ~, srcPath, ~, ~, ~] = DetermineLocalFolders %this returns mRNADyanmics/src, we need one level less
  cd(srcPath)
  cd('..')
  
  
  expectedFrameInfoPath = ['./test/integrationTest/expected/', stepName, filesep, testCase.Prefix, filesep, fileName];
  actualFrameInfoPath = ['../Data/DynamicsResults/', testCase.Prefix, filesep, fileName];
  testCase.assertTrue(exist(actualFrameInfoPath, 'file') == 2);
  
  ExpectedFile = load(expectedFrameInfoPath);
  ActualFile = load(actualFrameInfoPath);

  [common, d1, d2] = comp_struct(ExpectedFile, ActualFile)
  
  testCase.assertNotEmpty(common);
  testCase.assertEmpty(d1);
  testCase.assertEmpty(d2);
end
