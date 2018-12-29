try
  tic;
  fprintf('Test suite started running at %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  import matlab.unittest.TestSuite

  exportDataForLivemRNASuite = TestSuite.fromFolder('test/testExportDataForLivemRNA');
  segmentSpotsSuite = TestSuite.fromFolder('test/testSegmentSpots');
  segmentSpotsMLSuite = TestSuite.fromFolder('test/testSegmentSpotsML');
  trackmRNADynamicsSuite = TestSuite.fromFolder('test/testTrackmRNADynamics');
  
  completeTestSuite = [exportDataForLivemRNASuite, segmentSpotsSuite, segmentSpotsMLSuite, trackmRNADynamicsSuite];
  testResults = run(completeTestSuite);

  elapsedTime = toc;
  fprintf('Test suite run ended at %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  fprintf('Elapsed time for test suite was %d minutes and %f seconds\n', floor(elapsedTime/60), rem(elapsedTime,60));

  exit(any([testResults.Failed]));
catch ME
  disp('Cannot execute test suite.');
  disp(['Exception: ', ME.identifier]);
  disp(['Message: ', ME.message]);
  disp(['Cause: ', ME.cause]);
  exit(1);
end
