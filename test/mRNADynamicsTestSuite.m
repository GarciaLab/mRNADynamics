addpath(path, genpath('../testClasses'));

try
  tic;
  fprintf('Test suite started running at %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  import matlab.unittest.TestSuite

  exportDataForFISHSuite = TestSuite.fromFolder('testExportDataForFISH');
  segmentSpotsSuite = TestSuite.fromFolder('testSegmentSpots');
  segmentSpotsMLSuite = TestSuite.fromFolder('testSegmentSpotsML');
  trackmRNADynamicsSuite = TestSuite.fromFolder('testTrackmRNADynamics');
  
  completeTestSuite = [exportDataForFISHSuite, segmentSpotsSuite, segmentSpotsMLSuite, trackmRNADynamicsSuite];
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
  %exit(1);
end

%We can also do this, if we want to pick and choose the tests
%suite = TestSuite.fromFile('testKnownExperiments/test2015_07_25_P2P_75uW_bi.m');
%run(suite);
%suite = TestSuite.fromFile('testKnownExperiments/test2016_11_13_Hb_P2P_MS2V5_NB_MCP_mCherry.m');
%run(suite);
%suite = TestSuite.fromFile('testKnownExperiments/test2017_10_13_20170921_4.m');
%run(suite);