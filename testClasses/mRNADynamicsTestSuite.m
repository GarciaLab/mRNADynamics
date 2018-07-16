try
  import matlab.unittest.TestSuite

  exportDataForFISHSuite = TestSuite.fromFolder('testExportDataForFISH');
  segmentSpotsSuite = TestSuite.fromFolder('testSegmentSpots');
  
  completeTestSuite = [exportDataForFISHSuite, segmentSpotsSuite];
  testResults = run(completeTestSuite);

  exit(any([testResults.Failed]));
catch
  disp('Cannot execute test suite.');
  exit(1);
end

%We can also do this, if we want to pick and choose the tests
%suite = TestSuite.fromFile('testKnownExperiments/test2015_07_25_P2P_75uW_bi.m');
%run(suite);
%suite = TestSuite.fromFile('testKnownExperiments/test2016_11_13_Hb_P2P_MS2V5_NB_MCP_mCherry.m');
%run(suite);
%suite = TestSuite.fromFile('testKnownExperiments/test2017_10_13_20170921_4.m');
%run(suite);