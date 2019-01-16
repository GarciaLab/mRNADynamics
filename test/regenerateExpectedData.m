%exportTest = testExport_2015_07_25_P2P_75uW_bi_short;
%segmentSpotsMLTest = testSegmentSpotsML_2015_07_25_P2P_75uW_bi_short;
%trackmRNADynamicsTestCase = testTrackmRNADynamics_2015_07_25_P2P_75uW_bi_short;
%generateExpectedDataForPrefix(exportTest, segmentSpotsMLTest, trackmRNADynamicsTestCase);

% No histone path of traackmRNADynamics is processed separately
%trackmRNADynamicsNotHistoneTestCase = testTrackmRNADynamicsNotHistone_2015_07_25_P2P_75uW_bi_short;
%generateExpectedDataTrackmRNADynamicsNotHistone(trackmRNADynamicsNotHistoneTestCase);

% Segment spots only has one test so far
segmentSpotsTest = testSegmentSpots_2015_07_25_P2P_75uW_bi_short;
generateExpectedDataSegmentSpots(segmentSpotsTest);

% Other movies to test export data
exportTest = testExport_2016_11_13_Hb_P2P_MS2V5_NB_MCP_mCherry;
generateExpectedExportData(exportTest);

exportTest = testExport_2017_10_13_20170921_4;
generateExpectedExportData(exportTest);

exportTest = testExport_2018_06_05_A140P_MSE_30uW_550V;
generateExpectedExportData(exportTest);
