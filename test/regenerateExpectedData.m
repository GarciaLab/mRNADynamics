% ExportDataForFISH
run('testExportDataForFISH/generateAllExpectedDataForExport.m');

% Segment spots only has one test so far
tc = testSegmentSpots_2015_07_25_P2P_75uW_bi_short;
generateExpectedDataSegmentSpots(tc);