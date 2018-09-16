testCases{1} = testExport_2015_07_25_P2P_75uW_bi_short;
testCases{2} = testExport_2016_11_13_Hb_P2P_MS2V5_NB_MCP_mCherry;
testCases{3} = testExport_2017_10_13_20170921_4;
testCases{4} = testExport_2018_06_05_A140P_MSE_30uW_550V;

for i = 1:length(testCases)
   generateExpectedExportData(testCases{i});
end