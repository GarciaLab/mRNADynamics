classdef testExport_2016_11_13_Hb_P2P_MS2V5_NB_MCP_mCherry < matlab.unittest.TestCase
%runs the export data process and compares preprocessed data with a known result set

  properties
      %Hardcoded with the path of the experiment that the test will use
      Prefix = '2016-11-13-Hb-P2P-MS2V5-NB-MCP-mCherry';
  end

  methods(Test)

    function testRun(testCase)
      testCase = testExportDataForFISH(testCase);
    end

  end
end

