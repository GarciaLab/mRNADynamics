classdef test2015_07_25_P2P_75uW_bi < matlab.unittest.TestCase
%runs the export data process and compares preprocessed data with a known result set

  properties
      %Hardcoded with the path of the experiment that the test will use
      Prefix = '2015-07-25-P2P_75uW_bi';
      PreferredFileName = PreferredFileForTest('P2P_75uW_bi.lif');
  end

  methods(Test)

    function testRun(testCase)
      testCase = testExportDataForFISH(testCase);
    end

  end
end

