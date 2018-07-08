classdef test2018_06_05_A140P_MSE_30uW_550V < matlab.unittest.TestCase
%runs the export data process and compares preprocessed data with a known result set

  properties
      %Hardcoded with the path of the experiment that the test will use
      Prefix = '2018-06-05-A140P_MSE_30uW_550V';
  end

  methods(Test)

    function testRun(testCase)
      testCase = testExportDataForFISH(testCase);
    end

  end
end

