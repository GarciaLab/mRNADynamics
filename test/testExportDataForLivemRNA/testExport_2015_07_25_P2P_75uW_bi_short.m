classdef testExport_2015_07_25_P2P_75uW_bi_short < matlab.unittest.TestCase
%runs the export data process and
%compares preprocessed data with a known result set

  properties
      %Hardcoded with the path of the experiment that the test will use
      Prefix = '2015-07-25-P2P_75uW_bi_short';
      FileName = 'P2P_75uW_bi_short.lif';
      PreferredFileName;
  end

  methods(TestClassSetup)
    function instance = initializeTestCase(testCase)
        testCase.PreferredFileName = PreferredFileForTest(testCase.FileName);
    end
  end

  methods(Test)

    function testRun(testCase)
      testCase = testExportDataForLivemRNA(testCase);
    end

  end
end