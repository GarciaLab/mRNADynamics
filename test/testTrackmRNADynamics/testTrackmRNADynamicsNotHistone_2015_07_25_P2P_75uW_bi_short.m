classdef testTrackmRNADynamicsNotHistone_2015_07_25_P2P_75uW_bi_short < matlab.unittest.TestCase
%runs the Segment Spots data process and compares preprocessed data with a known result set

  properties
      %Hardcoded with the path of the experiment that the test will use
      Prefix = '2015-07-25-P2P_75uW_bi_short';
      Threshold1 = 0;
      Threshold2 = 0;
  end

  methods(Test)

    function testRun(testCase)
      testTrackmRNADynamicsNotHistone(testCase);
    end

  end
end

