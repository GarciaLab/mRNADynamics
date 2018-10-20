classdef testGenerateDogsWeka_2015_07_25_P2P_75uW_bi_short < matlab.unittest.TestCase
%runs the generate dogs process for Weka and compares preprocessed data with a known result set

  properties
      %Hardcoded with the path of the experiment that the test will use
      Prefix = '2015-07-25-P2P_75uW_bi_short';
      classifier = 'test_classifier.model';
  end

  methods(Test)

    function testRun(testCase)
      testGenerateDogsWeka(testCase);
    end

  end
end

