classdef integrationTest_2015_07_25_P2P_75uW_bi_short < matlab.unittest.TestCase

  properties
      Prefix = '2015-07-25-P2P_75uW_bi_short';
      DoG = 60;
      CheckParticleTrackingOption = 'ForCompileAll';
  end

  methods(Test)

    function testRun(testCase)
        integrationTest(testCase);
    end

  end
end

