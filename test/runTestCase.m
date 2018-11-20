function runTestCase(testCaseName)
  diary 'tests.log';
  
  try
    tc = testCaseName;
    testResult = tc.run;
    exit(testResult.Failed);
  catch ME
    disp('Cannot execute test suite.');
    disp(['Exception: ', ME.identifier]);
    disp(['Message: ', ME.message]);
    disp(['Cause: ', ME.cause]);
    exit(1);
  end

  diary off;
  
end
