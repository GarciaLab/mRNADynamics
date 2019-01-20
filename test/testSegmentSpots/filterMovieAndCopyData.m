function filterMovieAndCopyData(testCase, testPath, codePath, processedDataPath)
  expectedDataFolder = [testPath, '/filterMovie'];
  if 7 ~= exist(expectedDataFolder, 'dir')
    mkdir(expectedDataFolder);
  end

  processedDataExpectedPath = [expectedDataFolder, '/ProcessedData/', testCase.Prefix, '_'];

  % Classifier path must end with filesep
  classifiersPath = strcat(codePath, '/src/classifiers/');
  classifierForTest = ClassifierForTest(classifiersPath, testCase.classifier);

  % Generates expected data for DoGs creation
  filterMovie(testCase.Prefix, classifierForTest, 'ignoreMemoryCheck');

  % Then copy expected data for filterMovie pass
  copyExpectedDataFolder(processedDataPath, processedDataExpectedPath, 'filterMovie', testCase.Prefix)

  disp(['Test data for filterMovie completed for Prefix ', testCase.Prefix]);
end
