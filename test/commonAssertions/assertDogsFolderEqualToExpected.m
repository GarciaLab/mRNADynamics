% Given the path of the processed data folder and the expected data path asserts dogs a sub
% folder with expected data set.
function assertDogsFolderEqualToExpected(testCase, processedDataExperimentPath, testPath, expectedDataSubFolder);
  dogsPath = [processedDataExperimentPath, filesep, 'dogs'];
  expectedDataFolder = [testPath, filesep, expectedDataSubFolder, filesep, 'ProcessedData',...
    filesep, testCase.Prefix, '_', filesep, 'dogs']; 
  compareExpectedDataDir(testCase, dogsPath, expectedDataFolder);
end
