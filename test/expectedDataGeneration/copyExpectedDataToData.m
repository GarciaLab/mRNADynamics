function copyExpectedDataToData(Prefix, dynamicsResultsPath, dynamicsResultsExpectedPath, preprocessedDataPath, preProcessedDataExpectedPath, processedDataPath, processedDataExpectedPath)
  disp(['Copying expected data for Prefix ', Prefix]);
  
  copyfile([dynamicsResultsExpectedPath, filesep, '*'], dynamicsResultsPath);
  disp(['Expected data copied to folder ', dynamicsResultsPath]);

  copyfile([preProcessedDataExpectedPath, filesep, '*'], preprocessedDataPath);
  disp(['Expected data copied to folder ', preprocessedDataPath]);

  copyfile([processedDataExpectedPath, filesep, '*'], processedDataPath);
  disp(['Expected data copied to folder ', processedDataPath]);

end
