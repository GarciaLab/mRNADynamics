function copyExpectedDataFolder(sourceFolder, expectedDataFolder, pipelineStep, Prefix)
  disp(['Copying DynamicResults data for ', pipelineStep, ' and Prefix ', Prefix]);
  sourceFolder = [sourceFolder, '/*'];
  disp(['Copying ', sourceFolder, ' to ', expectedDataFolder]);
  copyfile(sourceFolder, expectedDataFolder);
  disp(['Expected data copied to folder ', expectedDataFolder]);
end