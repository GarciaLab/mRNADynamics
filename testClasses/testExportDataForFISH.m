function testCase = testExportDataForFISH(testCase)
    disp('Running with prefix ');
    disp(testCase.Prefix);

    if (~isprop(testCase, 'PreferredFileName')) 
      ExportDataForFISH(testCase.Prefix);
    else 
      ExportDataForFISH(testCase.Prefix, testCase.PreferredFileName);
    end
    compareExportResults(testCase);
end

function testCase = compareExportResults(testCase)
    %Figure out the initial folders. 
    CONFIG_CSV_PATH = ['ComputerFolders.csv'];

    configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');

    PreProcPath = getConfigValue(configValues, 'PreProcPath');
    testPath = getConfigValue(configValues, 'TestPath');

    expectedDataFolder = strcat(testPath, filesep, 'ExportDataForFISH', filesep, testCase.Prefix);

    %Get file names to compare in preprocessed data folder
    preprocessedDataFolder = strcat(PreProcPath, filesep, testCase.Prefix);
    testCase = compareExpectedDataDir(testCase, preprocessedDataFolder, expectedDataFolder);
end
