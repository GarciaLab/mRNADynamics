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

    expectedDataFolder = strcat(testPath,filesep,testCase.Prefix);

    %Get file names to compare in preprocessed data folder
    preprocessedDataFolder = strcat(PreProcPath,filesep,testCase.Prefix);
    disp(preprocessedDataFolder);
    cd(preprocessedDataFolder);
    sourceFiles = dir;

    numberOfFiles = length(sourceFiles)-2;
    filesToCompare = {numberOfFiles};
    filesToCompareIndex = 1;

    %Traverse the files in folder, but ignore the first two as they are . and ..
    for i = 3:length(sourceFiles)
        filesToCompare{filesToCompareIndex} = sourceFiles(i).name;
        filesToCompareIndex = filesToCompareIndex + 1;
    end

    %Run comparison
    expectedCompareResult = 'FC: no differences encountered';
    
    for i = 1:length(filesToCompare)
        sourceFile = strcat(preprocessedDataFolder,filesep,filesToCompare(i));
        targetFile = strcat(expectedDataFolder,filesep,filesToCompare(i));
        compareCommand = strcat("fc ",sourceFile," ",targetFile);
        [status,cmdout] = system(compareCommand, '-echo');
        disp(status);
        testCase.assertFalse(~contains(cmdout, expectedCompareResult))
    end
end
