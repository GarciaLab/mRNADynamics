function testCase = compareExpectedDataDir(testCase, dataFolder, expectedDataFolder)
    %Get file names to compare in preprocessed data folder
    disp(dataFolder);
    cd(dataFolder);
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
        sourceFile = strcat(dataFolder, filesep, filesToCompare(i));
        targetFile = strcat(expectedDataFolder, filesep, filesToCompare(i));
        compareCommand = strcat("fc ", sourceFile, " ", targetFile);
        [status,cmdout] = system(compareCommand, '-echo');
        disp(status);
        testCase.assertFalse(~contains(cmdout, expectedCompareResult))
    end
end