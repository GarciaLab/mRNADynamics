classdef exportDataForFISHTest < matlab.unittest.TestCase
    %runs the export data process and compares preprocessed data with a known result set
       
    properties (TestParameter)
        Prefix = '';
        expectedDataFolder = '';
    end
    
    methods(Test)
        function testCase = exportDataForFISHTest(dataFolder, pref)
            testCase.Prefix = pref;
            testCase.expectedDataFolder = strcat(dataFolder,filesep,pref);
        end           

        function testRun(testCase)
            %Figure out the initial folders. 
            [~,~,~,~, PreProcPath, ~, ~]=...
                DetermineLocalFolders;

            %Get file names to compare in preprocessed data folder
            preprocessedDataFolder = strcat(PreProcPath,filesep,testCase.Prefix)
            cd(preprocessedDataFolder)
            sourceFiles = dir;

            filesToCompare = {};
            filesToCompareIndex = 1;

            %Traverse the files in folder, but ignore the first two as they are . and ..
            for i = 3:length(sourceFiles)
                filesToCompare{filesToCompareIndex} = sourceFiles(i).name;
                filesToCompareIndex = filesToCompareIndex + 1;
            end

            %Run comparison
            compareCommand = '';
            sourceFile = '';
            targetFile = '';
            expectedCompareResult = 'FC: no differences encountered';
            
            for i = 1:length(filesToCompare)
                sourceFile = strcat(preprocessedDataFolder,filesep,filesToCompare(i));
                targetFile = strcat(testCase.expectedDataFolder,filesep,filesToCompare(i));
                compareCommand = strcat("fc ",sourceFile," ",targetFile);
                [status,cmdout] = system(compareCommand, '-echo');
                testCase.assertFalse(isempty(strfind(cmdout, expectedCompareResult)))
            end
        end
    end
end