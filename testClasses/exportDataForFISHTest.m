classdef exportDataForFISHTest < matlab.unittest.TestCase
    %runs the export data process and compares preprocessed data with a known result set
       
    properties
        Prefix = '';            
        dirName = '2015-07-25-P2P_75uW_bi'; 
        expectedDataFolder = 'd:/Documents/Berkeley/PreProcessedData/2015-07-25-P2P_75uW_bi';
    end
    
    methods(Test)
        function testCase = exportDataForFISHTest()
            %[~,~,DropboxFolder,~,~]= DetermineLocalFolders(pref);
            %testCase.dirKnownData = [DropboxFolder,filesep,pref];
            %testCase.dirPreprocessedData = [DropboxFolder,filesep,pref];


            %Figure out the initial folders. 
            [SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath, configValues, movieDatabasePath]=...
                DetermineLocalFolders;

            %Get file names to compare in preprocessed data folder
            preprocessedDataFolder = [PreProcPath,filesep,testCase.dirName]
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
        function testRun(testCase)
            cd(PreProcPath);
            sourceFiles = dir;
            hello = 'hello';
            disp(hello);
            disp(sourceFiles);

            testCase.assertTrue(1);
        end
    end
end