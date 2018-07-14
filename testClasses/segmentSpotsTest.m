classdef segmentSpotsTest < matlab.unittest.TestCase
    %segmentSpotsTest Runs unit testing on segmentSpots using the sample P2P
    %data set. Note that it does acccept other data sets, but certain tests
    %will fail. 
       
    properties
        Prefix = '';            
        dirSpots = '';  
        dirDog = '';
       % nFrames = 2;
    end
    
    methods(Test)
        function testCase = segmentSpotsTest(pref)
            testCase.Prefix = pref;
          %  testCase.nFrames = N;
          [~,FISHPath,DropboxFolder,~,~]= DetermineLocalFolders(testCase.Prefix);
            testCase.dirSpots = [DropboxFolder,filesep,testCase.Prefix];
            testCase.dirDog = [FISHPath, filesep, testCase.Prefix, '_'];
        end
       
        function testDoGGeneration(testCase)
            %Verifies that segmentSpots actually created dog files.
            nFrames = 2;
            delete(testCase.dirDog);
            segmentSpots(testCase.Prefix, [], 'Frames', nFrames);
            testCase.assertTrue(logical(exist(testCase.dirDog, 'dir')))
        end
        function testSpotsGeneration(testCase)
            %Verifies that segmentSpots makes Spots.mat
            nFrames = 2;
            spotsFile = [testCase.dirSpots,filesep,'Spots.mat'];
            delete(spotsFile);
            segmentSpots(testCase.Prefix, 10, 'Frames', nFrames);
            testCase.assertTrue(logical(exist(spotsFile, 'file')));
        end
        function testLogGeneration(testCase)
            %Verifies that segmentSpots makes log.mat
            nFrames = 2;
            logFile = [testCase.dirSpots,filesep, 'log.mat'];
            delete(logFile);
            segmentSpots(testCase.Prefix, 10, 'Frames', nFrames);
            testCase.assertTrue(logical(exist(logFile, 'file')));
        end
        function testSpotLength(testCase)
            %Verifies segmentSpots makes a Spots structure with the correct
            %length.
            nFrames = 2;
            spotsFile = [testCase.dirSpots,filesep,'Spots.mat'];
            delete(spotsFile);
            segmentSpots(testCase.Prefix, 10, 'Frames', nFrames);
            in = load([testCase.dirSpots, filesep, 'Spots.mat']);
            testCase.assertTrue(length(in.Spots) == nFrames);           
        end
        function testNumberDetectedCircles(testCase)
            %Verifies segmentSpots makes a Spots structure with the correct
            %number of detections (if the algorithm changes significantly,...
            %this test will fail.)
            nFrames = 2;
            logFile = [testCase.dirSpots,filesep, 'log.mat'];
            delete(logFile);
            segmentSpots(testCase.Prefix, 10, 'Frames', nFrames);
            in = load([testCase.dirSpots,filesep, 'log.mat']);
            testCase.assertTrue(in.log.totalCircles == 72);      
        end
    end
end