classdef trackmRNADynamicsTest < matlab.unittest.TestCase
    %trackmRNADynamicsTest runs unit testing on trackmRNADynamics using the sample P2P
    %data set. Note that it does acccept other data sets, but certain tests
    %will fail. Also note that this test fails as a whole for an unknown
    %reason, but each individual test passes (AR 6/19/17). 
       
    properties
        Prefix = '';            
        dirParticles = '';  
       % nFrames = 2;
    end
    
    methods(Test)
        function testCase = trackmRNADynamicsTest(pref)
            testCase.Prefix = pref;
          %  testCase.nFrames = N;
          [~,~,DropboxFolder,~,~]= DetermineLocalFolders(pref);
            testCase.dirParticles = [DropboxFolder,filesep,pref];
        end           
        function testParticlesGeneration(testCase)
            %Verifies that trackmRNADynamics makes Particles.mat
            nFrames = 2;
            ParticlesFile = [testCase.dirParticles,filesep,'Particles.mat'];
            delete(ParticlesFile);
            TrackmRNADynamics(testCase.Prefix, 1, 1);
            testCase.assertTrue(logical(exist(ParticlesFile, 'file')));
        end
        function testParticlesLength(testCase)
            %Verifies trackmRNADynamics makes a Particles structure with the correct
            %length.
            nFrames = 2;
            ParticlesFile = [testCase.dirParticles,filesep,'Particles.mat'];
            delete(ParticlesFile);
            TrackmRNADynamics(testCase.Prefix, 1, 1);
            in = load(ParticlesFile);
            testCase.assertTrue(length(in.Particles) == 12);           
        end
    end
end