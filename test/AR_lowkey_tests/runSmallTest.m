function runSmallTest()

%% Setup
nWorkers = 1;
Prefix = '2020-01-21-1Dg-8D_EfEfEf_9_small';

%% Exercise
tic;
ExportDataForLivemRNA(Prefix, 'nuclearGUI', false);
% TrackNuclei(Prefix, 'nWorkers', nWorkers);
filterMovie(Prefix, 'nWorkers', nWorkers);
segmentSpots(Prefix, 10018, 'nWorkers', nWorkers);
CompileParticles(Prefix, 'SkipAll', 'ApproveAll')
toc

%% Verify

%% Teardown

end