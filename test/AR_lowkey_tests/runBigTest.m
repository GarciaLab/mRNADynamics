function runBigTest()

%% Setup
nWorkers = 1;
Prefix = '2020-01-21-1Dg-8D_EfEfEf_9_sandbox';

%% Exercise

tic;
ExportDataForLivemRNA(Prefix, 'nuclearGUI', false);
TrackNuclei(Prefix, 'nWorkers', nWorkers);
filterMovie(Prefix, 'nWorkers', nWorkers);
segmentSpots(Prefix, 10030, 'nWorkers', nWorkers);
TrackmRNADynamics(Prefix);
CompileParticles(Prefix, 'SkipAll', 'ApproveAll')
toc

%redo it parallelized
nWorkers = 12;

tic;
ExportDataForLivemRNA(Prefix, 'nuclearGUI', false);
% filterMovie(Prefix, 'Tifs', 'nWorkers', nWorkers);
filterMovie(Prefix, 'nWorkers', nWorkers);
segmentSpots(Prefix, 10030, 'nWorkers', nWorkers);
CompileParticles(Prefix, 'SkipAll', 'ApproveAll')
toc

%% Verify

%% Teardown
end