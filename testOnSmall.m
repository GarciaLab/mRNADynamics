function runSmallTest()

nWorkers = 1;
Prefix = '2020-01-21-1Dg-8D_EfEfEf_9_small';

tic;
ExportDataForLivemRNA(Prefix, 'nuclearGUI', false);
filterMovie(Prefix, 'Tifs', 'nWorkers', nWorkers);
% TrackNuclei(Prefix, 'nWorkers', nWorkers);
filterMovie(Prefix, 10018, 'nWorkers', nWorkers);
segmentSpots(Prefix, 10023, 'nWorkers', nWorkers);
toc

end