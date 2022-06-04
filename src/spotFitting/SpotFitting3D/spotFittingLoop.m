function SpotsFr = spotFittingLoop(FitsFr, liveExperiment, imStack, nSpots)
%
% DESCRIPTION
% Sub-function called by fit3DGaussiansToAllSpots that loops over all
% segmented spots in a single frame and fits 3D Gaussian(s) to each.
%
% INPUT ARGUMENTS
% FitsFr: 2D fits for all spots detected in a single frame
% liveExperiment: LiveExperiment instance for this particular dataset
% imStack: 3D array containing the image data for this single frame
% nSpots: number of Gaussians to fit. Should be 2 for MS2 spots (to
%         account for sister chromatids) and 1 for transcription factor
%         clusters
% 
% OPTIONS
% N/A
%
% OUPUT
%
% Author (contact): Nicholas Lammers (nlammers@berkeley.edu)
% Created: 2019-2020ish
%
% Documented by: Meghan Turner, 2022-05-31
%
    SpotsFr.Fits = FitsFr;
    for i = 1:length(SpotsFr.Fits)
        SpotsFr = fitSnip3D(SpotsFr, i, liveExperiment, imStack, nSpots);    
    end    