function [Particles, schnitzcells] = TrackmRNADynamics(Prefix, varargin)
%
% DESCRIPTION
% %This function tracks transcription loci over time after
% segmentation and z-tracking have been performed. If nuclei have been
% tracked it uses that information for the particle tracking.
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
% [Options]: See below.
%
% OPTIONS
% 'app': This is used exclusively by the livemRNAApp to change the location
% of the display figures
%
% 'retrack': use this to retrack from previously generated
%              trackmRNADynamics results. The default is tracking from
%              scratch every time this function is called
%
% 'displayFigures': don't display figures while tracking
%
%
% OUTPUT
% Particles.mat : List of time traces found in the movie.
% *_lin.mat : Schnitzcell nuclear tracking is modified here.
%
% Author (contact): Hernan Garcia (hggarcia@berkeley.edu)
% Created: 01/01/2013 ish.
% Last Updated: 9/11/2018
%
% Documented by: Armando Reimer (areimer@berkeley.edu)

disp(['Running TrackmRNADynamics on ', Prefix, '...']);

liveExperiment = LiveExperiment(Prefix);
makeTrackingFigures = true;

% [app, retrack, optionalResults, displayFigures] =...
%     parseTrackmRNADynamicsArguments(varargin{:});


DropboxFolder = liveExperiment.userResultsFolder;
PreProcPath = liveExperiment.userPreFolder;

ExperimentType = liveExperiment.experimentType;

Channels = liveExperiment.Channels;
% Channel1 = Channels{1};
% Channel2 = Channels{2};
% Channel3 = Channels{3};

anaphaseFrames = liveExperiment.anaphaseFrames';
nc9 = anaphaseFrames(1);
nc10 = anaphaseFrames(2);
nc11 = anaphaseFrames(3);
nc12 = anaphaseFrames(4);
nc13 = anaphaseFrames(5);
nc14 = anaphaseFrames(6);

% Set the destination folder
OutputFolder = [DropboxFolder, filesep, Prefix];

% Load the information about this image
% Check if we have FrameInfo otherwise try to get the information straight
% from the file.
[FrameInfo, PixelSize_um] = loadFrameInfo(OutputFolder, PreProcPath, Prefix);

% Check if we have tracked the lineages of the nuclei
if exist([DropboxFolder, filesep, Prefix, filesep, Prefix, '_lin.mat'], 'file')
    useHistone = true;
else
    useHistone = false;
    warning('Warning: No nuclei lineage tracking found. Proceeding with tracking particles only.')
end

validateExperimentTypeSupported(ExperimentType);

% Particles = loadParticlesAndSelectForRetracking(OutputFolder, NCh, retrack);
% 
% Spots = loadSpotsAndCreateSpotFilter(DropboxFolder, Prefix, NCh);

% if displayFigures
%     [ParticlesFig, particlesAxes, NucleiFig, nucAxes] = generateTrackingFigures(app, useHistone);
% else
%     ParticlesFig = []; particlesAxes = []; NucleiFig = []; nucAxes = [];
% end

[Particles, SpotFilter] = performTracking(Prefix, useHistone,varargin{1:end}); 

mkdir([OutputFolder, filesep]);
save([OutputFolder, filesep, 'Particles.mat'], 'Particles','SpotFilter');

createFieldNCAndSaveFrameInfo(FrameInfo, OutputFolder, nc9, nc10, nc11, nc12, nc13, nc14);

end