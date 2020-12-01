function Particles = TrackmRNADynamics(Prefix, varargin)
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
% makeTrackingFigures = true;

DropboxFolder = liveExperiment.userResultsFolder;
ExperimentType = liveExperiment.experimentType;

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
FrameInfo= getFrameInfo(liveExperiment);

validateExperimentTypeSupported(ExperimentType);

mkdir([OutputFolder, filesep]);
createFieldNCAndSaveFrameInfo(FrameInfo, OutputFolder, nc9, nc10, nc11, nc12, nc13, nc14);

[Particles, ~] = performTracking(Prefix,varargin{1:end}); 


% save([OutputFolder, filesep, 'Particles.mat'], 'Particles','SpotFilter');



end