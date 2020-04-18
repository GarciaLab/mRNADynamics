%[Particles,schnitzcells]=TrackmRNADynamics(varargin)
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
% 'noRetracking': Use this to track from scratch instead of retracking if
%                trackmRNADynamics has been run before.
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
function [Particles, schnitzcells] = TrackmRNADynamics(varargin)

disp(['Running TrackmRNADynamics on ', varargin{1}, '...']);

Prefix = varargin{1};
thisExperiment = liveExperiment(Prefix);

[app, retrack, optionalResults, displayFigures] =...
    parseTrackmRNADynamicsArguments(varargin{:});


DropboxFolder = thisExperiment.userResultsFolder;
PreProcPath = thisExperiment.userPreFolder;

ExperimentType = thisExperiment.experimentType;

Channels = thisExperiment.Channels;
Channel1 = Channels{1};
Channel2 = Channels{2};
Channel3 = Channels{3};

anaphaseFrames = thisExperiment.anaphaseFrames';
nc9 = anaphaseFrames(1);
nc10 = anaphaseFrames(2);
nc11 = anaphaseFrames(3);
nc12 = anaphaseFrames(4);
nc13 = anaphaseFrames(5);
nc14 = anaphaseFrames(6);


spotChannels = getCoatChannel(Channel1, Channel2, Channel3);
NCh = length(spotChannels);

% Set the destination folder
OutputFolder = [DropboxFolder, filesep, Prefix];

% Determine the search radius based on the imaging conditions
SearchRadiusMicrons = 0.5; % Search radius in um

% Load the information about this image
% Check if we have FrameInfo otherwise try to get the information straight
% from the file.
[FrameInfo, PixelSize_um] = loadFrameInfo(OutputFolder, PreProcPath, Prefix);

SearchRadius = ceil(SearchRadiusMicrons / PixelSize_um);

% Check if we have tracked the lineages of the nuclei
if exist([DropboxFolder, filesep, Prefix, filesep, Prefix, '_lin.mat'], 'file')
    UseHistone = true;
    
    % Load the nuclei segmentation information
    load([DropboxFolder, filesep, Prefix, filesep, 'Ellipses.mat'], 'Ellipses')
    
    % Load the nuclear tracking information
    load([DropboxFolder, filesep, Prefix, filesep, Prefix, '_lin.mat'], 'schnitzcells')
    
else
    UseHistone = false;
    schnitzcells = [];
    Ellipses = [];
    warning('Warning: No nuclei lineage tracking found. Proceeding with tracking particles only.')
end

validateExperimentTypeSupported(ExperimentType);

Particles = loadParticlesAndSelectForRetracking(OutputFolder, NCh, retrack);

Spots = loadSpotsAndCreateSpotFilter(DropboxFolder, Prefix, NCh);

if displayFigures
    [ParticlesFig, particlesAxes, NucleiFig, nucAxes] = generateTrackingFigures(app, UseHistone);
else
    ParticlesFig = []; particlesAxes = []; NucleiFig = []; nucAxes = [];
end

[Particles, SpotFilter] = performTracking(Particles, schnitzcells,...
    NCh, Spots, app, SpotFilter, PreProcPath, ...
    Prefix, UseHistone, ParticlesFig, spotChannels,...
    NucleiFig, particlesAxes, nucAxes, Ellipses, ...
    PixelSize_um, SearchRadius, ExperimentType,...
    FrameInfo, retrack, displayFigures, thisExperiment);

mkdir([OutputFolder, filesep]);

try
    save([OutputFolder, filesep, 'Particles.mat'], 'Particles', 'SpotFilter', '-v6');
catch
    save([OutputFolder, filesep, 'Particles.mat'], 'Particles', 'SpotFilter', '-v7.3', '-nocompression');
end

createFieldNCAndSaveFrameInfo(FrameInfo, OutputFolder, nc9, nc10, nc11, nc12, nc13, nc14);

end