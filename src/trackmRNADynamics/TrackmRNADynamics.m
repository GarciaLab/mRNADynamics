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
% 'displayFigures': display figures while tracking
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

[~, ~, DefaultDropboxFolder, ~, ~] = DetermineLocalFolders;

[Prefix, app, retrack, optionalResults, displayFigures] = parseTrackmRNADynamicsArguments(DefaultDropboxFolder, varargin{:});

% Get the actual folder now that we have the Prefix


[~, ~, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(Prefix, optionalResults);


% What type of experiment are we dealing with? Get this out of MovieDatabase
[~, ExperimentType, ~, ~, ~, ~, Channel1, Channel2, ~, ~, ~, ~, ~, ...
    nc9, nc10, nc11, nc12, nc13, nc14, CF, Channel3,...
    prophase, metaphase, anaphase, DVResolution] = ...
    getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);

spotChannels = getCoatChannel(Channel1, Channel2, Channel3);
NCh = length(spotChannels);

% Set the destination folder
OutputFolder = [DropboxFolder, filesep, Prefix];

% Determine the search radius based on the imaging conditions
SearchRadiusMicrons = 3; % Search radius in um

% Load the information about this image
% Check if we have FrameInfo otherwise try to get the information straight
% from the file.
[FrameInfo, PixelSize] = loadFrameInfo(OutputFolder, PreProcPath, Prefix);

SearchRadius = ceil(SearchRadiusMicrons / PixelSize);

% Check if we have tracked the lineages of the nuclei
if exist([DropboxFolder, filesep, Prefix, filesep, Prefix, '_lin.mat'], 'file')
    UseHistone = 1;
    
    % Load the nuclei segmentation information
    load([DropboxFolder, filesep, Prefix, filesep, 'Ellipses.mat'], 'Ellipses')
    
    % Load the nuclear tracking information
    load([DropboxFolder, filesep, Prefix, filesep, Prefix, '_lin.mat'], 'schnitzcells')
    
else
    UseHistone = 0;
    schnitzcells = [];
    Ellipses = [];
    warning('Warning: No nuclei lineage tracking found. Proceeding with tracking particles only.')
end

validateExperimentTypeSupported(ExperimentType);

Particles = loadParticlesAndSelectForRetracking(OutputFolder, NCh, retrack);

[Spots, SpotFilter] = loadSpotsAndCreateSpotFilter(DropboxFolder, Prefix, NCh);

if displayFigures
    [ParticlesFig, particlesAxes, NucleiFig, nucAxes] = generateTrackingFigures(app, UseHistone);
else
    ParticlesFig = []; particlesAxes = []; NucleiFig = []; nucAxes = [];
end

[Particles, SpotFilter] = performTracking(Particles, schnitzcells, NCh, Spots, app, SpotFilter, PreProcPath, ...
    Prefix, UseHistone, ParticlesFig, spotChannels, NucleiFig, particlesAxes, nucAxes, Ellipses, ...
    PixelSize, SearchRadius, ExperimentType, FrameInfo, retrack, displayFigures);

mkdir([OutputFolder, filesep]);

save([OutputFolder, filesep, 'Particles.mat'], 'Particles', 'SpotFilter', '-v7.3');

createFieldNCAndSaveFrameInfo(FrameInfo, OutputFolder, nc9, nc10, nc11, nc12, nc13, nc14);

end