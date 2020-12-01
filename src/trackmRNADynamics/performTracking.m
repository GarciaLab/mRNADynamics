function [Particles, SpotFilter] = performTracking(Prefix,varargin)
close all force

% initialize structure to store tracking options
trackingOptions = struct;   

% Process user options
[trackingOptions.useHistone, trackingOptions.retrack,trackingOptions.noRetrack,...
 trackingOptions.displayFigures, trackingOptions.use3DInfo] = ...
                                determinePerformTrackingOptions(varargin);
                           
% trackingOptions.useHistone = useHistone; 

% Get all the required data for this Prefix
liveExperiment = LiveExperiment(Prefix);

% load Spots file
disp('loading Spots mat...')
Spots = getSpots(liveExperiment);
if isempty(Spots)
  error('No Spots file found. Have you run segmentSpots?')
end
if ~iscell(Spots)% NL: added for backwards compatibility
  Spots = {Spots};
end
trackingOptions.has3DInfo =  isfield(Spots{1}(1).Fits,'GaussPos3D');
trackingOptions.use3DInfo = trackingOptions.use3DInfo && trackingOptions.has3DInfo;

% handle tracking/retracking options that require user input
if trackingOptions.noRetrack && trackingOptions.retrack
  error('Conflicting retracking options specified')
elseif ~trackingOptions.noRetrack
  trackingOptions.retrack =  handleTrackingPrompts(liveExperiment,Spots,trackingOptions.retrack);
end

% reset histone option to 0 if we have no nucleus data
trackingOptions.useHistone = trackingOptions.useHistone && liveExperiment.hasSchnitzcellsFile;

% process key tracking options
% Set max spots per nucleus per frame, can be different between channels
trackingOptions = parseTrackingOptions(Spots, liveExperiment, trackingOptions);

% Perform main tracking
[Particles, trackingOptions, SpotFilter] = ParticleTrackingKalman(Spots, liveExperiment, trackingOptions);

% Add QC flags
[Particles, trackingOptions] = addQCFlags(Particles, liveExperiment, trackingOptions);


% save
disp('Saving...')
save([liveExperiment.resultsFolder, filesep, 'trackingOptions.mat'],'trackingOptions')
save([liveExperiment.resultsFolder, filesep, 'Particles.mat'],'Particles','SpotFilter');

disp('Done.')
end