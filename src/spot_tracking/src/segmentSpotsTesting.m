% Script to compare different spot fitting methods (2D, 1spot 3D, 2spot 3D)
% Relevant criteria are 1) localization accuracy and 2) intensity
% "accuracy" (difficult to evaluate this one quantitatively)
clear
close all

% specify Prefix corresponding to project we wish to examine 
Prefix = '2019-04-30-Dl_Venus_snaBAC_MCPmCherry_Leica_Zoom2_7uW14uW_04';

% Get paths
[~,ProcPath,DropboxFolder,~, PreProcPath,...
    ~, Prefix, ~,Channel1,Channel2,~, Channel3, spotChannels] = readMovieDatabase(Prefix);
% check for existence of truncat4ed Spots structure
sp_trunc_flag = exist([DropboxFolder '/' Prefix '/SpotsTrunc.mat' ]);
% if it hasn't been made, we have to load full spots struct and create
% truncated version (sigh)
if ~sp_trunc_flag
    % Load particles set
    load([DropboxFolder '/' Prefix '/Particles.mat' ])
    load([DropboxFolder '/' Prefix '/Spots.mat' ])
    % find and keep Particles that exist for long periods of time
    ParticlesTrunc = [];
    min_len = 90;
    for i = 1:numel(Particles)
        if numel(Particles(i).zPos) >= min_len
            ParticlesTrunc = [ParticlesTrunc Particles(i)];
        end
    end

    % now generate corresponding truncated spot structure
    frame_vec = [ParticlesTrunc.Frame];
    index_vec = [ParticlesTrunc.Index];
    frame_index = unique(frame_vec);
    SpotsTrunc = struct;
    for f = 1:numel(frame_index)
        frame = frame_index(f);
        frame_indices = frame_vec == frame;
        particle_indices = index_vec(frame_indices);
        FitsTemp = [];
        for i = 1:numel(particle_indices)
            fragment = Spots(frame).Fits(particle_indices(i));
            fragment.index_orig = particle_indices(i);
            FitsTemp = [FitsTemp fragment];
        end
        SpotsTrunc(frame).Fits = FitsTemp;
    end
    % save
    save([DropboxFolder '/' Prefix '/SpotsTrunc.mat'],'SpotsTrunc');
    save([DropboxFolder '/' Prefix '/ParticlesTrunc.mat'],'ParticlesTrunc');
else
    load([DropboxFolder '/' Prefix '/SpotsTrunc.mat'],'SpotsTrunc');
    load([DropboxFolder '/' Prefix '/ParticlesTrunc.mat'],'ParticlesTrunc');
end
%% 
% Use info from truncated structures to extract 3D snips for each spot and
% fit Gaussians
SpotsTrunc = fit3DGaussiansToAllSpots(Prefix, 1,'noSave','segmentSpots',SpotsTrunc);
%%
% compare 2D and 3D trajectories
% First Pull 3D positions into Particles (makes use of pre-existing
% tracking)
for i = 1:numel(ParticlesTrunc)
    frame_vec = ParticlesTrunc(i).Frame;
    index_vec = ParticlesTrunc(i).Index;
    xPos3D = NaN(size(frame_vec));
    yPos3D = NaN(size(frame_vec));
    zPos3D = NaN(size(frame_vec));
    Fluo3D = NaN(size(frame_vec));
    Fluo3DRaw = NaN(size(frame_vec));
    Fluo2D = NaN(size(frame_vec));
    % iterate through frames
    for f = 1:numel(frame_vec)
        frame = frame_vec(f);
        index_orig = index_vec(f);
        index_ft = [SpotsTrunc(frame).Fits.index_orig] == index_orig;
        xPos3D(f) = SpotsTrunc(frame).Fits(index_ft).GaussPos3D(1);
        yPos3D(f) = SpotsTrunc(frame).Fits(index_ft).GaussPos3D(2);
        zPos3D(f) = SpotsTrunc(frame).Fits(index_ft).GaussPos3D(3);
        Fluo3D(f) = SpotsTrunc(frame).Fits(index_ft).gauss3DIntensity;        
        Fluo2D(f) = SpotsTrunc(frame).Fits(index_ft).FixedAreaIntensity3;
    end
    ParticlesTrunc(i).xPos3D = xPos3D;
    ParticlesTrunc(i).yPos3D = yPos3D;
    ParticlesTrunc(i).zPos3D = zPos3D;
    ParticlesTrunc(i).Fluo3D = Fluo3D;    
    ParticlesTrunc(i).Fluo2D = Fluo2D;
end

%%
initialFrame = min([ParticlesTrunc.Frame]);
[displayFigures, numFrames, numShadows, intScale, keepPool, ...
    autoThresh, initialFrame, useIntegralCenter, Weka, keepProcessedData,...
    fit3D, skipChannel, optionalResults, filterMovieFlag, gpu, nWorkers, saveAsMat, saveType, nuclearMask]...
    = determineSegmentSpotsOptions([]);
Ellipses = {};
% extract additional info needed to run segmentTranscriptionalLoci
[~, ~, ~, ~, ~, ~, ~, ExperimentType, Channel1, Channel2,~, ~, spotChannels] = readMovieDatabase(Prefix);

[~, ProcPath, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(Prefix);

% laod FrameInfo
load([DropboxFolder, filesep, Prefix, filesep, 'FrameInfo.mat'], 'FrameInfo');

ProcessedDataFolder = [ProcPath, filesep, Prefix, '_'];
DogOutputFolder = [ProcessedDataFolder, filesep, 'dogs'];
mkdir(DogOutputFolder)

microscope = FrameInfo(1).FileMode;

zSize = 2;
for i = 1:size(FrameInfo,2)
    if (FrameInfo(i).NumberSlices+2)>zSize
        zSize = FrameInfo(i).NumberSlices + 2;
    end
end
numFrames = length(FrameInfo);
nCh = length(spotChannels);
[ffim, doFF] = loadSegmentSpotsFlatField(PreProcPath, Prefix, FrameInfo);

% The spot finding algorithm first segments the image into regions that are
% above the threshold. Then, it finds global maxima within these regions by searching in a region "neighborhood"
% within the regions.
pixelSize = FrameInfo(1).PixelSize * 1000; %nm
neighborhood = round(1300 / pixelSize); %nm
snippet_size = 2 * (floor(1300 / (2 * pixelSize))) + 1; % nm. note that this is forced to be odd
coatChannel = spotChannels;

% Call 2D spot segmentation algorithm 
channelIndex = 1;
[tempSpots, dogs] = segmentTranscriptionalLoci(nCh, coatChannel, channelIndex, initialFrame, numFrames, zSize, ...
        PreProcPath, Prefix, DogOutputFolder, displayFigures, doFF, ffim, Threshold(channelIndex), neighborhood, ...
        snippet_size, pixelSize, microscope, intScale, Weka,...
        useIntegralCenter, filterMovieFlag, optionalResults, gpu, saveAsMat, saveType, Ellipses);
