function [Particles, SpotFilter] = performTracking(Prefix,useHistone,varargin)
close all force
% Process user options
% searchRadiusMicrons = 5; by default
[~,maxSearchRadiusMicrons,retrack,noRetrack,displayFigures] = ...
            determinePerformTrackingOptions(varargin);

% Get all the required data for this Prefix
liveExperiment = LiveExperiment(Prefix);
matchCostMax = 1;
NCh = numel(liveExperiment.spotChannels);   %Only grabbing the spot channels - might cause issues if spot channels aren't the first n channels
pixelSize = liveExperiment.pixelSize_um;    %NL: pixel size is in um
channels = liveExperiment.Channels;
channelNames = cell(1,NCh);                 %Will fill this in the plotting loop
FrameInfo = getFrameInfo(liveExperiment);

% load Spots file
disp('loading Spots mat...')
Spots = getSpots(liveExperiment);
if isempty(Spots)
  error('No Spots file found. Have you run segmentSpots?')
end
if ~iscell(Spots)% NL: added for backwards compatibility
  Spots = {Spots};
end

% handle tracking/retracking options that require user input
if noRetrack && retrack
  error('conflicting retracking options specified')
elseif ~noRetrack
  retrack =  handleTrackingPrompts(liveExperiment,Spots,retrack);
end

schnitzCells = getSchnitzcells(liveExperiment);
dropboxFolder = liveExperiment.userResultsFolder;
resultsFolder = liveExperiment.resultsFolder;
  
tic
disp('Performing intitial particle linking...')
[RawParticles, SpotFilter, ParticleStitchInfo, ReviewedParticlesFull, FrameInfo] = track01ParticleProximity(FrameInfo, Spots, schnitzCells, ...
                liveExperiment, pixelSize, maxSearchRadiusMicrons, useHistone, retrack, ...
                displayFigures);
toc

% Iterate over all channels
h = waitbar(0, 'Performing tracking ...');
for Channel = 1:NCh
    
    % Iterate over all frames

    for CurrentFrame = 1:length(Spots{Channel})
        waitbar(CurrentFrame/length(Spots{Channel}),h,'Performing tracking ...');
        if isempty(app) && displayFigures
            figure(ParticlesFig)
            set(ParticlesFig, 'units', 'normalized', 'position', [0.01, .55, .33, .33]);
        end
        
        % Get the filter for this frame
        CurrentFrameFilter = logical(SpotFilter{Channel}(CurrentFrame,...
            ~isnan(SpotFilter{Channel}(CurrentFrame, :))));
        
        xPos = displayParticlesFigure(app, particlesAxes,...
            ParticlesFig, Spots, Channel, CurrentFrame, ...
            CurrentFrameFilter, PreProcPath, Prefix, SpotsChannel, FrameInfo, displayFigures);
        
        if UseHistone
            % we supply ZStep as Z axis resolution. this is unused by most people, unless you are dealing
            % with 3D nuclei tracking (for example, importing it from Imaris)
            [Particles, SpotFilter] = trackParticlesBasedOnNuclei(PreProcPath, Prefix,...
                CurrentFrame, NDigits, app, nucAxes, Ellipses, ...
                ExperimentType, Channel, schnitzcells, Particles, Spots,...
                SpotFilter, PixelSize, FrameInfo(1).ZStep, SearchRadius, retrack, displayFigures);
        else
            [Particles] = trackParticlesBasedOnProximity(Particles, Spots,...
                xPos, SpotFilter, Channel, CurrentFrame, PixelSize, SearchRadius,...
                retrack, displayFigures);
        end
        
    end
    
end
% warning('MT: Skipping dynamicsStitchBeta because it''s broken')
% Particles = ParticlesFull.FullParticles;

% Add QC flags
Particles = addQCFields(Particles,useHistone,FrameInfo,retrack,liveExperiment);


% save
disp('Saving...')
save([resultsFolder, filesep, 'ParticlesFull.mat'],'ParticlesFull')
save([resultsFolder, filesep, 'ParticleStitchInfo.mat'],'ParticleStitchInfo');


% make figures if desired
if displayFigures
    disp('Displaying figures...')
    for Channel = 1:NCh
          %Get names of channels for labeling the plots
          colonPos = strfind(channels{Channel},':');
          if isempty(colonPos)
              channelNames{Channel} = channels{Channel};
          else
              channelNames{Channel} = channels{Channel}(1:colonPos(1)-1);
          end
          
          %Make the plots
          trackFigFolder = [dropboxFolder filesep Prefix '\TrackingFigures\'];
          mkdir(trackFigFolder)
          xDim = FrameInfo(1).PixelsPerLine;
          yDim = FrameInfo(1).LinesPerFrame;
          % figure showing just detection points
          f0 = figure('Position',[0 0 1024 1024]);
          hold on
          scatter([ParticlesFull.RawParticles{Channel}.xPos],[ParticlesFull.RawParticles{Channel}.yPos],4,'k','filled');
          xlabel('x position (pixels)')
          ylabel('y position (pixels)')
          title(['Channel ' num2str(Channel) ': ' channelNames{Channel}])
          set(gca,'Fontsize',14)
          xlim([0 xDim])
          ylim([0 yDim])
          saveas(f0,[trackFigFolder 'detection_events_ch' num2str(Channel) '.png'])

          % figure showing initial linkages
          f1 = figure('Position',[0 0 1024 1024]);
          hold on  
          scatter([ParticlesFull.RawParticles{Channel}.xPos],[ParticlesFull.RawParticles{Channel}.yPos],4,'k','filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',0);
          for i = 1:length(ParticlesFull.RawParticles{Channel})
            plot(ParticlesFull.RawParticles{Channel}(i).xPos,ParticlesFull.RawParticles{Channel}(i).yPos,'LineWidth',1.25);
          end
          xlabel('x position (pixels)')
          ylabel('y position (pixels)')
          title(['Channel ' num2str(Channel) ': ' channelNames{Channel}])
          set(gca,'Fontsize',14)
          xlim([0 xDim])
          ylim([0 yDim])
          saveas(f1,[trackFigFolder 'initial_linkages_ch' num2str(Channel) '.png'])
    end

close(h)
end