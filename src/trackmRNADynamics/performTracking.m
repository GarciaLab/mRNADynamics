function [Particles, SpotFilter] = performTracking(Prefix,useHistone,varargin)
close all force
% Process user options
% searchRadiusMicrons = 5; by default
[~,maxSearchRadiusMicrons,retrack,noRetrack,displayFigures] = ...
            determinePerformTrackingOptions(varargin);

% Get all the required data for this Prefix
liveExperiment = LiveExperiment(Prefix);
matchCostMax = 3;
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

tic
disp('Inferring particle motion model...')
[HMMParticles, globalMotionModel] = track02TrainGHMM(RawParticles, [], displayFigures);
toc

% save motion model
save([resultsFolder, filesep, 'globalMotionModel.mat'],'globalMotionModel');

tic
disp('Simulating particle tracks...')
SimParticles = track03PredictParticlePaths(HMMParticles, FrameInfo, displayFigures);
toc

% generate master particles structure
ParticlesFull.RawParticles = RawParticles;
ParticlesFull.HMMParticles = HMMParticles;
ParticlesFull.SimParticles = SimParticles;

tic
disp('Stitching particle tracks...')
[ParticlesFull, ParticleStitchInfo, SpotFilter] = track04StitchTracks(ParticlesFull, SpotFilter, ReviewedParticlesFull, ParticleStitchInfo, Prefix,...
                                useHistone, retrack, displayFigures);
toc 

matchCostVec = determineMatchOptions(Prefix,useHistone,matchCostMax);
for Channel = 1:NCh
  Particles{Channel} = dynamicStitchBeta(ParticlesFull.FullParticles{Channel},ParticlesFull.SimParticles{Channel},...
    ParticleStitchInfo{Channel},Prefix,matchCostVec(Channel));
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

          % show simulated paths 
          f3 = figure('Position',[0 0 1024 1024]);
          nPlot = min([length(ParticlesFull.SimParticles{Channel}), 25]);
          rng(123);
          errIndices = randsample(1:length(ParticlesFull.SimParticles{Channel}),nPlot,false);  
          hold on
          scatter([ParticlesFull.RawParticles{Channel}.xPos],[ParticlesFull.RawParticles{Channel}.yPos],4,'k','filled','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',0);
          for i = errIndices 
            x = ParticlesFull.SimParticles{Channel}(i).hmmModel(1).pathVec;
            y = ParticlesFull.SimParticles{Channel}(i).hmmModel(2).pathVec;
            dx = ParticlesFull.SimParticles{Channel}(i).hmmModel(1).sigmaVec;
            dy = ParticlesFull.SimParticles{Channel}(i).hmmModel(2).sigmaVec;
            r = mean([dx' dy'],2) * ones(numel(x),1);
            viscircles([x,y],r,'Color',[0.3 0.3 0.3 0.05]);
          end
          for i = 1:length(ParticlesFull.RawParticles{Channel})
            plot(ParticlesFull.RawParticles{Channel}(i).xPos,ParticlesFull.RawParticles{Channel}(i).yPos,'k');
          end
          for i = errIndices
            plot(ParticlesFull.RawParticles{Channel}(i).xPos,ParticlesFull.RawParticles{Channel}(i).yPos,'LineWidth',1.5);
          end

          xlabel('x position (pixels)')
          ylabel('y position (pixels)')
          title(['Channel ' num2str(Channel) ': ' channelNames{Channel}])
          set(gca,'Fontsize',14)
          xlim([0 xDim])
          ylim([0 yDim])
          saveas(f3,[trackFigFolder 'projected_paths_ch' num2str(Channel) '.png'])

          f4 = figure('Position',[0 0 856 856]);
          hold on  
          for i = 1:length(Particles{Channel})
%             extantFilter = min(Particles{Channel}(i).Frame):max(Particles{Channel}(i).Frame);
            plot(Particles{Channel}(i).xPos,Particles{Channel}(i).yPos,'LineWidth',1.25);
          end
          scatter([ParticlesFull.RawParticles{Channel}.xPos],[ParticlesFull.RawParticles{Channel}.yPos],4,'k','filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0);
          % for i = 1:length(Particles{Channel})
          %   plot(Particles{Channel}(i).xPos,Particles{Channel}(i).yPos);
          % end
          xlabel('x position (pixels)')
          ylabel('y position (pixels)')
          title(['Channel ' num2str(Channel) ': ' channelNames{Channel}])
          set(gca,'Fontsize',14)
          xlim([0 xDim])
          ylim([0 yDim])
          saveas(f4,[trackFigFolder 'final_paths_ch' num2str(Channel) '.png'])
    end
  
end

disp('Done.')
end