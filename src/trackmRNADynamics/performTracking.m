function [Particles, SpotFilter] = performTracking(Prefix,useHistone,varargin)
close all force
% Process user options
% searchRadiusMicrons = 5; by default
[~,maxSearchRadiusMicrons,retrack,displayFigures] = ...
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
if ~iscell(Spots)% NL: added for backwards compatibility
  Spots = {Spots};
end

% create spotFilter if not retracking
if retrack && ~liveExperiment.hasParticlesFile
  error('No Particles structure found. Re-run without "retrack" option')  
end

schnitzCells = getSchnitzcells(liveExperiment);
dropboxFolder = liveExperiment.userResultsFolder;
resultsFolder = liveExperiment.resultsFolder;
  
tic
disp('Performing intitial particle linking...')
[RawParticles, SpotFilter, ParticleStitchInfo, ReviewedParticlesFull] = track01ParticleProximity(FrameInfo, Spots, schnitzCells, ...
                liveExperiment, pixelSize, maxSearchRadiusMicrons, useHistone, retrack, ...
                displayFigures);
toc

tic
disp('Inferring particle motion model...')
[HMMParticles, globalMotionModel] = track02TrainGHMM(RawParticles, retrack, displayFigures);
toc

tic
disp('Simulating particle tracks...')
SimParticles = track03PredictParticlePaths(HMMParticles, FrameInfo, retrack, displayFigures);
toc

tic
disp('Stitching particle tracks...')
[FullParticles,ParticleStitchInfo] = track04StitchTracks(SimParticles, SpotFilter, ReviewedParticlesFull, ParticleStitchInfo, Prefix,...
                                useHistone, retrack, displayFigures);
toc 

% % re-combine with previously approved particles that were excluded from
% % tracking process
% FullParticles = [ApprovedParticles FullParticles];

tic
matchCostVec = determineMatchOptions(Prefix,useHistone,matchCostMax);
for Channel = 1:NCh
  Particles = dynamicStitchBeta(FullParticles,SimParticles,ParticleStitchInfo,Prefix,matchCostVec,Channel);
end
toc

disp('Adding QC fields...')
% Iterate over all channels and generate additional QC flags

%%% flag long time gaps
Time = [FrameInfo.Time];
% dT = median(diff(Time));
timeThresh = 5*60;
%%% flag distance gaps 
distThresh = 2.5; % um
PixelSize = FrameInfo(1).PixelSize;
zSize = FrameInfo(1).ZStep;
%%% flag unlikely linkages 
if useHistone
  costThresh = repelem(0.9*matchCostMax,NCh);
else  
  costThresh = 0.95*matchCostVec;
end
%%% flag spots that are far from their assigned nuclei
ncDistPrctile = 99.5;

for Channel = 1:NCh
  if useHistone
    allDistanceVec = [Particles{Channel}.NucleusDist];  
    threshDist = prctile(allDistanceVec,ncDistPrctile);
  end
  for p = 1:length(Particles{Channel})
    % flag cases when particle is far away from nearest nucleus
    if useHistone
      Particles{Channel}(p).ncDistFlags = Particles{Channel}(p).NucleusDist>threshDist;
    else
      Particles{Channel}(p).ncDistFlags = false(size(Particles{Channel}(p).Frame));
    end
    %%% distance
    dx = diff(Particles{Channel}(p).xPos)*PixelSize;
    dy = diff(Particles{Channel}(p).yPos)*PixelSize;
    dz = diff(Particles{Channel}(p).zPosDetrended)*zSize;
    dr = sqrt(dx.^2+dy.^2+dz.^2)>distThresh;
    Particles{Channel}(p).distShiftFlags = [false dr] | [dr false];
    %%% time 
    dt = diff(Time(Particles{Channel}(p).Frame)) > timeThresh;
    Particles{Channel}(p).timeShiftFlags = [false dt] |  [dt false];
    % flag unlikely linkages
    Particles{Channel}(p).linkCostFlags = Particles{Channel}(p).linkCostCell>costThresh;    
    Particles{Channel}(p).costThresh = costThresh;
  end    
end

% make figures if desired
if displayFigures
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
          scatter([RawParticles{Channel}.xPos],[RawParticles{Channel}.yPos],4,'k','filled');
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
          scatter([RawParticles{Channel}.xPos],[RawParticles{Channel}.yPos],4,'k','filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',0);
          for i = 1:length(RawParticles{Channel})
            plot(RawParticles{Channel}(i).xPos,RawParticles{Channel}(i).yPos,'LineWidth',1.25);
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
          nPlot = min([length(SimParticles{Channel}), 25]);
          rng(123);
          errIndices = randsample(1:length(SimParticles{Channel}),nPlot,false);  
          hold on
          scatter([RawParticles{Channel}.xPos],[RawParticles{Channel}.yPos],4,'k','filled','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',0);
          for i = errIndices 
            x = SimParticles{Channel}(i).hmmModel(1).pathVec;
            y = SimParticles{Channel}(i).hmmModel(2).pathVec;
            dx = SimParticles{Channel}(i).hmmModel(1).sigmaVec;
            dy = SimParticles{Channel}(i).hmmModel(2).sigmaVec;
            r = mean([dx' dy'],2);
            viscircles([x',y'],r,'Color',[0.3 0.3 0.3 0.05]);
          end
          for i = 1:length(RawParticles{Channel})
            plot(RawParticles{Channel}(i).xPos,RawParticles{Channel}(i).yPos,'k');
          end
          for i = errIndices
            plot(RawParticles{Channel}(i).xPos,RawParticles{Channel}(i).yPos,'LineWidth',1.5);
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
            extantFilter = min(Particles{Channel}(i).Frame):max(Particles{Channel}(i).Frame);
            plot(Particles{Channel}(i).xPos,Particles{Channel}(i).yPos,'LineWidth',1.25);
          end
          scatter([RawParticles{Channel}.xPos],[RawParticles{Channel}.yPos],4,'k','filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0);
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

% save
save([resultsFolder, filesep, 'ParticlesFull.mat'],'RawParticles','HMMParticles', 'SimParticles','FullParticles')
save([resultsFolder, filesep, 'ParticleStitchInfo.mat'],'ParticleStitchInfo');
save([resultsFolder, filesep, 'globalMotionModel.mat'],'globalMotionModel');
end