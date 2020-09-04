% This is a standalone script to visualize Dorsal cluster traces and
% compare them to both the raw images and the MS2 spot traces

clear
close all

%% Load all the relevant info and data
Prefix = '2019-11-26-2xDl_Venus_snaBAC_MCPmCherry_Leica_Zoom45_21uW14uW_01';
liveExperiment = LiveExperiment(Prefix);

nCh = numel(liveExperiment.spotChannels);
channels = liveExperiment.Channels;
channelNames = cell(1,nCh);

xDim = liveExperiment.xDim;
yDim = liveExperiment.yDim;
zDim = liveExperiment.zDim;
nFrames = liveExperiment.nFrames;

Spots = getSpots(liveExperiment);
if ~iscell(Spots)% NL: added for backwards compatibility
  Spots = {Spots};
end

dropboxFolder = liveExperiment.userResultsFolder;
resultsFolder = liveExperiment.resultsFolder;
preProcFolder = liveExperiment.preFolder;
trackFigFolder = [dropboxFolder filesep Prefix '\TrackingFigures\'];
ParticlesFull = getParticlesFull(liveExperiment);


%% Grab and process the raw images

load([trackFigFolder, filesep, 'rawImagesZMaxProj.mat'])

% rawImDir{1} = dir([preProcFolder filesep '*ch01.tif']);
% rawImDir{2} = dir([preProcFolder filesep '*ch02.tif']);
% 
% zMaxProjArray{1} = zeros(xDim,yDim,nFrames);
% zMaxProjArray{2} = zeros(xDim,yDim,nFrames);
% 
% for channel = 1:nCh
%           
%     for n = 1:nFrames
%         tic
%         currImPath = [rawImDir{1,channel}(n).folder filesep rawImDir{1,channel}(n).name];
%         
%         %using bfopen is slow, but fits into 2 lines  - figure out a faster
%         %way to do this with the Tiff class
%         evalc('currIm = bfOpen3DVolume(currImPath);');   %using evalc to suppress fprint statement inside bfopen
%         imStack = currIm{1,1}{1,1};     %this is the xDim x yDim x zDim image matrix
%         
% %         imZSumProj = sum(imStack,3);
% %         figure(1)
% %         imshow(imZSumProj,[])
%     
%         imZMaxProj = max(imStack,[],3);
% %         figure(2)
% %         imshow(imZMaxProj,[])
%         zMaxProjArray{channel}(:,:,n) = imZMaxProj;
%         toc
%     end
% end
% 
% save([trackFigFolder, filesep, 'rawImagesZMaxProj.mat'],'zMaxProjArray')


%% Make the plots


% Plot raw detection events on top of max projected Dl cluster images
nCh = 1;
for Channel = 1:nCh
    %Get names of channels for labeling the plots
    colonPos = strfind(channels{Channel},':');
    if isempty(colonPos)
      channelNames{Channel} = channels{Channel};
    else
      channelNames{Channel} = channels{Channel}(1:colonPos(1)-1);
    end
      
    for n = 1:nFrames
        [spotsX, spotsY, ~] = SpotsXYZ(Spots{Channel}(n));
        
        rawDetectionFig = figure('Position',[0 0 1024 1024]);
        imshow(zMaxProjArray{1}(:,:,n),[])
        hold on
        scatter(spotsX,spotsY,4,'g','o');
        xlabel('x position (pixels)')
        ylabel('y position (pixels)')
        title(['Channel ' num2str(Channel) ': ' channelNames{Channel}])
        set(gca,'Fontsize',14)
        xlim([0 xDim])
        ylim([0 yDim])
        saveas(rawDetectionFig,[trackFigFolder 'detection_events_ch' num2str(Channel) '_t' num2str(n,'%02d') '.png'])
        hold off

        pause(2)
        close
    end
end