% Script to experiment with procedures for estimating optimal kalman filter
% parameters for particle tracking
clear 
close all

% specfiy Prefix to use
Prefix = '2019-03-04-Dl_Venus_snaBAC_MCPmCherry_Leica_Zoom2_7uW_14uW_02';

% extract key exp characteristics
thisExperiment = liveExperiment(Prefix);
DropboxFolder = thisExperiment.userResultsFolder;

% load particles
[Particles, ~] = getParticles(thisExperiment);

% frame info
FrameInfo = getFrameInfo(thisExperiment);

% find and extract lon tracts of continuos detections
minLen = 20;
longTrackCell = {};
trackLenVec = [];
iter = 1;
for p = 1:length(Particles)
  Frames = Particles(p).Frame;
  posArray = NaN(length(Frames),2);
  if length(Frames)>=minLen
    posArray(:,1) = Particles(p).xPos;
    posArray(:,2) = Particles(p).yPos;
    % look for jumps
    dF = diff(Frames);
    jumpIndices = find([2 dF(1:end-1) 2]>1);
    % find long runs
    longRunFT = diff(jumpIndices)>minLen;
    if any(longRunFT)
      lrStarts = jumpIndices([longRunFT false]);
      lrStops = jumpIndices([false longRunFT])-1;
      for j = 1:length(lrStarts)
        longTrackCell{iter} = posArray(lrStarts(j):lrStops(j),:);
        trackLenVec(iter) = lrStops(j) - lrStarts(j) + 1;
        iter = iter + 1;
      end
    end    
  end
end

%% 

close all
test_id = 27; 
trackLen = trackLenVec(test_id);
filterSize = 1;
% extract position data
posData = longTrackCell{test_id};
% set random ininital guess at param values
% n0 = rand(1,4);
smoothed_z = [imgaussfilt(posData(:,1),filterSize), imgaussfilt(posData(:,2),filterSize)];
diffs = posData-smoothed_z;

mean1 = mean(smoothed_z);
mean2 = mean(diffs);
% 
% Q = mean(sum((smoothed_z-mean1).^2) / (trackLen-1))/10;
R = mean(sum((diffs-mean2).^2) / (trackLen-1));
Q = 1e-5;
noiseVec = [[R R R] [Q Q Q] R];
timePoints = 1:trackLen;
[pdTrack, ctTrack,pdTrackSE,kF] = kFilterFwd(posData(1:end-5,:),noiseVec,timePoints);

figure;
hold on
errorbar(pdTrack(end-4:end,1),pdTrack(end-4:end,2),sqrt(pdTrackSE(end-4:end,1)),'both','-o')
% plot(pdTrack(end-4:end,1),pdTrack(end-4:end,2),'-o')
% plot(pdTrack(:,1),pdTrack(:,2),'-o')
plot(ctTrack(:,1),ctTrack(:,2))
plot(posData(:,1),posData(:,2),'-','Color','k')
plot(posData(end-5:end,1),posData(end-5:end,2),'-o')

% axis([0.95*min(pdTrack(:,1)) 1.05*max(pdTrack(:,1)) 0.95*min(pdTrack(:,2)) 1.05*max(pdTrack(:,2))])

%%

% define noise characteristics
sigma = 1;
mFactor = 1e-1;%1e0;
MeasurementNoise = sigma; 
InitialEstimateError = ones(1,3)*sigma;
MotionNoise = [sigma, sigma sigma]*mFactor;

tic
% initialize particle filter
kalmanFilter =  configureKalmanFilter('ConstantAcceleration', ...
          posData(1,:), InitialEstimateError , MotionNoise,  MeasurementNoise );


posCorBuf = [];
posCorBufSE = [];
posPredBuf = [];
posPredBufSE = [];

for k = 1:trackLen 
  if k < trackLen-9
    [CorrectedMeasurement,CorrectedState,CorrectedStateCovariance] = correct(kalmanFilter,posData(k,:));  
    posCorBuf(k,:) = CorrectedState([1 4])'; 
    posCorBufSE(k,:) = sqrt([CorrectedStateCovariance(1,1) CorrectedStateCovariance(4,4)]);
  else
    posCorBuf(k,:) = NaN;
    posCorBufSE(k,:) = NaN;
  end
  % make prediction
  [PredictedMeasurement,PredictedState,PredictedStateCovariance] = predict(kalmanFilter);
  
  % record
  posPredBuf(k,:) = PredictedState([1 4])';
  posPredBufSE(k,:) = sqrt([PredictedStateCovariance(1,1) PredictedStateCovariance(4,4)]);      
end

toc
figure;
hold on
errorbar(posPredBuf(:,1),posPredBuf(:,2),posPredBufSE(:,1),'both','-o')
plot(posCorBuf(:,1),posCorBuf(:,2))
plot(posData(:,1),posData(:,2),'-o','Color','k')

axis([0.99*min(posData(:,1)) 1.01*max(posData(:,1)) 0.99*min(posData(:,2)) 1.01*max(posData(:,2))])









%%
sigma = 1e-5;

InitialEstimateError = ones(1,3)*1;
MotionNoise = [sigma, sigma sigma];
MeasurementNoise = sigma*1e5; 

xTrack = linspace(1,10,100)+rand(1,100);%Particles(TrainIDs(tr)).xPos;
yTrack = linspace(1,10,100)+rand(1,100);%yTrack = Particles(TrainIDs(tr)).yPos;

% initialize particle filter
kalmanFilter =  configureKalmanFilter('ConstantAcceleration', ...
          [xTrack(1) yTrack(1)], InitialEstimateError , MotionNoise,  MeasurementNoise );   

posCorBuf = [];
posCorBufSE = [];
posPredBuf = [];
posPredBufSE = [];

for k = 1:length(xTrack)    
  if k < length(xTrack)-10
    [CorrectedMeasurement,CorrectedState,CorrectedStateCovariance] = correct(kalmanFilter,[xTrack(k) yTrack(k)]);  
    posCorBuf(k,:) = CorrectedState([1 4])'; 
    posCorBufSE(k,:) = sqrt([CorrectedStateCovariance(1,1) CorrectedStateCovariance(4,4)]);
  else
    posCorBuf(k,:) = NaN;
    posCorBufSE(k,:) = NaN;
  end
  % make prediction
  [PredictedMeasurement,PredictedState,PredictedStateCovariance] = predict(kalmanFilter);
  
  % record
  posPredBuf(k,:) = PredictedState([1 4])';
  posPredBufSE(k,:) = sqrt([PredictedStateCovariance(1,1) PredictedStateCovariance(4,4)]);      
 
end

close all
figure;
hold on
plot(xTrack,yTrack,'-o')

xp = posPredBuf(:,1) + posPredBufSE(:,1);
xn = posPredBuf(:,1) - posPredBufSE(:,1);
yp = posPredBuf(:,2) + posPredBufSE(:,2);
yn = posPredBuf(:,2) - posPredBufSE(:,2);

errorbar(posPredBuf(:,1),posPredBuf(:,2),posPredBufSE(:,1),'both','-o')
% plot(xPredBuf(:,1),xPredBuf(:,2),'-o')
plot(posCorBuf(:,1),posCorBuf(:,2),'-o')

