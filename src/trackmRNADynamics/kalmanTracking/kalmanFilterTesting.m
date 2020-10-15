% Script to experiment with procedures for using kalman filtering to:
% (1) predict and link particle fracgments
% (2) apply backward smoothing to initial trajectories to obtain estimate
% of particle position between observations
clear 
close all

% Load data
load("C:\Users\nlamm\Dropbox (Personal)\DynamicsResults\2019-03-04-Dl_Venus_snaBAC_MCPmCherry_Leica_Zoom2_7uW_14uW_02\ParticlesFull.mat")

SimParticles = ParticlesFull.SimParticles{1};
Particles = ParticlesFull.FullParticles{1};
%% 

close all
test_id = 128; 
filterSize = 1;

% extract position data 
frameVec = Particles(test_id).Frame;
framesFull = frameVec(1):frameVec(end);
posData = NaN(length(framesFull),2);
posData(ismember(framesFull,frameVec),:) = [Particles(test_id).xPos' Particles(test_id).yPos'];

% set noise parameters
MeasurementNoise = 1;
MotionNoise = MeasurementNoise*5e-3; % NL: this seems to work
InitNoise = MeasurementNoise;

% call forward filtering function
tic
KFTrack = kalmanFilterFwd(posData,InitNoise,MotionNoise,MeasurementNoise);
KFTrack = kalmanFilterBkd(KFTrack);
toc

figure;
hold on
errorbar(KFTrack.priorTrack(:,1),KFTrack.priorTrack(:,2),sqrt(KFTrack.priorTrackSE(:,1)),'both','o-')
errorbar(KFTrack.smoothedTrack(:,1),KFTrack.smoothedTrack(:,2),sqrt(KFTrack.smoothedTrackSE(:,1)),'both','-')

plot(posData(:,1),posData(:,2),'-x','Color','k')

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

