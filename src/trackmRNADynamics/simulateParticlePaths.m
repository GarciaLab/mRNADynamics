function [pathMean, pathSE, positionArray] = simulateParticlePaths(stateArray,Mu,Sigma)

% Function to gerate expected mean particle path, as well as estimated
% error in path projection. Optionally can return all simulated paths

% convert states to motion and motion variances
Sigma = sqrt(reshape(Sigma,[],1));
motionArray = Mu(stateArray);
sigmaArray = Sigma(stateArray);

% simulate jumps
jumpArray = normrnd(motionArray,sigmaArray);
positionArray = cumsum(jumpArray,2);

% find average and standard error
pathMean = mean(positionArray);
pathSE = std(positionArray);