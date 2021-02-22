function [tracks, nextId] = createNewTracks(tracks, measurements,...
    bboxes, unassignedDetections, nextId)

%ceny, cenx, smaj, smin, angle
nMeasurements = 5; %N
nStates = nMeasurements * 2; %M

measurements = measurements(...
    unassignedDetections, :);

bboxes = bboxes(unassignedDetections, :);


centroidError_pxSq = 4^2;
radiusError_pxSq = 5^2;
angleError_radSq = (pi/2)^2;
vErrorSq = 1;



InitialEstimateError = 1*1e3;
% MotionNoise = [centroidError_pxSq, centroidError_pxSq ,...
%     radiusError_pxSq];
MotionNoise = [centroidError_pxSq, 1]; %process noise
MeasurementNoise = radiusError_pxSq;
stateModel = [1 1;0 1]; 
measurementModel = [1 0];
A = [];
H = [];
for k = 1:nMeasurements
   A = blkdiag(A, stateModel); %state transition model
   H = blkdiag(H, measurementModel); %measurement model
end

MotionNoise = 0;
M = nMeasurements*2; %number of things to track along with their velocities
%state estimation covariance matrix
P = diag(repmat(InitialEstimateError, [1, M]));
%process noise covariance matrix
Q = diag(repmat(MotionNoise, [1, M]));
%measurement noise covariance matrix
R = diag([centroidError_pxSq, centroidError_pxSq,...
    radiusError_pxSq, radiusError_pxSq, angleError_radSq]);

numEllipses = size(measurements, 1);

for ellipse = 1:numEllipses
    
    bbox = bboxes(ellipse, :);
    measurement = measurements(ellipse, :);
    
    
    % Create a Kalman filter object.
    
    %     kalmanFilter = configureKalmanFilter(MotionModel,InitialLocation,InitialEstimateError,...
    %         MotionNoise,MeasurementNoise)
%     kalmanFilter = configureKalmanFilter('ConstantVelocity', ...
%         measurement, InitialEstimateError , MotionNoise,  MeasurementNoise );
%       
       kalmanFilter = vision.KalmanFilter(A,H,...
           'ProcessNoise',Q,'MeasurementNoise',R, 'StateCovariance',P,...
           'ControlModel', []);

%     
    % Create a new track.
    newTrack = struct(...
        'id', nextId, ...
        'bbox', bbox, ...
        'kalmanFilter', kalmanFilter, ...
        'age', 1, ...
        'idxHistory', nextId,...
        'totalVisibleCount', 1, ...
        'consecutiveInvisibleCount', 0);
    
    % Add it to the array of tracks.
    tracks(end + 1) = newTrack;
    
    % Increment the next id.
    nextId = nextId + 1;
end

end