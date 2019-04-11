function [ trackingStartingPoints ] = choseFramesToStartTracking(FrameInfo, indMitosis, numberOfFrames )
%CHOSEFRAMESTOSTARTTRACKING This function decides where the tracking should
% start. Tracking is done forward and backward through time, so any time
% between two mitosis could be chosen. A good choice corresponds to a frame
% where the nuclei are all circularly-shaped. A first approximation would
% be in the middle of the interphase.
%
% FIRST AND LAST NUCLEAR CYCLE
% For the case of the first nuclear cycle, the tracking is started at the
% midpoint between the first frame and the first mitosis. Since nuclei are
% usually not clearly visible at that time, this should be replaced by a
% simple back-tracking from the following nuclear cycle.
% For the case of the NC14, the tracking is started a defined amount of
% time after the mitosis. If there are not enough frames, the last one is
% picked.

% Get Parameters
time_resolution = getDefaultParameters(FrameInfo,'time resolution');
delta_time_tracking = getDefaultParameters(FrameInfo,'delta t','choseFramesToStartTracking')/time_resolution;

%indicesMitosis = unique([1;indMitosis]);

trackingStartingPoints = round(0.5*(indMitosis(2:end-1,1)+[1; indMitosis(2:end-2,2)]));

if size(indMitosis,1) > 1
    trackingStartingPoints(end+1) = min(numberOfFrames, round(indMitosis(end-1,2)+delta_time_tracking));
end
    
end

