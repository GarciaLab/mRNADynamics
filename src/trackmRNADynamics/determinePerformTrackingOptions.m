function [useHistone,searchRadiusMicrons,retrack,displayFigures] = ...
            determinePerformTrackingOptions(varargin)
        
% [useHistone,searchRadiusMicrons,retrack,displayFigures] = ...
%             determinePerformTrackingOptions(varargin)
% 
% DESCRIPTION
% Processes the user input options for performTracking.m.
%
% PARAMETERS
% varagin: user input options
%
% OPTIONS
% N/A
%
% OUTPUT
% useHistone
% searchRadiusMicrons
% retrack
% displayFigures
%
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 6/30/2020
% Last Updated: N/A

% Default options
useHistone = false;
searchRadiusMicrons = 5;
retrack = false;
displayFigures = false;

for i= 1:length(varargin)
    if strcmpi(varargin{i}, 'useHistone')
        useHistone = true;
    elseif strcmpi(varargin{i}, 'searchRadius')
        searchRadiusMicrons = varargin{i+1};
    elseif strcmpi(varargin{i}, 'retrack')
        retrack = true;    
    elseif strcmpi(varargin{i}, 'displayFigures')
        displayFigures = true;
    end
end