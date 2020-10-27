function [useHistone,retrack,noRetrack,displayFigures] = ...
                                  determinePerformTrackingOptions(options)
        
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
retrack = false;
displayFigures = false;
noRetrack = false;

for i= 1:length(options)
    if strcmpi(options{i}, 'useHistone')
        useHistone = true;
    elseif strcmpi(options{i}, 'retrack')
        retrack = true;
    elseif strcmpi(options{i}, 'noRetrack')
        noRetrack = true;
    elseif strcmpi(options{i}, 'displayFigures')
        displayFigures = true;
    end
end