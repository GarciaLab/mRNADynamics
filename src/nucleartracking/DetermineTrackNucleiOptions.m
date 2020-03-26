function [stitchSchnitz, ExpandedSpaceTolerance,...
    NoBulkShift, retrack, nWorkers, track,...
    noBreak, noStitch, fish,...
    markandfind, intFlag, chooseHis, segmentBetter,...
    min_rad_um, max_rad_um]...
    = DetermineTrackNucleiOptions(varargin)
%
%DETERMINETRACKNUCLEIOPTIONS Processes varargin for TrackNuclei,
%   See description of possible input arguments in TrackNuclei

min_rad_um = 2; %good for early fly embryos
max_rad_um = 6; %good for early fly embryos
stitchSchnitz=true;
ExpandedSpaceTolerance = 1.5;
NoBulkShift = true;
retrack = false;
nWorkers = 1;
track = true;
noBreak = false;
noStitch = false;
fish = false;
markandfind =  false;
intFlag = false;
chooseHis = false;
segmentBetter = true;

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'ExpandedSpaceTolerance')
        ExpandedSpaceTolerance = str2double(varargin{i+1});
    elseif strcmpi(varargin{i}, 'bulkShift')
        NoBulkShift = false;
    elseif strcmpi(varargin{i}, 'retrack')
        retrack = true;
    elseif strcmpi(varargin{i}, 'nWorkers')
        nWorkers = varargin{i+1};
    elseif strcmpi(varargin{i}, 'noTrack')
        track = false;
    elseif strcmpi(varargin{i}, 'min_rad_um')
        min_rad_um = varargin{i+1};
    elseif strcmpi(varargin{i}, 'max_rad_um')
        max_rad_um = varargin{i+1};
        
    elseif strcmpi(varargin{i}, 'noBreak')
        noBreak = true;
    elseif strcmpi(varargin{i}, 'chooseHis')
        chooseHis = true;
    elseif strcmpi(varargin{i}, 'segmentBetter')
        segmentBetter = varargin{i+1};
    elseif strcmpi(varargin{i}, 'noStitch')
        noStitch = true;
    elseif strcmpi(varargin{i}, 'markandfind') || strcmpi(varargin{i}, 'fish')
        noStitch = true;
        noBreak = true;
        track = true;
        fish = true;
        markandfind = true;
    elseif strcmpi(varargin{i}, 'integrate')
        intFlag = true;
    end
end

startParallelPool(nWorkers, 0,1);

end

