function [stitchSchnitz, ExpandedSpaceTolerance,...
    NoBulkShift, retrack, nWorkers, track...
    noBreak, noStitch, fish, markandfind, intFlag, useMultithresh]...
    = DetermineTrackNucleiOptions(varargin)
%
%DETERMINETRACKNUCLEIOPTIONS Processes varargin for TrackNuclei,
%   See description of possible input arguments in TrackNuclei
    

    stitchSchnitz=true;
    ExpandedSpaceTolerance = 1.5;
    NoBulkShift = true;
    retrack = false;
    nWorkers = 8;
    track = true;
    noBreak = false;
    noStitch = false;
    fish = false;
    markandfind =  false;
    intFlag = false;
    useMultithresh=false;
    
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
        elseif strcmpi(varargin{i}, 'noBreak')
            noBreak = true;
        elseif strcmpi(varargin{i}, 'noStitch')
            noStitch = true;
        elseif strcmpi(varargin{i}, 'markandfind') | strcmpi(varargin{i}, 'fish')
            noStitch = true;
            noBreak = true;
            track = true;
            fish = true;
            markandfind = true;
        elseif strcmpi(varargin{i}, 'integrate')
            intFlag = true;
        elseif strcmpi(varargin{i}, 'useMultithresh')
            useMultithresh=true;
        end
    end
    
    startParallelPool(nWorkers, 0,1);

end

