function [stitchSchnitz, ExpandedSpaceTolerance, NoBulkShift, retrack, nWorkers, track] = DetermineTrackNucleiOptions(varargin)
%DETERMINETRACKNUCLEIOPTIONS Processes varargin for TrackNuclei,
%   See description of possible input arguments in TrackNuclei
    

    stitchSchnitz=true;
    ExpandedSpaceTolerance = 1.5;
    NoBulkShift = true;
    retrack = false;
    nWorkers = 8;
    track = true;
    
    for i = 1:length(varargin)
        if strcmpi(varargin{i},'skipStitchSchnitz')
           stitchSchnitz=false;
        elseif strcmpi(varargin{i}, 'ExpandedSpaceTolerance')
            ExpandedSpaceTolerance = str2double(varargin{i+1});
        elseif strcmpi(varargin{i}, 'bulkShift')
            NoBulkShift = false;
        elseif strcmpi(varargin{i}, 'retrack')
            retrack = true;
        elseif strcmpi(varargin{i}, 'nWorkers')
            nWorkers = varargin{i+1};
        elseif strcmpi(varargin{i}, 'noTrack')
            track = false;
        end
    end
end

