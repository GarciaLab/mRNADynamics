function [SkipStitchSchnitz, ExpandedSpaceTolerance, NoBulkShift, retrack, nWorkers] = DetermineTrackNucleiOptions(varargin)
%DETERMINETRACKNUCLEIOPTIONS Processes varargin for TrackNuclei,
%   See description of possible input arguments in TrackNuclei
    

    SkipStitchSchnitz=true;
    ExpandedSpaceTolerance = 1;
    NoBulkShift = false;
    retrack = false;
    nWorkers = 8;
    
    for i = 1:length(varargin)
        if strcmpi(varargin{i},'stitchschnitz')
            SkipStitchSchnitz=0;
        elseif strcmpi(varargin{i}, 'ExpandedSpaceTolerance')
            ExpandedSpaceTolerance = str2double(varargin{i+1});
        elseif strcmpi(varargin{i}, 'NoBulkShift')
            NoBulkShift = true;
        elseif strcmpi(varargin{i}, 'retrack')
            retrack = true;
        elseif strcmpi(varargin{i}, 'nWorkers')
            nWorkers = varargin{i+1};
        end
    end
end

