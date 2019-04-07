function [SkipStitchSchnitz, ExpandedSpaceTolerance, NoBulkShift] = DetermineTrackNucleiOptions(varargin)
%DETERMINETRACKNUCLEIOPTIONS Processes varargin for TrackNuclei
%   See description of possible input arguments in TrackNuclei
    
    varargin = varargin{1};

    SkipStitchSchnitz=true;
    ExpandedSpaceTolerance = 1;
    NoBulkShift = false;
    
    i = 1;
    while i <= length(varargin)
        if strcmpi(varargin{i},'stitchschnitz')
            SkipStitchSchnitz=0;
            i = i + 1;
        elseif strcmpi(varargin{i}, 'ExpandedSpaceTolerance')
            ExpandedSpaceTolerance = str2double(varargin{i+1});
            i = i + 2;
        elseif strcmpi(varargin{i}, 'NoBulkShift')
            NoBulkShift = true;
            i = i + 1;
        else
            error('Input parameter not recognized')
        end
    end
end

