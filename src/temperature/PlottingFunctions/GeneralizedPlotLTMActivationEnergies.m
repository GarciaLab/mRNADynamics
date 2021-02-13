function GeneralizedPlotLTMActivationEnergies(this, parameter, outdir, varargin)
%%

% PlotTitle, PlottingColors, UseDifferentColors,
% UseDiffProfiles, UsePhysicalAPLength
SkipParamsVsAP = false;

subfn_varargin = {};
x = 1;
while x <= length(varargin)
    if strcmp(lower(varargin{x}), 'plottitle')
        PlotTitle = varargin{x+1};
        subfn_varargin = [subfn_varargin, 'PlotTitle', PlotTitle];
        x = x+1;

    elseif strcmpi(varargin{x}, 'SkipParamsVsAP')
        SkipParamsVsAP = true;
    elseif strcmp(lower(varargin{x}), 'tracetype')
        TraceType = lower(varargin{x+1});
        x = x+1;
    end
    x = x+1;
end


if ~exist('TraceType', 'var')
    TraceType = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'Fluo3D') | strcmpi(TraceType, 'Unaligned3D') 
    TraceType = 'Unaligned3D';
elseif strcmpi(TraceType, 'Fluo')| strcmpi(TraceType, 'Unaligned') 
    TraceType = 'Unaligned';
elseif strcmpi(TraceType, 'AnaphaseAligned')
    TraceType = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'AnaphaseAligned3D')
    TraceType = 'AnaphaseAligned3D';
elseif strcmpi(TraceType, 'Tbinned')
    TraceType = 'Tbinned';
elseif strcmpi(TraceType, 'Tbinned3D')
    TraceType = 'Tbinned3D';
else
    error('Invalid choice of trace type. Can use either "fluo", "fluo3d", "anaphasealigned", or "anaphasealigned3d".') % change to error
end


if strcmpi(parameter, 'TimeOns') | strcmpi(parameter, 'TimeOn')
    parameter = 'TimeOns';
elseif strcmpi(parameter, 'TimeOffs') | strcmpi(parameter, 'TimeOff')
    parameter = 'TimeOffs';
elseif strcmpi(parameter, 'TranscriptionWindows')| strcmpi(parameter, 'TranscriptionWindow')
    parameter = 'TranscriptionWindows';
elseif strcmpi(parameter, 'ElongationTimes') | strcmpi(parameter, 'ElongationTime')
    parameter = 'ElongationTimes';
elseif strcmpi(parameter, 'ElongationRates') | strcmpi(parameter, 'ElongationRate')
    parameter = 'ElongationRates';
elseif strcmpi(parameter, 'LoadingRates') | strcmpi(parameter, 'LoadingRate')
    parameter = 'LoadingRates';
elseif strcmpi(parameter, 'PostTranscriptionDurations') | strcmpi(parameter, 'PostTranscriptionDuration')
    parameter = 'PostTranscriptionDurations';
elseif strcmpi(parameter, 'PlateauHeights') | strcmpi(parameter, 'PlateauHeight')
    parameter = 'PlateauHeights';
elseif strcmpi(parameter, 'MaxFluos') | strcmpi(parameter, 'MaxFluo')
    parameter = 'MaxFluos';
end
subfn_varargin = [subfn_varargin, 'TraceType', TraceType];


if ~SkipParamsVsAP
    PlotLTMActivationEnergiesVsAP(this, parameter, outdir, subfn_varargin{:});
end