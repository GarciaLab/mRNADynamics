function GeneralizedPlotLTMParameters(this, parameter, outdir, varargin)
%%

% PlotTitle, PlottingColors, UseDifferentColors,
% UseDiffProfiles, UsePhysicalAPLength
UsePhysicalAPLength = false;
UseDifferentColors = true;
UseLines = true;
SkipParamsVsAP = false;
SkipSingleTempParamsVsAP = false;
SkipParamsVsTemp = false;
SkipBinnedParamsVsAP = false;
SkipBinnedParamsVsTemp = false;
SkipAPSubplots = false;
IncludeFits = true;
UseRescaledParamTiming = false;
UseRescaledFluo = false;
UsePerNucleusTraces = false;
UseBinnedTraces = false;
UseBinnedPerNucleusTraces = false;

subfn_varargin = {};
x = 1;
while x <= length(varargin)
    if strcmp(lower(varargin{x}), 'plottitle')
        PlotTitle = varargin{x+1};
        subfn_varargin = [subfn_varargin, 'PlotTitle', PlotTitle];
        x = x+1;
    elseif strcmp(lower(varargin{x}), 'plottingcolors')
        PlottingColors = varargin{x+1};
        
        x = x+1;
    elseif strcmp(lower(varargin{x}), 'usephysicalaplength')
        UsePhysicalAPLength = true;
        subfn_varargin = [subfn_varargin, 'UsePhysicalAPLength'];
    elseif strcmp(lower(varargin{x}), 'noline')
        UseLines = false;
        subfn_varargin = [subfn_varargin, 'UseLines'];
    elseif strcmpi(varargin{x}, 'SkipParamsVsAP')
        SkipParamsVsAP = true;
    elseif strcmpi(varargin{x}, 'SkipSingleTempParamsVsAP')
        SkipSingleTempParamsVsAP = true;
    elseif strcmpi(varargin{x}, 'SkipParamsVsTemp')
        SkipParamsVsTemp = true;
    elseif strcmpi(varargin{x}, 'SkipBinnedParamsVsAP')
        SkipBinnedParamsVsAP = true;
    elseif strcmpi(varargin{x}, 'SkipBinnedParamsVsTemp')
        SkipBinnedParamsVsTemp = true;
    elseif strcmpi(varargin{x}, 'SkipAPSubplots')
        SkipAPSubplots = true;
    elseif strcmpi(varargin{x}, 'ExcludeFits')
        IncludeFits = false;
        subfn_varargin = [subfn_varargin, 'ExcludeFits'];
    elseif strcmpi(varargin{x}, 'UsePerNucleusTraces')
        UsePerNucleusTraces = true;
    elseif strcmpi(varargin{x}, 'UseBinnedTraces')
        UseBinnedTraces = true;
        SkipBinnedParamsVsAP = true;
        SkipBinnedParamsVsTemp = true;
    elseif strcmpi(varargin{x}, 'UseBinnedPerNucleusTraces')
        UseBinnedPerNucleusTraces = true;
        SkipBinnedParamsVsAP = true;
        SkipBinnedParamsVsTemp = true;
        UseBinnedTraces = false;
        UsePerNucleusTraces = false;
    elseif strcmp(lower(varargin{x}), 'tracetype')
        TraceType = lower(varargin{x+1});
        x = x+1;
    elseif strcmpi(varargin{x}, 'userescaledtime') | strcmpi(varargin{x}, 'rescaletime') | ...
            strcmpi(varargin{x}, 'rescaletiming') | strcmpi(varargin{x}, 'userescaledtiming') | ...
            strcmpi(varargin{x}, 'userescaledparamtime') | strcmpi(varargin{x}, 'rescaleparamtime') | ...
            strcmpi(varargin{x}, 'rescaleparamtiming') | strcmpi(varargin{x}, 'userescaledparamtiming')
        UseRescaledParamTiming = true;
        subfn_varargin = [subfn_varargin, 'UseRescaledParamTiming'];
    elseif strcmpi(varargin{x}, 'rescalefluo') | strcmpi(varargin{x}, 'userescaledfluo')
        UseRescaledFluo = true;
        subfn_varargin = [subfn_varargin, 'UseRescaledFluo'];
    end
    x = x+1;
end

if UseBinnedPerNucleusTraces
    UsePerNucleusTraces = false;
    UseBinnedTraces = false;
    subfn_varargin = [subfn_varargin, 'UseBinnedPerNucleusTraces'];
elseif UsePerNucleusTraces & UseBinnedTraces
    UseBinnedPerNucleusTraces = true;
    subfn_varargin = [subfn_varargin, 'UseBinnedPerNucleusTraces'];
    UsePerNucleusTraces = false;
    UseBinnedTraces = false;
elseif UsePerNucleusTraces 
    subfn_varargin = [subfn_varargin, 'UsePerNucleusTraces'];
elseif UseBinnedTraces
    subfn_varargin = [subfn_varargin, 'UseBinnedTraces'];
end

if ~exist('PlottingColors', 'var')
    PlottingColors = 'default';
    subfn_varargin = [subfn_varargin, 'PlottingColors', PlottingColors];
elseif ~strcmpi(PlottingColors, 'gradient') &~strcmp(lower(PlottingColors), 'default')  & ~strcmp(lower(PlottingColors), 'pboc')
    error('Invalid choice of plotting colors. Can use either "default", "pboc", or "gradient".') % change to error
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
subfn_varargin = [subfn_varargin, 'TraceType', TraceType];
%%


if ~SkipParamsVsAP & ~UseBinnedTraces & ~UseBinnedPerNucleusTraces
    PlotLTMTrapParamsVsAP(this, parameter, outdir, subfn_varargin{:});
end
    
if ~SkipSingleTempParamsVsAP & ~UseBinnedTraces & ~UseBinnedPerNucleusTraces
    PlotLTMSingleTempTrapParamsVsAP(this, parameter, outdir, subfn_varargin{:});
end
    
    
if ~SkipParamsVsTemp
    PlotLTMTrapParamsVsTemp(this, parameter, outdir, subfn_varargin{:});
end

if ~SkipBinnedParamsVsAP
    PlotLTMBinnedTrapParamsVsAP(this, parameter, outdir, subfn_varargin{:});
end
if ~SkipBinnedParamsVsTemp
    PlotLTMBinnedTrapParamsVsTemp(this, parameter, outdir, subfn_varargin{:});
end

if ~SkipAPSubplots
    PlotLTMSubplotsTrapParamsVsTemp(this, parameter, outdir, subfn_varargin{:});
    PlotLTMLogSubplotsTrapParamsVsTemp(this, parameter, outdir, subfn_varargin{:});
end

close all   
    
    