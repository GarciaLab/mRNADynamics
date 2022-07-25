function [NCAPtrace, StdErrorTrace, NCTimes] = GetSpecifiedPerNucleusTracesForFitting(this, SetIndex, NC,APindex, TraceType)
% author: G. Martini
% date created: 1/10/20
% date last modified: 1/13/20

if ~(strcmp(lower(TraceType), 'fluo3d') | strcmp(lower(TraceType), 'fluo')  | strcmp(lower(TraceType), 'anaphasealigned')|...
        strcmp(lower(TraceType), 'anaphasealigned3d') | strcmp(lower(TraceType), 'tbinned') | strcmp(lower(TraceType), 'tbinned3d'))
    error('Invalid choice of trace type. Can use either "fluo", "fluo3d", "tbinned", "tbinned3d", "anaphasealigned", or "anaphasealigned3d".') % change to error
end
if strcmpi(TraceType, 'anaphasealigned')
    traceName = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'anaphasealigned3d')
    traceName = 'AnaphaseAligned3D';
elseif strcmpi(TraceType, 'fluo')
    traceName = 'Unaligned';
elseif strcmpi(TraceType, 'fluo3d')
    traceName = 'Unaligned3D';
elseif strcmpi(TraceType, 'tbinned')
    traceName = 'Tbinned';
elseif strcmpi(TraceType, 'tbinned3d')
    traceName = 'Tbinned3D';
end



MeanProfileTraces = this.MeanProfiles{SetIndex}.([traceName, 'CycleMeanPerNucleusTraces']);
StdErrorTraces =  this.MeanProfiles{SetIndex}.([traceName, 'CycleTraceStdErrors']);
NumNuclei = this.MeanProfiles{SetIndex}.([traceName, 'CycleNumOnNuclei']);
CycleTimes = this.MeanProfiles{SetIndex}.([traceName, 'CycleFrameTimes']);


NCTimes = CycleTimes{NC-8}/60;
LastIndex = min(length(NCTimes), size(MeanProfileTraces, 1));
NCAPtrace = squeeze(MeanProfileTraces(1:LastIndex,APindex,NC-8)).';
NumNucleiTrace = squeeze(NumNuclei(1:LastIndex,APindex,NC-8)).';
StdErrorTrace = squeeze(StdErrorTraces(1:LastIndex,APindex,NC-8)).';
NCTimes = NCTimes(1:LastIndex);
if isempty(NCTimes)
    NCTimes = [];
    StdErrorTrace = [];
    NCAPtrace = [];
else
    NCTimes = NCTimes(~isnan(NCAPtrace) & (NCAPtrace > 0) & (NumNucleiTrace >= this.MinimumTraceCount));
    StdErrorTrace = StdErrorTrace(~isnan(NCAPtrace) & (NCAPtrace > 0) & (NumNucleiTrace >= this.MinimumTraceCount));
    NCAPtrace = NCAPtrace(~isnan(NCAPtrace) & (NCAPtrace > 0) & (NumNucleiTrace >= this.MinimumTraceCount));
end



end