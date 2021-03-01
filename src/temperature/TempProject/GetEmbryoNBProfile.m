function [NCAPtrace, StdErrorTrace, NCTimes] = GetEmbryoNBProfile(this, SetIndex, NC,APindex, TraceType)
% author: G. Martini
% date created: 2/22/21
% date last modified: 2/22/21

if ~( strcmp(lower(TraceType), 'unaligned')  | strcmp(lower(TraceType), 'anaphasealigned')|...
         strcmp(lower(TraceType), 'tbinned'))
    error('Invalid choice of trace type. Can use either "unaligned", "tbinned", or "anaphasealigned".') % change to error
end
if strcmpi(TraceType, 'anaphasealigned')
    traceName = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'unaligned')
    traceName = 'Unaligned';
elseif strcmpi(TraceType, 'tbinned')
    traceName = 'Tbinned';
end



MeanProfileTraces = this.TFProfiles{SetIndex}.([traceName, 'CycleMeanTraces']);
StdErrorTraces =  this.TFProfiles{SetIndex}.([traceName, 'CycleTraceStdErrors']);
NumNuclei = this.TFProfiles{SetIndex}.([traceName, 'CycleNumOnNuclei']);
CycleTimes = this.TFProfiles{SetIndex}.([traceName, 'CycleFrameTimes']);


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
    NCTimes = NCTimes(~isnan(NCAPtrace) & (NCAPtrace > 0) & (NumNucleiTrace >= this.MinimumNuclearCount));
    StdErrorTrace = StdErrorTrace(~isnan(NCAPtrace) & (NCAPtrace > 0) & (NumNucleiTrace >= this.MinimumNuclearCount));
    NCAPtrace = NCAPtrace(~isnan(NCAPtrace) & (NCAPtrace > 0) & (NumNucleiTrace >= this.MinimumNuclearCount));
end



end