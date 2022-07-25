function [NCAPtrace, StdErrorTrace, NCTimes] = GetSpecifiedTracesForTbinnedFitting(this, TempIndex, NC,APindex, TraceType)
% author: G. Martini
% date created: 7/19/22
% date last modified: 7/19/22

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



MeanProfileTraces = this.BinnedMeanProfiles.([traceName, 'CycleMeanTraces'])(:,APindex,NC-8,TempIndex).';
StdErrorTraces =  this.BinnedMeanProfiles.([traceName, 'CycleStdErrors'])(:,APindex,NC-8,TempIndex).';
EmbryoCountTrace = this.BinnedMeanProfiles.([traceName, 'CycleNumEmbryos'])(:,APindex,NC-8,TempIndex).';
NumNuclei = this.BinnedMeanProfiles.([traceName, 'CycleTotalNuclei'])(:,APindex,NC-8,TempIndex).';
FractionOnNuclei = this.BinnedMeanProfiles.([traceName, 'CycleFractionOn'])(:,APindex,NC-8,TempIndex).';
NCTimes = this.BinnedMeanProfiles.([traceName, 'CycleTimes'])/60;


LastIndex = min(length(NCTimes), length(MeanProfileTraces));
NCAPtrace = MeanProfileTraces(1:LastIndex);
NumNucleiTrace = NumNuclei(1:LastIndex);
StdErrorTrace = StdErrorTraces(1:LastIndex);
EmbryoCountTrace = EmbryoCountTrace(1:LastIndex);
NCTimes = NCTimes(1:LastIndex);

if isempty(NCTimes)
    NCTimes = [];
    StdErrorTrace = [];
    NCAPtrace = [];
else
    NCTimes = NCTimes(~isnan(NCAPtrace) & (NCAPtrace > 0) & (EmbryoCountTrace >= this.MinimumEmbryos));
    StdErrorTrace = StdErrorTrace(~isnan(NCAPtrace) & (NCAPtrace > 0) & (EmbryoCountTrace >= this.MinimumEmbryos));
    NCAPtrace = NCAPtrace(~isnan(NCAPtrace) & (NCAPtrace > 0)& (EmbryoCountTrace >= this.MinimumEmbryos));
end



end