function MaxTimePoints = GetMaxTimePointsForAllNC(this, SetIndex, TraceType)
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


NCTimes = CycleTimes{end}/60;
MaxTimePoints = min( max(NCTimes)/(this.time_delta/60)+1, size(MeanProfileTraces, 1));
