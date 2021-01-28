function [output_mat, se_output_mat] = getMaxFluoMatForPlotting(this, TraceType)
if strcmp(lower(TraceType), 'anaphasealigned')
    fluoLabel = 'AnaphaseAlignedCycleMeanTraces';
    seLabel = 'AnaphaseAlignedCycleTraceStdErrors';
    countLabel = 'AnaphaseAlignedCycleNumOnNuclei';
elseif strcmp(lower(TraceType), 'anaphasealigned3d')
    fluoLabel = 'AnaphaseAligned3DCycleMeanTraces';
    seLabel = 'AnaphaseAligned3DCycleTraceStdErrors';
    countLabel = 'AnaphaseAligned3DCycleNumOnNuclei';
elseif strcmp(lower(TraceType), 'tbinned')
    fluoLabel = 'TbinnedCycleMeanTraces';
    seLabel = 'TbinnedCycleTraceStdErrors';
    countLabel = 'TbinnedCycleNumOnNuclei';
elseif strcmp(lower(TraceType), 'tbinned3d')
    fluoLabel = 'Tbinned3DCycleMeanTraces';
    seLabel = 'Tbinned3DCycleTraceStdErrors';
    countLabel = 'Tbinned3DCycleNumOnNuclei';
elseif strcmp(lower(TraceType), 'fluo')
    fluoLabel = 'UnalignedCycleMeanTraces';
    seLabel = 'UnalignedCycleTraceStdErrors';
    countLabel = 'UnalignedCycleNumOnNuclei';
elseif strcmp(lower(TraceType), 'fluo3d')
    fluoLabel = 'Unaligned3DCycleMeanTraces';
    seLabel = 'Unaligned3DCycleTraceStdErrors';
    countLabel = 'Unaligned3DCycleNumOnNuclei';
else
    error(['Invalid choice of TraceType: ', TraceType]);
end
Nsets = length(this.ExperimentPrefixes);
Nbins = 1/(this.Experiments{1}.APResolution)+1;
output_mat = NaN(Nsets, Nbins, 6);
se_output_mat = NaN(Nsets, Nbins, 6);
test_output_mat = NaN(Nsets, Nbins, 6);
d1dummy = ones(1, Nbins*6, 'uint8');
d2dummy = repmat(1:41,1,6);
d3dummy = repmat(1:6,Nbins, 1);
d3dummy = reshape(d3dummy, 1, 6*Nbins);
for i=1:Nsets
    if ismember(i, this.ProcessedExperiments)
        means = this.MeanProfiles{i}.(fluoLabel);
        ses = this.MeanProfiles{i}.(seLabel);
        counts = this.MeanProfiles{i}.(countLabel);
        means(counts < this.MinimumTraceCount) = NaN;
        ses(counts < this.MinimumTraceCount) = NaN;
        [output_mat(i,:,:), idx] = max(means, [], 1);
        idx = squeeze(idx);
        flatidx = reshape(idx, 1, size(idx, 1)*size(idx, 2));
        storage_I = sub2ind([Nsets, Nbins, 6],d1dummy*i, d2dummy, d3dummy);
        access_I = sub2ind(size(ses),flatidx, d2dummy, d3dummy);
        se_output_mat(storage_I) = ses(access_I);
        test_output_mat(storage_I) = means(access_I);
    end
end