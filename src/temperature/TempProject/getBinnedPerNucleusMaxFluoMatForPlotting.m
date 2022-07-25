function [output_mat, se_output_mat] = getBinnedPerNucleusMaxFluoMatForPlotting(this, TraceType)
if strcmp(lower(TraceType), 'anaphasealigned')
    fluoLabel = 'AnaphaseAlignedCycleNuclearMeanTraces';
    seLabel = 'AnaphaseAlignedCycleNuclearStdErrors';
    countLabel = 'AnaphaseAlignedCycleNumEmbryos';
elseif strcmp(lower(TraceType), 'anaphasealigned3d')
    fluoLabel = 'AnaphaseAligned3DCycleNuclearMeanTraces';
    seLabel = 'AnaphaseAligned3DCycleNuclearStdErrors';
    countLabel = 'AnaphaseAligned3DCycleNumEmbryos';
elseif strcmp(lower(TraceType), 'tbinned')
    fluoLabel = 'TbinnedCycleNuclearMeanTraces';
    seLabel = 'TbinnedCycleNuclearStdErrors';
    countLabel = 'TbinnedCycleNumEmbryos';
elseif strcmp(lower(TraceType), 'tbinned3d')
    fluoLabel = 'Tbinned3DCycleNuclearMeanTraces';
    seLabel = 'Tbinned3DCycleNuclearStdErrors';
    countLabel = 'Tbinned3DCycleNumEmbryos';
else
    error(['Invalid choice of TraceType: ', TraceType]);
end
Nsets = length(this.ExperimentPrefixes);
NTemperatures = length(unique(this.Temp_sps));
Nbins = 1/(this.Experiments{1}.APResolution)+1;
output_mat = NaN(NTemperatures, Nbins, 6);
se_output_mat = NaN(NTemperatures, Nbins, 6);
test_output_mat = NaN(NTemperatures, Nbins, 6);
d1dummy = ones(1, Nbins*6, 'uint8');
d2dummy = repmat(1:41,1,6);
d3dummy = repmat(1:6,Nbins, 1);
d3dummy = reshape(d3dummy, 1, 6*Nbins);
for i=1:NTemperatures
    means = this.BinnedMeanProfiles.(fluoLabel)(:,:,:,i);
    ses = this.BinnedMeanProfiles.(seLabel)(:,:,:,i);
    counts= this.BinnedMeanProfiles.(countLabel)(:,:,:,i);
    if ~all(size(counts) == size(means))
        counts = counts(1:size(means, 1), :,:);
    end
    means(counts< this.MinimumEmbryos) = NaN;
    ses(counts< this.MinimumEmbryos) = NaN;
    [output_mat(i,:,:), idx] = max(means, [], 1);
    idx = squeeze(idx);
    flatidx = reshape(idx, 1, size(idx, 1)*size(idx, 2));
    storage_I = sub2ind([NTemperatures, Nbins, 6],d1dummy*i, d2dummy, d3dummy);
    access_I = sub2ind(size(ses),flatidx, d2dummy, d3dummy);
    se_output_mat(storage_I) = ses(access_I);
    test_output_mat(storage_I) = means(access_I);
end
