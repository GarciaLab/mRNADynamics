function MeanProfiles = CalculateNBAPProfiles(Prefix, deltaTbinWidth, MinTimePoints, verbose)
if ~exist('deltaTbinWidth', 'var')
    deltaTbinWidth = 60;
end
if ~exist('MinTimePoints', 'var')
    MinTimePoints = 5;
end
if ~exist('verbose', 'var')
    verbose = true;
end

liveExperiment = LiveExperiment(Prefix);
HasHistone = false;
for i = 1:length(liveExperiment.Channels)
    if contains(lower(liveExperiment.Channels{i}), 'his')
       HasHistone = true;
    end
end

APbins = 0:liveExperiment.APResolution:1;
FrameInfo = getFrameInfo(liveExperiment);
FrameTimes = [FrameInfo(:).Time]; % in seconds
numFrames = length(FrameTimes);
nc_info = [liveExperiment.nc9, liveExperiment.nc10, liveExperiment.nc11,...
    liveExperiment.nc12, liveExperiment.nc13, liveExperiment.nc14, length(FrameInfo)];
PixelSize = liveExperiment.pixelSize_um;
load([liveExperiment.resultsFolder, 'CompiledNuclei.mat'], 'CompiledNuclei')
%%
% First make average profiles binning everything by first anaphase of the
% nuclear cycle
NChannels = length(liveExperiment.inputChannels);
if NChannels > 1
    error('Currently multiple input channels are not supported.')
elseif iscell(CompiledNuclei)
    CompiledNuclei = CompiledNuclei{1};
end

ApprovedCount = zeros(1, length(CompiledNuclei));
AllApproved = zeros(1, length(CompiledNuclei));
FirstApproved = zeros(1, length(CompiledNuclei));
PercentageFramesApproved = zeros(1, length(CompiledNuclei));
HasSpotsBeforeDisapproved = zeros(1, length(CompiledNuclei));
HasMinimumTimePoints = zeros(1, length(CompiledNuclei));
for cn=1:length(CompiledNuclei)
    if (CompiledNuclei(cn).Approved == 1)  & ...
            (CompiledNuclei(cn).Flag6 == 0)
        ApprovedCount(cn) = 1;
    end
    if CompiledNuclei(cn).VelocityInfo.SchnitzHasAllFrames
        AllApproved(cn) =  1;
    end
    if HasHistone
        if isfield('schnitzcells', 'containsFirstFrameOfCycle')
            if CompiledNuclei(cn).containsFirstFrameOfCycle
                FirstApproved(cn) =  1;
            end
        else
            if ~isempty(CompiledNuclei(cn).anaphaseFrame)
                if ~CompiledNuclei(cn).inferredAnaphaseFrame
                    FirstApproved(cn) =  1;
                end
            end
        end
    else
        FirstApproved(cn) =  1;
    end
    if all(CompiledNuclei(cn).FlaggingInfo.FrameApprovedFinal == 1)
        PercentageFramesApproved(cn) = 1;
    else
        PercentageFramesApproved(cn) = sum(CompiledNuclei(cn).FlaggingInfo.FrameApprovedFinal == 1)/length(CompiledNuclei(cn).FlaggingInfo.FrameApprovedFinal);
    end
    if ~isempty(find(CompiledNuclei(cn).FlaggingInfo.FrameApprovedFinal == 0, 1))
        if find(CompiledNuclei(cn).FlaggingInfo.FrameApprovedFinal == 0, 1) >...
                CompiledNuclei(cn).Frames(1)
            HasSpotsBeforeDisapproved(cn) = 1;
        end
    else
        HasSpotsBeforeDisapproved(cn)  = 1;
    end
    if length(CompiledNuclei(cn).Frames(CompiledNuclei(cn).FrameApproved == 1)) >= MinTimePoints
        HasMinimumTimePoints(cn) = 1;
    end
    
end

%%
CNcycles =  [CompiledNuclei(:).nc];
if length(CNcycles) ~= length(CompiledNuclei)
    error('must have a nuclear cycle assigned for all CompiledNuclei.')
end

UnalignedCycleMeanTraces = NaN(numFrames, length(APbins), 6);
UnalignedCycleTraceStdErrors =  NaN(numFrames, length(APbins), 6);
UnalignedCycleNumNuclei =  NaN(numFrames, length(APbins), 6);
UnalignedCycleNumOnNuclei =  NaN(numFrames, length(APbins), 6);
CycleTraces = cell(1, 6);
MeanAPs =  cell(1, 6);
ParticleAPbinIndices =  cell(1, 6);
UnalignedCycleFrameTimes = cell(1, 6);

numBinnedFrames = ceil((max(FrameTimes)-FrameTimes(nc_info(6)))/deltaTbinWidth)+1;
AnaphaseAlignedCycleMeanTraces = NaN(numBinnedFrames, length(APbins), 6);
AnaphaseAlignedCycleTraceStdErrors =  NaN(numBinnedFrames, length(APbins), 6);
AnaphaseAlignedCycleNumNuclei =  NaN(numBinnedFrames, length(APbins), 6);
AnaphaseAlignedCycleNumOnNuclei =  NaN(numBinnedFrames, length(APbins), 6);
CycleTracesWithAnaphaseAlignment = cell(1, 6);
AnaphaseAlignedCycleFrameTimes =  cell(1, 6);

TbinnedCycleMeanTraces = NaN(numBinnedFrames, length(APbins), 6);
TbinnedCycleTraceStdErrors =  NaN(numBinnedFrames, length(APbins), 6);
TbinnedCycleNumNuclei =  NaN(numBinnedFrames, length(APbins), 6);
TbinnedCycleNumOnNuclei =  NaN(numBinnedFrames, length(APbins), 6);
TbinnedCycleTraces =  cell(1, 6);
TbinnedCycleFrameTimes = cell(1, 6);


IncludedNCs = min(CNcycles):max(CNcycles);
%%
for k =1:length(IncludedNCs)
    NC = IncludedNCs(k);
    IncludedTraceIndices = find((CNcycles == NC) & (FirstApproved == 1) &...
        (ApprovedCount == 1) & (HasMinimumTimePoints == 1) & (PercentageFramesApproved >= .7));
    if isempty(IncludedTraceIndices)
        continue
    end
    
    % First calculate means for all traces with no time binning or
    % anaphase alignment.
    
    
    CycleTraces{NC-8} = NaN(numFrames, length(IncludedTraceIndices));
    %CycleTraces3Derror =  zeros(numFrames, length(IncludedTraceIndices));
    MeanAPs{NC-8} = zeros(1, length(IncludedTraceIndices));
    ParticleAPbinIndices{NC-8} = zeros(1, length(IncludedTraceIndices), 'uint16');
    for i = 1:length(IncludedTraceIndices)
        trace_index = IncludedTraceIndices(i);
        CurrentCompiledParticle = CompiledNuclei(trace_index);
        CycleTraces{NC-8}(CurrentCompiledParticle.FlaggingInfo.TrueFrames(CurrentCompiledParticle.FlaggingInfo.FrameApprovedFinal == 1), i) = 0;
        CycleTraces{NC-8}(CurrentCompiledParticle.Frames(CurrentCompiledParticle.FrameApproved == 1), i) = CurrentCompiledParticle.FluoTimeTrace(CurrentCompiledParticle.FrameApproved == 1);
        CycleTraces{NC-8}(CurrentCompiledParticle.Frames(CurrentCompiledParticle.FrameApproved == 0), i) = NaN;
        %CycleTraces3Derror(CurrentCompiledParticle.Frame, i) = CurrentCompiledParticle.Fluo3DGaussz;
        MeanAPs{NC-8}(i) = mean(CurrentCompiledParticle.APpos);
        BinIndex = find((MeanAPs{NC-8}(i) < APbins(2:end)) & (MeanAPs{NC-8}(i) >= APbins(1:end-1)));
        if ~isempty(BinIndex)
            ParticleAPbinIndices{NC-8}(i) = BinIndex;
        else
            ParticleAPbinIndices{NC-8}(i) = NaN;
        end
        
    end
    
    APbinsInMovie = nanmin(ParticleAPbinIndices{NC-8}):nanmax(ParticleAPbinIndices{NC-8});
    NCCycleTraces = CycleTraces{NC-8};
    IncludedRowsInBin = find(sum(~isnan(NCCycleTraces),2).' > 0);
    UnalignedCycleFrameTimes{NC-8}= FrameTimes(IncludedRowsInBin)-min(FrameTimes(IncludedRowsInBin));
    CycleTraces{NC-8} = NCCycleTraces(IncludedRowsInBin,:);
    NCCycleTraces = NCCycleTraces(IncludedRowsInBin,:);
    NCCycleTraces(NCCycleTraces == 0) = NaN;

    
    for APbinIndex = APbinsInMovie
        IncludedColumnsInBin = find(ParticleAPbinIndices{NC-8} == APbinIndex);
        IncludedRowsInBin = find(sum(~isnan(CycleTraces{NC-8}(:,IncludedColumnsInBin)),2).' > 0);
        UnalignedCycleMeanTraces(IncludedRowsInBin, APbinIndex, NC-8) = 0;
        UnalignedCycleNumNuclei(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(CycleTraces{NC-8}(IncludedRowsInBin, IncludedColumnsInBin)), 2);
        UnalignedCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8) = 0;
        
        IncludedRowsInBin = find(sum(~isnan(NCCycleTraces(:,IncludedColumnsInBin)),2).' > 0);
        UnalignedCycleMeanTraces(IncludedRowsInBin, APbinIndex, NC-8) = nanmean(NCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin), 2);
        UnalignedCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8) = nanstd(NCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin),0, 2)./...
            sqrt(length(find(~isnan(NCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin)))));
        UnalignedCycleNumOnNuclei(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(NCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin)), 2);
    end
    
    
    
    
    % Do anaphase alignment and bin
    CurrentNCTraces = CycleTraces{NC-8};
    IncludedRowsInBin = find(sum(~isnan(CurrentNCTraces),2).' > 0);
    NCFrameTimes = FrameTimes(IncludedRowsInBin);
    NumBins = ceil((max(NCFrameTimes)-min(NCFrameTimes))/deltaTbinWidth)+1;
    AnaphaseAlignedCycleFrameTimes{NC-8} = (0:NumBins-1)*deltaTbinWidth;
    if length(AnaphaseAlignedCycleFrameTimes{NC-8}) > size(AnaphaseAlignedCycleMeanTraces, 1)
        AnaphaseAlignedCycleMeanTraces = [AnaphaseAlignedCycleMeanTraces; NaN(length(AnaphaseAlignedCycleFrameTimes{NC-8})-size(AnaphaseAlignedCycleMeanTraces, 1),length(APbins),6)];
        AnaphaseAlignedCycleTraceStdErrors = [AnaphaseAlignedCycleTraceStdErrors; NaN(length(AnaphaseAlignedCycleFrameTimes{NC-8})-size(AnaphaseAlignedCycleTraceStdErrors, 1),length(APbins),6)];
        AnaphaseAlignedCycleNumNuclei = [AnaphaseAlignedCycleNumNuclei; NaN(length(AnaphaseAlignedCycleFrameTimes{NC-8})-size(AnaphaseAlignedCycleNumNuclei, 1),length(APbins),6)];
        AnaphaseAlignedCycleNumOnNuclei = [AnaphaseAlignedCycleNumOnNuclei; NaN(length(AnaphaseAlignedCycleFrameTimes{NC-8})-size(AnaphaseAlignedCycleNumOnNuclei, 1),length(APbins),6)];
    end
    CycleTracesWithAnaphaseAlignment{NC-8} = NaN(NumBins, size(CurrentNCTraces, 2));
    for trace_index = 1:size(CurrentNCTraces, 2)
        CurrentTrace = CurrentNCTraces(IncludedRowsInBin,trace_index);
        TraceIncludedRowsInBin = find(~isnan(CurrentTrace));
        TraceFrameTimes = NCFrameTimes(TraceIncludedRowsInBin)-min(NCFrameTimes(TraceIncludedRowsInBin));
        CurrentTrace = CurrentTrace(TraceIncludedRowsInBin);
        IncludedInterpolatedFrames = ones(1, length(AnaphaseAlignedCycleFrameTimes{NC-8}));
        IncludedInterpolatedFrames(AnaphaseAlignedCycleFrameTimes{NC-8} > max(TraceFrameTimes)) = 0;
        for bin_index = 1:NumBins
            if ~isempty(find(AnaphaseAlignedCycleFrameTimes{NC-8}(bin_index) == TraceFrameTimes, 1))
                continue
            end
            if IncludedInterpolatedFrames(bin_index) == 0
                continue
            end
            TimeInterval = TraceFrameTimes(find(TraceFrameTimes > AnaphaseAlignedCycleFrameTimes{NC-8}(bin_index), 1))-TraceFrameTimes(find(TraceFrameTimes < AnaphaseAlignedCycleFrameTimes{NC-8}(bin_index), 1, 'last'));
            if TimeInterval > deltaTbinWidth*3
                IncludedInterpolatedFrames(bin_index) = 0;
            end
        end
        if length(TraceFrameTimes) >= MinTimePoints
            InterpolatedCurrentTrace = interp1(TraceFrameTimes,CurrentTrace,AnaphaseAlignedCycleFrameTimes{NC-8}(find(IncludedInterpolatedFrames == 1)));
            CycleTracesWithAnaphaseAlignment{NC-8}(find(IncludedInterpolatedFrames == 1), trace_index) =InterpolatedCurrentTrace;
        end
    end

   
    APbinsInMovie = min(ParticleAPbinIndices{NC-8}):max(ParticleAPbinIndices{NC-8});
    AnaphaseAlignedNCCycleTraces = CycleTracesWithAnaphaseAlignment{NC-8};
    AnaphaseAlignedNCCycleTraces(AnaphaseAlignedNCCycleTraces == 0) = NaN;
    for APbinIndex = APbinsInMovie
        IncludedColumnsInBin = find(ParticleAPbinIndices{NC-8} == APbinIndex);
        IncludedRowsInBin = find(sum(~isnan(CycleTracesWithAnaphaseAlignment{NC-8}(:,IncludedColumnsInBin)),2).' > 0);
        AnaphaseAlignedCycleMeanTraces(IncludedRowsInBin, APbinIndex, NC-8) = 0;
        AnaphaseAlignedCycleNumNuclei(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(CycleTracesWithAnaphaseAlignment{NC-8}(IncludedRowsInBin, IncludedColumnsInBin)), 2);
        AnaphaseAlignedCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8) = 0;
    
        IncludedRowsInBin = find(sum(~isnan(AnaphaseAlignedNCCycleTraces(:,IncludedColumnsInBin)),2).' > 0);
        AnaphaseAlignedCycleMeanTraces(IncludedRowsInBin, APbinIndex, NC-8) = nanmean(AnaphaseAlignedNCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin), 2);
        AnaphaseAlignedCycleNumOnNuclei(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(AnaphaseAlignedNCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin)), 2);
        AnaphaseAlignedCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8) = nanstd(AnaphaseAlignedNCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin),0, 2)./...
            sqrt(AnaphaseAlignedCycleNumOnNuclei(IncludedRowsInBin, APbinIndex, NC-8));
     
    end
    
    % Bin without anaphase alignment
    CurrentNCTraces = CycleTraces{NC-8};
    IncludedRowsInBin = find(sum(~isnan(CurrentNCTraces),2).' > 0);
    NCFrameTimes = FrameTimes(IncludedRowsInBin);
    NumBins = ceil((max(NCFrameTimes)-min(NCFrameTimes))/deltaTbinWidth)+1;
    TbinnedCycleFrameTimes{NC-8} = (0:NumBins-1)*deltaTbinWidth;
    if length(TbinnedCycleFrameTimes{NC-8}) > size(TbinnedCycleMeanTraces, 1)
        TbinnedCycleMeanTraces = [TbinnedCycleMeanTraces; NaN(length(TbinnedCycleFrameTimes{NC-8})-size(TbinnedCycleMeanTraces, 1),length(APbins),6)];
        TbinnedCycleTraceStdErrors = [TbinnedCycleTraceStdErrors; NaN(length(TbinnedCycleFrameTimes{NC-8})-size(TbinnedCycleTraceStdErrors, 1),length(APbins),6)];
        TbinnedCycleNumNuclei = [TbinnedCycleNumNuclei; NaN(length(TbinnedCycleFrameTimes{NC-8})-size(TbinnedCycleNumNuclei, 1),length(APbins),6)];
        TbinnedCycleNumOnNuclei = [TbinnedCycleNumOnNuclei; NaN(length(TbinnedCycleFrameTimes{NC-8})-size(TbinnedCycleNumOnNuclei, 1),length(APbins),6)];
    end
    TbinnedCycleTraces{NC-8} = NaN(NumBins, size(CurrentNCTraces, 2));
    
    for trace_index = 1:size(CurrentNCTraces, 2)
        CurrentTrace = CurrentNCTraces(IncludedRowsInBin,trace_index);
        TraceIncludedRowsInBin = find(~isnan(CurrentTrace));
        TraceFrameTimes = NCFrameTimes(TraceIncludedRowsInBin)-min(NCFrameTimes);
        CurrentTrace = CurrentTrace(TraceIncludedRowsInBin);
        IncludedInterpolatedFrames = ones(1, length(TbinnedCycleFrameTimes{NC-8}));
        IncludedInterpolatedFrames(TbinnedCycleFrameTimes{NC-8} > max(TraceFrameTimes)) = 0;
        IncludedInterpolatedFrames(TbinnedCycleFrameTimes{NC-8} < min(TraceFrameTimes)) = 0;
        for bin_index = 1:NumBins
            if ~isempty(find(TbinnedCycleFrameTimes{NC-8}(bin_index) == TraceFrameTimes, 1))
                continue
            end
            if IncludedInterpolatedFrames(bin_index) == 0
                continue
            end
            TimeInterval = TraceFrameTimes(find(TraceFrameTimes > TbinnedCycleFrameTimes{NC-8}(bin_index), 1))-TraceFrameTimes(find(TraceFrameTimes < TbinnedCycleFrameTimes{NC-8}(bin_index), 1, 'last'));
            if TimeInterval > deltaTbinWidth*3
                IncludedInterpolatedFrames(bin_index) = 0;
            end
        end
        if length(TraceFrameTimes) >= MinTimePoints
            InterpolatedCurrentTrace = interp1(TraceFrameTimes,CurrentTrace,TbinnedCycleFrameTimes{NC-8}(find(IncludedInterpolatedFrames == 1)));
            TbinnedCycleTraces{NC-8}(find(IncludedInterpolatedFrames == 1), trace_index) =InterpolatedCurrentTrace;
        end
    end

    
    APbinsInMovie = min(ParticleAPbinIndices{NC-8}):max(ParticleAPbinIndices{NC-8});
    TbinnedNCCycleTraces = TbinnedCycleTraces{NC-8};
    TbinnedNCCycleTraces(TbinnedNCCycleTraces == 0) = NaN;
    for APbinIndex = APbinsInMovie
        IncludedColumnsInBin = find(ParticleAPbinIndices{NC-8} == APbinIndex);
        IncludedRowsInBin = find(sum(~isnan(TbinnedCycleTraces{NC-8}(:,IncludedColumnsInBin)),2).' > 0);
     
        TbinnedCycleMeanTraces(IncludedRowsInBin, APbinIndex, NC-8) = 0;
        TbinnedCycleNumNuclei(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(TbinnedCycleTraces{NC-8}(IncludedRowsInBin, IncludedColumnsInBin)), 2);
        TbinnedCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8) = 0;
       
        IncludedRowsInBin = find(sum(~isnan(TbinnedNCCycleTraces(:,IncludedColumnsInBin)),2).' > 0);
     
        TbinnedCycleMeanTraces(IncludedRowsInBin, APbinIndex, NC-8) = nanmean(TbinnedNCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin), 2);
        
        TbinnedCycleNumOnNuclei(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(TbinnedNCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin)), 2);
        TbinnedCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8) = nanstd(TbinnedNCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin),0, 2)./...
            sqrt(TbinnedCycleNumOnNuclei(IncludedRowsInBin, APbinIndex, NC-8));
    end
    
    % Check to make sure calculations make sense
    TestAllTracesMat = squeeze(UnalignedCycleMeanTraces(:,:,NC-8));
    TestAllNumNucMat = squeeze(UnalignedCycleNumNuclei(:,:,NC-8));
    if ~all(TestAllNumNucMat(TestAllTracesMat > 0) > 0)
        disp(['Inconsistencies between Mean Traces and Nuclear count in NC ', num2str(NC), ' for standard mean traces'])
    end

    TestAnaphaseAlignedTracesMat = squeeze(AnaphaseAlignedCycleMeanTraces(:,:,NC-8));
    TestAnaphaseAlignedNumNucMat = squeeze(AnaphaseAlignedCycleNumNuclei(:,:,NC-8));
    if ~all(TestAnaphaseAlignedNumNucMat(TestAnaphaseAlignedTracesMat > 0) > 0)
        disp(['Inconsistencies between Mean Traces and Nuclear count in NC ', num2str(NC), ' for aligned mean traces'])
    end
   
    TestTbinnedTracesMat = squeeze(TbinnedCycleMeanTraces(:,:,NC-8));
    TestTbinnedNumNucMat = squeeze(TbinnedCycleNumNuclei(:,:,NC-8));
    if ~all(TestTbinnedNumNucMat(TestTbinnedTracesMat > 0) > 0)
        disp(['Inconsistencies between Mean Traces and Nuclear count in NC ', num2str(NC), ' for t-binned mean traces'])
    end

    
    
end

TestLengths1 = [size(UnalignedCycleMeanTraces, 1), size(UnalignedCycleNumNuclei, 1), size(UnalignedCycleNumOnNuclei, 1), size(UnalignedCycleTraceStdErrors, 1)];
if ~all(TestLengths1 == TestLengths1(1))
    disp(['Inconsistence size for fluo calculations.'])
end



TestLengths3 = [size(AnaphaseAlignedCycleMeanTraces, 1), size(AnaphaseAlignedCycleNumNuclei, 1), size(AnaphaseAlignedCycleNumOnNuclei, 1), size(AnaphaseAlignedCycleTraceStdErrors, 1)];
if ~all(TestLengths3 == TestLengths3(1))
    disp(['Inconsistence size for anaphase aligned calculations.'])
end



MeanProfiles = {};
MeanProfiles.UnalignedCycleMeanTraces = UnalignedCycleMeanTraces;
MeanProfiles.UnalignedCycleNumNuclei = UnalignedCycleNumNuclei;
MeanProfiles.UnalignedCycleNumOnNuclei = UnalignedCycleNumOnNuclei;
MeanProfiles.UnalignedCycleTraceStdErrors = UnalignedCycleTraceStdErrors;
MeanProfiles.AnaphaseAlignedCycleMeanTraces = AnaphaseAlignedCycleMeanTraces;
MeanProfiles.AnaphaseAlignedCycleNumNuclei = AnaphaseAlignedCycleNumNuclei;
MeanProfiles.AnaphaseAlignedCycleNumOnNuclei = AnaphaseAlignedCycleNumOnNuclei;
MeanProfiles.AnaphaseAlignedCycleTraceStdErrors = AnaphaseAlignedCycleTraceStdErrors;
MeanProfiles.TbinnedCycleMeanTraces = TbinnedCycleMeanTraces;
MeanProfiles.TbinnedCycleNumNuclei = TbinnedCycleNumNuclei;
MeanProfiles.TbinnedCycleNumOnNuclei = TbinnedCycleNumOnNuclei;
MeanProfiles.TbinnedCycleTraceStdErrors = TbinnedCycleTraceStdErrors;
MeanProfiles.AnaphaseAlignedCycleFrameTimes = AnaphaseAlignedCycleFrameTimes;
MeanProfiles.UnalignedCycleFrameTimes = UnalignedCycleFrameTimes;
MeanProfiles.TbinnedCycleFrameTimes = TbinnedCycleFrameTimes;

savedVariables = {'UnalignedCycleMeanTraces', 'UnalignedCycleNumNuclei',  'UnalignedCycleNumOnNuclei', 'UnalignedCycleTraceStdErrors',...
    'AnaphaseAlignedCycleMeanTraces', 'AnaphaseAlignedCycleNumNuclei', 'AnaphaseAlignedCycleNumOnNuclei','AnaphaseAlignedCycleTraceStdErrors',...
    'TbinnedCycleMeanTraces', 'TbinnedCycleNumNuclei','TbinnedCycleNumOnNuclei', 'TbinnedCycleTraceStdErrors',...
    'UnalignedCycleFrameTimes', 'AnaphaseAlignedCycleFrameTimes',  'TbinnedCycleFrameTimes'};

if NChannels > 1
    pathname = [liveExperiment.resultsFolder,'MeanProfiles_', num2str(ChN),'.mat'];
else
    pathname = [liveExperiment.resultsFolder,'MeanProfiles.mat'];
end
try
    save(pathname, savedVariables{:});
catch
    save(pathname, savedVariables{:}, '-nocompression');
end









