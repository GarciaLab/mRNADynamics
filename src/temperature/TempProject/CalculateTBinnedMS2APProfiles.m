function ltm = CalculateTBinnedMS2APProfiles(ltm, deltaTbinWidth, MinTimePoints, verbose)
if ~exist('deltaTbinWidth', 'var')
    deltaTbinWidth = 30;
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
load([liveExperiment.resultsFolder, 'CompiledParticles.mat'], 'CompiledParticles')
%%
% First make average profiles binning everything by first anaphase of the
% nuclear cycle
NChannels = length(liveExperiment.spotChannels);
for ChN=1:NChannels
    ApprovedCount = zeros(1, length(CompiledParticles{ChN}));
    AllApproved = zeros(1, length(CompiledParticles{ChN}));
    FirstApproved = zeros(1, length(CompiledParticles{ChN}));
    PercentageFramesApproved = zeros(1, length(CompiledParticles{ChN}));
    HasSpotsBeforeDisapproved = zeros(1, length(CompiledParticles{ChN}));
    HasMinimumTimePoints = zeros(1, length(CompiledParticles{ChN}));
    for cp=1:length(CompiledParticles{ChN})
        if (CompiledParticles{ChN}(cp).Approved == 1) & (CompiledParticles{ChN}(cp).schnitzcell.Approved > 0) & ...
                (CompiledParticles{ChN}(cp).schnitzcell.Flag ~= 6)
            ApprovedCount(cp) = 1;
        end
        if CompiledParticles{ChN}(cp).schnitzcell.VelocityInfo.SchnitzHasAllFrames
            AllApproved(cp) =  1;
        end
        if HasHistone
            if isfield('schnitzcells', 'containsFirstFrameOfCycle')
                if CompiledParticles{ChN}(cp).schnitzcell.containsFirstFrameOfCycle
                    FirstApproved(cp) =  1;
                end
            else
                if ~isempty(CompiledParticles{ChN}(cp).schnitzcell.anaphaseFrame)
                    if ~CompiledParticles{ChN}(cp).schnitzcell.inferredAnaphaseFrame
                        FirstApproved(cp) =  1;
                    end
                end
            end
        else
            FirstApproved(cp) =  1;
        end
        if all(CompiledParticles{ChN}(cp).FlaggingInfo.FrameApprovedFinal == 1)
            PercentageFramesApproved(cp) = 1;
        else
            PercentageFramesApproved(cp) = sum(CompiledParticles{ChN}(cp).FlaggingInfo.FrameApprovedFinal == 1)/length(CompiledParticles{ChN}(cp).FlaggingInfo.FrameApprovedFinal);
        end
        if ~isempty(find(CompiledParticles{ChN}(cp).FlaggingInfo.FrameApprovedFinal == 0, 1))
            if find(CompiledParticles{ChN}(cp).FlaggingInfo.FrameApprovedFinal == 0, 1) >...
                    CompiledParticles{ChN}(cp).Frame(1)
                HasSpotsBeforeDisapproved(cp) = 1;
            end
        else
            HasSpotsBeforeDisapproved(cp)  = 1;
        end
        if length(CompiledParticles{ChN}(cp).Frame(CompiledParticles{ChN}(cp).FrameApproved == 1)) >= MinTimePoints
            HasMinimumTimePoints(cp) = 1;
        end
        
    end
    
    CPcycles =  [CompiledParticles{ChN}(:).cycle];
    
    UnalignedCycleMeanTraces = NaN(numFrames, length(APbins), 6);
    UnalignedCycleTraceStdErrors =  NaN(numFrames, length(APbins), 6);
    UnalignedCycleNumNuclei =  NaN(numFrames, length(APbins), 6);
    UnalignedCycleNumOnNuclei =  NaN(numFrames, length(APbins), 6);
    Unaligned3DCycleMeanTraces = NaN(numFrames, length(APbins), 6);
    Unaligned3DCycleTraceStdErrors =  NaN(numFrames, length(APbins), 6);
    Unaligned3DCycleNumNuclei =  NaN(numFrames, length(APbins), 6);
    Unaligned3DCycleNumOnNuclei =  NaN(numFrames, length(APbins), 6);
    CycleTraces = cell(1, 6);
    CycleTraces3D = cell(1, 6);
    MeanAPs = cell(1, 6);
    ParticleAPbinIndices = cell(1, 6);
    UnalignedCycleFrameTimes = cell(1, 6);
    Unaligned3DCycleFrameTimes =cell(1, 6);
    
    numBinnedFrames = ceil((max(FrameTimes)-FrameTimes(nc_info(6)))/deltaTbinWidth)+1;
    AnaphaseAlignedCycleMeanTraces = NaN(numBinnedFrames, length(APbins), 6);
    AnaphaseAlignedCycleTraceStdErrors =  NaN(numBinnedFrames, length(APbins), 6);
    AnaphaseAlignedCycleNumNuclei =  NaN(numBinnedFrames, length(APbins), 6);
    AnaphaseAlignedCycleNumOnNuclei =  NaN(numBinnedFrames, length(APbins), 6);
    AnaphaseAligned3DCycleMeanTraces = NaN(numBinnedFrames, length(APbins), 6);
    AnaphaseAligned3DCycleTraceStdErrors =  NaN(numBinnedFrames, length(APbins), 6);
    AnaphaseAligned3DCycleNumNuclei =  NaN(numBinnedFrames, length(APbins), 6);
    AnaphaseAligned3DCycleNumOnNuclei =  NaN(numBinnedFrames, length(APbins), 6);
    CycleTracesWithAnaphaseAlignment = cell(1, 6);
    Cycle3DTracesWithAnaphaseAlignment = cell(1, 6);
    AnaphaseAlignedCycleFrameTimes = cell(1, 6);
    AnaphaseAligned3DCycleFrameTimes = cell(1, 6);
    
    TbinnedCycleMeanTraces = NaN(numBinnedFrames, length(APbins), 6);
    TbinnedCycleTraceStdErrors =  NaN(numBinnedFrames, length(APbins), 6);
    TbinnedCycleNumNuclei =  NaN(numBinnedFrames, length(APbins), 6);
    TbinnedCycleNumOnNuclei =  NaN(numBinnedFrames, length(APbins), 6);
    Tbinned3DCycleMeanTraces = NaN(numBinnedFrames, length(APbins), 6);
    Tbinned3DCycleTraceStdErrors =  NaN(numBinnedFrames, length(APbins), 6);
    Tbinned3DCycleNumNuclei =  NaN(numBinnedFrames, length(APbins), 6);
    Tbinned3DCycleNumOnNuclei =  NaN(numBinnedFrames, length(APbins), 6);
    TbinnedCycleTraces = cell(1, 6);
    Tbinned3DCycleTraces = cell(1, 6);
    TbinnedCycleFrameTimes = cell(1, 6);
    Tbinned3DCycleFrameTimes = cell(1, 6);
    
    
    IncludedNCs = min(CPcycles):max(CPcycles);
    for k = 1:length(IncludedNCs)
        NC = IncludedNCs(k);
        IncludedTraceIndices = find((CPcycles == NC) & (FirstApproved == 1) &...
            (ApprovedCount == 1) & (HasMinimumTimePoints == 1) & (PercentageFramesApproved >= .7));
        if isempty(IncludedTraceIndices)
            continue
        end
        
        % First calculate means for all traces with no time binning or
        % anaphase alignment.
        
        
        CycleTraces{NC-8} = NaN(numFrames, length(IncludedTraceIndices));
        CycleTraces3D{NC-8} =  NaN(numFrames, length(IncludedTraceIndices));
        %CycleTraces3Derror =  zeros(numFrames, length(IncludedTraceIndices));
        MeanAPs{NC-8} = zeros(1, length(IncludedTraceIndices));
        ParticleAPbinIndices{NC-8} = zeros(1, length(IncludedTraceIndices), 'uint16');
        for i = 1:length(IncludedTraceIndices)
            trace_index = IncludedTraceIndices(i);
            CurrentCompiledParticle = CompiledParticles{ChN}(trace_index);
            CycleTraces{NC-8}(CurrentCompiledParticle.FlaggingInfo.TrueFrames(CurrentCompiledParticle.FlaggingInfo.FrameApprovedFinal == 1), i) = 0;
            CycleTraces{NC-8}(CurrentCompiledParticle.Frame(CurrentCompiledParticle.FrameApproved == 1), i) = CurrentCompiledParticle.Fluo(CurrentCompiledParticle.FrameApproved == 1);
            CycleTraces{NC-8}(CurrentCompiledParticle.Frame(CurrentCompiledParticle.FrameApproved == 0), i) = NaN;
            CycleTraces3D{NC-8}(CurrentCompiledParticle.FlaggingInfo.TrueFrames(CurrentCompiledParticle.FlaggingInfo.FrameApprovedFinal == 1), i) = 0;
            CycleTraces3D{NC-8}(CurrentCompiledParticle.Frame(CurrentCompiledParticle.FrameApproved == 1), i) = CurrentCompiledParticle.Fluo3DGauss(CurrentCompiledParticle.FrameApproved == 1);
            CycleTraces3D{NC-8}(CurrentCompiledParticle.Frame(CurrentCompiledParticle.FrameApproved == 0), i) = NaN;
            %CycleTraces3Derror(CurrentCompiledParticle.Frame, i) = CurrentCompiledParticle.Fluo3DGaussz;
            MeanAPs{NC-8}(i) = mean(CurrentCompiledParticle.schnitzcell.APpos);
            ParticleAPbinIndices{NC-8}(i) = find((MeanAPs{NC-8}(i) < APbins(2:end)) & (MeanAPs{NC-8}(i) >= APbins(1:end-1)));
            
        end
        
        APbinsInMovie = min(ParticleAPbinIndices{NC-8}):max(ParticleAPbinIndices{NC-8});
        NCCycleTraces = CycleTraces{NC-8};
        NC3DCycleTraces = CycleTraces3D{NC-8};
        IncludedRowsInBin = find(sum(~isnan(NCCycleTraces),2).' > 0);
        UnalignedCycleFrameTimes{NC-8}= FrameTimes(IncludedRowsInBin)-min(FrameTimes(IncludedRowsInBin));
        CycleTraces{NC-8} = NCCycleTraces(IncludedRowsInBin,:);
        NCCycleTraces = NCCycleTraces(IncludedRowsInBin,:);
        
        IncludedRowsInBin3D = find(sum(~isnan(NC3DCycleTraces),2).' > 0);
        Unaligned3DCycleFrameTimes{NC-8}= FrameTimes(IncludedRowsInBin3D)-min(FrameTimes(IncludedRowsInBin3D));
        CycleTraces3D{NC-8} = NC3DCycleTraces(IncludedRowsInBin3D,:);
        NC3DCycleTraces = NC3DCycleTraces(IncludedRowsInBin3D,:);
        
        NCCycleTraces(NCCycleTraces == 0) = NaN;
        NC3DCycleTraces(NC3DCycleTraces == 0) = NaN;
        
        for APbinIndex = APbinsInMovie
            IncludedColumnsInBin = find(ParticleAPbinIndices{NC-8} == APbinIndex);
            IncludedRowsInBin = find(sum(~isnan(CycleTraces{NC-8}(:,IncludedColumnsInBin)),2).' > 0);
            IncludedRowsInBin3D = find(sum(~isnan(CycleTraces3D{NC-8}(:,IncludedColumnsInBin)),2).' > 0);
            UnalignedCycleMeanTraces(IncludedRowsInBin, APbinIndex, NC-8) = 0;
            UnalignedCycleNumNuclei(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(CycleTraces{NC-8}(IncludedRowsInBin, IncludedColumnsInBin)), 2);
            UnalignedCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8) = 0;
            Unaligned3DCycleMeanTraces(IncludedRowsInBin3D, APbinIndex, NC-8) = 0;
            Unaligned3DCycleNumNuclei(IncludedRowsInBin3D, APbinIndex, NC-8) = sum(~isnan(CycleTraces3D{NC-8}(IncludedRowsInBin3D, IncludedColumnsInBin)), 2);
            Unaligned3DCycleTraceStdErrors(IncludedRowsInBin3D, APbinIndex, NC-8) = 0;
            
            
            
            IncludedRowsInBin = find(sum(~isnan(NCCycleTraces(:,IncludedColumnsInBin)),2).' > 0);
            IncludedRowsInBin3D = find(sum(~isnan(NC3DCycleTraces(:,IncludedColumnsInBin)),2).' > 0);
            UnalignedCycleMeanTraces(IncludedRowsInBin, APbinIndex, NC-8) = nanmean(NCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin), 2);
            if size(NCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin), 2) > 1
                for r = IncludedRowsInBin
                    UnalignedCycleTraceStdErrors(r, APbinIndex, NC-8) = std(NCCycleTraces(r, IncludedColumnsInBin),'omitnan');
                end
            else
                UnalignedCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8)  = NaN;
            end
                %./...sqrt(length(find(~isnan(NCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin)))));
            UnalignedCycleNumOnNuclei(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(NCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin)), 2);
            if ~isempty(IncludedRowsInBin3D)
                Unaligned3DCycleMeanTraces(IncludedRowsInBin3D, APbinIndex, NC-8) = nanmean(NC3DCycleTraces(IncludedRowsInBin3D, IncludedColumnsInBin), 2);
                
                if size(NC3DCycleTraces(IncludedRowsInBin3D, IncludedColumnsInBin), 2) > 1
                    for r = IncludedRowsInBin3D
                        Unaligned3DCycleTraceStdErrors(r, APbinIndex, NC-8) = std(NC3DCycleTraces(r, IncludedColumnsInBin),'omitnan');
                    end
                 
                else
                    Unaligned3DCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8)  = NaN;
                end
                Unaligned3DCycleNumOnNuclei(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(NC3DCycleTraces(IncludedRowsInBin, IncludedColumnsInBin)), 2);
            else
                Unaligned3DCycleMeanTraces(IncludedRowsInBin3D, APbinIndex, NC-8) =NaN;
                Unaligned3DCycleTraceStdErrors(IncludedRowsInBin3D, APbinIndex, NC-8) = NaN;
                Unaligned3DCycleNumOnNuclei(IncludedRowsInBin, APbinIndex, NC-8) = NaN;
            end
            
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
        % Do anaphase alignment and bin with 3D fluo info
        CurrentNCTraces3D = CycleTraces3D{NC-8};
        IncludedRowsInBin3D = find(sum(~isnan(CurrentNCTraces3D),2).' > 0);
        NCFrameTimes3D = FrameTimes(IncludedRowsInBin3D);
        
        NumBins3D = ceil((max(NCFrameTimes3D)-min(NCFrameTimes3D))/deltaTbinWidth)+1;
        AnaphaseAligned3DCycleFrameTimes{NC-8} = (0:NumBins3D-1)*deltaTbinWidth;
        if length(AnaphaseAligned3DCycleFrameTimes{NC-8}) > size(AnaphaseAligned3DCycleMeanTraces, 1)
            AnaphaseAligned3DCycleMeanTraces = [AnaphaseAligned3DCycleMeanTraces; NaN(length(AnaphaseAligned3DCycleFrameTimes{NC-8})-size(AnaphaseAligned3DCycleMeanTraces, 1),length(APbins),6)];
            AnaphaseAligned3DCycleTraceStdErrors = [AnaphaseAligned3DCycleTraceStdErrors; NaN(length(AnaphaseAligned3DCycleFrameTimes{NC-8})-size(AnaphaseAligned3DCycleTraceStdErrors, 1),length(APbins),6)];
            AnaphaseAligned3DCycleNumNuclei = [AnaphaseAligned3DCycleNumNuclei; NaN(length(AnaphaseAligned3DCycleFrameTimes{NC-8})-size(AnaphaseAligned3DCycleNumNuclei, 1),length(APbins),6)];
            AnaphaseAligned3DCycleNumOnNuclei = [AnaphaseAligned3DCycleNumOnNuclei; NaN(length(AnaphaseAligned3DCycleFrameTimes{NC-8})-size(AnaphaseAligned3DCycleNumOnNuclei, 1),length(APbins),6)];
        end
        Cycle3DTracesWithAnaphaseAlignment{NC-8} = NaN(NumBins3D, size(CurrentNCTraces3D, 2));
        for trace_index = 1:size(CurrentNCTraces3D, 2)
            CurrentTrace = CurrentNCTraces3D(IncludedRowsInBin3D,trace_index);
            TraceIncludedRowsInBin = find(~isnan(CurrentTrace));
            TraceFrameTimes = NCFrameTimes3D(TraceIncludedRowsInBin)-min(NCFrameTimes3D(TraceIncludedRowsInBin));
            CurrentTrace = CurrentTrace(TraceIncludedRowsInBin);
            IncludedInterpolatedFrames = ones(1, length(AnaphaseAligned3DCycleFrameTimes{NC-8}));
            IncludedInterpolatedFrames(AnaphaseAligned3DCycleFrameTimes{NC-8} > max(TraceFrameTimes)) = 0;
            for bin_index = 1:NumBins3D
                if ~isempty(find(AnaphaseAligned3DCycleFrameTimes{NC-8}(bin_index) == TraceFrameTimes, 1))
                    continue
                end
                if IncludedInterpolatedFrames(bin_index) == 0
                    continue
                end
                TimeInterval = TraceFrameTimes(find(TraceFrameTimes > AnaphaseAligned3DCycleFrameTimes{NC-8}(bin_index), 1))-TraceFrameTimes(find(TraceFrameTimes < AnaphaseAligned3DCycleFrameTimes{NC-8}(bin_index), 1, 'last'));
                if TimeInterval > deltaTbinWidth*3
                    IncludedInterpolatedFrames(bin_index) = 0;
                end
            end
            if length(TraceFrameTimes) >= MinTimePoints
                InterpolatedCurrentTrace = interp1(TraceFrameTimes,CurrentTrace,AnaphaseAligned3DCycleFrameTimes{NC-8}(find(IncludedInterpolatedFrames == 1)));
                Cycle3DTracesWithAnaphaseAlignment{NC-8}(find(IncludedInterpolatedFrames == 1), trace_index) =InterpolatedCurrentTrace;
            end
        end
        
        APbinsInMovie = min(ParticleAPbinIndices{NC-8}):max(ParticleAPbinIndices{NC-8});
        AnaphaseAlignedNCCycleTraces = CycleTracesWithAnaphaseAlignment{NC-8};
        AnaphaseAlignedNCCycleTraces(AnaphaseAlignedNCCycleTraces == 0) = NaN;
        AnaphaseAlignedNC3DCycleTraces = Cycle3DTracesWithAnaphaseAlignment{NC-8};
        AnaphaseAlignedNC3DCycleTraces(AnaphaseAlignedNC3DCycleTraces == 0) = NaN;
        for APbinIndex = APbinsInMovie
            IncludedColumnsInBin = find(ParticleAPbinIndices{NC-8} == APbinIndex);
            IncludedRowsInBin = find(sum(~isnan(CycleTracesWithAnaphaseAlignment{NC-8}(:,IncludedColumnsInBin)),2).' > 0);
            IncludedRowsInBin3D = find(sum(~isnan(Cycle3DTracesWithAnaphaseAlignment{NC-8}(:,IncludedColumnsInBin)),2).' > 0);
            AnaphaseAlignedCycleMeanTraces(IncludedRowsInBin, APbinIndex, NC-8) = 0;
            AnaphaseAlignedCycleNumNuclei(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(CycleTracesWithAnaphaseAlignment{NC-8}(IncludedRowsInBin, IncludedColumnsInBin)), 2);
            AnaphaseAlignedCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8) = 0;
            AnaphaseAligned3DCycleMeanTraces(IncludedRowsInBin3D, APbinIndex, NC-8) = 0;
            AnaphaseAligned3DCycleNumNuclei(IncludedRowsInBin3D, APbinIndex, NC-8) = sum(~isnan(Cycle3DTracesWithAnaphaseAlignment{NC-8}(IncludedRowsInBin3D, IncludedColumnsInBin)), 2);
            AnaphaseAligned3DCycleTraceStdErrors(IncludedRowsInBin3D, APbinIndex, NC-8) = 0;
            
            IncludedRowsInBin = find(sum(~isnan(AnaphaseAlignedNCCycleTraces(:,IncludedColumnsInBin)),2).' > 0);
            IncludedRowsInBin3D = find(sum(~isnan(AnaphaseAlignedNC3DCycleTraces(:,IncludedColumnsInBin)),2).' > 0);
            AnaphaseAlignedCycleMeanTraces(IncludedRowsInBin, APbinIndex, NC-8) = nanmean(AnaphaseAlignedNCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin), 2);
            AnaphaseAlignedCycleNumOnNuclei(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(AnaphaseAlignedNCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin)), 2);
            if size(AnaphaseAlignedNCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin), 2) > 1
                for r = IncludedRowsInBin3D
                    AnaphaseAlignedCycleTraceStdErrors(r, APbinIndex, NC-8) = std(AnaphaseAlignedNCCycleTraces(r, IncludedColumnsInBin),'omitnan');
                end
            else
                AnaphaseAlignedCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8)  = NaN;
            end
            
            
            if ~isempty(IncludedRowsInBin3D)
                AnaphaseAligned3DCycleMeanTraces(IncludedRowsInBin3D, APbinIndex, NC-8) = nanmean(AnaphaseAlignedNC3DCycleTraces(IncludedRowsInBin3D, IncludedColumnsInBin), 2);
                AnaphaseAligned3DCycleNumOnNuclei(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(AnaphaseAlignedNC3DCycleTraces(IncludedRowsInBin, IncludedColumnsInBin)), 2);
                if size(AnaphaseAlignedNC3DCycleTraces(IncludedRowsInBin3D, IncludedColumnsInBin), 2) > 1
                    for r = IncludedRowsInBin3D
                        AnaphaseAligned3DCycleTraceStdErrors(r, APbinIndex, NC-8) = std(AnaphaseAlignedNC3DCycleTraces(r, IncludedColumnsInBin),'omitnan');
                    end
                else
                    AnaphaseAligned3DCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8)  = NaN;
                end
            else
                AnaphaseAligned3DCycleMeanTraces(IncludedRowsInBin3D, APbinIndex, NC-8) = NaN;
                AnaphaseAligned3DCycleNumOnNuclei(IncludedRowsInBin, APbinIndex, NC-8) = NaN;
                AnaphaseAligned3DCycleTraceStdErrors(IncludedRowsInBin3D, APbinIndex, NC-8) =NaN;
            end
            
            
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
        % Do anaphase alignment and bin with 3D fluo info
        CurrentNCTraces3D = CycleTraces3D{NC-8};
        IncludedRowsInBin3D = find(sum(~isnan(CurrentNCTraces3D),2).' > 0);
        NCFrameTimes3D = FrameTimes(IncludedRowsInBin3D);
        
        NumBins3D = ceil((max(NCFrameTimes3D)-min(NCFrameTimes3D))/deltaTbinWidth)+1;
        Tbinned3DCycleFrameTimes{NC-8} = (0:NumBins3D-1)*deltaTbinWidth;
        if length(Tbinned3DCycleFrameTimes{NC-8}) > size(Tbinned3DCycleMeanTraces, 1)
            Tbinned3DCycleMeanTraces = [Tbinned3DCycleMeanTraces; NaN(length(Tbinned3DCycleFrameTimes{NC-8}) -size(Tbinned3DCycleMeanTraces, 1),length(APbins),6)];
            Tbinned3DCycleTraceStdErrors = [Tbinned3DCycleTraceStdErrors; NaN(length(Tbinned3DCycleFrameTimes{NC-8}) -size(Tbinned3DCycleTraceStdErrors, 1),length(APbins),6)];
            Tbinned3DCycleNumNuclei = [Tbinned3DCycleNumNuclei; NaN(length(Tbinned3DCycleFrameTimes{NC-8}) -size(Tbinned3DCycleNumNuclei, 1),length(APbins),6)];
            Tbinned3DCycleNumOnNuclei = [Tbinned3DCycleNumOnNuclei; NaN(length(Tbinned3DCycleFrameTimes{NC-8}) -size(Tbinned3DCycleNumOnNuclei, 1),length(APbins),6)];
        end
        Tbinned3DCycleTraces{NC-8} = NaN(NumBins3D, size(CurrentNCTraces3D, 2));
        for trace_index = 1:size(CurrentNCTraces3D, 2)
            CurrentTrace = CurrentNCTraces3D(IncludedRowsInBin3D,trace_index);
            TraceIncludedRowsInBin = find(~isnan(CurrentTrace));
            TraceFrameTimes = NCFrameTimes3D(TraceIncludedRowsInBin)-min(NCFrameTimes3D);
            CurrentTrace = CurrentTrace(TraceIncludedRowsInBin);
            IncludedInterpolatedFrames = ones(1, length(AnaphaseAligned3DCycleFrameTimes{NC-8}));
            IncludedInterpolatedFrames(Tbinned3DCycleFrameTimes{NC-8} > max(TraceFrameTimes)) = 0;
            IncludedInterpolatedFrames(Tbinned3DCycleFrameTimes{NC-8} < min(TraceFrameTimes)) = 0;
            for bin_index = 1:NumBins3D
                if ~isempty(find(Tbinned3DCycleFrameTimes{NC-8}(bin_index) == TraceFrameTimes, 1))
                    continue
                end
                if IncludedInterpolatedFrames(bin_index) == 0
                    continue
                end
                TimeInterval = TraceFrameTimes(find(TraceFrameTimes > Tbinned3DCycleFrameTimes{NC-8}(bin_index), 1))-TraceFrameTimes(find(TraceFrameTimes < Tbinned3DCycleFrameTimes{NC-8}(bin_index), 1, 'last'));
                if TimeInterval > deltaTbinWidth*3
                    IncludedInterpolatedFrames(bin_index) = 0;
                end
            end
            if length(TraceFrameTimes)>= MinTimePoints
                InterpolatedCurrentTrace = interp1(TraceFrameTimes,CurrentTrace,Tbinned3DCycleFrameTimes{NC-8}(find(IncludedInterpolatedFrames == 1)));
                Tbinned3DCycleTraces{NC-8}(find(IncludedInterpolatedFrames == 1), trace_index) =InterpolatedCurrentTrace;
            end
        end
        
        APbinsInMovie = min(ParticleAPbinIndices{NC-8}):max(ParticleAPbinIndices{NC-8});
        TbinnedNCCycleTraces = TbinnedCycleTraces{NC-8};
        TbinnedNCCycleTraces(TbinnedNCCycleTraces == 0) = NaN;
        TbinnedNC3DCycleTraces = Tbinned3DCycleTraces{NC-8};
        TbinnedNC3DCycleTraces(TbinnedNC3DCycleTraces == 0) = NaN;
        for APbinIndex = APbinsInMovie
            IncludedColumnsInBin = find(ParticleAPbinIndices{NC-8} == APbinIndex);
            IncludedRowsInBin = find(sum(~isnan(TbinnedCycleTraces{NC-8}(:,IncludedColumnsInBin)),2).' > 0);
            IncludedRowsInBin3D = find(sum(~isnan(Tbinned3DCycleTraces{NC-8}(:,IncludedColumnsInBin)),2).' > 0);
            TbinnedCycleMeanTraces(IncludedRowsInBin, APbinIndex, NC-8) = 0;
            TbinnedCycleNumNuclei(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(TbinnedCycleTraces{NC-8}(IncludedRowsInBin, IncludedColumnsInBin)), 2);
            TbinnedCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8) = 0;
            Tbinned3DCycleMeanTraces(IncludedRowsInBin3D, APbinIndex, NC-8) = 0;
            Tbinned3DCycleNumNuclei(IncludedRowsInBin3D, APbinIndex, NC-8) = sum(~isnan(Tbinned3DCycleTraces{NC-8}(IncludedRowsInBin3D, IncludedColumnsInBin)), 2);
            Tbinned3DCycleTraceStdErrors(IncludedRowsInBin3D, APbinIndex, NC-8) = 0;
            
            IncludedRowsInBin = find(sum(~isnan(TbinnedNCCycleTraces(:,IncludedColumnsInBin)),2).' > 0);
            IncludedRowsInBin3D = find(sum(~isnan(TbinnedNC3DCycleTraces(:,IncludedColumnsInBin)),2).' > 0);
            TbinnedCycleMeanTraces(IncludedRowsInBin, APbinIndex, NC-8) = nanmean(TbinnedNCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin), 2);
            
            TbinnedCycleNumOnNuclei(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(TbinnedNCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin)), 2);
            if size(TbinnedNCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin), 2) > 1
                for r = IncludedRowsInBin3D
                    TbinnedCycleTraceStdErrors(r, APbinIndex, NC-8) = std(TbinnedNCCycleTraces(r, IncludedColumnsInBin),'omitnan');
                end
            else
                TbinnedCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8)  = NaN;
            end
            
            
            if ~isempty(IncludedRowsInBin3D)
                Tbinned3DCycleMeanTraces(IncludedRowsInBin3D, APbinIndex, NC-8) = nanmean(TbinnedNC3DCycleTraces(IncludedRowsInBin3D, IncludedColumnsInBin), 2);
                
                Tbinned3DCycleNumOnNuclei(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(TbinnedNC3DCycleTraces(IncludedRowsInBin, IncludedColumnsInBin)), 2);
               
                if size(TbinnedNC3DCycleTraces(IncludedRowsInBin3D, IncludedColumnsInBin), 2) > 1
                    for r = IncludedRowsInBin3D
                        Tbinned3DCycleTraceStdErrors(r, APbinIndex, NC-8) = std(TbinnedNC3DCycleTraces(r, IncludedColumnsInBin),'omitnan');
                    end
                else
                    Tbinned3DCycleTraceStdErrors(IncludedRowsInBin3D, APbinIndex, NC-8)  = NaN;
                end
            else
                Tbinned3DCycleMeanTraces(IncludedRowsInBin3D, APbinIndex, NC-8) = NaN;
                
                Tbinned3DCycleNumOnNuclei(IncludedRowsInBin, APbinIndex, NC-8) = NaN;
                Tbinned3DCycleTraceStdErrors(IncludedRowsInBin3D, APbinIndex, NC-8) =NaN;
            end
        end
        
        % Check to make sure calculations make sense
        TestAllTracesMat = squeeze(UnalignedCycleMeanTraces(:,:,NC-8));
        TestAllNumNucMat = squeeze(UnalignedCycleNumNuclei(:,:,NC-8));
        if ~all(TestAllNumNucMat(TestAllTracesMat > 0) > 0)
            disp(['Inconsistencies between Mean Traces and Nuclear count in NC ', num2str(NC), ' for standard mean traces'])
        end
        TestAll3DTracesMat = squeeze(Unaligned3DCycleMeanTraces(:,:,NC-8));
        TestAll3DNumNucMat = squeeze(Unaligned3DCycleNumNuclei(:,:,NC-8));
        if ~all(TestAll3DNumNucMat(TestAll3DTracesMat > 0) > 0)
            disp(['Inconsistencies between Mean Traces and Nuclear count in NC ', num2str(NC), ' for standard 3D mean traces'])
        end
        TestAnaphaseAlignedTracesMat = squeeze(AnaphaseAlignedCycleMeanTraces(:,:,NC-8));
        TestAnaphaseAlignedNumNucMat = squeeze(AnaphaseAlignedCycleNumNuclei(:,:,NC-8));
        if ~all(TestAnaphaseAlignedNumNucMat(TestAnaphaseAlignedTracesMat > 0) > 0)
            disp(['Inconsistencies between Mean Traces and Nuclear count in NC ', num2str(NC), ' for aligned mean traces'])
        end
        TestAnaphaseAligned3DTracesMat = squeeze(AnaphaseAligned3DCycleMeanTraces(:,:,NC-8));
        TestAnaphaseAligned3DNumNucMat = squeeze(AnaphaseAligned3DCycleNumNuclei(:,:,NC-8));
        if ~all(TestAnaphaseAligned3DNumNucMat(TestAnaphaseAligned3DTracesMat > 0) > 0)
            disp(['Inconsistencies between Mean Traces and Nuclear count in NC ', num2str(NC), ' for aligned 3D mean traces'])
        end
        TestTbinnedTracesMat = squeeze(TbinnedCycleMeanTraces(:,:,NC-8));
        TestTbinnedNumNucMat = squeeze(TbinnedCycleNumNuclei(:,:,NC-8));
        if ~all(TestTbinnedNumNucMat(TestTbinnedTracesMat > 0) > 0)
            disp(['Inconsistencies between Mean Traces and Nuclear count in NC ', num2str(NC), ' for t-binned mean traces'])
        end
        TestTbinned3DTracesMat = squeeze(Tbinned3DCycleMeanTraces(:,:,NC-8));
        TestTbinned3DNumNucMat = squeeze(Tbinned3DCycleNumNuclei(:,:,NC-8));
        if ~all(TestTbinned3DNumNucMat(TestTbinned3DTracesMat > 0) > 0)
            disp(['Inconsistencies between Mean Traces and Nuclear count in NC ', num2str(NC), ' for t-binned 3D mean traces'])
        end
        
        
    end
    
    TestLengths1 = [size(UnalignedCycleMeanTraces, 1), size(UnalignedCycleNumNuclei, 1), size(UnalignedCycleNumOnNuclei, 1), size(UnalignedCycleTraceStdErrors, 1)];
    if ~all(TestLengths1 == TestLengths1(1))
        disp(['Inconsistence size for fluo calculations.'])
    end
    
    TestLengths2 = [size(Unaligned3DCycleMeanTraces, 1), size(Unaligned3DCycleNumNuclei, 1), size(Unaligned3DCycleNumOnNuclei, 1), size(Unaligned3DCycleTraceStdErrors, 1)];
    if ~all(TestLengths2 == TestLengths2(1))
        disp(['Inconsistence size for 3D fluo calculations.'])
    end
    
    TestLengths3 = [size(AnaphaseAlignedCycleMeanTraces, 1), size(AnaphaseAlignedCycleNumNuclei, 1), size(AnaphaseAlignedCycleNumOnNuclei, 1), size(AnaphaseAlignedCycleTraceStdErrors, 1)];
    if ~all(TestLengths3 == TestLengths3(1))
        disp(['Inconsistence size for anaphase aligned calculations.'])
    end
    
    TestLengths4 = [size(AnaphaseAligned3DCycleMeanTraces, 1), size(AnaphaseAligned3DCycleNumNuclei, 1), size(AnaphaseAligned3DCycleNumOnNuclei, 1), size(AnaphaseAligned3DCycleTraceStdErrors, 1)];
    if ~all(TestLengths4 == TestLengths4(1))
        disp(['Inconsistence size for 3D anaphase aligned calculations.'])
    end
    
    MeanProfiles = {};
    MeanProfiles.UnalignedCycleMeanTraces = UnalignedCycleMeanTraces;
    MeanProfiles.UnalignedCycleNumNuclei = UnalignedCycleNumNuclei;
    MeanProfiles.UnalignedCycleNumOnNuclei = UnalignedCycleNumOnNuclei;
    MeanProfiles.UnalignedCycleTraceStdErrors = UnalignedCycleTraceStdErrors;
    MeanProfiles.Unaligned3DCycleMeanTraces = Unaligned3DCycleMeanTraces;
    MeanProfiles.Unaligned3DCycleNumNuclei = Unaligned3DCycleNumNuclei;
    MeanProfiles.Unaligned3DCycleNumOnNuclei = Unaligned3DCycleNumOnNuclei;
    MeanProfiles.Unaligned3DCycleTraceStdErrors = Unaligned3DCycleTraceStdErrors;
    MeanProfiles.AnaphaseAlignedCycleMeanTraces = AnaphaseAlignedCycleMeanTraces;
    MeanProfiles.AnaphaseAlignedCycleNumNuclei = AnaphaseAlignedCycleNumNuclei;
    MeanProfiles.AnaphaseAlignedCycleNumOnNuclei = AnaphaseAlignedCycleNumOnNuclei;
    MeanProfiles.AnaphaseAlignedCycleTraceStdErrors = AnaphaseAlignedCycleTraceStdErrors;
    MeanProfiles.AnaphaseAligned3DCycleMeanTraces = AnaphaseAligned3DCycleMeanTraces;
    MeanProfiles.AnaphaseAligned3DCycleNumNuclei = AnaphaseAligned3DCycleNumNuclei;
    MeanProfiles.AnaphaseAligned3DCycleNumOnNuclei = AnaphaseAligned3DCycleNumOnNuclei;
    MeanProfiles.AnaphaseAligned3DCycleTraceStdErrors = AnaphaseAligned3DCycleTraceStdErrors;
    MeanProfiles.TbinnedCycleMeanTraces = TbinnedCycleMeanTraces;
    MeanProfiles.TbinnedCycleNumNuclei = TbinnedCycleNumNuclei;
    MeanProfiles.TbinnedCycleNumOnNuclei = TbinnedCycleNumOnNuclei;
    MeanProfiles.TbinnedCycleTraceStdErrors = TbinnedCycleTraceStdErrors;
    MeanProfiles.Tbinned3DCycleMeanTraces = Tbinned3DCycleMeanTraces;
    MeanProfiles.Tbinned3DCycleNumNuclei = Tbinned3DCycleNumNuclei;
    MeanProfiles.Tbinned3DCycleNumOnNuclei = Tbinned3DCycleNumOnNuclei;
    MeanProfiles.Tbinned3DCycleTraceStdErrors = Tbinned3DCycleTraceStdErrors;
    MeanProfiles.AnaphaseAlignedCycleFrameTimes = AnaphaseAlignedCycleFrameTimes;
    MeanProfiles.AnaphaseAligned3DCycleFrameTimes = AnaphaseAligned3DCycleFrameTimes;
    MeanProfiles.UnalignedCycleFrameTimes = UnalignedCycleFrameTimes;
    MeanProfiles.Unaligned3DCycleFrameTimes = Unaligned3DCycleFrameTimes;
    MeanProfiles.TbinnedCycleFrameTimes = TbinnedCycleFrameTimes;
    MeanProfiles.Tbinned3DCycleFrameTimes = Tbinned3DCycleFrameTimes;
    
    savedVariables = {'UnalignedCycleMeanTraces', 'UnalignedCycleNumNuclei',  'UnalignedCycleNumOnNuclei', 'UnalignedCycleTraceStdErrors',...
        'Unaligned3DCycleMeanTraces','Unaligned3DCycleNumNuclei', 'Unaligned3DCycleNumOnNuclei', 'Unaligned3DCycleTraceStdErrors',...
        'AnaphaseAlignedCycleMeanTraces', 'AnaphaseAlignedCycleNumNuclei', 'AnaphaseAlignedCycleNumOnNuclei','AnaphaseAlignedCycleTraceStdErrors',...
        'AnaphaseAligned3DCycleMeanTraces','AnaphaseAligned3DCycleNumNuclei', 'AnaphaseAligned3DCycleNumOnNuclei','AnaphaseAligned3DCycleTraceStdErrors',...
        'TbinnedCycleMeanTraces', 'TbinnedCycleNumNuclei','TbinnedCycleNumOnNuclei', 'TbinnedCycleTraceStdErrors',...
        'Tbinned3DCycleMeanTraces','Tbinned3DCycleNumNuclei','Tbinned3DCycleNumOnNuclei',  'Tbinned3DCycleTraceStdErrors',...
        'UnalignedCycleFrameTimes', 'Unaligned3DCycleFrameTimes', 'AnaphaseAlignedCycleFrameTimes', 'AnaphaseAligned3DCycleFrameTimes',...
        'TbinnedCycleFrameTimes', 'Tbinned3DCycleFrameTimes'};
    
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
end









