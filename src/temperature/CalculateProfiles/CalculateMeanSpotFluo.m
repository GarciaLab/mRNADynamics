function this = CalculateMeanSpotFluo(this, verbose)
deltaTbinWidth = this.time_delta;
MinTimePoints = this.MinimumTimePoints;


if ~exist('verbose', 'var')
    verbose = true;
end
%%
NumSets = length(this.Experiments);
APResolution = this.Experiments{1}.APResolution;
NumAPbins = uint16(1/APResolution) + 1;
APbins = 0:APResolution:1;
this.MeanSpotFluo = {};

traceNames = {'AnaphaseAligned', 'AnaphaseAligned3D', 'Unaligned', 'Unaligned3D', 'Tbinned', 'Tbinned3D'};
for trIndex = 1:length(traceNames)
    this.MeanSpotFluo.(traceNames{trIndex}) = NaN(NumSets, NumAPbins, 6);
    this.MeanSpotFluo.([traceNames{trIndex}, 'StdError']) = NaN(NumSets, NumAPbins, 6);
    this.MeanSpotFluo.([traceNames{trIndex}, 'NumOnNuclei']) = NaN(NumSets, NumAPbins, 6);
end
%%
for SetIndex = 1:NumSets
    if ~ismember(SetIndex, this.ProcessedExperiments)
        continue
    end
    Prefix = this.ExperimentPrefixes{SetIndex};
    liveExperiment = this.Experiments{SetIndex};
    
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
        
        
        CycleTraces = {};
        CycleTraces3D = {};
        MeanAPs = {};
        ParticleAPbinIndices = {};
        UnalignedCycleFrameTimes = {};
        Unaligned3DCycleFrameTimes = {};
        
        numBinnedFrames = ceil((max(FrameTimes)-FrameTimes(nc_info(6)))/deltaTbinWidth)+1;
        CycleTracesWithAnaphaseAlignment = {};
        Cycle3DTracesWithAnaphaseAlignment = {};
        AnaphaseAlignedCycleFrameTimes = {};
        AnaphaseAligned3DCycleFrameTimes = {};
        
        TbinnedCycleTraces = {};
        Tbinned3DCycleTraces = {};
        TbinnedCycleFrameTimes = {};
        Tbinned3DCycleFrameTimes = {};
        
        
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
            
            UnalignedMinTimes = this.TimeOns.Unaligned(SetIndex,:,NC-8)-this.TimeOns.UnalignedStdError(SetIndex,:,NC-8);
            Unaligned3DMinTimes = this.TimeOns.Unaligned3D(SetIndex,:,NC-8)-this.TimeOns.Unaligned3DStdError(SetIndex,:,NC-8);
            
            for APbinIndex = APbinsInMovie
                IncludedColumnsInBin = find(ParticleAPbinIndices{NC-8} == APbinIndex);
                NCCycleTraces_AP = NCCycleTraces(:,IncludedColumnsInBin);
                ExcludedRows = UnalignedCycleFrameTimes{NC-8}/60 < UnalignedMinTimes(APbinIndex);
                NCCycleTraces_AP(ExcludedRows, :) = NaN;
                NCCycleTraces_AP = reshape(NCCycleTraces_AP, [1, size(NCCycleTraces_AP, 1)*size(NCCycleTraces_AP, 2)]);
                NCCycleTraces_AP = NCCycleTraces_AP(~isnan(NCCycleTraces_AP));
                if ~isempty(NCCycleTraces_AP)
                    this.MeanSpotFluo.Unaligned(SetIndex, APbinIndex, NC-8) = mean(NCCycleTraces_AP);
                    this.MeanSpotFluo.UnalignedStdError(SetIndex, APbinIndex, NC-8) = std(NCCycleTraces_AP)/sqrt(length(NCCycleTraces_AP));
                    this.MeanSpotFluo.UnalignedNumOnNuclei(SetIndex, APbinIndex, NC-8) = length(NCCycleTraces_AP);
                end
                
                NC3DCycleTraces_AP = NC3DCycleTraces(:,IncludedColumnsInBin);
                ExcludedRows = Unaligned3DCycleFrameTimes{NC-8}/60 < Unaligned3DMinTimes(APbinIndex);
                NC3DCycleTraces_AP(ExcludedRows, :) = NaN;
                NC3DCycleTraces_AP = reshape(NC3DCycleTraces_AP, [1, size(NC3DCycleTraces_AP, 1)*size(NC3DCycleTraces_AP, 2)]);
                NC3DCycleTraces_AP = NC3DCycleTraces_AP(~isnan(NC3DCycleTraces_AP));
                if ~isempty(NC3DCycleTraces_AP)
                    this.MeanSpotFluo.Unaligned3D(SetIndex, APbinIndex, NC-8) = mean(NC3DCycleTraces_AP);
                    this.MeanSpotFluo.Unaligned3DStdError(SetIndex, APbinIndex, NC-8) = std(NC3DCycleTraces_AP)/sqrt(length(NC3DCycleTraces_AP));
                    this.MeanSpotFluo.Unaligned3DNumOnNuclei(SetIndex, APbinIndex, NC-8) = length(NC3DCycleTraces_AP);
                end
                
                
            end
            
            
            
            
            % Do anaphase alignment and bin
            CurrentNCTraces = CycleTraces{NC-8};
            IncludedRowsInBin = find(sum(~isnan(CurrentNCTraces),2).' > 0);
            NCFrameTimes = FrameTimes(IncludedRowsInBin);
            NumBins = ceil((max(NCFrameTimes)-min(NCFrameTimes))/deltaTbinWidth)+1;
            AnaphaseAlignedCycleFrameTimes{NC-8} = (0:NumBins-1)*deltaTbinWidth;
       
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
            
            AnaphaseAlignedMinTimes = this.TimeOns.AnaphaseAligned(SetIndex,:,NC-8)-this.TimeOns.AnaphaseAlignedStdError(SetIndex,:,NC-8);
            AnaphaseAligned3DMinTimes = this.TimeOns.AnaphaseAligned3D(SetIndex,:,NC-8)-this.TimeOns.AnaphaseAligned3DStdError(SetIndex,:,NC-8);
            
            
            
            for APbinIndex = APbinsInMovie
                IncludedColumnsInBin = find(ParticleAPbinIndices{NC-8} == APbinIndex);
                NCCycleTraces_AP = AnaphaseAlignedNCCycleTraces(:,IncludedColumnsInBin);
                ExcludedRows = AnaphaseAlignedCycleFrameTimes{NC-8}/60 < AnaphaseAlignedMinTimes(APbinIndex);
                NCCycleTraces_AP(ExcludedRows, :) = NaN;
                NCCycleTraces_AP = reshape(NCCycleTraces_AP, [1, size(NCCycleTraces_AP, 1)*size(NCCycleTraces_AP, 2)]);
                NCCycleTraces_AP = NCCycleTraces_AP(~isnan(NCCycleTraces_AP));
                if ~isempty(NCCycleTraces_AP)
                    this.MeanSpotFluo.AnaphaseAligned(SetIndex, APbinIndex, NC-8) = mean(NCCycleTraces_AP);
                    this.MeanSpotFluo.AnaphaseAlignedStdError(SetIndex, APbinIndex, NC-8) = std(NCCycleTraces_AP)/sqrt(length(NCCycleTraces_AP));
                    this.MeanSpotFluo.AnaphaseAlignedNumOnNuclei(SetIndex, APbinIndex, NC-8) = length(NCCycleTraces_AP);
                end
                
                NC3DCycleTraces_AP = AnaphaseAlignedNC3DCycleTraces(:,IncludedColumnsInBin);
                ExcludedRows = AnaphaseAligned3DCycleFrameTimes{NC-8}/60 < AnaphaseAligned3DMinTimes(APbinIndex);
                NC3DCycleTraces_AP(ExcludedRows, :) = NaN;
                NC3DCycleTraces_AP = reshape(NC3DCycleTraces_AP, [1, size(NC3DCycleTraces_AP, 1)*size(NC3DCycleTraces_AP, 2)]);
                NC3DCycleTraces_AP = NC3DCycleTraces_AP(~isnan(NC3DCycleTraces_AP));
                if ~isempty(NC3DCycleTraces_AP)
                    this.MeanSpotFluo.AnaphaseAligned3D(SetIndex, APbinIndex, NC-8) = mean(NC3DCycleTraces_AP);
                    this.MeanSpotFluo.AnaphaseAligned3DStdError(SetIndex, APbinIndex, NC-8) = std(NC3DCycleTraces_AP)/sqrt(length(NC3DCycleTraces_AP));
                    this.MeanSpotFluo.AnaphaseAligned3DNumOnNuclei(SetIndex, APbinIndex, NC-8) = length(NC3DCycleTraces_AP);
                end
                
                
            end
            
            % Bin without anaphase alignment
            CurrentNCTraces = CycleTraces{NC-8};
            IncludedRowsInBin = find(sum(~isnan(CurrentNCTraces),2).' > 0);
            NCFrameTimes = FrameTimes(IncludedRowsInBin);
            NumBins = ceil((max(NCFrameTimes)-min(NCFrameTimes))/deltaTbinWidth)+1;
            TbinnedCycleFrameTimes{NC-8} = (0:NumBins-1)*deltaTbinWidth;
 
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
            
            TbinnedMinTimes = this.TimeOns.Tbinned(SetIndex,:,NC-8)-this.TimeOns.TbinnedStdError(SetIndex,:,NC-8);
            Tbinned3DMinTimes = this.TimeOns.Tbinned3D(SetIndex,:,NC-8)-this.TimeOns.Tbinned3DStdError(SetIndex,:,NC-8);
            
            
            
            for APbinIndex = APbinsInMovie
                IncludedColumnsInBin = find(ParticleAPbinIndices{NC-8} == APbinIndex);
                NCCycleTraces_AP = TbinnedNCCycleTraces(:,IncludedColumnsInBin);
                ExcludedRows = TbinnedCycleFrameTimes{NC-8}/60 < TbinnedMinTimes(APbinIndex);
                NCCycleTraces_AP(ExcludedRows, :) = NaN;
                NCCycleTraces_AP = reshape(NCCycleTraces_AP, [1, size(NCCycleTraces_AP, 1)*size(NCCycleTraces_AP, 2)]);
                NCCycleTraces_AP = NCCycleTraces_AP(~isnan(NCCycleTraces_AP));
                if ~isempty(NCCycleTraces_AP)
                    this.MeanSpotFluo.Tbinned(SetIndex, APbinIndex, NC-8) = mean(NCCycleTraces_AP);
                    this.MeanSpotFluo.TbinnedStdError(SetIndex, APbinIndex, NC-8) = std(NCCycleTraces_AP)/sqrt(length(NCCycleTraces_AP));
                    this.MeanSpotFluo.TbinnedNumOnNuclei(SetIndex, APbinIndex, NC-8) = length(NCCycleTraces_AP);
                end
                
                NC3DCycleTraces_AP = TbinnedNC3DCycleTraces(:,IncludedColumnsInBin);
                ExcludedRows = Tbinned3DCycleFrameTimes{NC-8}/60 < Tbinned3DMinTimes(APbinIndex);
                NC3DCycleTraces_AP(ExcludedRows, :) = NaN;
                NC3DCycleTraces_AP = reshape(NC3DCycleTraces_AP, [1, size(NC3DCycleTraces_AP, 1)*size(NC3DCycleTraces_AP, 2)]);
                NC3DCycleTraces_AP = NC3DCycleTraces_AP(~isnan(NC3DCycleTraces_AP));
                if ~isempty(NC3DCycleTraces_AP)
                    this.MeanSpotFluo.Tbinned3D(SetIndex, APbinIndex, NC-8) = mean(NC3DCycleTraces_AP);
                    this.MeanSpotFluo.Tbinned3DStdError(SetIndex, APbinIndex, NC-8) = std(NC3DCycleTraces_AP)/sqrt(length(NC3DCycleTraces_AP));
                    this.MeanSpotFluo.Tbinned3DNumOnNuclei(SetIndex, APbinIndex, NC-8) = length(NC3DCycleTraces_AP);
                end
                
                
            end
            
            
            
            
        end
    end
    
end







