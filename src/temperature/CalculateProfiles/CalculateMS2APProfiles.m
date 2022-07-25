function MeanProfiles = CalculateMS2APProfiles(Prefix, deltaTbinWidth, MinTimePoints,UseManualApproval, verbose)
if ~exist('deltaTbinWidth', 'var')
    deltaTbinWidth = 30;
elseif isempty(deltaTbinWidth)
    deltaTbinWidth = 30;
end
if ~exist('MinTimePoints', 'var')
    MinTimePoints = 4;
elseif isempty(MinTimePoints)
    MinTimePoints = 4;
end
if ~exist('UseManualApproval', 'var')
    UseManualApproval = false;
end
if ~exist('verbose', 'var')
    verbose = true;
end

liveExperiment = LiveExperiment(Prefix);
schnitzcells = getSchnitzcells(liveExperiment);
HasHistone = false;
for i = 1:length(liveExperiment.Channels)
    if contains(lower(liveExperiment.Channels{i}), 'his')
        HasHistone = true;
    end
end
APbins = 0:liveExperiment.APResolution:1;
APResolution = liveExperiment.APResolution;
FrameInfo = getFrameInfo(liveExperiment);
anaphaseFrames = liveExperiment.anaphaseFrames.';
anaphaseFrames = [anaphaseFrames length([FrameInfo(:).Time])];
NCFrames = cell(1, 6);
for i = 1:6
    if anaphaseFrames(i) > 0
        NCFrames{i} = anaphaseFrames(i):anaphaseFrames(i+1);
    end
end

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
    % Store schnitz info
    schnitzCycles = NaN(1, length(schnitzcells));
    schnitzFlags = NaN(1, length(schnitzcells));
    schnitzApproved = NaN(1, length(schnitzcells));
    schnitzAnaphaseFrames = NaN(1, length(schnitzcells));
    schnitzInferredAnaphaseFrames = zeros(1, length(schnitzcells), 'logical');
    schnitzFirstFrameContainedInTrace =zeros(1, length(schnitzcells), 'logical');
    schnitzLastFrameContainedInTrace = zeros(1, length(schnitzcells), 'logical');
    schnitzAPpos = NaN(1, length(schnitzcells));
    schnitzFractionFramesApproved =NaN(1, length(schnitzcells));
    schnitzAllFramesApproved = zeros(1, length(schnitzcells), 'logical');
    NumSchnitz = length(schnitzcells);
    AllSchnitzIndices = 1:NumSchnitz;
    for sc_idx = 1:length(schnitzcells)
        if ~isempty(schnitzcells(sc_idx).cycle)
            schnitzCycles(sc_idx) = schnitzcells(sc_idx).cycle;
        end
        if ~isempty(schnitzcells(sc_idx).Flag)
            schnitzFlags(sc_idx) = schnitzcells(sc_idx).Flag;
        end
        if ~isempty(schnitzcells(sc_idx).Approved)
            schnitzApproved(sc_idx) = schnitzcells(sc_idx).Approved;
        end
        if ~isempty(schnitzcells(sc_idx).anaphaseFrame)
            schnitzAnaphaseFrames(sc_idx) = schnitzcells(sc_idx).anaphaseFrame;
        end
        if ~isempty(schnitzcells(sc_idx).inferredAnaphaseFrame)
            schnitzInferredAnaphaseFrames(sc_idx) = schnitzcells(sc_idx).inferredAnaphaseFrame;
        end
        if ~isempty(schnitzcells(sc_idx).containsFirstFrameOfCycle)
            schnitzFirstFrameContainedInTrace(sc_idx) = schnitzcells(sc_idx).containsFirstFrameOfCycle;
        end
        if ~isempty(schnitzcells(sc_idx).containsLastFrameOfCycle)
            schnitzLastFrameContainedInTrace(sc_idx) = schnitzcells(sc_idx).containsLastFrameOfCycle;
        end
        if ~isempty(schnitzcells(sc_idx).APpos)
            schnitzAPpos(sc_idx) = mean(schnitzcells(sc_idx).APpos);
        end
        if ~isempty(schnitzcells(sc_idx).VelocityInfo.SchnitzHasAllFrames)
            schnitzAllFramesApproved(sc_idx) = [schnitzcells(sc_idx).VelocityInfo.SchnitzHasAllFrames];
        end
        if ~isempty(schnitzcells(sc_idx).FrameApproved)
            schnitzFractionFramesApproved(sc_idx) = sum(schnitzcells(sc_idx).FrameApproved)/length(schnitzcells(sc_idx).FrameApproved);
        end
    end
    
    ApprovedCount = zeros(1, length(CompiledParticles{ChN}));
    AllApproved = zeros(1, length(CompiledParticles{ChN}));
    FirstApproved = zeros(1, length(CompiledParticles{ChN}));
    PercentageFramesApproved = zeros(1, length(CompiledParticles{ChN}));
    HasSpotsBeforeDisapproved = zeros(1, length(CompiledParticles{ChN}));
    HasMinimumTimePoints = zeros(1, length(CompiledParticles{ChN}));
    ApprovedVector = zeros(1, length(CompiledParticles{ChN}));
    CPcycles = NaN(1, length(CompiledParticles{ChN}));
    ParticleSchnitzIndices = NaN(1, length(CompiledParticles{ChN}));
    ConstantlyPresentNucleiVector = zeros(1,length(CompiledParticles{ChN}), 'logical');
    for cp=1:length(CompiledParticles{ChN})
        if (CompiledParticles{ChN}(cp).ManualApproved >= 1) & (CompiledParticles{ChN}(cp).schnitzcell.Approved > 0) & ...
                (CompiledParticles{ChN}(cp).schnitzcell.Flag ~= 6)
            ApprovedVector(cp) = 1;
        end
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
        
        if (AllApproved(cp) == 1) & (ApprovedCount(cp) == 1) &  (ApprovedVector(cp) == 1) & ...
                (CompiledParticles{ChN}(cp).schnitzcell.containsFirstFrameOfCycle) & ...
                (PercentageFramesApproved(cp) ==1)
            ConstantlyPresentNucleiVector(cp) = true;
        end
        if ~isempty(CompiledParticles{ChN}(cp).cycle)
            CPcycles(cp) = CompiledParticles{ChN}(cp).cycle;
        end
        if ~isempty(CompiledParticles{ChN}(cp).schnitz)
            ParticleSchnitzIndices(cp) = CompiledParticles{ChN}(cp).schnitz;
        end
    end
    
    % Not supported for NChannels > 1
    SchnitzUnalignedTotalNuclei = zeros(numFrames, length(APbins), 6, NumSchnitz);
    SchnitzUnalignedOffNuclei = zeros(numFrames, length(APbins), 6, NumSchnitz);
    SchnitzUnalignedOnNuclei = zeros(numFrames, length(APbins), 6, NumSchnitz);
    SchnitzUnalignedQuiescentNuclei = zeros(numFrames, length(APbins), 6, NumSchnitz);
    SchnitzUnalignedFinishedTranscribingNuclei =  zeros(numFrames, length(APbins), 6, NumSchnitz);
    SchnitzUnalignedMeanTraces = NaN(numFrames, length(APbins), 6, NumSchnitz);
    SchnitzUnalignedTraceCount = zeros(numFrames, length(APbins), 6, NumSchnitz);
    SchnitzUnaligned3DMeanTraces = NaN(numFrames, length(APbins), 6, NumSchnitz);
    SchnitzUnaligned3DTraceCount = zeros(numFrames, length(APbins), 6, NumSchnitz);
    
    sc_approved_p_not = zeros(1, length(schnitzcells), 'logical');
    for sc_idx = 1:length(schnitzcells)%496:496
        
        if (schnitzFlags(sc_idx) ~= 0) | (schnitzApproved(sc_idx) ~= 1) | ...
                ~schnitzFirstFrameContainedInTrace(sc_idx) | ~schnitzLastFrameContainedInTrace(sc_idx) | ...
                (schnitzFractionFramesApproved(sc_idx) < 1)
            continue
        end
        scAP = schnitzAPpos(sc_idx);
        scAPvector = schnitzcells(sc_idx).APpos;
        scAPbin = find((scAP < APbins+APResolution/2) & (scAP >= APbins-APResolution/2));
        scNC = schnitzCycles(sc_idx);
        SchnitzFrames = schnitzcells(sc_idx).frames.';
        CandidateAnaphaseFrames = schnitzAnaphaseFrames(schnitzCycles == scNC & abs(schnitzAnaphaseFrames-anaphaseFrames(scNC-8)) <= 3 & schnitzAnaphaseFrames > 0);
        if ~isempty(CandidateAnaphaseFrames)
            SchnitzFrameMin = min([CandidateAnaphaseFrames anaphaseFrames(scNC-8)]);
        else
            SchnitzFrameMin = anaphaseFrames(scNC-8);
        end
        SchnitzHasParticle = ismember(sc_idx, ParticleSchnitzIndices);
        if SchnitzHasParticle
            p_idx = find(ParticleSchnitzIndices == sc_idx);
            if length(p_idx) > 1
                continue
            elseif length(p_idx) == 0
                if (scNC >= 9) & (scNC <= 14) & ~isempty(scAPvector)
                    for fr_idx = 1:length(SchnitzFrames)
                        scAPbin = find((scAPvector(fr_idx) < APbins+APResolution/2) & (scAPvector(fr_idx) >= APbins-APResolution/2));
                        SchnitzUnalignedTotalNuclei(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                        SchnitzUnalignedOffNuclei(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                        SchnitzUnalignedMeanTraces(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 0;
                        SchnitzUnaligned3DMeanTraces(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 0;
                        SchnitzUnalignedTraceCount(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 0;
                        SchnitzUnaligned3DTraceCount(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 0;
                    end
                end
            elseif length(p_idx) == 1
                if CompiledParticles{ChN}(p_idx).ManualApproved
                    if (~all(CompiledParticles{ChN}(p_idx).FlaggingInfo.PositionApproved) | ~all(CompiledParticles{ChN}(p_idx).FlaggingInfo.FluoApproved))
                        sc_approved_p_not(sc_idx) = true;
                    else
                        pFluoVector = CompiledParticles{ChN}(p_idx).FlaggingInfo.MaxSpotFluoLevel;
                        pFrames = CompiledParticles{ChN}(p_idx).FlaggingInfo.TrueFrames;
                        pOrigFluoVector = CompiledParticles{ChN}(p_idx).Fluo;
                        pOrigFrames = CompiledParticles{ChN}(p_idx).Frame;
                        pOrigFluo3DVector = CompiledParticles{ChN}(p_idx).Fluo3DGauss;
                        if sum(~isnan(pFluoVector)) < MinTimePoints
                            for fr_idx = 1:length(SchnitzFrames)
                                scAPbin = find((scAPvector(fr_idx) < APbins+APResolution/2) & (scAPvector(fr_idx) >= APbins-APResolution/2));
                                SchnitzUnalignedTotalNuclei(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                                SchnitzUnalignedOffNuclei(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                                SchnitzUnalignedMeanTraces(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 0;
                                SchnitzUnaligned3DMeanTraces(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 0;
                                SchnitzUnalignedTraceCount(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 0;
                                SchnitzUnaligned3DTraceCount(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 0;
                            end
                            continue
                        end
                        
                        
                        for fr_idx = 1:length(SchnitzFrames)
                            scAPbin = find((scAPvector(fr_idx) < APbins+APResolution/2) & (scAPvector(fr_idx) >= APbins-APResolution/2));
                            SchnitzUnalignedTotalNuclei(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                            pfr_idx = find(pFrames == SchnitzFrames(fr_idx));
                            orig_pfr_idx = find(pOrigFrames == SchnitzFrames(fr_idx));
                            if ~isempty(pfr_idx)
                                if ~isnan(pFluoVector(pfr_idx))
                                    SchnitzUnalignedOnNuclei(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                                elseif sum(pFrames > SchnitzFrames(fr_idx) & ~isnan(pFluoVector)) > 0
                                    SchnitzUnalignedQuiescentNuclei(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                                elseif sum(pFrames < SchnitzFrames(fr_idx) & ~isnan(pFluoVector)) > 0
                                    SchnitzUnalignedFinishedTranscribingNuclei(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                                else
                                    SchnitzUnalignedOffNuclei(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                                end
                            else
                                if sum(pFrames > SchnitzFrames(fr_idx) & ~isnan(pFluoVector)) > 0
                                    SchnitzUnalignedQuiescentNuclei(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                                elseif sum(pFrames < SchnitzFrames(fr_idx) & ~isnan(pFluoVector)) > 0
                                    SchnitzUnalignedFinishedTranscribingNuclei(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                                else
                                    SchnitzUnalignedOffNuclei(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-6, sc_idx) = 1;
                                end
                            end
                            
                            if ~isempty(orig_pfr_idx)
                                if ~isnan(pOrigFluoVector(orig_pfr_idx))
                                    SchnitzUnalignedMeanTraces(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = pOrigFluoVector(orig_pfr_idx);
                                    SchnitzUnalignedTraceCount(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                                else
                                    SchnitzUnalignedMeanTraces(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 0;
                                end
                                
                            else
                                SchnitzUnalignedMeanTraces(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 0;
                            end
                            
                            if ~isempty(orig_pfr_idx)
                                if ~isnan(pOrigFluo3DVector(orig_pfr_idx))
                                    SchnitzUnaligned3DMeanTraces(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = pOrigFluo3DVector(orig_pfr_idx);
                                    SchnitzUnaligned3DTraceCount(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                                else
                                    SchnitzUnaligned3DMeanTraces(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 0;
                                end
                            else
                                SchnitzUnaligned3DMeanTraces(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 0;
                            end
                        end
                        
                    end
                else % Fill in what to do if particle is there but is not approved!
                    if length(CompiledParticles{ChN}(p_idx).Frame) < MinTimePoints | length(CompiledParticles{ChN}(p_idx).Frame)/length(CompiledParticles{ChN}(p_idx).TrueFrames) < .1
                        if (scNC >= 9) & (scNC <= 14) & ~isempty(scAPvector)
                            for fr_idx = 1:length(SchnitzFrames)
                                scAPbin = find((scAPvector(fr_idx) < APbins+APResolution/2) & (scAPvector(fr_idx) >= APbins-APResolution/2));
                                SchnitzUnalignedTotalNuclei(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                                SchnitzUnalignedOffNuclei(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                                SchnitzUnalignedMeanTraces(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 0;
                                SchnitzUnaligned3DMeanTraces(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 0;
                            end
                        end
                    end
                end
                
                
            end
        else
            if (scNC >= 9) & (scNC <= 14) & ~isempty(scAPvector)
                for fr_idx = 1:length(SchnitzFrames)
                    scAPbin = find((scAPvector(fr_idx) < APbins+APResolution/2) & (scAPvector(fr_idx) >= APbins-APResolution/2));
                    SchnitzUnalignedTotalNuclei(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                    SchnitzUnalignedOffNuclei(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                    SchnitzUnalignedMeanTraces(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 0;
                    SchnitzUnaligned3DMeanTraces(SchnitzFrames(fr_idx)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 0;
                end
            end
        end
        
        
    end
    
    CandidateAnaphaseFrames = schnitzAnaphaseFrames(schnitzCycles == 14 & abs(schnitzAnaphaseFrames-anaphaseFrames(14-8)) <= 3 & schnitzAnaphaseFrames > 0);
    if ~isempty(CandidateAnaphaseFrames)
        SchnitzFrameMin = min([CandidateAnaphaseFrames anaphaseFrames(14-8)]);
    else
        SchnitzFrameMin = anaphaseFrames(14-8);
    end
    numBinnedFrames = ceil((max(FrameTimes)-FrameTimes(SchnitzFrameMin))/deltaTbinWidth)+1;
    SchnitzTbinnedTotalNuclei = zeros(numBinnedFrames, length(APbins), 6, NumSchnitz);
    SchnitzTbinnedOffNuclei = zeros(numBinnedFrames, length(APbins), 6, NumSchnitz);
    SchnitzTbinnedOnNuclei = zeros(numBinnedFrames, length(APbins), 6, NumSchnitz);
    SchnitzTbinnedQuiescentNuclei = zeros(numBinnedFrames, length(APbins), 6, NumSchnitz);
    SchnitzTbinnedFinishedTranscribingNuclei =  zeros(numBinnedFrames, length(APbins), 6, NumSchnitz);
    SchnitzTbinnedMeanTraces = NaN(numBinnedFrames, length(APbins), 6, NumSchnitz);
    SchnitzTbinnedTraceCount = zeros(numBinnedFrames, length(APbins), 6, NumSchnitz);
    SchnitzTbinned3DMeanTraces = NaN(numBinnedFrames, length(APbins), 6, NumSchnitz);
    SchnitzTbinned3DTraceCount = zeros(numBinnedFrames, length(APbins), 6, NumSchnitz);
    
    TbinnedTimeVector = 0:deltaTbinWidth:(max(FrameTimes)-FrameTimes(nc_info(6)));
    if length(TbinnedTimeVector) < numBinnedFrames
        TbinnedTimeVector(end+1) = max(TbinnedTimeVector)+deltaTbinWidth;
    end
    TbinnedTimeVecIndex = 1:length(TbinnedTimeVector);
    for sc_idx = 1:length(schnitzcells)%496:496
        
        if  (schnitzFlags(sc_idx) ~= 0) | (schnitzApproved(sc_idx) ~= 1) | ...
                ~schnitzFirstFrameContainedInTrace(sc_idx) | ~schnitzLastFrameContainedInTrace(sc_idx) | ...
                (schnitzFractionFramesApproved(sc_idx) < 0.9)
            continue
        end
        scAP = schnitzAPpos(sc_idx);
        scAPvector = schnitzcells(sc_idx).APpos;
        scAPbin = find((scAP < APbins+APResolution/2) & (scAP >= APbins-APResolution/2));
        scNC = schnitzCycles(sc_idx);
        SchnitzFrames = schnitzcells(sc_idx).frames.';
        CandidateAnaphaseFrames = schnitzAnaphaseFrames(schnitzCycles == scNC & abs(schnitzAnaphaseFrames-anaphaseFrames(scNC-8)) <= 3 & schnitzAnaphaseFrames > 0);
        if ~isempty(CandidateAnaphaseFrames)
            SchnitzFrameMin = min([CandidateAnaphaseFrames anaphaseFrames(scNC-8)]);
        else
            SchnitzFrameMin = anaphaseFrames(scNC-8);
        end
        SchnitzTimes = [FrameInfo(SchnitzFrames).Time]-FrameInfo(SchnitzFrameMin).Time;
        
        SchnitzHasParticle = ismember(sc_idx, ParticleSchnitzIndices);
        if SchnitzHasParticle
            p_idx = find(ParticleSchnitzIndices == sc_idx);
            if length(p_idx) > 1
                continue
            elseif length(p_idx) == 0
                if (scNC >= 9) & (scNC <= 14) & ~isempty(scAPvector)
                    InterpAPvec = interp1(SchnitzTimes,scAPvector,TbinnedTimeVector);
                    for fr_idx = 1:length(TbinnedTimeVector)
                        if ~isnan(InterpAPvec(fr_idx))
                            scAPbin = find((InterpAPvec(fr_idx) < APbins+APResolution/2) & (InterpAPvec(fr_idx) >= APbins-APResolution/2));
                            SchnitzTbinnedTotalNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                            SchnitzTbinnedOffNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                            SchnitzTbinnedMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = 0;
                            SchnitzTbinned3DMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = 0;
                        end
                    end
                end
            elseif length(p_idx) == 1
                if CompiledParticles{ChN}(p_idx).ManualApproved
                    if (~all(CompiledParticles{ChN}(p_idx).FlaggingInfo.PositionApproved) | ~all(CompiledParticles{ChN}(p_idx).FlaggingInfo.FluoApproved))
                        continue
                    else
                        pFluoVector = CompiledParticles{ChN}(p_idx).FlaggingInfo.MaxSpotFluoLevel;
                        pFrames = CompiledParticles{ChN}(p_idx).FlaggingInfo.TrueFrames;
                        pOrigFluoVector = CompiledParticles{ChN}(p_idx).Fluo;
                        pOrigFrames = CompiledParticles{ChN}(p_idx).Frame;
                        pOrigFluo3DVector = CompiledParticles{ChN}(p_idx).Fluo3DGauss;
                        if sum(~isnan(pFluoVector)) < MinTimePoints
                            if length(scAPvector) <= 1
                                continue
                            end
                            InterpAPvec = interp1(SchnitzTimes,scAPvector,TbinnedTimeVector);
                            for fr_idx = 1:length(TbinnedTimeVector)
                                if ~isnan(InterpAPvec(fr_idx))
                                    scAPbin = find((InterpAPvec(fr_idx) < APbins+APResolution/2) & (InterpAPvec(fr_idx) >= APbins-APResolution/2));
                                    SchnitzTbinnedTotalNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                                    SchnitzTbinnedOffNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                                    SchnitzTbinnedMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = 0;
                                    SchnitzTbinned3DMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = 0;
                                end
                            end
                            continue
                        end
                        
                        FullFluoVector = NaN(1, length(min(SchnitzFrames):max(SchnitzFrames)));
                        FullOrigFluoVector = NaN(1, length(min(SchnitzFrames):max(SchnitzFrames)));
                        FullOrigFluo3DVector = NaN(1, length(min(SchnitzFrames):max(SchnitzFrames)));
                        if min(pFrames) < min(SchnitzFrames) | max(pFrames) > max(SchnitzFrames)
                            continue
                        end
                        if min(pOrigFrames) < min(SchnitzFrames) | max(pOrigFrames) > max(SchnitzFrames)
                            continue
                        end
                        FullFluoVector(find(ismember(min(SchnitzFrames):max(SchnitzFrames), pFrames))) = pFluoVector;
                        FullOrigFluoVector(find(ismember(min(SchnitzFrames):max(SchnitzFrames), pOrigFrames))) = pOrigFluoVector;
                        FullOrigFluo3DVector(find(ismember(min(SchnitzFrames):max(SchnitzFrames), pOrigFrames))) = pOrigFluo3DVector;
                        if length(scAPvector) <= 1
                            continue
                        end
                        InterpAPvec = interp1(SchnitzTimes,scAPvector,TbinnedTimeVector);
                        SchnitzPatchedTime = [FrameInfo(min(SchnitzFrames):max(SchnitzFrames)).Time]-FrameInfo(anaphaseFrames(scNC-8)).Time;
                        InterpFluoVec = interp1(SchnitzPatchedTime,FullFluoVector,TbinnedTimeVector);
                        InterpOrigFluoVec = interp1(SchnitzPatchedTime,FullOrigFluoVector,TbinnedTimeVector);
                        InterpOrigFluo3DVec = interp1(SchnitzPatchedTime,FullOrigFluo3DVector,TbinnedTimeVector);
                        for fr_idx = 1:length(TbinnedTimeVector)
                            if ~isnan(InterpAPvec(fr_idx))
                                scAPbin = find((InterpAPvec(fr_idx) < APbins+APResolution/2) & (InterpAPvec(fr_idx) >= APbins-APResolution/2));
                                SchnitzTbinnedTotalNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                                
                                if ~isnan(InterpFluoVec(fr_idx))
                                    SchnitzTbinnedOnNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                                elseif sum(TbinnedTimeVector > TbinnedTimeVector(fr_idx) & ~isnan(InterpFluoVec)) > 0
                                    SchnitzTbinnedQuiescentNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                                elseif sum(TbinnedTimeVector < TbinnedTimeVector(fr_idx) & ~isnan(InterpFluoVec)) > 0
                                    SchnitzTbinnedFinishedTranscribingNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                                else
                                    SchnitzTbinnedOffNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                                end
                                
                                if ~isnan(InterpOrigFluoVec(fr_idx))
                                    SchnitzTbinnedMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = InterpOrigFluoVec(fr_idx);
                                    SchnitzTbinnedTraceCount(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                                else
                                    SchnitzTbinnedMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = 0;
                                end
                                
                                if ~isnan(InterpOrigFluo3DVec(fr_idx))
                                    SchnitzTbinned3DMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = InterpOrigFluo3DVec(fr_idx);
                                    SchnitzTbinned3DTraceCount(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                                else
                                    SchnitzTbinned3DMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = 0;
                                end
                                
                            end
                        end
                        
                    end
                else % Fill in what to do if particle is there but is not approved!
                    if length(CompiledParticles{ChN}(p_idx).Frame) < MinTimePoints | length(CompiledParticles{ChN}(p_idx).Frame)/length(CompiledParticles{ChN}(p_idx).TrueFrames) < .1
                        if (scNC >= 9) & (scNC <= 14) & ~isempty(scAPvector)
                            InterpAPvec = interp1(SchnitzTimes,scAPvector,TbinnedTimeVector);
                            for fr_idx = 1:length(TbinnedTimeVector)
                                if ~isnan(InterpAPvec(fr_idx))
                                    scAPbin = find((InterpAPvec(fr_idx) < APbins+APResolution/2) & (InterpAPvec(fr_idx) >= APbins-APResolution/2));
                                    SchnitzTbinnedTotalNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                                    SchnitzTbinnedOffNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                                    SchnitzTbinnedMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = 0;
                                    SchnitzTbinned3DMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = 0;
                                end
                            end
                        end
                    end
                end
                
                
            end
        else
            if (scNC >= 9) & (scNC <= 14) & ~isempty(scAPvector) & length(scAPvector) > 1
                InterpAPvec = interp1(SchnitzTimes,scAPvector,TbinnedTimeVector);
                for fr_idx = 1:length(TbinnedTimeVector)
                    if ~isnan(InterpAPvec(fr_idx))
                        scAPbin = find((InterpAPvec(fr_idx) < APbins+APResolution/2) & (InterpAPvec(fr_idx) >= APbins-APResolution/2));
                        SchnitzTbinnedTotalNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                        SchnitzTbinnedOffNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                        SchnitzTbinnedMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = 0;
                        SchnitzTbinned3DMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = 0;
                    end
                end
            end
            
            
        end
    end
    
    SchnitzAnaphaseAlignedTotalNuclei = zeros(numBinnedFrames, length(APbins), 6, NumSchnitz);
    SchnitzAnaphaseAlignedOffNuclei = zeros(numBinnedFrames, length(APbins), 6, NumSchnitz);
    SchnitzAnaphaseAlignedOnNuclei = zeros(numBinnedFrames, length(APbins), 6, NumSchnitz);
    SchnitzAnaphaseAlignedQuiescentNuclei = zeros(numBinnedFrames, length(APbins), 6, NumSchnitz);
    SchnitzAnaphaseAlignedFinishedTranscribingNuclei =  zeros(numBinnedFrames, length(APbins), 6, NumSchnitz);
    SchnitzAnaphaseAlignedMeanTraces = NaN(numBinnedFrames, length(APbins), 6, NumSchnitz);
    SchnitzAnaphaseAlignedTraceCount = zeros(numBinnedFrames, length(APbins), 6, NumSchnitz);
    SchnitzAnaphaseAligned3DMeanTraces = NaN(numBinnedFrames, length(APbins), 6, NumSchnitz);
    SchnitzAnaphaseAligned3DTraceCount = zeros(numBinnedFrames, length(APbins), 6, NumSchnitz);
    
    for sc_idx = 1:length(schnitzcells)%496:496
        
        if  (schnitzFlags(sc_idx) ~= 0) | (schnitzApproved(sc_idx) ~= 1) | ...
                ~schnitzFirstFrameContainedInTrace(sc_idx) | ~schnitzLastFrameContainedInTrace(sc_idx) | ...
                (schnitzFractionFramesApproved(sc_idx) < 0.9)
            continue
        end
        scAP = schnitzAPpos(sc_idx);
        scAPvector = schnitzcells(sc_idx).APpos;
        scAPbin = find((scAP < APbins+APResolution/2) & (scAP >= APbins-APResolution/2));
        scNC = schnitzCycles(sc_idx);
        SchnitzFrames = schnitzcells(sc_idx).frames.';
        if isempty(schnitzcells(sc_idx).anaphaseFrame) | isnan(schnitzcells(sc_idx).anaphaseFrame) | (schnitzcells(sc_idx).anaphaseFrame == 0)
            continue
        end
        SchnitzTimes = [FrameInfo(SchnitzFrames).Time]-FrameInfo(schnitzcells(sc_idx).anaphaseFrame).Time;
        %
        SchnitzHasParticle = ismember(sc_idx, ParticleSchnitzIndices);
        if SchnitzHasParticle
            p_idx = find(ParticleSchnitzIndices == sc_idx);
            if length(p_idx) > 1
                continue
            elseif length(p_idx) == 0
                if (scNC >= 9) & (scNC <= 14) & ~isempty(scAPvector)
                    InterpAPvec = interp1(SchnitzTimes,scAPvector,TbinnedTimeVector);
                    for fr_idx = 1:length(TbinnedTimeVector)
                        if ~isnan(InterpAPvec(fr_idx))
                            scAPbin = find((InterpAPvec(fr_idx) < APbins+APResolution/2) & (InterpAPvec(fr_idx) >= APbins-APResolution/2));
                            SchnitzAnaphaseAlignedTotalNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                            SchnitzAnaphaseAlignedOffNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                            SchnitzAnaphaseAlignedMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = 0;
                            SchnitzAnaphaseAligned3DMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = 0;
                        end
                    end
                end
            elseif length(p_idx) == 1
                if CompiledParticles{ChN}(p_idx).ManualApproved
                    if (~all(CompiledParticles{ChN}(p_idx).FlaggingInfo.PositionApproved) | ~all(CompiledParticles{ChN}(p_idx).FlaggingInfo.FluoApproved))
                        continue
                    else
                        pFluoVector = CompiledParticles{ChN}(p_idx).FlaggingInfo.MaxSpotFluoLevel;
                        pFrames = CompiledParticles{ChN}(p_idx).FlaggingInfo.TrueFrames;
                        pOrigFluoVector = CompiledParticles{ChN}(p_idx).Fluo;
                        pOrigFrames = CompiledParticles{ChN}(p_idx).Frame;
                        pOrigFluo3DVector = CompiledParticles{ChN}(p_idx).Fluo3DGauss;
                        if sum(~isnan(pFluoVector)) < MinTimePoints
                            if length(scAPvector) <= 1
                                continue
                            end
                            InterpAPvec = interp1(SchnitzTimes,scAPvector,TbinnedTimeVector);
                            for fr_idx = 1:length(TbinnedTimeVector)
                                if ~isnan(InterpAPvec(fr_idx))
                                    scAPbin = find((InterpAPvec(fr_idx) < APbins+APResolution/2) & (InterpAPvec(fr_idx) >= APbins-APResolution/2));
                                    SchnitzAnaphaseAlignedTotalNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                                    SchnitzAnaphaseAlignedOffNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                                    SchnitzAnaphaseAlignedMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = 0;
                                    SchnitzAnaphaseAligned3DMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = 0;
                                end
                            end
                        end
                        
                        FullFluoVector = NaN(1, length(min(SchnitzFrames):max(SchnitzFrames)));
                        FullOrigFluoVector = NaN(1, length(min(SchnitzFrames):max(SchnitzFrames)));
                        FullOrigFluo3DVector = NaN(1, length(min(SchnitzFrames):max(SchnitzFrames)));
                        if min(pFrames) < min(SchnitzFrames) | max(pFrames) > max(SchnitzFrames)
                            continue
                        end
                        if min(pOrigFrames) < min(SchnitzFrames) | max(pOrigFrames) > max(SchnitzFrames)
                            continue
                        end
                        FullFluoVector(find(ismember(min(SchnitzFrames):max(SchnitzFrames), pFrames))) = pFluoVector;
                        FullOrigFluoVector(find(ismember(min(SchnitzFrames):max(SchnitzFrames), pOrigFrames))) = pOrigFluoVector;
                        FullOrigFluo3DVector(find(ismember(min(SchnitzFrames):max(SchnitzFrames), pOrigFrames))) = pOrigFluo3DVector;
                        if length(scAPvector) <= 1
                            continue
                        end
                        InterpAPvec = interp1(SchnitzTimes,scAPvector,TbinnedTimeVector);
                        SchnitzPatchedTime = [FrameInfo(min(SchnitzFrames):max(SchnitzFrames)).Time]-FrameInfo(schnitzcells(sc_idx).anaphaseFrame).Time;
                        InterpFluoVec = interp1(SchnitzPatchedTime,FullFluoVector,TbinnedTimeVector);
                        InterpOrigFluoVec = interp1(SchnitzPatchedTime,FullOrigFluoVector,TbinnedTimeVector);
                        InterpOrigFluo3DVec = interp1(SchnitzPatchedTime,FullOrigFluo3DVector,TbinnedTimeVector);
                        for fr_idx = 1:length(TbinnedTimeVector)
                            if ~isnan(InterpAPvec(fr_idx))
                                scAPbin = find((InterpAPvec(fr_idx) < APbins+APResolution/2) & (InterpAPvec(fr_idx) >= APbins-APResolution/2));
                                SchnitzAnaphaseAlignedTotalNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                                
                                if ~isnan(InterpFluoVec(fr_idx))
                                    SchnitzAnaphaseAlignedOnNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                                elseif sum(TbinnedTimeVector > TbinnedTimeVector(fr_idx) & ~isnan(InterpFluoVec)) > 0
                                    SchnitzAnaphaseAlignedQuiescentNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                                elseif sum(TbinnedTimeVector < TbinnedTimeVector(fr_idx) & ~isnan(InterpFluoVec)) > 0
                                    SchnitzAnaphaseAlignedFinishedTranscribingNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                                else
                                    SchnitzAnaphaseAlignedOffNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                                end
                                
                                
                                if ~isnan(InterpOrigFluoVec(fr_idx))
                                    SchnitzAnaphaseAlignedMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx)= InterpOrigFluoVec(fr_idx);
                                    SchnitzAnaphaseAlignedTraceCount(fr_idx, scAPbin, scNC-8, sc_idx)= 1;
                                else
                                    SchnitzAnaphaseAlignedMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = 0;
                                end
                                
                                if ~isnan(InterpOrigFluo3DVec(fr_idx))
                                    SchnitzAnaphaseAligned3DMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx)= InterpOrigFluo3DVec(fr_idx);
                                    SchnitzAnaphaseAligned3DTraceCount(fr_idx, scAPbin, scNC-8, sc_idx)= 1;
                                else
                                    SchnitzAnaphaseAligned3DMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = 0;
                                end
                            end
                            
                        end
                    end
                else % Fill in what to do if particle is there but is not approved!
                    if length(CompiledParticles{ChN}(p_idx).Frame) <= 3 & length(CompiledParticles{ChN}(p_idx).Frame)/length(CompiledParticles{ChN}(p_idx).TrueFrames) < .1
                        if (scNC >= 9) & (scNC <= 14) & ~isempty(scAPvector)
                            InterpAPvec = interp1(SchnitzTimes,scAPvector,TbinnedTimeVector);
                            for fr_idx = 1:length(TbinnedTimeVector)
                                if ~isnan(InterpAPvec(fr_idx))
                                    scAPbin = find((InterpAPvec(fr_idx) < APbins+APResolution/2) & (InterpAPvec(fr_idx) >= APbins-APResolution/2));
                                    SchnitzAnaphaseAlignedTotalNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                                    SchnitzAnaphaseAlignedOffNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                                    SchnitzAnaphaseAlignedMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = 0;
                                    SchnitzAnaphaseAligned3DMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = 0;
                                end
                            end
                        end
                    end
                end
                
                
            end
        else
            if (scNC >= 9) & (scNC <= 14) & ~isempty(scAPvector) & length(scAPvector) > 1
                InterpAPvec = interp1(SchnitzTimes,scAPvector,TbinnedTimeVector);
                for fr_idx = 1:length(TbinnedTimeVector)
                    if ~isnan(InterpAPvec(fr_idx))
                        scAPbin = find((InterpAPvec(fr_idx) < APbins+APResolution/2) & (InterpAPvec(fr_idx) >= APbins-APResolution/2));
                        SchnitzAnaphaseAlignedTotalNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                        SchnitzAnaphaseAlignedOffNuclei(fr_idx, scAPbin, scNC-8, sc_idx) = 1;
                        SchnitzAnaphaseAlignedMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = 0;
                        SchnitzAnaphaseAligned3DMeanTraces(fr_idx, scAPbin, scNC-8, sc_idx) = 0;
                    end
                end
            end
            
            
        end
    end
    
    
    
    
    
    
    
    
    UnalignedCycleMeanTraces = NaN(numFrames, length(APbins), 6);
    UnalignedCycleTraceStdErrors =  NaN(numFrames, length(APbins), 6);
    UnalignedCycleTraceCount=  NaN(numFrames, length(APbins), 6);
    UnalignedCycleNumNuclei =  sum(SchnitzUnalignedTotalNuclei, 4);
    UnalignedCycleNumOnNuclei =  sum(SchnitzUnalignedOnNuclei, 4);
    UnalignedCycleFractionOn=  UnalignedCycleNumOnNuclei./UnalignedCycleNumNuclei;
    UnalignedCycleNumOffNuclei =  sum(SchnitzUnalignedOffNuclei, 4);
    UnalignedCycleNumQuiescentNuclei =  sum(SchnitzUnalignedQuiescentNuclei, 4);
    UnalignedCycleNumFinishedTranscribingNuclei =  sum(SchnitzUnalignedFinishedTranscribingNuclei, 4);
    Unaligned3DCycleMeanTraces = NaN(numFrames, length(APbins), 6);
    Unaligned3DCycleTraceStdErrors =  NaN(numFrames, length(APbins), 6);
    Unaligned3DCycleTraceCount=  NaN(numFrames, length(APbins), 6);
    Unaligned3DCycleNumNuclei =  sum(SchnitzUnalignedOnNuclei, 4);
    Unaligned3DCycleNumOnNuclei =  sum(SchnitzUnalignedOnNuclei, 4);
    Unaligned3DCycleNumOffNuclei =  sum(SchnitzUnalignedOffNuclei, 4);
    Unaligned3DCycleNumQuiescentNuclei =  sum(SchnitzUnalignedQuiescentNuclei, 4);
    Unaligned3DCycleNumFinishedTranscribingNuclei =  sum(SchnitzUnalignedFinishedTranscribingNuclei, 4);
    Unaligned3DCycleFractionOn=  Unaligned3DCycleNumOnNuclei./Unaligned3DCycleNumNuclei;
    CycleTraces = cell(1, 6);
    CycleTraces3D = cell(1, 6);
    MeanAPs = cell(1, 6);
    ParticleAPbinIndices = cell(1, 6);
    UnalignedCycleFrameTimes = cell(1, 6);
    Unaligned3DCycleFrameTimes =cell(1, 6);
    
    numBinnedFrames = ceil((max(FrameTimes)-FrameTimes(nc_info(6)))/deltaTbinWidth)+1;
    AnaphaseAlignedCycleMeanTraces = NaN(numBinnedFrames, length(APbins), 6);
    AnaphaseAlignedCycleTraceStdErrors =  NaN(numBinnedFrames, length(APbins), 6);
    AnaphaseAlignedCycleTraceCount=  NaN(numBinnedFrames, length(APbins), 6);
    AnaphaseAlignedCycleNumNuclei =  sum(SchnitzAnaphaseAlignedTotalNuclei, 4);
    AnaphaseAlignedCycleNumOnNuclei = sum(SchnitzAnaphaseAlignedOnNuclei, 4);
    AnaphaseAlignedCycleFractionOn=  AnaphaseAlignedCycleNumOnNuclei./AnaphaseAlignedCycleNumNuclei;
    AnaphaseAlignedCycleNumOffNuclei =  sum(SchnitzAnaphaseAlignedOffNuclei, 4);
    AnaphaseAlignedCycleNumQuiescentNuclei =  sum(SchnitzAnaphaseAlignedQuiescentNuclei, 4);
    AnaphaseAlignedCycleNumFinishedTranscribingNuclei =  sum(SchnitzAnaphaseAlignedFinishedTranscribingNuclei, 4);
    AnaphaseAligned3DCycleMeanTraces = NaN(numBinnedFrames, length(APbins), 6);
    AnaphaseAligned3DCycleTraceStdErrors =  NaN(numBinnedFrames, length(APbins), 6);
    AnaphaseAligned3DCycleTraceCount=  NaN(numBinnedFrames, length(APbins), 6);
    AnaphaseAligned3DCycleNumNuclei =  sum(SchnitzAnaphaseAlignedTotalNuclei, 4);
    AnaphaseAligned3DCycleNumOnNuclei = sum(SchnitzAnaphaseAlignedOnNuclei, 4);
    AnaphaseAligned3DCycleFractionOn=  AnaphaseAligned3DCycleNumOnNuclei./AnaphaseAligned3DCycleNumNuclei;
    AnaphaseAligned3DCycleNumOffNuclei =  sum(SchnitzAnaphaseAlignedOffNuclei, 4);
    AnaphaseAligned3DCycleNumQuiescentNuclei =  sum(SchnitzAnaphaseAlignedQuiescentNuclei, 4);
    AnaphaseAligned3DCycleNumFinishedTranscribingNuclei =  sum(SchnitzAnaphaseAlignedFinishedTranscribingNuclei, 4);
    CycleTracesWithAnaphaseAlignment = cell(1, 6);
    Cycle3DTracesWithAnaphaseAlignment = cell(1, 6);
    AnaphaseAlignedCycleFrameTimes = cell(1, 6);
    AnaphaseAligned3DCycleFrameTimes = cell(1, 6);
    
    TbinnedCycleMeanTraces = NaN(numBinnedFrames, length(APbins), 6);
    TbinnedCycleTraceStdErrors =  NaN(numBinnedFrames, length(APbins), 6);
    TbinnedCycleTraceCount=  NaN(numBinnedFrames, length(APbins), 6);
    TbinnedCycleNumNuclei =  sum(SchnitzTbinnedTotalNuclei, 4);
    TbinnedCycleNumOnNuclei =  sum(SchnitzTbinnedOnNuclei, 4);
    TbinnedCycleFractionOn=  TbinnedCycleNumOnNuclei./TbinnedCycleNumNuclei;
    TbinnedCycleNumOffNuclei =  sum(SchnitzTbinnedOffNuclei, 4);
    TbinnedCycleNumQuiescentNuclei =  sum(SchnitzTbinnedQuiescentNuclei, 4);
    TbinnedCycleNumFinishedTranscribingNuclei =  sum(SchnitzTbinnedFinishedTranscribingNuclei, 4);
    Tbinned3DCycleMeanTraces = NaN(numBinnedFrames, length(APbins), 6);
    Tbinned3DCycleTraceStdErrors =  NaN(numBinnedFrames, length(APbins), 6);
    Tbinned3DCycleTraceCount=  NaN(numBinnedFrames, length(APbins), 6);
    Tbinned3DCycleNumNuclei =  sum(SchnitzTbinnedTotalNuclei, 4);
    Tbinned3DCycleNumOnNuclei =  sum(SchnitzTbinnedOnNuclei, 4);
    Tbinned3DCycleFractionOn=  Tbinned3DCycleNumOnNuclei./Tbinned3DCycleNumNuclei;
    Tbinned3DCycleNumOffNuclei =  sum(SchnitzTbinnedOffNuclei, 4);
    Tbinned3DCycleNumQuiescentNuclei =  sum(SchnitzTbinnedQuiescentNuclei, 4);
    Tbinned3DCycleNumFinishedTranscribingNuclei =  sum(SchnitzTbinnedFinishedTranscribingNuclei, 4);
    TbinnedCycleTraces = cell(1, 6);
    Tbinned3DCycleTraces = cell(1, 6);
    TbinnedCycleFrameTimes = cell(1, 6);
    Tbinned3DCycleFrameTimes = cell(1, 6);
    
    
    IncludedNCs = min(CPcycles):max(CPcycles);
    for k = 1:length(IncludedNCs)
        NC = IncludedNCs(k);
        if UseManualApproval
            IncludedTraceIndices = find((CPcycles == NC) & (ApprovedVector == 1));
        else
            IncludedTraceIndices = find((CPcycles == NC) & (FirstApproved == 1) &...
                (ApprovedCount == 1) & (HasMinimumTimePoints == 1) & (PercentageFramesApproved == 1));
        end
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
        NCMinTraceFrames = NaN(1, length(IncludedTraceIndices));
        for i = 1:length(IncludedTraceIndices)
            trace_index = IncludedTraceIndices(i);
            CurrentCompiledParticle = CompiledParticles{ChN}(trace_index);
            NCMinTraceFrames(i) = min(CurrentCompiledParticle.FlaggingInfo.TrueFrames(CurrentCompiledParticle.FlaggingInfo.FrameApprovedFinal == 1));
            CycleTraces{NC-8}(CurrentCompiledParticle.FlaggingInfo.TrueFrames(CurrentCompiledParticle.FlaggingInfo.FrameApprovedFinal == 1), i) = 0;
            CycleTraces{NC-8}(CurrentCompiledParticle.Frame(CurrentCompiledParticle.FrameApproved == 1), i) = CurrentCompiledParticle.Fluo(CurrentCompiledParticle.FrameApproved == 1);
            CycleTraces{NC-8}(CurrentCompiledParticle.Frame(CurrentCompiledParticle.FrameApproved == 0), i) = NaN;
            CycleTraces3D{NC-8}(CurrentCompiledParticle.FlaggingInfo.TrueFrames(CurrentCompiledParticle.FlaggingInfo.FrameApprovedFinal == 1), i) = 0;
            CycleTraces3D{NC-8}(CurrentCompiledParticle.Frame(CurrentCompiledParticle.FrameApproved == 1), i) = CurrentCompiledParticle.Fluo3DGauss(CurrentCompiledParticle.FrameApproved == 1);
            CycleTraces3D{NC-8}(CurrentCompiledParticle.Frame(CurrentCompiledParticle.FrameApproved == 0), i) = NaN;
            %CycleTraces3Derror(CurrentCompiledParticle.Frame, i) = CurrentCompiledParticle.Fluo3DGaussz;
            MeanAPs{NC-8}(i) = mean(CurrentCompiledParticle.schnitzcell.APpos);
            ParticleAPbinIndices{NC-8}(i) = find((MeanAPs{NC-8}(i) < APbins+APResolution/2) & (MeanAPs{NC-8}(i) >= APbins-APResolution/2));
            
        end
        SchnitzNCFrameMins = schnitzAnaphaseFrames(schnitzCycles == NC & abs(schnitzAnaphaseFrames - anaphaseFrames(NC-8)) <= 2 & schnitzAnaphaseFrames ~= 0);
        
        APbinsInMovie = min(ParticleAPbinIndices{NC-8}):max(ParticleAPbinIndices{NC-8});
        NCCycleTraces = CycleTraces{NC-8};
        NC3DCycleTraces = CycleTraces3D{NC-8};
        IncludedRowsInBin = find(sum(~isnan(NCCycleTraces),2).' > 0);
        if isempty(SchnitzNCFrameMins)
            SchnitzNCFrameMin = min(IncludedRowsInBin);
        else
            SchnitzNCFrameMin = min(SchnitzNCFrameMins);
        end
        if min(IncludedRowsInBin) > SchnitzNCFrameMin
            IncludedRowsInBin = [SchnitzNCFrameMin:(min(IncludedRowsInBin)-1) IncludedRowsInBin];
        end
        NCRealMinFrame = min(IncludedRowsInBin,[], 'omitnan');
        UnalignedCycleFrameTimes{NC-8}= FrameTimes(IncludedRowsInBin)-min(FrameTimes(IncludedRowsInBin));
        CycleTraces{NC-8} = NCCycleTraces(IncludedRowsInBin,:);
        NCCycleTraces = NCCycleTraces(IncludedRowsInBin,:);
        
        IncludedRowsInBin3D = find(sum(~isnan(NC3DCycleTraces),2).' > 0);
        if isempty(SchnitzNCFrameMins)
            SchnitzNCFrameMin = min(IncludedRowsInBin3D);
        else
            SchnitzNCFrameMin = min(SchnitzNCFrameMins);
        end
        if min(IncludedRowsInBin3D) >SchnitzNCFrameMin
            IncludedRowsInBin3D = [SchnitzNCFrameMin:(min(IncludedRowsInBin3D)-1) IncludedRowsInBin3D];
        end
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
            UnalignedCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8) = 0;
            Unaligned3DCycleMeanTraces(IncludedRowsInBin3D, APbinIndex, NC-8) = 0;
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
            UnalignedCycleTraceCount(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(NCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin)), 2);
            if ~isempty(IncludedRowsInBin3D)
                Unaligned3DCycleMeanTraces(IncludedRowsInBin3D, APbinIndex, NC-8) = nanmean(NC3DCycleTraces(IncludedRowsInBin3D, IncludedColumnsInBin), 2);
                
                if size(NC3DCycleTraces(IncludedRowsInBin3D, IncludedColumnsInBin), 2) > 1
                    for r = IncludedRowsInBin3D
                        Unaligned3DCycleTraceStdErrors(r, APbinIndex, NC-8) = std(NC3DCycleTraces(r, IncludedColumnsInBin),'omitnan');
                    end
                    
                else
                    Unaligned3DCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8)  = NaN;
                end
                Unaligned3DCycleTraceCount(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(NC3DCycleTraces(IncludedRowsInBin, IncludedColumnsInBin)), 2);
            else
                Unaligned3DCycleMeanTraces(IncludedRowsInBin3D, APbinIndex, NC-8) =NaN;
                Unaligned3DCycleTraceStdErrors(IncludedRowsInBin3D, APbinIndex, NC-8) = NaN;
                Unaligned3DCycleTraceCount(IncludedRowsInBin, APbinIndex, NC-8) = NaN;
            end
            
        end
        
        
        
        
        % Do anaphase alignment and bin
        CurrentNCTraces = CycleTraces{NC-8};
        IncludedRowsInBin = find(sum(~isnan(CurrentNCTraces),2).' > 0);
    
        NCFrameTimes = UnalignedCycleFrameTimes{NC-8};
        NumBins = ceil((max(NCFrameTimes)-min(NCFrameTimes))/deltaTbinWidth)+1;
        AnaphaseAlignedCycleFrameTimes{NC-8} = (0:NumBins-1)*deltaTbinWidth;
        if length(AnaphaseAlignedCycleFrameTimes{NC-8}) > size(AnaphaseAlignedCycleMeanTraces, 1)
            AnaphaseAlignedCycleMeanTraces = [AnaphaseAlignedCycleMeanTraces; NaN(length(AnaphaseAlignedCycleFrameTimes{NC-8})-size(AnaphaseAlignedCycleMeanTraces, 1),length(APbins),6)];
            AnaphaseAlignedCycleTraceStdErrors = [AnaphaseAlignedCycleTraceStdErrors; NaN(length(AnaphaseAlignedCycleFrameTimes{NC-8})-size(AnaphaseAlignedCycleTraceStdErrors, 1),length(APbins),6)];
            AnaphaseAlignedCycleTraceCount = [AnaphaseAlignedCycleTraceCount; NaN(length(AnaphaseAlignedCycleFrameTimes{NC-8})-size(AnaphaseAlignedCycleTraceCount, 1),length(APbins),6)];
        end
        CycleTracesWithAnaphaseAlignment{NC-8} = NaN(NumBins, size(CurrentNCTraces, 2));
        for trace_index = 1:size(CurrentNCTraces, 2)
            CurrentTrace = CurrentNCTraces(:,trace_index);
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
        NCFrameTimes3D =  Unaligned3DCycleFrameTimes{NC-8};
        
        NumBins3D = ceil((max(NCFrameTimes3D)-min(NCFrameTimes3D))/deltaTbinWidth)+1;
        AnaphaseAligned3DCycleFrameTimes{NC-8} = (0:NumBins3D-1)*deltaTbinWidth;
        if length(AnaphaseAligned3DCycleFrameTimes{NC-8}) > size(AnaphaseAligned3DCycleMeanTraces, 1)
            AnaphaseAligned3DCycleMeanTraces = [AnaphaseAligned3DCycleMeanTraces; NaN(length(AnaphaseAligned3DCycleFrameTimes{NC-8})-size(AnaphaseAligned3DCycleMeanTraces, 1),length(APbins),6)];
            AnaphaseAligned3DCycleTraceStdErrors = [AnaphaseAligned3DCycleTraceStdErrors; NaN(length(AnaphaseAligned3DCycleFrameTimes{NC-8})-size(AnaphaseAligned3DCycleTraceStdErrors, 1),length(APbins),6)];
            AnaphaseAligned3DCycleTraceCount = [AnaphaseAligned3DCycleTraceCount; NaN(length(AnaphaseAligned3DCycleFrameTimes{NC-8})-size(AnaphaseAligned3DCycleTraceCount, 1),length(APbins),6)];
        end
        Cycle3DTracesWithAnaphaseAlignment{NC-8} = NaN(NumBins3D, size(CurrentNCTraces3D, 2));
        for trace_index = 1:size(CurrentNCTraces3D, 2)
            CurrentTrace = CurrentNCTraces3D(:,trace_index);
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
            AnaphaseAlignedCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8) = 0;
            AnaphaseAligned3DCycleMeanTraces(IncludedRowsInBin3D, APbinIndex, NC-8) = 0;
            AnaphaseAligned3DCycleTraceStdErrors(IncludedRowsInBin3D, APbinIndex, NC-8) = 0;
            
            IncludedRowsInBin = find(sum(~isnan(AnaphaseAlignedNCCycleTraces(:,IncludedColumnsInBin)),2).' > 0);
            IncludedRowsInBin3D = find(sum(~isnan(AnaphaseAlignedNC3DCycleTraces(:,IncludedColumnsInBin)),2).' > 0);
            AnaphaseAlignedCycleMeanTraces(IncludedRowsInBin, APbinIndex, NC-8) = nanmean(AnaphaseAlignedNCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin), 2);
            AnaphaseAlignedCycleTraceCount(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(AnaphaseAlignedNCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin)), 2);
            if size(AnaphaseAlignedNCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin), 2) > 1
                for r = IncludedRowsInBin3D
                    AnaphaseAlignedCycleTraceStdErrors(r, APbinIndex, NC-8) = std(AnaphaseAlignedNCCycleTraces(r, IncludedColumnsInBin),'omitnan');
                end
            else
                AnaphaseAlignedCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8)  = NaN;
            end
            
            
            if ~isempty(IncludedRowsInBin3D)
                AnaphaseAligned3DCycleMeanTraces(IncludedRowsInBin3D, APbinIndex, NC-8) = nanmean(AnaphaseAlignedNC3DCycleTraces(IncludedRowsInBin3D, IncludedColumnsInBin), 2);
                AnaphaseAligned3DCycleTraceCount(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(AnaphaseAlignedNC3DCycleTraces(IncludedRowsInBin, IncludedColumnsInBin)), 2);
                if size(AnaphaseAlignedNC3DCycleTraces(IncludedRowsInBin3D, IncludedColumnsInBin), 2) > 1
                    for r = IncludedRowsInBin3D
                        AnaphaseAligned3DCycleTraceStdErrors(r, APbinIndex, NC-8) = std(AnaphaseAlignedNC3DCycleTraces(r, IncludedColumnsInBin),'omitnan');
                    end
                else
                    AnaphaseAligned3DCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8)  = NaN;
                end
            else
                AnaphaseAligned3DCycleMeanTraces(IncludedRowsInBin3D, APbinIndex, NC-8) = NaN;
                AnaphaseAligned3DCycleTraceCount(IncludedRowsInBin, APbinIndex, NC-8) = NaN;
                AnaphaseAligned3DCycleTraceStdErrors(IncludedRowsInBin3D, APbinIndex, NC-8) =NaN;
            end
            
            
        end
        
        % Bin without anaphase alignment
        CurrentNCTraces = CycleTraces{NC-8};
        IncludedRowsInBin = find(sum(~isnan(CurrentNCTraces),2).' > 0);
        NCFrameTimes = UnalignedCycleFrameTimes{NC-8};
        NumBins = ceil((max(NCFrameTimes)-min(NCFrameTimes))/deltaTbinWidth)+1;
        TbinnedCycleFrameTimes{NC-8} = (0:NumBins-1)*deltaTbinWidth;
        if length(TbinnedCycleFrameTimes{NC-8}) > size(TbinnedCycleMeanTraces, 1)
            TbinnedCycleMeanTraces = [TbinnedCycleMeanTraces; NaN(length(TbinnedCycleFrameTimes{NC-8})-size(TbinnedCycleMeanTraces, 1),length(APbins),6)];
            TbinnedCycleTraceStdErrors = [TbinnedCycleTraceStdErrors; NaN(length(TbinnedCycleFrameTimes{NC-8})-size(TbinnedCycleTraceStdErrors, 1),length(APbins),6)];
            TbinnedCycleTraceCount = [TbinnedCycleTraceCount; NaN(length(TbinnedCycleFrameTimes{NC-8})-size(TbinnedCycleTraceCount, 1),length(APbins),6)];
        end
        TbinnedCycleTraces{NC-8} = NaN(NumBins, size(CurrentNCTraces, 2));
        for trace_index = 1:size(CurrentNCTraces, 2)
            CurrentTrace = CurrentNCTraces(:,trace_index);
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
        NCFrameTimes3D = Unaligned3DCycleFrameTimes{NC-8};
        
        NumBins3D = ceil((max(NCFrameTimes3D)-min(NCFrameTimes3D))/deltaTbinWidth)+1;
        Tbinned3DCycleFrameTimes{NC-8} = (0:NumBins3D-1)*deltaTbinWidth;
        if length(Tbinned3DCycleFrameTimes{NC-8}) > size(Tbinned3DCycleMeanTraces, 1)
            Tbinned3DCycleMeanTraces = [Tbinned3DCycleMeanTraces; NaN(length(Tbinned3DCycleFrameTimes{NC-8}) -size(Tbinned3DCycleMeanTraces, 1),length(APbins),6)];
            Tbinned3DCycleTraceStdErrors = [Tbinned3DCycleTraceStdErrors; NaN(length(Tbinned3DCycleFrameTimes{NC-8}) -size(Tbinned3DCycleTraceStdErrors, 1),length(APbins),6)];
            Tbinned3DCycleTraceCount = [Tbinned3DCycleTraceCount; NaN(length(Tbinned3DCycleFrameTimes{NC-8}) -size(Tbinned3DCycleTraceCount, 1),length(APbins),6)];
        end
        Tbinned3DCycleTraces{NC-8} = NaN(NumBins3D, size(CurrentNCTraces3D, 2));
        for trace_index = 1:size(CurrentNCTraces3D, 2)
            CurrentTrace = CurrentNCTraces3D(:,trace_index);
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
            TbinnedCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8) = 0;
            Tbinned3DCycleMeanTraces(IncludedRowsInBin3D, APbinIndex, NC-8) = 0;
            Tbinned3DCycleTraceStdErrors(IncludedRowsInBin3D, APbinIndex, NC-8) = 0;
            
            IncludedRowsInBin = find(sum(~isnan(TbinnedNCCycleTraces(:,IncludedColumnsInBin)),2).' > 0);
            IncludedRowsInBin3D = find(sum(~isnan(TbinnedNC3DCycleTraces(:,IncludedColumnsInBin)),2).' > 0);
            TbinnedCycleMeanTraces(IncludedRowsInBin, APbinIndex, NC-8) = nanmean(TbinnedNCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin), 2);
            
            TbinnedCycleTraceCount(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(TbinnedNCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin)), 2);
            if size(TbinnedNCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin), 2) > 1
                for r = IncludedRowsInBin3D
                    TbinnedCycleTraceStdErrors(r, APbinIndex, NC-8) = std(TbinnedNCCycleTraces(r, IncludedColumnsInBin),'omitnan');
                end
            else
                TbinnedCycleTraceStdErrors(IncludedRowsInBin, APbinIndex, NC-8)  = NaN;
            end
            
            
            if ~isempty(IncludedRowsInBin3D)
                Tbinned3DCycleMeanTraces(IncludedRowsInBin3D, APbinIndex, NC-8) = nanmean(TbinnedNC3DCycleTraces(IncludedRowsInBin3D, IncludedColumnsInBin), 2);
                
                Tbinned3DCycleTraceCount(IncludedRowsInBin, APbinIndex, NC-8) = sum(~isnan(TbinnedNC3DCycleTraces(IncludedRowsInBin, IncludedColumnsInBin)), 2);
                
                if size(TbinnedNC3DCycleTraces(IncludedRowsInBin3D, IncludedColumnsInBin), 2) > 1
                    for r = IncludedRowsInBin3D
                        Tbinned3DCycleTraceStdErrors(r, APbinIndex, NC-8) = std(TbinnedNC3DCycleTraces(r, IncludedColumnsInBin),'omitnan');
                    end
                else
                    Tbinned3DCycleTraceStdErrors(IncludedRowsInBin3D, APbinIndex, NC-8)  = NaN;
                end
            else
                Tbinned3DCycleMeanTraces(IncludedRowsInBin3D, APbinIndex, NC-8) = NaN;
                
                Tbinned3DCycleTraceCount(IncludedRowsInBin, APbinIndex, NC-8) = NaN;
                Tbinned3DCycleTraceStdErrors(IncludedRowsInBin3D, APbinIndex, NC-8) =NaN;
            end
        end
        
        % Check to make sure calculations make sense
        TestAllTracesMat = squeeze(UnalignedCycleMeanTraces(:,:,NC-8));
        TestAll3DTracesMat = squeeze(Unaligned3DCycleMeanTraces(:,:,NC-8));
        
        TestAnaphaseAlignedTracesMat = squeeze(AnaphaseAlignedCycleMeanTraces(:,:,NC-8));
        TestAnaphaseAligned3DTracesMat = squeeze(AnaphaseAligned3DCycleMeanTraces(:,:,NC-8));
        TestTbinnedTracesMat = squeeze(TbinnedCycleMeanTraces(:,:,NC-8));
        TestTbinned3DTracesMat = squeeze(Tbinned3DCycleMeanTraces(:,:,NC-8));
        
        
    end
    
    TestLengths1 = [size(UnalignedCycleMeanTraces, 1),  size(UnalignedCycleTraceCount, 1), size(UnalignedCycleTraceStdErrors, 1)];
    if ~all(TestLengths1 == TestLengths1(1))
        disp(['Inconsistence size for fluo calculations.'])
    end
    
    TestLengths2 = [size(Unaligned3DCycleMeanTraces, 1),size(Unaligned3DCycleTraceCount, 1), size(Unaligned3DCycleTraceStdErrors, 1)];
    if ~all(TestLengths2 == TestLengths2(1))
        disp(['Inconsistence size for 3D fluo calculations.'])
    end
    
    TestLengths3 = [size(AnaphaseAlignedCycleMeanTraces, 1), size(AnaphaseAlignedCycleTraceCount, 1), size(AnaphaseAlignedCycleTraceStdErrors, 1)];
    if ~all(TestLengths3 == TestLengths3(1))
        disp(['Inconsistence size for anaphase aligned calculations.'])
    end
    
    TestLengths4 = [size(AnaphaseAligned3DCycleMeanTraces, 1),size(AnaphaseAligned3DCycleTraceCount, 1), size(AnaphaseAligned3DCycleTraceStdErrors, 1)];
    if ~all(TestLengths4 == TestLengths4(1))
        disp(['Inconsistence size for 3D anaphase aligned calculations.'])
    end
    
    %%
    MeanProfiles = {};
    MeanProfiles.UnalignedCycleMeanTraces = UnalignedCycleMeanTraces;
    MeanProfiles.UnalignedCycleNumNuclei = UnalignedCycleNumNuclei;
    MeanProfiles.UnalignedCycleNumOnNuclei = UnalignedCycleNumOnNuclei;
    MeanProfiles.UnalignedCycleFractionOn = UnalignedCycleFractionOn;
    MeanProfiles.UnalignedCycleNumOffNuclei = UnalignedCycleNumOffNuclei;
    MeanProfiles.UnalignedCycleNumQuiescentNuclei = UnalignedCycleNumQuiescentNuclei;
    MeanProfiles.UnalignedCycleNumFinishedTranscribingNuclei = UnalignedCycleNumFinishedTranscribingNuclei;
    MeanProfiles.UnalignedCycleTraceStdErrors = UnalignedCycleTraceStdErrors;
    MeanProfiles.UnalignedCycleTraceCount = UnalignedCycleTraceCount;
    UnalignedCycleMeanPerNucleusTraces = UnalignedCycleMeanTraces.*UnalignedCycleFractionOn;
    MeanProfiles.UnalignedCycleMeanPerNucleusTraces = UnalignedCycleMeanPerNucleusTraces;
    
    MeanProfiles.Unaligned3DCycleMeanTraces = Unaligned3DCycleMeanTraces;
    MeanProfiles.Unaligned3DCycleNumNuclei = Unaligned3DCycleNumNuclei;
    MeanProfiles.Unaligned3DCycleNumOnNuclei = Unaligned3DCycleNumOnNuclei;
    MeanProfiles.Unaligned3DCycleFractionOn = Unaligned3DCycleFractionOn;
    MeanProfiles.Unaligned3DCycleNumOffNuclei = Unaligned3DCycleNumOffNuclei;
    MeanProfiles.Unaligned3DCycleNumQuiescentNuclei = Unaligned3DCycleNumQuiescentNuclei;
    MeanProfiles.Unaligned3DCycleNumFinishedTranscribingNuclei = Unaligned3DCycleNumFinishedTranscribingNuclei;
    MeanProfiles.Unaligned3DCycleTraceStdErrors = Unaligned3DCycleTraceStdErrors;
    MeanProfiles.Unaligned3DCycleTraceCount = Unaligned3DCycleTraceCount;
    Unaligned3DCycleMeanPerNucleusTraces = Unaligned3DCycleMeanTraces.*Unaligned3DCycleFractionOn;
    MeanProfiles.Unaligned3DCycleMeanPerNucleusTraces = Unaligned3DCycleMeanPerNucleusTraces;
    
    MeanProfiles.AnaphaseAlignedCycleMeanTraces = AnaphaseAlignedCycleMeanTraces;
    MeanProfiles.AnaphaseAlignedCycleNumNuclei = AnaphaseAlignedCycleNumNuclei;
    MeanProfiles.AnaphaseAlignedCycleNumOnNuclei = AnaphaseAlignedCycleNumOnNuclei;
    MeanProfiles.AnaphaseAlignedCycleFractionOn = AnaphaseAlignedCycleFractionOn;
    MeanProfiles.AnaphaseAlignedCycleNumOffNuclei = AnaphaseAlignedCycleNumOffNuclei;
    MeanProfiles.AnaphaseAlignedCycleNumQuiescentNuclei = AnaphaseAlignedCycleNumQuiescentNuclei;
    MeanProfiles.AnaphaseAlignedCycleNumFinishedTranscribingNuclei = AnaphaseAlignedCycleNumFinishedTranscribingNuclei;
    MeanProfiles.AnaphaseAlignedCycleTraceStdErrors = AnaphaseAlignedCycleTraceStdErrors;
    MeanProfiles.AnaphaseAlignedCycleTraceCount = AnaphaseAlignedCycleTraceCount;
    if size(AnaphaseAlignedCycleMeanTraces, 1) == size(AnaphaseAlignedCycleFractionOn, 1)
        AnaphaseAlignedCycleMeanPerNucleusTraces = AnaphaseAlignedCycleMeanTraces.*AnaphaseAlignedCycleFractionOn;
    end
    while size(AnaphaseAlignedCycleMeanTraces, 1) ~= size(AnaphaseAlignedCycleFractionOn, 1)
    if size(AnaphaseAlignedCycleMeanTraces, 1) == size(AnaphaseAlignedCycleFractionOn, 1)
        AnaphaseAlignedCycleMeanPerNucleusTraces = AnaphaseAlignedCycleMeanTraces.*AnaphaseAlignedCycleFractionOn;
    elseif size(AnaphaseAlignedCycleMeanTraces, 1)  > size(AnaphaseAlignedCycleFractionOn, 1)
        FirstMeanTraces = find(AnaphaseAlignedCycleMeanTraces(:,13,6).' > 0, 1);
        FirstFractionOns = find(AnaphaseAlignedCycleFractionOn(:,13,6).' > 0, 1);
        if FirstMeanTraces == FirstFractionOns
            AnaphaseAlignedCycleFractionOn = [AnaphaseAlignedCycleFractionOn ; AnaphaseAlignedCycleFractionOn(end,:,:)];
            MeanProfiles.AnaphaseAlignedCycleFractionOn = AnaphaseAlignedCycleFractionOn;
            
            AnaphaseAlignedCycleNumNuclei = [AnaphaseAlignedCycleNumNuclei ; AnaphaseAlignedCycleNumNuclei(end,:,:)];
            MeanProfiles.AnaphaseAlignedCycleNumNuclei = AnaphaseAlignedCycleNumNuclei;
            
            AnaphaseAlignedCycleNumOnNuclei = [AnaphaseAlignedCycleNumOnNuclei ; AnaphaseAlignedCycleNumOnNuclei(end,:,:)];
            MeanProfiles.AnaphaseAlignedCycleNumOnNuclei = AnaphaseAlignedCycleNumOnNuclei;
            
            AnaphaseAlignedCycleNumOffNuclei = [AnaphaseAlignedCycleNumOffNuclei ; AnaphaseAlignedCycleNumOffNuclei(end,:,:)];
            MeanProfiles.AnaphaseAlignedCycleNumOffNuclei = AnaphaseAlignedCycleNumOffNuclei;
            
            AnaphaseAlignedCycleNumQuiescentNuclei = [AnaphaseAlignedCycleNumQuiescentNuclei ; AnaphaseAlignedCycleNumQuiescentNuclei(end,:,:)];
            MeanProfiles.AnaphaseAlignedCycleNumQuiescentNuclei = AnaphaseAlignedCycleNumQuiescentNuclei;
            
            AnaphaseAlignedCycleNumFinishedTranscribingNuclei = [AnaphaseAlignedCycleNumFinishedTranscribingNuclei ; AnaphaseAlignedCycleNumFinishedTranscribingNuclei(end,:,:)];
            MeanProfiles.AnaphaseAlignedCycleNumFinishedTranscribingNuclei = AnaphaseAlignedCycleNumFinishedTranscribingNuclei;
        else
            AnaphaseAlignedCycleFractionOn = [AnaphaseAlignedCycleFractionOn(1,:,:) ; AnaphaseAlignedCycleFractionOn];
            MeanProfiles.AnaphaseAlignedCycleFractionOn = AnaphaseAlignedCycleFractionOn;
            
            AnaphaseAlignedCycleNumNuclei = [AnaphaseAlignedCycleNumNuclei(1,:,:) ; AnaphaseAlignedCycleNumNuclei];
            MeanProfiles.AnaphaseAlignedCycleNumNuclei = AnaphaseAlignedCycleNumNuclei;
            
            AnaphaseAlignedCycleNumOnNuclei = [AnaphaseAlignedCycleNumOnNuclei(1,:,:) ; AnaphaseAlignedCycleNumOnNuclei];
            MeanProfiles.AnaphaseAlignedCycleNumOnNuclei = AnaphaseAlignedCycleNumOnNuclei;
            
            AnaphaseAlignedCycleNumOffNuclei = [AnaphaseAlignedCycleNumOffNuclei(1,:,:) ; AnaphaseAlignedCycleNumOffNuclei];
            MeanProfiles.AnaphaseAlignedCycleNumOffNuclei = AnaphaseAlignedCycleNumOffNuclei;
            
            AnaphaseAlignedCycleNumQuiescentNuclei = [AnaphaseAlignedCycleNumQuiescentNuclei(1,:,:) ; AnaphaseAlignedCycleNumQuiescentNuclei];
            MeanProfiles.AnaphaseAlignedCycleNumQuiescentNuclei = AnaphaseAlignedCycleNumQuiescentNuclei;
            
            AnaphaseAlignedCycleNumFinishedTranscribingNuclei = [AnaphaseAlignedCycleNumFinishedTranscribingNuclei(1,:,:) ; AnaphaseAlignedCycleNumFinishedTranscribingNuclei];
            MeanProfiles.AnaphaseAlignedCycleNumFinishedTranscribingNuclei = AnaphaseAlignedCycleNumFinishedTranscribingNuclei;
        end
    else
        FirstMeanTraces = find(AnaphaseAlignedCycleMeanTraces(:,13,6).' > 0, 1);
        FirstFractionOns = find(AnaphaseAlignedCycleFractionOn(:,13,6).' > 0, 1);
        if FirstMeanTraces == FirstFractionOns
            AnaphaseAlignedCycleFractionOn = AnaphaseAlignedCycleFractionOn(1:end-1,:,:);
            MeanProfiles.AnaphaseAlignedCycleFractionOn = AnaphaseAlignedCycleFractionOn;
            
            AnaphaseAlignedCycleNumNuclei = AnaphaseAlignedCycleNumNuclei(1:end-1,:,:);
            MeanProfiles.AnaphaseAlignedCycleNumNuclei = AnaphaseAlignedCycleNumNuclei;
            
            AnaphaseAlignedCycleNumOnNuclei = AnaphaseAlignedCycleNumOnNuclei(1:end-1,:,:);
            MeanProfiles.AnaphaseAlignedCycleNumOnNuclei = AnaphaseAlignedCycleNumOnNuclei;
            
            AnaphaseAlignedCycleNumOffNuclei = AnaphaseAlignedCycleNumOffNuclei(1:end-1,:,:);
            MeanProfiles.AnaphaseAlignedCycleNumOffNuclei = AnaphaseAlignedCycleNumOffNuclei;
            
            AnaphaseAlignedCycleNumQuiescentNuclei = AnaphaseAlignedCycleNumQuiescentNuclei(1:end-1,:,:);
            MeanProfiles.AnaphaseAlignedCycleNumQuiescentNuclei = AnaphaseAlignedCycleNumQuiescentNuclei;
            
            AnaphaseAlignedCycleNumFinishedTranscribingNuclei = AnaphaseAlignedCycleNumFinishedTranscribingNuclei(1:end-1,:,:);
            MeanProfiles.AnaphaseAlignedCycleNumFinishedTranscribingNuclei = AnaphaseAlignedCycleNumFinishedTranscribingNuclei;
        else
            AnaphaseAlignedCycleFractionOn = AnaphaseAlignedCycleFractionOn(2:end,:,:);
            MeanProfiles.AnaphaseAlignedCycleFractionOn = AnaphaseAlignedCycleFractionOn;
            
            AnaphaseAlignedCycleNumNuclei = AnaphaseAlignedCycleNumNuclei(2:end,:,:);
            MeanProfiles.AnaphaseAlignedCycleNumNuclei = AnaphaseAlignedCycleNumNuclei;
            
            AnaphaseAlignedCycleNumOnNuclei = AnaphaseAlignedCycleNumOnNuclei(2:end,:,:);
            MeanProfiles.AnaphaseAlignedCycleNumOnNuclei = AnaphaseAlignedCycleNumOnNuclei;
            
            AnaphaseAlignedCycleNumOffNuclei = AnaphaseAlignedCycleNumOffNuclei(2:end,:,:);
            MeanProfiles.AnaphaseAlignedCycleNumOffNuclei = AnaphaseAlignedCycleNumOffNuclei;
            
            AnaphaseAlignedCycleNumQuiescentNuclei = AnaphaseAlignedCycleNumQuiescentNuclei(2:end,:,:);
            MeanProfiles.AnaphaseAlignedCycleNumQuiescentNuclei = AnaphaseAlignedCycleNumQuiescentNuclei;
            
            AnaphaseAlignedCycleNumFinishedTranscribingNuclei = AnaphaseAlignedCycleNumFinishedTranscribingNuclei(2:end,:,:);
            MeanProfiles.AnaphaseAlignedCycleNumFinishedTranscribingNuclei = AnaphaseAlignedCycleNumFinishedTranscribingNuclei;
        end
    end
    end
    AnaphaseAlignedCycleMeanPerNucleusTraces = AnaphaseAlignedCycleMeanTraces.*AnaphaseAlignedCycleFractionOn;
    MeanProfiles.AnaphaseAlignedCycleMeanPerNucleusTraces = AnaphaseAlignedCycleMeanPerNucleusTraces;
    
    MeanProfiles.AnaphaseAligned3DCycleMeanTraces = AnaphaseAligned3DCycleMeanTraces;
    MeanProfiles.AnaphaseAligned3DCycleNumNuclei = AnaphaseAligned3DCycleNumNuclei;
    MeanProfiles.AnaphaseAligned3DCycleNumOnNuclei = AnaphaseAligned3DCycleNumOnNuclei;
    MeanProfiles.AnaphaseAligned3DCycleFractionOn = AnaphaseAligned3DCycleFractionOn;
    MeanProfiles.AnaphaseAligned3DCycleNumOffNuclei = AnaphaseAligned3DCycleNumOffNuclei;
    MeanProfiles.AnaphaseAligned3DCycleNumQuiescentNuclei = AnaphaseAligned3DCycleNumQuiescentNuclei;
    MeanProfiles.AnaphaseAligned3DCycleNumFinishedTranscribingNuclei = AnaphaseAligned3DCycleNumFinishedTranscribingNuclei;
    MeanProfiles.AnaphaseAligned3DCycleTraceStdErrors = AnaphaseAligned3DCycleTraceStdErrors;
    MeanProfiles.AnaphaseAligned3DCycleTraceCount = AnaphaseAligned3DCycleTraceCount;
    if size(AnaphaseAligned3DCycleMeanTraces, 1) == size(AnaphaseAligned3DCycleFractionOn, 1)
        AnaphaseAligned3DCycleMeanPerNucleusTraces = AnaphaseAligned3DCycleMeanTraces.*AnaphaseAligned3DCycleFractionOn;
    end
    while size(AnaphaseAligned3DCycleMeanTraces, 1) ~= size(AnaphaseAligned3DCycleFractionOn, 1)
    if size(AnaphaseAligned3DCycleMeanTraces, 1) == size(AnaphaseAligned3DCycleFractionOn, 1)
        AnaphaseAligned3DCycleMeanPerNucleusTraces = AnaphaseAligned3DCycleMeanTraces.*AnaphaseAligned3DCycleFractionOn;
    elseif size(AnaphaseAligned3DCycleMeanTraces, 1)  > size(AnaphaseAligned3DCycleFractionOn, 1)
        FirstMeanTraces = find(AnaphaseAligned3DCycleMeanTraces(:,13,6).' > 0, 1);
        FirstFractionOns = find(AnaphaseAligned3DCycleFractionOn(:,13,6).' > 0, 1);
        if FirstMeanTraces == FirstFractionOns
            AnaphaseAligned3DCycleFractionOn = [AnaphaseAligned3DCycleFractionOn ; AnaphaseAligned3DCycleFractionOn(end,:,:)];
            MeanProfiles.AnaphaseAligned3DCycleFractionOn = AnaphaseAligned3DCycleFractionOn;
            
            AnaphaseAligned3DCycleNumNuclei = [AnaphaseAligned3DCycleNumNuclei ; AnaphaseAligned3DCycleNumNuclei(end,:,:)];
            MeanProfiles.AnaphaseAligned3DCycleNumNuclei = AnaphaseAligned3DCycleNumNuclei;
            
            AnaphaseAligned3DCycleNumOnNuclei = [AnaphaseAligned3DCycleNumOnNuclei ; AnaphaseAligned3DCycleNumOnNuclei(end,:,:)];
            MeanProfiles.AnaphaseAligned3DCycleNumOnNuclei = AnaphaseAligned3DCycleNumOnNuclei;
            
            AnaphaseAligned3DCycleNumOffNuclei = [AnaphaseAligned3DCycleNumOffNuclei ; AnaphaseAligned3DCycleNumOffNuclei(end,:,:)];
            MeanProfiles.AnaphaseAligned3DCycleNumOffNuclei = AnaphaseAligned3DCycleNumOffNuclei;
            
            AnaphaseAligned3DCycleNumQuiescentNuclei = [AnaphaseAligned3DCycleNumQuiescentNuclei ; AnaphaseAligned3DCycleNumQuiescentNuclei(end,:,:)];
            MeanProfiles.AnaphaseAligned3DCycleNumQuiescentNuclei = AnaphaseAligned3DCycleNumQuiescentNuclei;
            
            AnaphaseAligned3DCycleNumFinishedTranscribingNuclei = [AnaphaseAligned3DCycleNumFinishedTranscribingNuclei ; AnaphaseAligned3DCycleNumFinishedTranscribingNuclei(end,:,:)];
            MeanProfiles.AnaphaseAligned3DCycleNumFinishedTranscribingNuclei = AnaphaseAligned3DCycleNumFinishedTranscribingNuclei;
        else
            AnaphaseAligned3DCycleFractionOn = [AnaphaseAligned3DCycleFractionOn(1,:,:) ; AnaphaseAligned3DCycleFractionOn];
            MeanProfiles.AnaphaseAligned3DCycleFractionOn = AnaphaseAligned3DCycleFractionOn;
            
            AnaphaseAligned3DCycleNumNuclei = [AnaphaseAligned3DCycleNumNuclei(1,:,:) ; AnaphaseAligned3DCycleNumNuclei];
            MeanProfiles.AnaphaseAligned3DCycleNumNuclei = AnaphaseAligned3DCycleNumNuclei;
            
            AnaphaseAligned3DCycleNumOnNuclei = [AnaphaseAligned3DCycleNumOnNuclei(1,:,:) ; AnaphaseAligned3DCycleNumOnNuclei];
            MeanProfiles.AnaphaseAligned3DCycleNumOnNuclei = AnaphaseAligned3DCycleNumOnNuclei;
            
            AnaphaseAligned3DCycleNumOffNuclei = [AnaphaseAligned3DCycleNumOffNuclei(1,:,:) ; AnaphaseAligned3DCycleNumOffNuclei];
            MeanProfiles.AnaphaseAligned3DCycleNumOffNuclei = AnaphaseAligned3DCycleNumOffNuclei;
            
            AnaphaseAligned3DCycleNumQuiescentNuclei = [AnaphaseAligned3DCycleNumQuiescentNuclei(1,:,:) ; AnaphaseAligned3DCycleNumQuiescentNuclei];
            MeanProfiles.AnaphaseAligned3DCycleNumQuiescentNuclei = AnaphaseAligned3DCycleNumQuiescentNuclei;
            
            AnaphaseAligned3DCycleNumFinishedTranscribingNuclei = [AnaphaseAligned3DCycleNumFinishedTranscribingNuclei(1,:,:) ; AnaphaseAligned3DCycleNumFinishedTranscribingNuclei];
            MeanProfiles.AnaphaseAligned3DCycleNumFinishedTranscribingNuclei = AnaphaseAligned3DCycleNumFinishedTranscribingNuclei;
        end
    else
        FirstMeanTraces = find(AnaphaseAligned3DCycleMeanTraces(:,13,6).' > 0, 1);
        FirstFractionOns = find(AnaphaseAligned3DCycleFractionOn(:,13,6).' > 0, 1);
        if FirstMeanTraces == FirstFractionOns
            AnaphaseAligned3DCycleFractionOn = AnaphaseAligned3DCycleFractionOn(1:end-1,:,:);
            MeanProfiles.AnaphaseAligned3DCycleFractionOn = AnaphaseAligned3DCycleFractionOn;
            
            AnaphaseAligned3DCycleNumNuclei = AnaphaseAligned3DCycleNumNuclei(1:end-1,:,:);
            MeanProfiles.AnaphaseAligned3DCycleNumNuclei = AnaphaseAligned3DCycleNumNuclei;
            
            AnaphaseAligned3DCycleNumOnNuclei = AnaphaseAligned3DCycleNumOnNuclei(1:end-1,:,:);
            MeanProfiles.AnaphaseAligned3DCycleNumOnNuclei = AnaphaseAligned3DCycleNumOnNuclei;
            
            AnaphaseAligned3DCycleNumOffNuclei = AnaphaseAligned3DCycleNumOffNuclei(1:end-1,:,:);
            MeanProfiles.AnaphaseAligned3DCycleNumOffNuclei = AnaphaseAligned3DCycleNumOffNuclei;
            
            AnaphaseAligned3DCycleNumQuiescentNuclei = AnaphaseAligned3DCycleNumQuiescentNuclei(1:end-1,:,:);
            MeanProfiles.AnaphaseAligned3DCycleNumQuiescentNuclei = AnaphaseAligned3DCycleNumQuiescentNuclei;
            
            AnaphaseAligned3DCycleNumFinishedTranscribingNuclei = AnaphaseAligned3DCycleNumFinishedTranscribingNuclei(1:end-1,:,:);
            MeanProfiles.AnaphaseAligned3DCycleNumFinishedTranscribingNuclei = AnaphaseAligned3DCycleNumFinishedTranscribingNuclei;
        else
            AnaphaseAligned3DCycleFractionOn = AnaphaseAligned3DCycleFractionOn(2:end,:,:);
            MeanProfiles.AnaphaseAligned3DCycleFractionOn = AnaphaseAligned3DCycleFractionOn;
            
            AnaphaseAligned3DCycleNumNuclei = AnaphaseAligned3DCycleNumNuclei(2:end,:,:);
            MeanProfiles.AnaphaseAligned3DCycleNumNuclei = AnaphaseAligned3DCycleNumNuclei;
            
            AnaphaseAligned3DCycleNumOnNuclei = AnaphaseAligned3DCycleNumOnNuclei(2:end,:,:);
            MeanProfiles.AnaphaseAligned3DCycleNumOnNuclei = AnaphaseAligned3DCycleNumOnNuclei;
            
            AnaphaseAligned3DCycleNumOffNuclei = AnaphaseAligned3DCycleNumOffNuclei(2:end,:,:);
            MeanProfiles.AnaphaseAligned3DCycleNumOffNuclei = AnaphaseAligned3DCycleNumOffNuclei;
            
            AnaphaseAligned3DCycleNumQuiescentNuclei = AnaphaseAligned3DCycleNumQuiescentNuclei(2:end,:,:);
            MeanProfiles.AnaphaseAligned3DCycleNumQuiescentNuclei = AnaphaseAligned3DCycleNumQuiescentNuclei;
            
            AnaphaseAligned3DCycleNumFinishedTranscribingNuclei = AnaphaseAligned3DCycleNumFinishedTranscribingNuclei(2:end,:,:);
            MeanProfiles.AnaphaseAligned3DCycleNumFinishedTranscribingNuclei = AnaphaseAligned3DCycleNumFinishedTranscribingNuclei;
        end
    end
    end
    
    AnaphaseAligned3DCycleMeanPerNucleusTraces = AnaphaseAligned3DCycleMeanTraces.*AnaphaseAligned3DCycleFractionOn;
    MeanProfiles.AnaphaseAligned3DCycleMeanPerNucleusTraces = AnaphaseAligned3DCycleMeanPerNucleusTraces;
    
    
    MeanProfiles.TbinnedCycleMeanTraces = TbinnedCycleMeanTraces;
    MeanProfiles.TbinnedCycleNumNuclei = TbinnedCycleNumNuclei;
    MeanProfiles.TbinnedCycleNumOnNuclei = TbinnedCycleNumOnNuclei;
    MeanProfiles.TbinnedCycleFractionOn = TbinnedCycleFractionOn;
    MeanProfiles.TbinnedCycleNumOffNuclei = TbinnedCycleNumOffNuclei;
    MeanProfiles.TbinnedCycleNumQuiescentNuclei = TbinnedCycleNumQuiescentNuclei;
    MeanProfiles.TbinnedCycleNumFinishedTranscribingNuclei = TbinnedCycleNumFinishedTranscribingNuclei;
    MeanProfiles.TbinnedCycleTraceStdErrors = TbinnedCycleTraceStdErrors;
    MeanProfiles.TbinnedCycleTraceCount = TbinnedCycleTraceCount;
    
    if size(TbinnedCycleMeanTraces, 1) == size(TbinnedCycleFractionOn, 1)
        TbinnedCycleMeanPerNucleusTraces = TbinnedCycleMeanTraces.*TbinnedCycleFractionOn;
    end
    while size(TbinnedCycleMeanTraces, 1) ~= size(TbinnedCycleFractionOn, 1)
    if size(TbinnedCycleMeanTraces, 1) == size(TbinnedCycleFractionOn, 1)
        TbinnedCycleMeanPerNucleusTraces = TbinnedCycleMeanTraces.*TbinnedCycleFractionOn;
    elseif size(TbinnedCycleMeanTraces, 1)  > size(TbinnedCycleFractionOn, 1)
        FirstMeanTraces = find(TbinnedCycleMeanTraces(:,13,6).' > 0, 1);
        FirstFractionOns = find(TbinnedCycleFractionOn(:,13,6).' > 0, 1);
        if FirstMeanTraces == FirstFractionOns
            TbinnedCycleFractionOn = [TbinnedCycleFractionOn ; TbinnedCycleFractionOn(end,:,:)];
            MeanProfiles.TbinnedCycleFractionOn = TbinnedCycleFractionOn;
            
            TbinnedCycleNumNuclei = [TbinnedCycleNumNuclei ; TbinnedCycleNumNuclei(end,:,:)];
            MeanProfiles.TbinnedCycleNumNuclei = TbinnedCycleNumNuclei;
            
            TbinnedCycleNumOnNuclei = [TbinnedCycleNumOnNuclei ; TbinnedCycleNumOnNuclei(end,:,:)];
            MeanProfiles.TbinnedCycleNumOnNuclei = TbinnedCycleNumOnNuclei;
            
            TbinnedCycleNumOffNuclei = [TbinnedCycleNumOffNuclei ; TbinnedCycleNumOffNuclei(end,:,:)];
            MeanProfiles.TbinnedCycleNumOffNuclei = TbinnedCycleNumOffNuclei;
            
            TbinnedCycleNumQuiescentNuclei = [TbinnedCycleNumQuiescentNuclei ; TbinnedCycleNumQuiescentNuclei(end,:,:)];
            MeanProfiles.TbinnedCycleNumQuiescentNuclei = TbinnedCycleNumQuiescentNuclei;
            
            TbinnedCycleNumFinishedTranscribingNuclei = [TbinnedCycleNumFinishedTranscribingNuclei ; TbinnedCycleNumFinishedTranscribingNuclei(end,:,:)];
            MeanProfiles.TbinnedCycleNumFinishedTranscribingNuclei = TbinnedCycleNumFinishedTranscribingNuclei;
        else
            TbinnedCycleFractionOn = [TbinnedCycleFractionOn(1,:,:) ; TbinnedCycleFractionOn];
            MeanProfiles.TbinnedCycleFractionOn = TbinnedCycleFractionOn;
            
            TbinnedCycleNumNuclei = [TbinnedCycleNumNuclei(1,:,:) ; TbinnedCycleNumNuclei];
            MeanProfiles.TbinnedCycleNumNuclei = TbinnedCycleNumNuclei;
            
            TbinnedCycleNumOnNuclei = [TbinnedCycleNumOnNuclei(1,:,:) ; TbinnedCycleNumOnNuclei];
            MeanProfiles.TbinnedCycleNumOnNuclei = TbinnedCycleNumOnNuclei;
            
            TbinnedCycleNumOffNuclei = [TbinnedCycleNumOffNuclei(1,:,:) ; TbinnedCycleNumOffNuclei];
            MeanProfiles.TbinnedCycleNumOffNuclei = TbinnedCycleNumOffNuclei;
            
            TbinnedCycleNumQuiescentNuclei = [TbinnedCycleNumQuiescentNuclei(1,:,:) ; TbinnedCycleNumQuiescentNuclei];
            MeanProfiles.TbinnedCycleNumQuiescentNuclei = TbinnedCycleNumQuiescentNuclei;
            
            TbinnedCycleNumFinishedTranscribingNuclei = [TbinnedCycleNumFinishedTranscribingNuclei(1,:,:) ; TbinnedCycleNumFinishedTranscribingNuclei];
            MeanProfiles.TbinnedCycleNumFinishedTranscribingNuclei = TbinnedCycleNumFinishedTranscribingNuclei;
        end
    else
        FirstMeanTraces = find(TbinnedCycleMeanTraces(:,13,6).' > 0, 1);
        FirstFractionOns = find(TbinnedCycleFractionOn(:,13,6).' > 0, 1);
        if FirstMeanTraces == FirstFractionOns
            TbinnedCycleFractionOn = TbinnedCycleFractionOn(1:end-1,:,:);
            MeanProfiles.TbinnedCycleFractionOn = TbinnedCycleFractionOn;
            
            TbinnedCycleNumNuclei = TbinnedCycleNumNuclei(1:end-1,:,:);
            MeanProfiles.TbinnedCycleNumNuclei = TbinnedCycleNumNuclei;
            
            TbinnedCycleNumOnNuclei = TbinnedCycleNumOnNuclei(1:end-1,:,:);
            MeanProfiles.TbinnedCycleNumOnNuclei = TbinnedCycleNumOnNuclei;
            
            TbinnedCycleNumOffNuclei = TbinnedCycleNumOffNuclei(1:end-1,:,:);
            MeanProfiles.TbinnedCycleNumOffNuclei = TbinnedCycleNumOffNuclei;
            
            TbinnedCycleNumQuiescentNuclei = TbinnedCycleNumQuiescentNuclei(1:end-1,:,:);
            MeanProfiles.TbinnedCycleNumQuiescentNuclei = TbinnedCycleNumQuiescentNuclei;
            
            TbinnedCycleNumFinishedTranscribingNuclei = TbinnedCycleNumFinishedTranscribingNuclei(1:end-1,:,:);
            MeanProfiles.TbinnedCycleNumFinishedTranscribingNuclei = TbinnedCycleNumFinishedTranscribingNuclei;
        else
            TbinnedCycleFractionOn = TbinnedCycleFractionOn(2:end,:,:);
            MeanProfiles.TbinnedCycleFractionOn = TbinnedCycleFractionOn;
            
            TbinnedCycleNumNuclei = TbinnedCycleNumNuclei(2:end,:,:);
            MeanProfiles.TbinnedCycleNumNuclei = TbinnedCycleNumNuclei;
            
            TbinnedCycleNumOnNuclei = TbinnedCycleNumOnNuclei(2:end,:,:);
            MeanProfiles.TbinnedCycleNumOnNuclei = TbinnedCycleNumOnNuclei;
            
            TbinnedCycleNumOffNuclei = TbinnedCycleNumOffNuclei(2:end,:,:);
            MeanProfiles.TbinnedCycleNumOffNuclei = TbinnedCycleNumOffNuclei;
            
            TbinnedCycleNumQuiescentNuclei = TbinnedCycleNumQuiescentNuclei(2:end,:,:);
            MeanProfiles.TbinnedCycleNumQuiescentNuclei = TbinnedCycleNumQuiescentNuclei;
            
            TbinnedCycleNumFinishedTranscribingNuclei = TbinnedCycleNumFinishedTranscribingNuclei(2:end,:,:);
            MeanProfiles.TbinnedCycleNumFinishedTranscribingNuclei = TbinnedCycleNumFinishedTranscribingNuclei;
        end
    end
    end
    TbinnedCycleMeanPerNucleusTraces = TbinnedCycleMeanTraces.*TbinnedCycleFractionOn;
    MeanProfiles.TbinnedCycleMeanPerNucleusTraces = TbinnedCycleMeanPerNucleusTraces;
    
    MeanProfiles.Tbinned3DCycleMeanTraces = Tbinned3DCycleMeanTraces;
    MeanProfiles.Tbinned3DCycleNumNuclei = Tbinned3DCycleNumNuclei;
    MeanProfiles.Tbinned3DCycleNumOnNuclei = Tbinned3DCycleNumOnNuclei;
    MeanProfiles.Tbinned3DCycleFractionOn = Tbinned3DCycleFractionOn;
    MeanProfiles.Tbinned3DCycleNumOffNuclei = Tbinned3DCycleNumOffNuclei;
    MeanProfiles.Tbinned3DCycleNumQuiescentNuclei = Tbinned3DCycleNumQuiescentNuclei;
    MeanProfiles.Tbinned3DCycleNumFinishedTranscribingNuclei = Tbinned3DCycleNumFinishedTranscribingNuclei;
    MeanProfiles.Tbinned3DCycleTraceStdErrors = Tbinned3DCycleTraceStdErrors;
    MeanProfiles.Tbinned3DCycleTraceCount = Tbinned3DCycleTraceCount;
    
    MeanProfiles.AnaphaseAlignedCycleFrameTimes = AnaphaseAlignedCycleFrameTimes;
    MeanProfiles.AnaphaseAligned3DCycleFrameTimes = AnaphaseAligned3DCycleFrameTimes;
    MeanProfiles.UnalignedCycleFrameTimes = UnalignedCycleFrameTimes;
    MeanProfiles.Unaligned3DCycleFrameTimes = Unaligned3DCycleFrameTimes;
    MeanProfiles.TbinnedCycleFrameTimes = TbinnedCycleFrameTimes;
    MeanProfiles.Tbinned3DCycleFrameTimes = Tbinned3DCycleFrameTimes;
    if size(Tbinned3DCycleMeanTraces, 1) == size(Tbinned3DCycleFractionOn, 1)
        Tbinned3DCycleMeanPerNucleusTraces = Tbinned3DCycleMeanTraces.*Tbinned3DCycleFractionOn;
    end
    while size(Tbinned3DCycleMeanTraces, 1) ~= size(Tbinned3DCycleFractionOn, 1)
    if size(Tbinned3DCycleMeanTraces, 1) == size(Tbinned3DCycleFractionOn, 1)
        Tbinned3DCycleMeanPerNucleusTraces = Tbinned3DCycleMeanTraces.*Tbinned3DCycleFractionOn;
    elseif size(Tbinned3DCycleMeanTraces, 1)  > size(Tbinned3DCycleFractionOn, 1)
        FirstMeanTraces = find(Tbinned3DCycleMeanTraces(:,13,6).' > 0, 1);
        FirstFractionOns = find(Tbinned3DCycleFractionOn(:,13,6).' > 0, 1);
        if FirstMeanTraces == FirstFractionOns
            Tbinned3DCycleFractionOn = [Tbinned3DCycleFractionOn ; Tbinned3DCycleFractionOn(end,:,:)];
            MeanProfiles.Tbinned3DCycleFractionOn = Tbinned3DCycleFractionOn;
            
            Tbinned3DCycleNumNuclei = [Tbinned3DCycleNumNuclei ; Tbinned3DCycleNumNuclei(end,:,:)];
            MeanProfiles.Tbinned3DCycleNumNuclei = Tbinned3DCycleNumNuclei;
            
            Tbinned3DCycleNumOnNuclei = [Tbinned3DCycleNumOnNuclei ; Tbinned3DCycleNumOnNuclei(end,:,:)];
            MeanProfiles.Tbinned3DCycleNumOnNuclei = Tbinned3DCycleNumOnNuclei;
            
            Tbinned3DCycleNumOffNuclei = [Tbinned3DCycleNumOffNuclei ; Tbinned3DCycleNumOffNuclei(end,:,:)];
            MeanProfiles.Tbinned3DCycleNumOffNuclei = Tbinned3DCycleNumOffNuclei;
            
            Tbinned3DCycleNumQuiescentNuclei = [Tbinned3DCycleNumQuiescentNuclei ; Tbinned3DCycleNumQuiescentNuclei(end,:,:)];
            MeanProfiles.Tbinned3DCycleNumQuiescentNuclei = Tbinned3DCycleNumQuiescentNuclei;
            
            Tbinned3DCycleNumFinishedTranscribingNuclei = [Tbinned3DCycleNumFinishedTranscribingNuclei ; Tbinned3DCycleNumFinishedTranscribingNuclei(end,:,:)];
            MeanProfiles.Tbinned3DCycleNumFinishedTranscribingNuclei = Tbinned3DCycleNumFinishedTranscribingNuclei;
        else
            Tbinned3DCycleFractionOn = [Tbinned3DCycleFractionOn(1,:,:) ; Tbinned3DCycleFractionOn];
            MeanProfiles.Tbinned3DCycleFractionOn = Tbinned3DCycleFractionOn;
            
            Tbinned3DCycleNumNuclei = [Tbinned3DCycleNumNuclei(1,:,:) ; Tbinned3DCycleNumNuclei];
            MeanProfiles.Tbinned3DCycleNumNuclei = Tbinned3DCycleNumNuclei;
            
            Tbinned3DCycleNumOnNuclei = [Tbinned3DCycleNumOnNuclei(1,:,:) ; Tbinned3DCycleNumOnNuclei];
            MeanProfiles.Tbinned3DCycleNumOnNuclei = Tbinned3DCycleNumOnNuclei;
            
            Tbinned3DCycleNumOffNuclei = [Tbinned3DCycleNumOffNuclei(1,:,:) ; Tbinned3DCycleNumOffNuclei];
            MeanProfiles.Tbinned3DCycleNumOffNuclei = Tbinned3DCycleNumOffNuclei;
            
            Tbinned3DCycleNumQuiescentNuclei = [Tbinned3DCycleNumQuiescentNuclei(1,:,:) ; Tbinned3DCycleNumQuiescentNuclei];
            MeanProfiles.Tbinned3DCycleNumQuiescentNuclei = Tbinned3DCycleNumQuiescentNuclei;
            
            Tbinned3DCycleNumFinishedTranscribingNuclei = [Tbinned3DCycleNumFinishedTranscribingNuclei(1,:,:) ; Tbinned3DCycleNumFinishedTranscribingNuclei];
            MeanProfiles.Tbinned3DCycleNumFinishedTranscribingNuclei = Tbinned3DCycleNumFinishedTranscribingNuclei;
        end
    else
        FirstMeanTraces = find(Tbinned3DCycleMeanTraces(:,13,6).' > 0, 1);
        FirstFractionOns = find(Tbinned3DCycleFractionOn(:,13,6).' > 0, 1);
        if FirstMeanTraces == FirstFractionOns
            Tbinned3DCycleFractionOn = Tbinned3DCycleFractionOn(1:end-1,:,:);
            MeanProfiles.Tbinned3DCycleFractionOn = Tbinned3DCycleFractionOn;
            
            Tbinned3DCycleNumNuclei = Tbinned3DCycleNumNuclei(1:end-1,:,:);
            MeanProfiles.Tbinned3DCycleNumNuclei = Tbinned3DCycleNumNuclei;
            
            Tbinned3DCycleNumOnNuclei = Tbinned3DCycleNumOnNuclei(1:end-1,:,:);
            MeanProfiles.Tbinned3DCycleNumOnNuclei = Tbinned3DCycleNumOnNuclei;
            
            Tbinned3DCycleNumOffNuclei = Tbinned3DCycleNumOffNuclei(1:end-1,:,:);
            MeanProfiles.Tbinned3DCycleNumOffNuclei = Tbinned3DCycleNumOffNuclei;
            
            Tbinned3DCycleNumQuiescentNuclei = Tbinned3DCycleNumQuiescentNuclei(1:end-1,:,:);
            MeanProfiles.Tbinned3DCycleNumQuiescentNuclei = Tbinned3DCycleNumQuiescentNuclei;
            
            Tbinned3DCycleNumFinishedTranscribingNuclei = Tbinned3DCycleNumFinishedTranscribingNuclei(1:end-1,:,:);
            MeanProfiles.Tbinned3DCycleNumFinishedTranscribingNuclei = Tbinned3DCycleNumFinishedTranscribingNuclei;
        else
            Tbinned3DCycleFractionOn = Tbinned3DCycleFractionOn(2:end,:,:);
            MeanProfiles.Tbinned3DCycleFractionOn = Tbinned3DCycleFractionOn;
            
            Tbinned3DCycleNumNuclei = Tbinned3DCycleNumNuclei(2:end,:,:);
            MeanProfiles.Tbinned3DCycleNumNuclei = Tbinned3DCycleNumNuclei;
            
            Tbinned3DCycleNumOnNuclei = Tbinned3DCycleNumOnNuclei(2:end,:,:);
            MeanProfiles.Tbinned3DCycleNumOnNuclei = Tbinned3DCycleNumOnNuclei;
            
            Tbinned3DCycleNumOffNuclei = Tbinned3DCycleNumOffNuclei(2:end,:,:);
            MeanProfiles.Tbinned3DCycleNumOffNuclei = Tbinned3DCycleNumOffNuclei;
            
            Tbinned3DCycleNumQuiescentNuclei = Tbinned3DCycleNumQuiescentNuclei(2:end,:,:);
            MeanProfiles.Tbinned3DCycleNumQuiescentNuclei = Tbinned3DCycleNumQuiescentNuclei;
            
            Tbinned3DCycleNumFinishedTranscribingNuclei = Tbinned3DCycleNumFinishedTranscribingNuclei(2:end,:,:);
            MeanProfiles.Tbinned3DCycleNumFinishedTranscribingNuclei = Tbinned3DCycleNumFinishedTranscribingNuclei;
        end
    end
    end
    TbinnedCycleMeanPerNucleusTraces = TbinnedCycleMeanTraces.*TbinnedCycleFractionOn;
    MeanProfiles.TbinnedCycleMeanPerNucleusTraces = TbinnedCycleMeanPerNucleusTraces;
    Tbinned3DCycleMeanPerNucleusTraces = Tbinned3DCycleMeanTraces.*Tbinned3DCycleFractionOn;
    MeanProfiles.Tbinned3DCycleMeanPerNucleusTraces = Tbinned3DCycleMeanPerNucleusTraces;
    
    savedVariables = {'UnalignedCycleMeanTraces', 'UnalignedCycleNumNuclei',  'UnalignedCycleNumOnNuclei',...
        'UnalignedCycleNumOffNuclei', 'UnalignedCycleFractionOn',  'UnalignedCycleNumQuiescentNuclei',...
        'UnalignedCycleNumFinishedTranscribingNuclei','UnalignedCycleTraceStdErrors','UnalignedCycleTraceCount',...
        'Unaligned3DCycleMeanTraces','Unaligned3DCycleNumNuclei', 'Unaligned3DCycleNumOnNuclei',...
        'Unaligned3DCycleNumOffNuclei', 'Unaligned3DCycleFractionOn',  'Unaligned3DCycleNumQuiescentNuclei',...
        'Unaligned3DCycleNumFinishedTranscribingNuclei','Unaligned3DCycleTraceStdErrors','Unaligned3DCycleTraceCount',...
        'AnaphaseAlignedCycleMeanTraces', 'AnaphaseAlignedCycleNumNuclei', 'AnaphaseAlignedCycleNumOnNuclei',...
        'AnaphaseAlignedCycleNumOffNuclei', 'AnaphaseAlignedCycleFractionOn',  'AnaphaseAlignedCycleNumQuiescentNuclei',...
        'AnaphaseAlignedCycleNumFinishedTranscribingNuclei','AnaphaseAlignedCycleTraceStdErrors','AnaphaseAlignedCycleTraceCount',...
        'AnaphaseAligned3DCycleMeanTraces','AnaphaseAligned3DCycleNumNuclei', 'AnaphaseAligned3DCycleNumOnNuclei',...
        'AnaphaseAligned3DCycleNumOffNuclei', 'AnaphaseAligned3DCycleFractionOn',  'AnaphaseAligned3DCycleNumQuiescentNuclei',...
        'AnaphaseAligned3DCycleNumFinishedTranscribingNuclei','AnaphaseAligned3DCycleTraceStdErrors','AnaphaseAligned3DCycleTraceCount',...
        'TbinnedCycleMeanTraces', 'TbinnedCycleNumNuclei','TbinnedCycleNumOnNuclei',...
        'TbinnedCycleNumOffNuclei', 'TbinnedCycleFractionOn',  'TbinnedCycleNumQuiescentNuclei',...
        'TbinnedCycleNumFinishedTranscribingNuclei','TbinnedCycleTraceStdErrors','TbinnedCycleTraceCount',...
        'Tbinned3DCycleMeanTraces','Tbinned3DCycleNumNuclei','Tbinned3DCycleNumOnNuclei',...
        'Tbinned3DCycleNumOffNuclei', 'Tbinned3DCycleFractionOn',  'Tbinned3DCycleNumQuiescentNuclei',...
        'Tbinned3DCycleNumFinishedTranscribingNuclei','Tbinned3DCycleTraceStdErrors','Tbinned3DCycleTraceCount',...
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








