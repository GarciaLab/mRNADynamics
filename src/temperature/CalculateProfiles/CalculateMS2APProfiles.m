function MeanProfiles = CalculateMS2APProfiles(Prefix, deltaTbinWidth, MinTimePoints,UseManualApproval, verbose)
if ~exist('deltaTbinWidth', 'var')
    deltaTbinWidth = 30;
elseif isempty(deltaTbinWidth)
    deltaTbinWidth = 30;
end
if ~exist('MinTimePoints', 'var')
    MinTimePoints = 5;
elseif isempty(MinTimePoints)
    MinTimePoints = 5;
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
ParticleFluos = GetParticleFluoDist(CompiledParticles);
MaxParticleFluo = 6000;%max(ParticleFluos);
MaxTraceLengths = GetParticleTraceLengthDist(CompiledParticles);
pOnCutoffs = 5:5:50;
ChN = 1;
EarliestTurnOnTimeCells = cell(1, 6);
EarliestTurnOnTimes = NaN(1, 6);
for NC = 9:14
    EarliestTurnOnTimeCells{NC-8} = [];
    for i =1:length(CompiledParticles{ChN})
        if schnitzcells(CompiledParticles{ChN}(i).schnitz).cycle == NC & ...
                min(CompiledParticles{ChN}(i).Frame) >= nc_info(NC-8)-1
            if sum(CompiledParticles{ChN}(i).FrameApproved) > 0
                EarliestTurnOnTimeCells{NC-8}(end+1) = min(CompiledParticles{ChN}(i).FlaggingInfo.AllSchnitzFrames(CompiledParticles{ChN}(i).FlaggingInfo.UseTraceFluo));
            end
        end
    end
    if ~isempty(EarliestTurnOnTimeCells{NC-8})
        EarliestTurnOnTimes(NC-8) = prctile(EarliestTurnOnTimeCells{NC-8}, 20);
        EarliestTurnOnTimes(NC-8) = max([EarliestTurnOnTimes(NC-8) nc_info(NC-8)+0.2*(nc_info(NC-7)-nc_info(NC-8))]);
    end
end
NumSchnitz = length(schnitzcells);
%%
% First make average profiles binning everything by first anaphase of the
% nuclear cycle
NChannels = length(liveExperiment.spotChannels);
for ChN=1:NChannels
    % First Do Fraction On Calculations
    
    ParticlesAreGood = ones(1, length(CompiledParticles{ChN}), 'logical');
    ParticleCycles = NaN(1, length(CompiledParticles{ChN}));
    ParticleSchnitzIndex = NaN(1, length(CompiledParticles{ChN}));
    for p_idx = 1:length(CompiledParticles{ChN})
        ParticlesAreGood(p_idx) = true;
        if ~isempty(CompiledParticles{ChN}(p_idx).Approved)
            if CompiledParticles{ChN}(p_idx).Approved ~= 1
                ParticlesAreGood(p_idx) = false;
            end
        end
        if ~CompiledParticles{ChN}(p_idx).FlaggingInfo.SpotStateCanBeCalled
            ParticlesAreGood(p_idx) = false;
        end
        if ~isempty(CompiledParticles{ChN}(p_idx).cycle)
            ParticleCycles(p_idx) = CompiledParticles{ChN}(p_idx).cycle;
        end
        if ~isempty(CompiledParticles{ChN}(p_idx).schnitz)
            ParticleSchnitzIndex(p_idx) = CompiledParticles{ChN}(p_idx).schnitz;
            CompiledParticles{ChN}(p_idx).schnitzcell = schnitzcells(ParticleSchnitzIndex(p_idx));
            if ~all(CompiledParticles{ChN}(p_idx).schnitzcell.FlaggingInfo.SchnitzAwayFromBoundary(CompiledParticles{ChN}(p_idx).schnitzcell.FlaggingInfo.AllSchnitzFrames >= EarliestTurnOnTimes(CompiledParticles{ChN}(p_idx).schnitzcell.cycle-8) ))
                if sum(CompiledParticles{ChN}(p_idx).schnitzcell.FlaggingInfo.SchnitzAwayFromBoundary(CompiledParticles{ChN}(p_idx).schnitzcell.FlaggingInfo.AllSchnitzFrames >= EarliestTurnOnTimes(CompiledParticles{ChN}(p_idx).schnitzcell.cycle-8) ))/...
                        length(CompiledParticles{ChN}(p_idx).schnitzcell.FlaggingInfo.AllSchnitzFrames(CompiledParticles{ChN}(p_idx).schnitzcell.FlaggingInfo.AllSchnitzFrames >= EarliestTurnOnTimes(CompiledParticles{ChN}(p_idx).schnitzcell.cycle-8) )) < 0.95
                    ParticlesAreGood(p_idx) = false;
                end
            end
            if ~isempty(CompiledParticles{ChN}(p_idx).schnitzcell.anaphaseFrame)
                if CompiledParticles{ChN}(p_idx).schnitzcell.anaphaseFrame == 0 | isnan(CompiledParticles{ChN}(p_idx).schnitzcell.anaphaseFrame)
                    ParticlesAreGood(p_idx)= false;
                end
            end
            if CompiledParticles{ChN}(p_idx).schnitzcell.Approved ~= 1 | CompiledParticles{ChN}(p_idx).schnitzcell.Flag ~= 0
                ParticlesAreGood(p_idx)= false;
            end
        else
            ParticlesAreGood(p_idx) = false;
        end
        
    end
    
    SchnitzesAreGood = ones(1, length(schnitzcells), 'logical');
    SchnitzCycles = NaN(1, length(schnitzcells));
    SchnitzAnaphaseFrames =  NaN(1, length(schnitzcells));
    SchnitzParticleIndex = NaN(1, length(schnitzcells));
    for sc_idx = 1:length(schnitzcells)
        SchnitzesAreGood(sc_idx) = true;
        if ismember(sc_idx, ParticleSchnitzIndex)
            SchnitzParticleIndex(sc_idx) = find(ParticleSchnitzIndex == sc_idx);
        end
        if ~isempty(schnitzcells(sc_idx).Flag)
            if schnitzcells(sc_idx).Flag ~= 0
                SchnitzesAreGood(sc_idx) = false;
            end
        end
        if ~isempty(schnitzcells(sc_idx).Approved)
            if schnitzcells(sc_idx).Approved ~= 1
                SchnitzesAreGood(sc_idx) = false;
            end
        end
        if ~isempty(schnitzcells(sc_idx).cycle)
            SchnitzCycles(sc_idx) = schnitzcells(sc_idx).cycle;
            if EarliestTurnOnTimes(schnitzcells(sc_idx).cycle-8) < min(schnitzcells(sc_idx).FlaggingInfo.AllSchnitzFrames(schnitzcells(sc_idx).FlaggingInfo.SchnitzAwayFromBoundary))
                SchnitzesAreGood(sc_idx) = false;
            end
        end
        if ~all(schnitzcells(sc_idx).FlaggingInfo.SchnitzPresent(schnitzcells(sc_idx).FlaggingInfo.AllSchnitzFrames >= EarliestTurnOnTimes(schnitzcells(sc_idx).cycle-8)  ))
            if sum(schnitzcells(sc_idx).FlaggingInfo.SchnitzPresent(schnitzcells(sc_idx).FlaggingInfo.AllSchnitzFrames >= EarliestTurnOnTimes(schnitzcells(sc_idx).cycle-8) ))/...
                    length(schnitzcells(sc_idx).FlaggingInfo.AllSchnitzFrames(schnitzcells(sc_idx).FlaggingInfo.AllSchnitzFrames >= EarliestTurnOnTimes(schnitzcells(sc_idx).cycle-8) )) < 0.95
                SchnitzesAreGood(sc_idx) = false;
            end
        end
        if ~isempty(schnitzcells(sc_idx).anaphaseFrame)
            if schnitzcells(sc_idx).anaphaseFrame == 0 | isnan(schnitzcells(sc_idx).anaphaseFrame)
                SchnitzesAreGood(sc_idx) = false;
            else
                SchnitzAnaphaseFrames(sc_idx) = schnitzcells(sc_idx).anaphaseFrame;
            end
        else
            SchnitzesAreGood(sc_idx) = false;
        end
        
    end
    
    
    
    
    
    SchnitzUnalignedTotalNuclei = zeros(numFrames, length(APbins), 6, NumSchnitz);
    SchnitzUnalignedOffNuclei = zeros(numFrames, length(APbins), 6, NumSchnitz);
    SchnitzUnalignedOnNuclei = zeros(numFrames, length(APbins), 6, NumSchnitz);
    SchnitzUnalignedQuiescentNuclei = zeros(numFrames, length(APbins), 6, NumSchnitz);
    SchnitzUnalignedFinishedTranscribingNuclei =  zeros(numFrames, length(APbins), 6, NumSchnitz);
    SchnitzUnalignedMeanTraces = NaN(numFrames, length(APbins), 6, NumSchnitz);
    SchnitzUnalignedTraceCount = zeros(numFrames, length(APbins), 6, NumSchnitz);
    SchnitzUnalignedActiveNuclei = zeros(length(APbins), 6, NumSchnitz);
    SchnitzUnalignedInactiveNuclei = zeros(length(APbins), 6, NumSchnitz);
    SchnitzUnalignedNucleiCount = zeros(length(APbins), 6, NumSchnitz);
    SchnitzUnalignedGradedOnNuclei = zeros(numFrames, length(APbins), 6,10, NumSchnitz);
    SchnitzUnalignedGradedOffNuclei = zeros(numFrames, length(APbins), 6,10, NumSchnitz);
    SchnitzUnalignedGradedTotalNuclei = zeros(numFrames, length(APbins), 6,10, NumSchnitz);
    SchnitzUnalignedGradedActiveNuclei = zeros(length(APbins), 6,10, NumSchnitz);
    SchnitzUnalignedGradedInactiveNuclei = zeros(length(APbins), 6,10, NumSchnitz);
    SchnitzUnalignedGradedNucleiCount = zeros(length(APbins), 6,10, NumSchnitz);
    
    for sc_idx =1:length(schnitzcells)
        if ~SchnitzesAreGood(sc_idx)
            continue
        end
        
        scAP = mean(schnitzcells(sc_idx).APpos);
        scAPvector = schnitzcells(sc_idx).APpos;
        Embryo_scAPbin = find((scAP < APbins+APResolution/2) & (scAP >= APbins-APResolution/2));
        scNC = SchnitzCycles(sc_idx);
        if (scNC < 9) | (scNC > 14) | isempty(scAPvector)
            continue
        end
        
        
        
        
        CandidateAnaphaseFrames = SchnitzAnaphaseFrames(SchnitzCycles == scNC & abs(SchnitzAnaphaseFrames-anaphaseFrames(scNC-8)) <= 3 & SchnitzAnaphaseFrames > 0);
        if ~isempty(CandidateAnaphaseFrames)
            SchnitzFrameMin = min([CandidateAnaphaseFrames anaphaseFrames(scNC-8)]);
        else
            SchnitzFrameMin = anaphaseFrames(scNC-8);
        end
        SchnitzHasParticle = ~isnan(SchnitzParticleIndex(sc_idx));
        if ~SchnitzHasParticle
            AllSchnitzFrames = schnitzcells(sc_idx).FlaggingInfo.AllSchnitzFrames;
            
            for frame_index = 1:length(AllSchnitzFrames)
                if schnitzcells(sc_idx).FlaggingInfo.SpotStateDefinitive(frame_index)
                    %scAPbin = find((scAPvector(frame_index) < APbins+APResolution/2) & (scAPvector(frame_index) >= APbins-APResolution/2));
                    scAPbin = Embryo_scAPbin;
                    SchnitzUnalignedTotalNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                    SchnitzUnalignedOffNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                    
                    for grade_index = 1:length(pOnCutoffs)
                        SchnitzUnalignedGradedOffNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                        SchnitzUnalignedGradedTotalNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                        
                    end
                end
            end
            SchnitzUnalignedNucleiCount(Embryo_scAPbin, scNC-8, sc_idx) = 1;
            SchnitzUnalignedInactiveNuclei(Embryo_scAPbin, scNC-8, sc_idx) = 1;
            
            for grade_index = 1:length(pOnCutoffs)
                SchnitzUnalignedGradedInactiveNuclei(scAPbin, scNC-8, grade_index, sc_idx) = 1;
                SchnitzUnalignedGradedNucleiCount(scAPbin, scNC-8, grade_index, sc_idx) = 1;
                
            end
            
            
            
            
        else
            p_idx = SchnitzParticleIndex(sc_idx);
            CurrentParticle = CompiledParticles{ChN}(p_idx);
            AllSchnitzFrames = CurrentParticle.FlaggingInfo.AllSchnitzFrames;
            if ParticlesAreGood(p_idx)
                if (sum(CurrentParticle.FlaggingInfo.SchnitzOn) >= MinTimePoints) & (scNC >= 9) & (scNC <= 14) & ~isempty(scAPvector)
                    for frame_index = 1:length(AllSchnitzFrames)
                        if CurrentParticle.FlaggingInfo.SpotStateDefinitive(frame_index)
                            %scAPbin = find((scAPvector(frame_index) < APbins+APResolution/2) & (scAPvector(frame_index) >= APbins-APResolution/2));
                            scAPbin = Embryo_scAPbin;
                            
                            if CurrentParticle.FlaggingInfo.SchnitzOn(frame_index)
                                SchnitzUnalignedOnNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                                SchnitzUnalignedTotalNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                                
                            elseif CurrentParticle.FlaggingInfo.SchnitzOff(frame_index)
                                SchnitzUnalignedOffNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                                SchnitzUnalignedTotalNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                            elseif CurrentParticle.FlaggingInfo.SchnitzTemporarilyOff(frame_index)
                                SchnitzUnalignedQuiescentNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                                SchnitzUnalignedTotalNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                            elseif CurrentParticle.FlaggingInfo.SchnitzFinishedTranscribing(frame_index)
                                SchnitzUnalignedFinishedTranscribingNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                                SchnitzUnalignedTotalNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                            end
                            
                            
                            if CurrentParticle.FlaggingInfo.SchnitzOn(frame_index)
                                for grade_index = 1:length(pOnCutoffs)
                                    if CurrentParticle.FlaggingInfo.UseTraceFluo(frame_index)
                                        SchnitzUnalignedGradedTotalNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                                        if CurrentParticle.FlaggingInfo.FluoVector(frame_index) > pOnCutoffs(grade_index)*MaxParticleFluo
                                            SchnitzUnalignedGradedOnNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                                        else
                                            SchnitzUnalignedGradedOffNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                                        end
                                    end
                                end
                                
                            elseif CurrentParticle.FlaggingInfo.SchnitzOff(frame_index)| ...
                                    CurrentParticle.FlaggingInfo.SchnitzTemporarilyOff(frame_index) | ...
                                    CurrentParticle.FlaggingInfo.SchnitzFinishedTranscribing(frame_index)
                                for grade_index = 1:length(pOnCutoffs)
                                    SchnitzUnalignedGradedOffNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                                    SchnitzUnalignedGradedTotalNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8,grade_index, sc_idx) = 1;
                                end
                            end
                            
                            
                            
                            
                        end
                    end
                    
                    if all(CurrentParticle.FlaggingInfo.SchnitzOff(CurrentParticle.FlaggingInfo.SpotStateDefinitive))
                        SchnitzUnalignedInactiveNuclei(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                    else
                        SchnitzUnalignedActiveNuclei(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                    end
                    SchnitzUnalignedNucleiCount(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                    
                    for grade_index = 1:length(pOnCutoffs)
                        if all(SchnitzUnalignedGradedTotalNuclei(:, scAPbin, scNC-8,grade_index, sc_idx) ==0 | ...
                                (SchnitzUnalignedGradedTotalNuclei(:, scAPbin, scNC-8,grade_index, sc_idx) ==1 & SchnitzUnalignedGradedOffNuclei(:, scAPbin, scNC-8,grade_index, sc_idx) ==1 ))
                            SchnitzUnalignedGradedInactiveNuclei(Embryo_scAPbin, scNC-8, grade_index, sc_idx) = 1;
                        else
                            SchnitzUnalignedGradedActiveNuclei(Embryo_scAPbin, scNC-8, grade_index,sc_idx) = 1;
                        end
                        SchnitzUnalignedGradedNucleiCount(Embryo_scAPbin, scNC-8, grade_index, sc_idx) = 1;
                    end
                    
                else
                    for frame_index = 1:length(AllSchnitzFrames)
                        if CurrentParticle.FlaggingInfo.SpotStateDefinitive(frame_index)
                            %scAPbin = find((scAPvector(frame_index) < APbins+APResolution/2) & (scAPvector(frame_index) >= APbins-APResolution/2));
                            scAPbin = Embryo_scAPbin;
                            SchnitzUnalignedTotalNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                            SchnitzUnalignedOffNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                            for grade_index = 1:length(pOnCutoffs)
                                SchnitzUnalignedGradedOffNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                                SchnitzUnalignedGradedTotalNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                            end
                        end
                    end
                    SchnitzUnalignedInactiveNuclei(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                    SchnitzUnalignedNucleiCount(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                    
                    for grade_index = 1:length(pOnCutoffs)
                        SchnitzUnalignedGradedInactiveNuclei(scAPbin, scNC-8, grade_index, sc_idx) = 1;
                        SchnitzUnalignedGradedNucleiCount(scAPbin, scNC-8, grade_index, sc_idx) = 1;
                    end
                end
            elseif CurrentParticle.Approved == -1 & CurrentParticle.FlaggingInfo.SpotStateCanBeCalled
                for frame_index = 1:length(AllSchnitzFrames)
                    if CurrentParticle.FlaggingInfo.SpotStateDefinitive(frame_index)
                        scAPbin = Embryo_scAPbin;
                        for grade_index = 1:length(pOnCutoffs)
                            SchnitzUnalignedGradedOffNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                            SchnitzUnalignedGradedTotalNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                        end
                        SchnitzUnalignedTotalNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                        SchnitzUnalignedOffNuclei(AllSchnitzFrames(frame_index)-SchnitzFrameMin+1, scAPbin, scNC-8, sc_idx) = 1;
                        
                    end
                end
                SchnitzUnalignedInactiveNuclei(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                SchnitzUnalignedNucleiCount(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                for grade_index = 1:length(pOnCutoffs)
                    SchnitzUnalignedGradedInactiveNuclei(scAPbin, scNC-8, grade_index, sc_idx) = 1;
                    SchnitzUnalignedGradedNucleiCount(scAPbin, scNC-8, grade_index, sc_idx) = 1;
                end
            end
        end
    end
    
    CandidateAnaphaseFrames = SchnitzAnaphaseFrames(SchnitzCycles == 14 & abs(SchnitzAnaphaseFrames-anaphaseFrames(14-8)) <= 3 & SchnitzAnaphaseFrames > 0);
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
    SchnitzTbinnedActiveNuclei = zeros(length(APbins), 6, NumSchnitz);
    SchnitzTbinnedInactiveNuclei = zeros(length(APbins), 6, NumSchnitz);
    SchnitzTbinnedNucleiCount = zeros(length(APbins), 6, NumSchnitz);
    SchnitzTbinnedGradedOnNuclei = zeros(numFrames, length(APbins), 6,10, NumSchnitz);
    SchnitzTbinnedGradedOffNuclei = zeros(numFrames, length(APbins), 6,10, NumSchnitz);
    SchnitzTbinnedGradedTotalNuclei = zeros(numFrames, length(APbins), 6,10, NumSchnitz);
    SchnitzTbinnedGradedActiveNuclei = zeros(length(APbins), 6,10, NumSchnitz);
    SchnitzTbinnedGradedInactiveNuclei = zeros(length(APbins), 6,10, NumSchnitz);
    SchnitzTbinnedGradedNucleiCount = zeros(length(APbins), 6,10, NumSchnitz);
    
    TbinnedTimeVector = 0:deltaTbinWidth:(max(FrameTimes)-FrameTimes(nc_info(6)));
    if length(TbinnedTimeVector) < numBinnedFrames
        TbinnedTimeVector(end+1) = max(TbinnedTimeVector)+deltaTbinWidth;
    end
    TbinnedTimeVecIndex = 1:length(TbinnedTimeVector);
    
    for sc_idx =1:length(schnitzcells)
        if ~SchnitzesAreGood(sc_idx)
            continue
        end
        
        scAP = mean(schnitzcells(sc_idx).APpos);
        scAPvector = schnitzcells(sc_idx).APpos;
        Embryo_scAPbin = find((scAP < APbins+APResolution/2) & (scAP >= APbins-APResolution/2));
        scNC = SchnitzCycles(sc_idx);
        if (scNC < 9) | (scNC > 14) | isempty(scAPvector)
            continue
        end
        
        
        
        
        CandidateAnaphaseFrames = SchnitzAnaphaseFrames(SchnitzCycles == scNC & abs(SchnitzAnaphaseFrames-anaphaseFrames(scNC-8)) <= 3 & SchnitzAnaphaseFrames > 0);
        if ~isempty(CandidateAnaphaseFrames)
            SchnitzFrameMin = min([CandidateAnaphaseFrames anaphaseFrames(scNC-8)]);
        else
            SchnitzFrameMin = anaphaseFrames(scNC-8);
        end
        SchnitzHasParticle = ~isnan(SchnitzParticleIndex(sc_idx));
        
        if ~SchnitzHasParticle
            AllSchnitzFrames = schnitzcells(sc_idx).FlaggingInfo.AllSchnitzFrames;
            AllSchnitzTimes = [FrameInfo(AllSchnitzFrames).Time]-FrameInfo(SchnitzFrameMin).Time;
            %AllscAPvector = NaN(1, length(AllSchnitzFrames));
            %AllscAPvector(find(ismember(schnitzcells(sc_idx).frames', schnitzcells(sc_idx).FlaggingInfo.AllSchnitzFrames))) = scAPvector;
            %InterpAPvec = interp1(AllSchnitzTimes,AllscAPvector,TbinnedTimeVector);
            
            DefinitiveSchnitzDouble = double(schnitzcells(sc_idx).FlaggingInfo.SpotStateDefinitive);
            InterpDefinitiveSchnitz = interp1(AllSchnitzTimes,DefinitiveSchnitzDouble,TbinnedTimeVector);
            
            
            for frame_index = 1:length(TbinnedTimeVector)
                if InterpDefinitiveSchnitz(frame_index) == 1 %&  ~isnan(InterpAPvec(frame_index))
                    %scAPbin = find((InterpAPvec(frame_index) < APbins+APResolution/2) & (InterpAPvec(frame_index) >= APbins-APResolution/2));
                    scAPbin = Embryo_scAPbin;
                    SchnitzTbinnedTotalNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                    SchnitzTbinnedOffNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                    
                    for grade_index = 1:length(pOnCutoffs)
                        SchnitzTbinnedGradedOffNuclei(frame_index, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                        SchnitzTbinnedGradedTotalNuclei(frame_index, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                        
                    end
                end
            end
            SchnitzTbinnedNucleiCount(Embryo_scAPbin, scNC-8, sc_idx) = 1;
            SchnitzTbinnedInactiveNuclei(Embryo_scAPbin, scNC-8, sc_idx) = 1;
            
            for grade_index = 1:length(pOnCutoffs)
                SchnitzTbinnedGradedInactiveNuclei(Embryo_scAPbin, scNC-8, grade_index, sc_idx) = 1;
                SchnitzTbinnedGradedNucleiCount(Embryo_scAPbin, scNC-8, grade_index, sc_idx) = 1;
            end
            
        else
            p_idx = SchnitzParticleIndex(sc_idx);
            CurrentParticle = CompiledParticles{ChN}(p_idx);
            CurrentSchnitz = schnitzcells(sc_idx);
            AllSchnitzFrames = CurrentParticle.FlaggingInfo.AllSchnitzFrames;
            AllSchnitzTimes = [FrameInfo(AllSchnitzFrames).Time]-FrameInfo(SchnitzFrameMin).Time;
            if ParticlesAreGood(p_idx)
                if (sum(CurrentParticle.FlaggingInfo.SchnitzOn) >= MinTimePoints) & (scNC >= 9) & (scNC <= 14) & ~isempty(scAPvector)
                    %AllscAPvectorSchnitz = NaN(1, length(CurrentSchnitz.FlaggingInfo.AllSchnitzFrames));
                    %AllscAPvectorSchnitz(find(ismember(CurrentSchnitz.frames', CurrentSchnitz.FlaggingInfo.AllSchnitzFrames))) = scAPvector;
                    %AllscAPvector =  NaN(1, length(CurrentParticle.FlaggingInfo.AllSchnitzFrames));
                    %AllscAPvector(find(ismember(CurrentSchnitz.FlaggingInfo.AllSchnitzFrames, CurrentParticle.FlaggingInfo.AllSchnitzFrames))) = AllscAPvectorSchnitz;
                    
                    %InterpAPvec = interp1(AllSchnitzTimes,AllscAPvector,TbinnedTimeVector);
                    
                    DefinitiveSchnitzDouble = double(CurrentParticle.FlaggingInfo.SpotStateDefinitive);
                    InterpDefinitiveSchnitz = interp1(AllSchnitzTimes,DefinitiveSchnitzDouble,TbinnedTimeVector);
                    
                    DefinitiveSchnitzOff = double(CurrentParticle.FlaggingInfo.SchnitzOff);
                    InterpSchnitzOff = interp1(AllSchnitzTimes,DefinitiveSchnitzOff,TbinnedTimeVector);
                    DefinitiveSchnitzOn = double(CurrentParticle.FlaggingInfo.SchnitzOn);
                    InterpSchnitzOn = interp1(AllSchnitzTimes,DefinitiveSchnitzOn,TbinnedTimeVector);
                    DefinitiveSchnitzQuiescent = double(CurrentParticle.FlaggingInfo.SchnitzTemporarilyOff);
                    InterpSchnitzQuiescent = interp1(AllSchnitzTimes,DefinitiveSchnitzQuiescent,TbinnedTimeVector);
                    DefinitiveSchnitzFinishedTranscribing = double(CurrentParticle.FlaggingInfo.SchnitzFinishedTranscribing);
                    InterpSchnitzFinishedTranscribing = interp1(AllSchnitzTimes,DefinitiveSchnitzFinishedTranscribing,TbinnedTimeVector);
                    DefinitiveUseTraceFluo = double(CurrentParticle.FlaggingInfo.UseTraceFluo);
                    InterpUseTraceFluo = interp1(AllSchnitzTimes,DefinitiveUseTraceFluo,TbinnedTimeVector);
                    DefinitiveFluoVector = double(CurrentParticle.FlaggingInfo.FluoVector);
                    InterpUseTraceFluo = interp1(AllSchnitzTimes,DefinitiveFluoVector,TbinnedTimeVector);
                    
                    for frame_index = 1:length(TbinnedTimeVector)
                        if (InterpDefinitiveSchnitz(frame_index) == 1) %% &  ~isnan(InterpAPvec(frame_index))
                            %scAPbin = find((InterpAPvec(frame_index) < APbins+APResolution/2) & (InterpAPvec(frame_index) >= APbins-APResolution/2));
                            scAPbin = Embryo_scAPbin;
                            if (InterpSchnitzOn(frame_index) == 1)
                                SchnitzTbinnedOnNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                                SchnitzTbinnedTotalNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                            elseif (InterpSchnitzOff(frame_index) == 1)
                                SchnitzTbinnedOffNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                                SchnitzTbinnedTotalNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                            elseif (InterpSchnitzQuiescent(frame_index) == 1)
                                SchnitzTbinnedQuiescentNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                                SchnitzTbinnedTotalNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                            elseif (InterpSchnitzFinishedTranscribing(frame_index) == 1)
                                SchnitzTbinnedFinishedTranscribingNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                                SchnitzTbinnedTotalNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                            end
                            if (InterpSchnitzOn(frame_index) == 1)
                                for grade_index = 1:length(pOnCutoffs)
                                    if InterpUseTraceFluo(frame_index) == 1
                                        SchnitzTbinnedGradedTotalNuclei(frame_index, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                                        if DefinitiveFluoVector(frame_index) > pOnCutoffs(grade_index)*MaxParticleFluo
                                            SchnitzTbinnedGradedOnNuclei(frame_index, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                                        else
                                            SchnitzTbinnedGradedOffNuclei(frame_index, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                                        end
                                    end
                                end
                            elseif (InterpSchnitzOff(frame_index) == 1)| ...
                                    (InterpSchnitzQuiescent(frame_index) == 1) | ...
                                    (InterpSchnitzFinishedTranscribing(frame_index) == 1)
                                for grade_index = 1:length(pOnCutoffs)
                                    SchnitzTbinnedGradedOffNuclei(frame_index, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                                    SchnitzTbinnedGradedTotalNuclei(frame_index, scAPbin, scNC-8,grade_index, sc_idx) = 1;
                                end
                            end
                            
                            
                        end
                        
                    end
                    
                    if all(InterpSchnitzOff(InterpDefinitiveSchnitz == 1))
                        SchnitzTbinnedInactiveNuclei(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                    else
                        SchnitzTbinnedActiveNuclei(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                    end
                    SchnitzTbinnedNucleiCount(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                    
                    for grade_index = 1:length(pOnCutoffs)
                        if all(SchnitzTbinnedGradedTotalNuclei(:, scAPbin, scNC-8,grade_index, sc_idx) ==0 | ...
                                (SchnitzTbinnedGradedTotalNuclei(:, scAPbin, scNC-8,grade_index, sc_idx) ==1 & SchnitzTbinnedGradedOffNuclei(:, scAPbin, scNC-8,grade_index, sc_idx) ==1 ))
                            SchnitzTbinnedGradedInactiveNuclei(Embryo_scAPbin, scNC-8, grade_index, sc_idx) = 1;
                        else
                            SchnitzTbinnedGradedActiveNuclei(Embryo_scAPbin, scNC-8, grade_index,sc_idx) = 1;
                        end
                        SchnitzTbinnedGradedNucleiCount(Embryo_scAPbin, scNC-8, grade_index, sc_idx) = 1;
                    end
                    
                else
                    %AllscAPvectorSchnitz = NaN(1, length(CurrentSchnitz.FlaggingInfo.AllSchnitzFrames));
                    %AllscAPvectorSchnitz(find(ismember(CurrentSchnitz.frames', CurrentSchnitz.FlaggingInfo.AllSchnitzFrames))) = scAPvector;
                    %AllscAPvector =  NaN(1, length(CurrentParticle.FlaggingInfo.AllSchnitzFrames));
                    %AllscAPvector(find(ismember(CurrentSchnitz.FlaggingInfo.AllSchnitzFrames, CurrentParticle.FlaggingInfo.AllSchnitzFrames))) = AllscAPvectorSchnitz;
                    
                    %InterpAPvec = interp1(AllSchnitzTimes,AllscAPvector,TbinnedTimeVector);
                    DefinitiveSchnitzDouble = double(CurrentParticle.FlaggingInfo.SpotStateDefinitive);
                    InterpDefinitiveSchnitz = interp1(AllSchnitzTimes,DefinitiveSchnitzDouble,TbinnedTimeVector);
                    
                    DefinitiveSchnitzOff = double(CurrentParticle.FlaggingInfo.SchnitzOff);
                    
                    for frame_index = 1:length(TbinnedTimeVector)
                        if (InterpDefinitiveSchnitz(frame_index) == 1)% &  ~isnan(InterpAPvec(frame_index))
                            %scAPbin = find((InterpAPvec(frame_index) < APbins+APResolution/2) & (InterpAPvec(frame_index) >= APbins-APResolution/2));
                            scAPbin = Embryo_scAPbin;
                            SchnitzTbinnedOffNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                            SchnitzTbinnedTotalNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                            for grade_index = 1:length(pOnCutoffs)
                                SchnitzTbinnedGradedOffNuclei(frame_index, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                                SchnitzTbinnedGradedTotalNuclei(frame_index, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                            end
                        end
                    end
                    SchnitzTbinnedInactiveNuclei(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                    SchnitzTbinnedNucleiCount(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                    for grade_index = 1:length(pOnCutoffs)
                        SchnitzTbinnedGradedInactiveNuclei(scAPbin, scNC-8, grade_index, sc_idx) = 1;
                        SchnitzTbinnedGradedNucleiCount(scAPbin, scNC-8, grade_index, sc_idx) = 1;
                    end
                end
            elseif CurrentParticle.Approved == -1 & CurrentParticle.FlaggingInfo.SpotStateCanBeCalled
                %AllscAPvectorSchnitz = NaN(1, length(CurrentSchnitz.FlaggingInfo.AllSchnitzFrames));
                %AllscAPvectorSchnitz(find(ismember(CurrentSchnitz.frames', CurrentSchnitz.FlaggingInfo.AllSchnitzFrames))) = scAPvector;
                %AllscAPvector =  NaN(1, length(CurrentParticle.FlaggingInfo.AllSchnitzFrames));
                %AllscAPvector(find(ismember(CurrentSchnitz.FlaggingInfo.AllSchnitzFrames, CurrentParticle.FlaggingInfo.AllSchnitzFrames))) = AllscAPvectorSchnitz;
                
                %InterpAPvec = interp1(AllSchnitzTimes,AllscAPvector,TbinnedTimeVector);
                DefinitiveSchnitzDouble = double(CurrentParticle.FlaggingInfo.SpotStateDefinitive);
                InterpDefinitiveSchnitz = interp1(AllSchnitzTimes,DefinitiveSchnitzDouble,TbinnedTimeVector);
                
                for frame_index = 1:length(AllSchnitzFrames)
                    if (InterpDefinitiveSchnitz(frame_index) == 1) %&  ~isnan(InterpAPvec(frame_index))
                        %scAPbin = find((InterpAPvec(frame_index) < APbins+APResolution/2) & (InterpAPvec(frame_index) >= APbins-APResolution/2));
                        scAPbin = Embryo_scAPbin;
                        SchnitzTbinnedTotalNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                        SchnitzTbinnedOffNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                        for grade_index = 1:length(pOnCutoffs)
                            SchnitzTbinnedGradedOffNuclei(frame_index, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                            SchnitzTbinnedGradedTotalNuclei(frame_index, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                        end
                    end
                end
                
                SchnitzTbinnedInactiveNuclei(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                SchnitzTbinnedNucleiCount(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                for grade_index = 1:length(pOnCutoffs)
                    SchnitzTbinnedGradedInactiveNuclei(scAPbin, scNC-8, grade_index, sc_idx) = 1;
                    SchnitzTbinnedGradedNucleiCount(scAPbin, scNC-8, grade_index, sc_idx) = 1;
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
    SchnitzAnaphaseAlignedActiveNuclei = zeros(length(APbins), 6, NumSchnitz);
    SchnitzAnaphaseAlignedInactiveNuclei = zeros(length(APbins), 6, NumSchnitz);
    SchnitzAnaphaseAlignedNucleiCount = zeros(length(APbins), 6, NumSchnitz);
    SchnitzAnaphaseAlignedGradedOnNuclei = zeros(numFrames, length(APbins), 6,10, NumSchnitz);
    SchnitzAnaphaseAlignedGradedOffNuclei = zeros(numFrames, length(APbins), 6,10, NumSchnitz);
    SchnitzAnaphaseAlignedGradedTotalNuclei = zeros(numFrames, length(APbins), 6,10, NumSchnitz);
    SchnitzAnaphaseAlignedGradedActiveNuclei = zeros(length(APbins), 6,10, NumSchnitz);
    SchnitzAnaphaseAlignedGradedInactiveNuclei = zeros(length(APbins), 6,10, NumSchnitz);
    SchnitzAnaphaseAlignedGradedNucleiCount = zeros(length(APbins), 6,10, NumSchnitz);
    
    for sc_idx =1:length(schnitzcells)
        if ~SchnitzesAreGood(sc_idx)
            continue
        end
        
        scAP = mean(schnitzcells(sc_idx).APpos);
        scAPvector = schnitzcells(sc_idx).APpos;
        Embryo_scAPbin = find((scAP < APbins+APResolution/2) & (scAP >= APbins-APResolution/2));
        scNC = SchnitzCycles(sc_idx);
        if (scNC < 9) | (scNC > 14) | isempty(scAPvector)
            continue
        end
        
        
        
        
        CandidateAnaphaseFrames = SchnitzAnaphaseFrames(SchnitzCycles == scNC & abs(SchnitzAnaphaseFrames-anaphaseFrames(scNC-8)) <= 3 & SchnitzAnaphaseFrames > 0);
        if ~isempty(CandidateAnaphaseFrames)
            SchnitzFrameMin = min([CandidateAnaphaseFrames anaphaseFrames(scNC-8)]);
        else
            SchnitzFrameMin = anaphaseFrames(scNC-8);
        end
        SchnitzHasParticle = ~isnan(SchnitzParticleIndex(sc_idx));
        
        if ~SchnitzHasParticle
            AllSchnitzFrames = schnitzcells(sc_idx).FlaggingInfo.AllSchnitzFrames;
            AllSchnitzTimes = [FrameInfo(AllSchnitzFrames).Time]-FrameInfo(schnitzcells(sc_idx).anaphaseFrame).Time;
            %AllscAPvector = NaN(1, length(AllSchnitzFrames));
            %AllscAPvector(find(ismember(schnitzcells(sc_idx).frames', schnitzcells(sc_idx).FlaggingInfo.AllSchnitzFrames))) = scAPvector;
            %InterpAPvec = interp1(AllSchnitzTimes,AllscAPvector,TbinnedTimeVector);
            DefinitiveSchnitzDouble = double(schnitzcells(sc_idx).FlaggingInfo.SpotStateDefinitive);
            InterpDefinitiveSchnitz = interp1(AllSchnitzTimes,DefinitiveSchnitzDouble,TbinnedTimeVector);
            
            
            for frame_index = 1:length(TbinnedTimeVector)
                if (InterpDefinitiveSchnitz(frame_index) == 1) %&  ~isnan(InterpAPvec(frame_index))
                    %scAPbin = find((InterpAPvec(frame_index) < APbins+APResolution/2) & (InterpAPvec(frame_index) >= APbins-APResolution/2));
                    scAPbin = Embryo_scAPbin;
                    SchnitzAnaphaseAlignedTotalNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                    SchnitzAnaphaseAlignedOffNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                end
                for grade_index = 1:length(pOnCutoffs)
                    SchnitzAnaphaseAlignedGradedOffNuclei(frame_index, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                    SchnitzAnaphaseAlignedGradedTotalNuclei(frame_index, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                    
                end
            end
            SchnitzAnaphaseAlignedNucleiCount(Embryo_scAPbin, scNC-8, sc_idx) = 1;
            SchnitzAnaphaseAlignedInactiveNuclei(Embryo_scAPbin, scNC-8, sc_idx) = 1;
            for grade_index = 1:length(pOnCutoffs)
                SchnitzAnaphaseAlignedGradedInactiveNuclei(Embryo_scAPbin, scNC-8, grade_index, sc_idx) = 1;
                SchnitzAnaphaseAlignedGradedNucleiCount(Embryo_scAPbin, scNC-8, grade_index, sc_idx) = 1;
            end
            
        else
            p_idx = SchnitzParticleIndex(sc_idx);
            CurrentParticle = CompiledParticles{ChN}(p_idx);
            CurrentSchnitz = schnitzcells(sc_idx);
            AllSchnitzFrames = CurrentParticle.FlaggingInfo.AllSchnitzFrames;
            AllSchnitzTimes = [FrameInfo(AllSchnitzFrames).Time]-FrameInfo(CurrentParticle.schnitzcell.anaphaseFrame).Time;
            if ParticlesAreGood(p_idx)
                if (sum(CurrentParticle.FlaggingInfo.SchnitzOn) >= MinTimePoints) & (scNC >= 9) & (scNC <= 14) & ~isempty(scAPvector)
                    %AllscAPvectorSchnitz = NaN(1, length(CurrentSchnitz.FlaggingInfo.AllSchnitzFrames));
                    %AllscAPvectorSchnitz(find(ismember(CurrentSchnitz.frames', CurrentSchnitz.FlaggingInfo.AllSchnitzFrames))) = scAPvector;
                    %AllscAPvector =  NaN(1, length(CurrentParticle.FlaggingInfo.AllSchnitzFrames));
                    %AllscAPvector(find(ismember(CurrentSchnitz.FlaggingInfo.AllSchnitzFrames, CurrentParticle.FlaggingInfo.AllSchnitzFrames))) = AllscAPvectorSchnitz;
                    
                    %InterpAPvec = interp1(AllSchnitzTimes,AllscAPvector,TbinnedTimeVector);
                    DefinitiveSchnitzDouble = double(CurrentParticle.FlaggingInfo.SpotStateDefinitive);
                    InterpDefinitiveSchnitz = interp1(AllSchnitzTimes,DefinitiveSchnitzDouble,TbinnedTimeVector);
                    
                    DefinitiveSchnitzOff = double(CurrentParticle.FlaggingInfo.SchnitzOff);
                    InterpSchnitzOff = interp1(AllSchnitzTimes,DefinitiveSchnitzOff,TbinnedTimeVector);
                    DefinitiveSchnitzOn = double(CurrentParticle.FlaggingInfo.SchnitzOn);
                    InterpSchnitzOn = interp1(AllSchnitzTimes,DefinitiveSchnitzOn,TbinnedTimeVector);
                    DefinitiveSchnitzQuiescent = double(CurrentParticle.FlaggingInfo.SchnitzTemporarilyOff);
                    InterpSchnitzQuiescent = interp1(AllSchnitzTimes,DefinitiveSchnitzQuiescent,TbinnedTimeVector);
                    DefinitiveSchnitzFinishedTranscribing = double(CurrentParticle.FlaggingInfo.SchnitzFinishedTranscribing);
                    InterpSchnitzFinishedTranscribing = interp1(AllSchnitzTimes,DefinitiveSchnitzFinishedTranscribing,TbinnedTimeVector);
                    DefinitiveUseTraceFluo = double(CurrentParticle.FlaggingInfo.UseTraceFluo);
                    InterpUseTraceFluo = interp1(AllSchnitzTimes,DefinitiveUseTraceFluo,TbinnedTimeVector);
                    DefinitiveFluoVector = double(CurrentParticle.FlaggingInfo.FluoVector);
                    InterpUseTraceFluo = interp1(AllSchnitzTimes,DefinitiveFluoVector,TbinnedTimeVector);
                    
                    for frame_index = 1:length(TbinnedTimeVector)
                        if (InterpDefinitiveSchnitz(frame_index) == 1)% &  ~isnan(InterpAPvec(frame_index))
                            %scAPbin = find((InterpAPvec(frame_index) < APbins+APResolution/2) & (InterpAPvec(frame_index) >= APbins-APResolution/2));
                            scAPbin = Embryo_scAPbin;
                            if (InterpSchnitzOn(frame_index) == 1)
                                SchnitzAnaphaseAlignedOnNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                                SchnitzAnaphaseAlignedTotalNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                            elseif (InterpSchnitzOff(frame_index) == 1)
                                SchnitzAnaphaseAlignedOffNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                                SchnitzAnaphaseAlignedTotalNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                            elseif (InterpSchnitzQuiescent(frame_index) == 1)
                                SchnitzAnaphaseAlignedQuiescentNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                                SchnitzAnaphaseAlignedTotalNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                            elseif (InterpSchnitzFinishedTranscribing(frame_index) == 1)
                                SchnitzAnaphaseAlignedFinishedTranscribingNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                                SchnitzAnaphaseAlignedTotalNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                            end
                            
                            if (InterpSchnitzOn(frame_index) == 1)
                                for grade_index = 1:length(pOnCutoffs)
                                    if InterpUseTraceFluo(frame_index) == 1
                                        SchnitzAnaphaseAlignedGradedTotalNuclei(frame_index, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                                        if DefinitiveFluoVector(frame_index) > pOnCutoffs(grade_index)*MaxParticleFluo
                                            SchnitzAnaphaseAlignedGradedOnNuclei(frame_index, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                                        else
                                            SchnitzAnaphaseAlignedGradedOffNuclei(frame_index, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                                        end
                                    end
                                end
                            elseif (InterpSchnitzOff(frame_index) == 1)| ...
                                    (InterpSchnitzQuiescent(frame_index) == 1) | ...
                                    (InterpSchnitzFinishedTranscribing(frame_index) == 1)
                                for grade_index = 1:length(pOnCutoffs)
                                    SchnitzAnaphaseAlignedGradedOffNuclei(frame_index, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                                    SchnitzAnaphaseAlignedGradedTotalNuclei(frame_index, scAPbin, scNC-8,grade_index, sc_idx) = 1;
                                end
                            end
                            
                        end
                        
                        
                    end
                    
                    if all(InterpSchnitzOff(InterpDefinitiveSchnitz == 1))
                        SchnitzAnaphaseAlignedInactiveNuclei(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                    else
                        SchnitzAnaphaseAlignedActiveNuclei(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                    end
                    SchnitzAnaphaseAlignedNucleiCount(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                    
                    for grade_index = 1:length(pOnCutoffs)
                        if all(SchnitzAnaphaseAlignedGradedTotalNuclei(:, scAPbin, scNC-8,grade_index, sc_idx) ==0 | ...
                                (SchnitzAnaphaseAlignedGradedTotalNuclei(:, scAPbin, scNC-8,grade_index, sc_idx) ==1 & SchnitzAnaphaseAlignedGradedOffNuclei(:, scAPbin, scNC-8,grade_index, sc_idx) ==1 ))
                            SchnitzAnaphaseAlignedGradedInactiveNuclei(Embryo_scAPbin, scNC-8, grade_index, sc_idx) = 1;
                        else
                            SchnitzAnaphaseAlignedGradedActiveNuclei(Embryo_scAPbin, scNC-8, grade_index,sc_idx) = 1;
                        end
                        SchnitzAnaphaseAlignedGradedNucleiCount(Embryo_scAPbin, scNC-8, grade_index, sc_idx) = 1;
                    end
                    
                else
                    %AllscAPvectorSchnitz = NaN(1, length(CurrentSchnitz.FlaggingInfo.AllSchnitzFrames));
                    %AllscAPvectorSchnitz(find(ismember(CurrentSchnitz.frames', CurrentSchnitz.FlaggingInfo.AllSchnitzFrames))) = scAPvector;
                    %AllscAPvector =  NaN(1, length(CurrentParticle.FlaggingInfo.AllSchnitzFrames));
                    %AllscAPvector(find(ismember(CurrentSchnitz.FlaggingInfo.AllSchnitzFrames, CurrentParticle.FlaggingInfo.AllSchnitzFrames))) = AllscAPvectorSchnitz;
                    
                    %InterpAPvec = interp1(AllSchnitzTimes,AllscAPvector,TbinnedTimeVector);
                    DefinitiveSchnitzDouble = double(CurrentParticle.FlaggingInfo.SpotStateDefinitive);
                    InterpDefinitiveSchnitz = interp1(AllSchnitzTimes,DefinitiveSchnitzDouble,TbinnedTimeVector);
                    
                    DefinitiveSchnitzOff = double(CurrentParticle.FlaggingInfo.SchnitzOff);
                    
                    for frame_index = 1:length(TbinnedTimeVector)
                        if (InterpDefinitiveSchnitz(frame_index) == 1) % &  ~isnan(InterpAPvec(frame_index))
                            %scAPbin = find((InterpAPvec(frame_index) < APbins+APResolution/2) & (InterpAPvec(frame_index) >= APbins-APResolution/2));
                            scAPbin = Embryo_scAPbin;
                            SchnitzAnaphaseAlignedOffNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                            SchnitzAnaphaseAlignedTotalNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                            for grade_index = 1:length(pOnCutoffs)
                                SchnitzAnaphaseAlignedGradedOffNuclei(frame_index, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                                SchnitzAnaphaseAlignedGradedTotalNuclei(frame_index, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                            end
                        end
                        
                    end
                    SchnitzAnaphaseAlignedInactiveNuclei(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                    SchnitzAnaphaseAlignedNucleiCount(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                    for grade_index = 1:length(pOnCutoffs)
                        SchnitzAnaphaseAlignedGradedInactiveNuclei(scAPbin, scNC-8, grade_index, sc_idx) = 1;
                        SchnitzAnaphaseAlignedGradedNucleiCount(scAPbin, scNC-8, grade_index, sc_idx) = 1;
                    end
                end
            elseif CurrentParticle.Approved == -1 & CurrentParticle.FlaggingInfo.SpotStateCanBeCalled
                %AllscAPvectorSchnitz = NaN(1, length(CurrentSchnitz.FlaggingInfo.AllSchnitzFrames));
                %AllscAPvectorSchnitz(find(ismember(CurrentSchnitz.frames', CurrentSchnitz.FlaggingInfo.AllSchnitzFrames))) = scAPvector;
                %AllscAPvector =  NaN(1, length(CurrentParticle.FlaggingInfo.AllSchnitzFrames));
                %AllscAPvector(find(ismember(CurrentSchnitz.FlaggingInfo.AllSchnitzFrames, CurrentParticle.FlaggingInfo.AllSchnitzFrames))) = AllscAPvectorSchnitz;
                
                %InterpAPvec = interp1(AllSchnitzTimes,AllscAPvector,TbinnedTimeVector);
                DefinitiveSchnitzDouble = double(CurrentParticle.FlaggingInfo.SpotStateDefinitive);
                InterpDefinitiveSchnitz = interp1(AllSchnitzTimes,DefinitiveSchnitzDouble,TbinnedTimeVector);
                
                for frame_index = 1:length(AllSchnitzFrames)
                    if (InterpDefinitiveSchnitz(frame_index) == 1)% &  ~isnan(InterpAPvec(frame_index))
                        %scAPbin = find((InterpAPvec(frame_index) < APbins+APResolution/2) & (InterpAPvec(frame_index) >= APbins-APResolution/2));
                        scAPbin = Embryo_scAPbin;
                        SchnitzAnaphaseAlignedTotalNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                        SchnitzAnaphaseAlignedOffNuclei(frame_index, scAPbin, scNC-8, sc_idx) = 1;
                        for grade_index = 1:length(pOnCutoffs)
                            SchnitzAnaphaseAlignedGradedOffNuclei(frame_index, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                            SchnitzAnaphaseAlignedGradedTotalNuclei(frame_index, scAPbin, scNC-8, grade_index, sc_idx) = 1;
                        end
                    end
                end
                
                SchnitzAnaphaseAlignedInactiveNuclei(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                SchnitzAnaphaseAlignedNucleiCount(Embryo_scAPbin, scNC-8, sc_idx) = 1;
                for grade_index = 1:length(pOnCutoffs)
                    SchnitzAnaphaseAlignedGradedInactiveNuclei(scAPbin, scNC-8, grade_index, sc_idx) = 1;
                    SchnitzAnaphaseAlignedGradedNucleiCount(scAPbin, scNC-8, grade_index, sc_idx) = 1;
                end
            end
        end
    end
    
    
    
    %%
    
    
    UnalignedCycleMeanTraces = NaN(numFrames, length(APbins), 6);
    UnalignedCycleTraceStdErrors =  NaN(numFrames, length(APbins), 6);
    UnalignedCycleTraceCount=  NaN(numFrames, length(APbins), 6);
    UnalignedCycleNumNuclei =  sum(SchnitzUnalignedTotalNuclei, 4);
    UnalignedCycleNumOnNuclei =  sum(SchnitzUnalignedOnNuclei, 4);
    UnalignedCycleGradedOnNuclei = sum(SchnitzUnalignedGradedOnNuclei, 5);
    UnalignedCycleGradedOffNuclei = sum(SchnitzUnalignedGradedOffNuclei, 5);
    UnalignedCycleGradedTotalNuclei = sum(SchnitzUnalignedGradedTotalNuclei, 5);
    UnalignedCycleGradedFractionOn = UnalignedCycleGradedOnNuclei./UnalignedCycleGradedTotalNuclei;
    
    UnalignedCycleFractionOn=  UnalignedCycleNumOnNuclei./UnalignedCycleNumNuclei;
    UnalignedCycleNumOffNuclei =  sum(SchnitzUnalignedOffNuclei, 4);
    
    UnalignedCycleNumQuiescentNuclei =  sum(SchnitzUnalignedQuiescentNuclei, 4);
    UnalignedCycleNumFinishedTranscribingNuclei =  sum(SchnitzUnalignedFinishedTranscribingNuclei, 4);
    UnalignedCycleActiveNuclei = sum(SchnitzUnalignedActiveNuclei, 3);
    UnalignedCycleInactiveNuclei = sum(SchnitzUnalignedInactiveNuclei, 3);
    UnalignedCycleNucleiCount = sum(SchnitzUnalignedNucleiCount, 3);
    UnalignedCycleGradedActiveNuclei = sum(SchnitzUnalignedGradedActiveNuclei, 4);
    UnalignedCycleGradedInactiveNuclei = sum(SchnitzUnalignedGradedInactiveNuclei, 4);
    UnalignedCycleGradedNucleiCount = sum(SchnitzUnalignedGradedNucleiCount, 4);
    
    Unaligned3DCycleMeanTraces = NaN(numFrames, length(APbins), 6);
    Unaligned3DCycleTraceStdErrors =  NaN(numFrames, length(APbins), 6);
    Unaligned3DCycleTraceCount=  NaN(numFrames, length(APbins), 6);
    Unaligned3DCycleNumNuclei =  sum(SchnitzUnalignedOnNuclei, 4);
    Unaligned3DCycleNumOnNuclei =  sum(SchnitzUnalignedOnNuclei, 4);
    Unaligned3DCycleNumOffNuclei =  sum(SchnitzUnalignedOffNuclei, 4);
    Unaligned3DCycleNumQuiescentNuclei =  sum(SchnitzUnalignedQuiescentNuclei, 4);
    Unaligned3DCycleNumFinishedTranscribingNuclei =  sum(SchnitzUnalignedFinishedTranscribingNuclei, 4);
    Unaligned3DCycleFractionOn=  Unaligned3DCycleNumOnNuclei./Unaligned3DCycleNumNuclei;
    Unaligned3DCycleActiveNuclei = sum(SchnitzUnalignedActiveNuclei, 3);
    Unaligned3DCycleInactiveNuclei = sum(SchnitzUnalignedInactiveNuclei, 3);
    Unaligned3DCycleNucleiCount = sum(SchnitzUnalignedNucleiCount, 3);
    Unaligned3DCycleGradedOnNuclei = sum(SchnitzUnalignedGradedOnNuclei, 5);
    Unaligned3DCycleGradedOffNuclei = sum(SchnitzUnalignedGradedOffNuclei, 5);
    Unaligned3DCycleGradedTotalNuclei = sum(SchnitzUnalignedGradedTotalNuclei, 5);
    Unaligned3DCycleGradedFractionOn = UnalignedCycleGradedOnNuclei./UnalignedCycleGradedTotalNuclei;
    Unaligned3DCycleGradedActiveNuclei = sum(SchnitzUnalignedGradedActiveNuclei, 4);
    Unaligned3DCycleGradedInactiveNuclei = sum(SchnitzUnalignedGradedInactiveNuclei, 4);
    Unaligned3DCycleGradedNucleiCount = sum(SchnitzUnalignedGradedNucleiCount, 4);
    
    
    
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
    AnaphaseAlignedCycleActiveNuclei = sum(SchnitzAnaphaseAlignedActiveNuclei, 3);
    AnaphaseAlignedCycleInactiveNuclei = sum(SchnitzAnaphaseAlignedInactiveNuclei, 3);
    AnaphaseAlignedCycleNucleiCount = sum(SchnitzAnaphaseAlignedNucleiCount, 3);
    AnaphaseAlignedCycleGradedOnNuclei = sum(SchnitzAnaphaseAlignedGradedOnNuclei, 5);
    AnaphaseAlignedCycleGradedOffNuclei = sum(SchnitzAnaphaseAlignedGradedOffNuclei, 5);
    AnaphaseAlignedCycleGradedTotalNuclei = sum(SchnitzAnaphaseAlignedGradedTotalNuclei, 5);
    AnaphaseAlignedCycleGradedFractionOn = AnaphaseAlignedCycleGradedOnNuclei./AnaphaseAlignedCycleGradedTotalNuclei;
    AnaphaseAlignedCycleGradedActiveNuclei = sum(SchnitzAnaphaseAlignedGradedActiveNuclei, 4);
    AnaphaseAlignedCycleGradedInactiveNuclei = sum(SchnitzAnaphaseAlignedGradedInactiveNuclei, 4);
    AnaphaseAlignedCycleGradedNucleiCount = sum(SchnitzAnaphaseAlignedGradedNucleiCount, 4);
    
    AnaphaseAligned3DCycleMeanTraces = NaN(numBinnedFrames, length(APbins), 6);
    AnaphaseAligned3DCycleTraceStdErrors =  NaN(numBinnedFrames, length(APbins), 6);
    AnaphaseAligned3DCycleTraceCount=  NaN(numBinnedFrames, length(APbins), 6);
    AnaphaseAligned3DCycleNumNuclei =  sum(SchnitzAnaphaseAlignedTotalNuclei, 4);
    AnaphaseAligned3DCycleNumOnNuclei = sum(SchnitzAnaphaseAlignedOnNuclei, 4);
    AnaphaseAligned3DCycleFractionOn=  AnaphaseAligned3DCycleNumOnNuclei./AnaphaseAligned3DCycleNumNuclei;
    AnaphaseAligned3DCycleNumOffNuclei =  sum(SchnitzAnaphaseAlignedOffNuclei, 4);
    AnaphaseAligned3DCycleNumQuiescentNuclei =  sum(SchnitzAnaphaseAlignedQuiescentNuclei, 4);
    AnaphaseAligned3DCycleNumFinishedTranscribingNuclei =  sum(SchnitzAnaphaseAlignedFinishedTranscribingNuclei, 4);
    AnaphaseAligned3DCycleNumFinishedTranscribingNuclei =  sum(SchnitzAnaphaseAlignedFinishedTranscribingNuclei, 4);
    AnaphaseAligned3DCycleActiveNuclei = sum(SchnitzAnaphaseAlignedActiveNuclei, 3);
    AnaphaseAligned3DCycleInactiveNuclei = sum(SchnitzAnaphaseAlignedInactiveNuclei, 3);
    AnaphaseAligned3DCycleNucleiCount = sum(SchnitzAnaphaseAlignedNucleiCount, 3);
    AnaphaseAligned3DCycleGradedOnNuclei = sum(SchnitzAnaphaseAlignedGradedOnNuclei, 5);
    AnaphaseAligned3DCycleGradedOffNuclei = sum(SchnitzAnaphaseAlignedGradedOffNuclei, 5);
    AnaphaseAligned3DCycleGradedTotalNuclei = sum(SchnitzAnaphaseAlignedGradedTotalNuclei, 5);
    AnaphaseAligned3DCycleGradedFractionOn = AnaphaseAlignedCycleGradedOnNuclei./AnaphaseAlignedCycleGradedTotalNuclei;
    AnaphaseAligned3DCycleGradedActiveNuclei = sum(SchnitzAnaphaseAlignedGradedActiveNuclei, 4);
    AnaphaseAligned3DCycleGradedInactiveNuclei = sum(SchnitzAnaphaseAlignedGradedInactiveNuclei, 4);
    AnaphaseAligned3DCycleGradedNucleiCount = sum(SchnitzAnaphaseAlignedGradedNucleiCount, 4);
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
    TbinnedCycleActiveNuclei = sum(SchnitzTbinnedActiveNuclei, 3);
    TbinnedCycleInactiveNuclei = sum(SchnitzTbinnedInactiveNuclei, 3);
    TbinnedCycleNucleiCount = sum(SchnitzTbinnedNucleiCount, 3);
    TbinnedCycleGradedOnNuclei = sum(SchnitzTbinnedGradedOnNuclei, 5);
    TbinnedCycleGradedOffNuclei = sum(SchnitzTbinnedGradedOffNuclei, 5);
    TbinnedCycleGradedTotalNuclei = sum(SchnitzTbinnedGradedTotalNuclei, 5);
    TbinnedCycleGradedFractionOn = TbinnedCycleGradedOnNuclei./TbinnedCycleGradedTotalNuclei;
    TbinnedCycleGradedActiveNuclei = sum(SchnitzTbinnedGradedActiveNuclei, 4);
    TbinnedCycleGradedInactiveNuclei = sum(SchnitzTbinnedGradedInactiveNuclei, 4);
    TbinnedCycleGradedNucleiCount = sum(SchnitzTbinnedGradedNucleiCount, 4);
    
    Tbinned3DCycleMeanTraces = NaN(numBinnedFrames, length(APbins), 6);
    Tbinned3DCycleTraceStdErrors =  NaN(numBinnedFrames, length(APbins), 6);
    Tbinned3DCycleTraceCount=  NaN(numBinnedFrames, length(APbins), 6);
    Tbinned3DCycleNumNuclei =  sum(SchnitzTbinnedTotalNuclei, 4);
    Tbinned3DCycleNumOnNuclei =  sum(SchnitzTbinnedOnNuclei, 4);
    Tbinned3DCycleFractionOn=  Tbinned3DCycleNumOnNuclei./Tbinned3DCycleNumNuclei;
    Tbinned3DCycleNumOffNuclei =  sum(SchnitzTbinnedOffNuclei, 4);
    Tbinned3DCycleNumQuiescentNuclei =  sum(SchnitzTbinnedQuiescentNuclei, 4);
    Tbinned3DCycleNumFinishedTranscribingNuclei =  sum(SchnitzTbinnedFinishedTranscribingNuclei, 4);
    Tbinned3DCycleActiveNuclei = sum(SchnitzTbinnedActiveNuclei, 3);
    Tbinned3DCycleInactiveNuclei = sum(SchnitzTbinnedInactiveNuclei, 3);
    Tbinned3DCycleNucleiCount = sum(SchnitzTbinnedNucleiCount, 3);
    Tbinned3DCycleGradedOnNuclei = sum(SchnitzTbinnedGradedOnNuclei, 5);
    Tbinned3DCycleGradedOffNuclei = sum(SchnitzTbinnedGradedOffNuclei, 5);
    Tbinned3DCycleGradedTotalNuclei = sum(SchnitzTbinnedGradedTotalNuclei, 5);
    Tbinned3DCycleGradedFractionOn = TbinnedCycleGradedOnNuclei./TbinnedCycleGradedTotalNuclei;
    Tbinned3DCycleGradedActiveNuclei = sum(SchnitzTbinnedGradedActiveNuclei, 4);
    Tbinned3DCycleGradedInactiveNuclei = sum(SchnitzTbinnedGradedInactiveNuclei, 4);
    Tbinned3DCycleGradedNucleiCount = sum(SchnitzTbinnedGradedNucleiCount, 4);
    TbinnedCycleTraces = cell(1, 6);
    Tbinned3DCycleTraces = cell(1, 6);
    TbinnedCycleFrameTimes = cell(1, 6);
    Tbinned3DCycleFrameTimes = cell(1, 6);
    %%
    ParticlesAreGood = ones(1, length(CompiledParticles{ChN}), 'logical');
    ParticleCycles = NaN(1, length(CompiledParticles{ChN}));
    ParticleSchnitzIndex = NaN(1, length(CompiledParticles{ChN}));
    ParticleHasAnaphaseFrame = zeros(1, length(CompiledParticles{ChN}), 'logical');
    for p_idx = 1:length(CompiledParticles{ChN})
        if ~isempty(CompiledParticles{ChN}(p_idx).Approved)
            if CompiledParticles{ChN}(p_idx).Approved < 1
                ParticlesAreGood(p_idx) = false;
            elseif CompiledParticles{ChN}(p_idx).Approved  == 2 & ~isempty(CompiledParticles{ChN}(p_idx).cycle)
                if sum(CompiledParticles{ChN}(p_idx).FlaggingInfo.SpotStateDefinitive | CompiledParticles{ChN}(p_idx).FlaggingInfo.UseTraceFluo)/MaxTraceLengths(CompiledParticles{ChN}(p_idx).cycle) < 0.5
                    ParticlesAreGood(p_idx) = false;
                end
            end
        end
        
        if ~isempty(CompiledParticles{ChN}(p_idx).cycle)
            ParticleCycles(p_idx) = CompiledParticles{ChN}(p_idx).cycle;
        else
            ParticlesAreGood(p_idx) = false;
        end
        if ~isempty(CompiledParticles{ChN}(p_idx).schnitz)
            if ~SchnitzesAreGood(CompiledParticles{ChN}(p_idx).schnitz)
                ParticlesAreGood(p_idx) = false;
            end
            ParticleSchnitzIndex(p_idx) = CompiledParticles{ChN}(p_idx).schnitz;
            CompiledParticles{ChN}(p_idx).schnitzcell = schnitzcells(ParticleSchnitzIndex(p_idx));
            %             if ~all(CompiledParticles{ChN}(p_idx).schnitzcell.FlaggingInfo.SchnitzAwayFromBoundary(CompiledParticles{ChN}(p_idx).schnitzcell.FlaggingInfo.AllSchnitzFrames >= EarliestTurnOnTimes(CompiledParticles{ChN}(p_idx).schnitzcell.cycle-8) ))
            %                 if sum(CompiledParticles{ChN}(p_idx).schnitzcell.FlaggingInfo.SchnitzAwayFromBoundary(CompiledParticles{ChN}(p_idx).schnitzcell.FlaggingInfo.AllSchnitzFrames >= EarliestTurnOnTimes(CompiledParticles{ChN}(p_idx).schnitzcell.cycle-8) ))/...
            %                         length(CompiledParticles{ChN}(p_idx).schnitzcell.FlaggingInfo.AllSchnitzFrames(CompiledParticles{ChN}(p_idx).schnitzcell.FlaggingInfo.AllSchnitzFrames >= EarliestTurnOnTimes(CompiledParticles{ChN}(p_idx).schnitzcell.cycle-8) )) < 0.8
            %                     ParticlesAreGood(p_idx) = false;
            %                 end
            %             end
            if ~isempty(CompiledParticles{ChN}(p_idx).schnitzcell.anaphaseFrame)
                if CompiledParticles{ChN}(p_idx).schnitzcell.inferredAnaphaseFrame | CompiledParticles{ChN}(p_idx).schnitzcell.anaphaseFrame == 0 | isnan(CompiledParticles{ChN}(p_idx).schnitzcell.anaphaseFrame)
                    ParticleHasAnaphaseFrame(p_idx) = false;
                elseif CompiledParticles{ChN}(p_idx).schnitzcell.anaphaseFrame >= nc_info(ParticleCycles(p_idx)-8)-2
                    ParticleHasAnaphaseFrame(p_idx) = true;
                end
            else
                ParticlesAreGood(p_idx)= false;
            end
            if CompiledParticles{ChN}(p_idx).schnitzcell.Approved ~= 1 | CompiledParticles{ChN}(p_idx).schnitzcell.Flag ~= 0
                ParticlesAreGood(p_idx)= false;
            end
            
        else
            ParticlesAreGood(p_idx) = false;
        end
        
    end
    %%
    
    
    IncludedNCs = min(ParticleCycles):max(ParticleCycles);
    for k = 1:length(IncludedNCs)
        NC = IncludedNCs(k);
        
        IncludedTraceIndices = find((ParticleCycles == NC) & ParticlesAreGood);
        
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
        NCParticleAnaphases = NaN(1, length(IncludedTraceIndices));
        for i = 1:length(IncludedTraceIndices)
            trace_index = IncludedTraceIndices(i);
            
            CurrentCompiledParticle = CompiledParticles{ChN}(trace_index);
            if ParticleHasAnaphaseFrame(trace_index)
                NCParticleAnaphases(i) = CurrentCompiledParticle.schnitzcell.anaphaseFrame;
            end
            NCMinTraceFrames(i) = min(CurrentCompiledParticle.FlaggingInfo.AllSchnitzFrames);
            CycleTraces{NC-8}(CurrentCompiledParticle.FlaggingInfo.AllSchnitzFrames, i) = 0;
            CycleTraces{NC-8}(CurrentCompiledParticle.FlaggingInfo.AllSchnitzFrames(CurrentCompiledParticle.FlaggingInfo.UseTraceFluo), i) = ...
                CurrentCompiledParticle.Fluo(CurrentCompiledParticle.FrameApproved);
            CycleTraces3D{NC-8}(CurrentCompiledParticle.FlaggingInfo.AllSchnitzFrames, i) = 0;
            CycleTraces3D{NC-8}(CurrentCompiledParticle.FlaggingInfo.AllSchnitzFrames(CurrentCompiledParticle.FlaggingInfo.UseTraceFluo), i) = ...
                CurrentCompiledParticle.Fluo3DGauss(CurrentCompiledParticle.FrameApproved);
            MeanAPs{NC-8}(i) = mean(CurrentCompiledParticle.schnitzcell.APpos);
            ParticleAPbinIndices{NC-8}(i) = find((MeanAPs{NC-8}(i) < APbins+APResolution/2) & (MeanAPs{NC-8}(i) >= APbins-APResolution/2));
            
        end
        SchnitzNCFrameMins = SchnitzAnaphaseFrames(SchnitzCycles == NC & abs(SchnitzAnaphaseFrames - anaphaseFrames(NC-8)) <= 2 & SchnitzAnaphaseFrames ~= 0);
        
        APbinsInMovie = min(ParticleAPbinIndices{NC-8}):max(ParticleAPbinIndices{NC-8});
        NCCycleTraces = CycleTraces{NC-8};
        NC3DCycleTraces = CycleTraces3D{NC-8};
        IncludedRowsInBin = find(sum(~isnan(NCCycleTraces),2).' > 0);
        if isempty(SchnitzNCFrameMins)
            SchnitzNCFrameMin = min(IncludedRowsInBin);
        else
            SchnitzNCFrameMin = min(SchnitzNCFrameMins);
        end
        if min(IncludedRowsInBin) > min([NCParticleAnaphases SchnitzNCFrameMins] )
            IncludedRowsInBin = [min([NCParticleAnaphases SchnitzNCFrameMins] ):(min(IncludedRowsInBin)-1) IncludedRowsInBin];
        end
        NCRealMinFrame = min(IncludedRowsInBin,[], 'omitnan');
        UnalignedCycleFrameTimes{NC-8}= FrameTimes(IncludedRowsInBin)-min(FrameTimes(IncludedRowsInBin));
        NCParticleAnaphases = NCParticleAnaphases-min(IncludedRowsInBin)+1;
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
            
            UnalignedCycleMeanTraces(IncludedRowsInBin, APbinIndex, NC-8) = mean(NCCycleTraces(IncludedRowsInBin, IncludedColumnsInBin), 2, 'omitnan');
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
            if ~isnan(NCParticleAnaphases(trace_index))
                CurrentTrace = CurrentNCTraces(:,trace_index);
                TraceIncludedRowsInBin = find(~isnan(CurrentTrace));
                TraceFrameTimes = NCFrameTimes(TraceIncludedRowsInBin)-NCFrameTimes(NCParticleAnaphases(trace_index));
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
    MeanProfiles.UnalignedCycleNumActiveNuclei = UnalignedCycleActiveNuclei;
    MeanProfiles.UnalignedCycleNumInactiveNuclei = UnalignedCycleInactiveNuclei;
    MeanProfiles.UnalignedCycleNucleiCount = UnalignedCycleNucleiCount;
    MeanProfiles.UnalignedCycleGradedOnNuclei = UnalignedCycleGradedOnNuclei ;
    MeanProfiles.UnalignedCycleGradedOffNuclei =  UnalignedCycleGradedOffNuclei ;
    MeanProfiles.UnalignedCycleGradedTotalNuclei = UnalignedCycleGradedTotalNuclei;
    MeanProfiles.UnalignedCycleGradedFractionOn = UnalignedCycleGradedFractionOn ;
    MeanProfiles.UnalignedCycleGradedActiveNuclei = UnalignedCycleGradedActiveNuclei ;
    MeanProfiles.UnalignedCycleGradedInactiveNuclei = UnalignedCycleGradedInactiveNuclei ;
    MeanProfiles.UnalignedCycleGradedNucleiCount =  UnalignedCycleGradedNucleiCount ;
    
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
    MeanProfiles.Unaligned3DCycleNumActiveNuclei = Unaligned3DCycleActiveNuclei;
    MeanProfiles.Unaligned3DCycleNumInactiveNuclei = Unaligned3DCycleInactiveNuclei;
    MeanProfiles.Unaligned3DCycleNucleiCount = Unaligned3DCycleNucleiCount;
    MeanProfiles.Unaligned3DCycleGradedOnNuclei = Unaligned3DCycleGradedOnNuclei ;
    MeanProfiles.Unaligned3DCycleGradedOffNuclei =  Unaligned3DCycleGradedOffNuclei ;
    MeanProfiles.Unaligned3DCycleGradedTotalNuclei = Unaligned3DCycleGradedTotalNuclei;
    MeanProfiles.Unaligned3DCycleGradedFractionOn = Unaligned3DCycleGradedFractionOn ;
    MeanProfiles.Unaligned3DCycleGradedActiveNuclei = Unaligned3DCycleGradedActiveNuclei ;
    MeanProfiles.Unaligned3DCycleGradedInactiveNuclei = Unaligned3DCycleGradedInactiveNuclei ;
    MeanProfiles.Unaligned3DCycleGradedNucleiCount =  Unaligned3DCycleGradedNucleiCount ;
    
    MeanProfiles.AnaphaseAlignedCycleMeanTraces = AnaphaseAlignedCycleMeanTraces;
    MeanProfiles.AnaphaseAlignedCycleNumNuclei = AnaphaseAlignedCycleNumNuclei;
    MeanProfiles.AnaphaseAlignedCycleNumOnNuclei = AnaphaseAlignedCycleNumOnNuclei;
    MeanProfiles.AnaphaseAlignedCycleFractionOn = AnaphaseAlignedCycleFractionOn;
    MeanProfiles.AnaphaseAlignedCycleNumOffNuclei = AnaphaseAlignedCycleNumOffNuclei;
    MeanProfiles.AnaphaseAlignedCycleNumQuiescentNuclei = AnaphaseAlignedCycleNumQuiescentNuclei;
    MeanProfiles.AnaphaseAlignedCycleNumFinishedTranscribingNuclei = AnaphaseAlignedCycleNumFinishedTranscribingNuclei;
    MeanProfiles.AnaphaseAlignedCycleTraceStdErrors = AnaphaseAlignedCycleTraceStdErrors;
    MeanProfiles.AnaphaseAlignedCycleTraceCount = AnaphaseAlignedCycleTraceCount;
    
    AnaphaseAlignedCycleMeanPerNucleusTraces = AnaphaseAlignedCycleMeanTraces.*AnaphaseAlignedCycleFractionOn;
    MeanProfiles.AnaphaseAlignedCycleMeanPerNucleusTraces = AnaphaseAlignedCycleMeanPerNucleusTraces;
    MeanProfiles.AnaphaseAlignedCycleNumActiveNuclei = AnaphaseAlignedCycleActiveNuclei;
    MeanProfiles.AnaphaseAlignedCycleNumInactiveNuclei = AnaphaseAlignedCycleInactiveNuclei;
    MeanProfiles.AnaphaseAlignedCycleNucleiCount = AnaphaseAlignedCycleNucleiCount;
    MeanProfiles.AnaphaseAlignedCycleGradedOnNuclei = AnaphaseAlignedCycleGradedOnNuclei ;
    MeanProfiles.AnaphaseAlignedCycleGradedOffNuclei =  AnaphaseAlignedCycleGradedOffNuclei ;
    MeanProfiles.AnaphaseAlignedCycleGradedTotalNuclei = AnaphaseAlignedCycleGradedTotalNuclei;
    MeanProfiles.AnaphaseAlignedCycleGradedFractionOn = AnaphaseAlignedCycleGradedFractionOn ;
    MeanProfiles.AnaphaseAlignedCycleGradedActiveNuclei = AnaphaseAlignedCycleGradedActiveNuclei ;
    MeanProfiles.AnaphaseAlignedCycleGradedInactiveNuclei = AnaphaseAlignedCycleGradedInactiveNuclei ;
    MeanProfiles.AnaphaseAlignedCycleGradedNucleiCount =  AnaphaseAlignedCycleGradedNucleiCount ;
    
    MeanProfiles.AnaphaseAligned3DCycleMeanTraces = AnaphaseAligned3DCycleMeanTraces;
    MeanProfiles.AnaphaseAligned3DCycleNumNuclei = AnaphaseAligned3DCycleNumNuclei;
    MeanProfiles.AnaphaseAligned3DCycleNumOnNuclei = AnaphaseAligned3DCycleNumOnNuclei;
    MeanProfiles.AnaphaseAligned3DCycleFractionOn = AnaphaseAligned3DCycleFractionOn;
    MeanProfiles.AnaphaseAligned3DCycleNumOffNuclei = AnaphaseAligned3DCycleNumOffNuclei;
    MeanProfiles.AnaphaseAligned3DCycleNumQuiescentNuclei = AnaphaseAligned3DCycleNumQuiescentNuclei;
    MeanProfiles.AnaphaseAligned3DCycleNumFinishedTranscribingNuclei = AnaphaseAligned3DCycleNumFinishedTranscribingNuclei;
    MeanProfiles.AnaphaseAligned3DCycleTraceStdErrors = AnaphaseAligned3DCycleTraceStdErrors;
    MeanProfiles.AnaphaseAligned3DCycleTraceCount = AnaphaseAligned3DCycleTraceCount;
    
    AnaphaseAligned3DCycleMeanPerNucleusTraces = AnaphaseAligned3DCycleMeanTraces.*AnaphaseAligned3DCycleFractionOn;
    MeanProfiles.AnaphaseAligned3DCycleMeanPerNucleusTraces = AnaphaseAligned3DCycleMeanPerNucleusTraces;
    MeanProfiles.AnaphaseAligned3DCycleNumActiveNuclei = AnaphaseAligned3DCycleActiveNuclei;
    MeanProfiles.AnaphaseAligned3DCycleNumInactiveNuclei = AnaphaseAligned3DCycleInactiveNuclei;
    MeanProfiles.AnaphaseAligned3DCycleNucleiCount = AnaphaseAligned3DCycleNucleiCount;
    MeanProfiles.AnaphaseAligned3DCycleGradedOnNuclei = AnaphaseAligned3DCycleGradedOnNuclei ;
    MeanProfiles.AnaphaseAligned3DCycleGradedOffNuclei =  AnaphaseAligned3DCycleGradedOffNuclei ;
    MeanProfiles.AnaphaseAligned3DCycleGradedTotalNuclei = AnaphaseAligned3DCycleGradedTotalNuclei;
    MeanProfiles.AnaphaseAligned3DCycleGradedFractionOn = AnaphaseAligned3DCycleGradedFractionOn ;
    MeanProfiles.AnaphaseAligned3DCycleGradedActiveNuclei = AnaphaseAligned3DCycleGradedActiveNuclei ;
    MeanProfiles.AnaphaseAligned3DCycleGradedInactiveNuclei = AnaphaseAligned3DCycleGradedInactiveNuclei ;
    MeanProfiles.AnaphaseAligned3DCycleGradedNucleiCount =  AnaphaseAligned3DCycleGradedNucleiCount ;
    
    
    MeanProfiles.TbinnedCycleMeanTraces = TbinnedCycleMeanTraces;
    MeanProfiles.TbinnedCycleNumNuclei = TbinnedCycleNumNuclei;
    MeanProfiles.TbinnedCycleNumOnNuclei = TbinnedCycleNumOnNuclei;
    MeanProfiles.TbinnedCycleFractionOn = TbinnedCycleFractionOn;
    MeanProfiles.TbinnedCycleNumOffNuclei = TbinnedCycleNumOffNuclei;
    MeanProfiles.TbinnedCycleNumQuiescentNuclei = TbinnedCycleNumQuiescentNuclei;
    MeanProfiles.TbinnedCycleNumFinishedTranscribingNuclei = TbinnedCycleNumFinishedTranscribingNuclei;
    MeanProfiles.TbinnedCycleTraceStdErrors = TbinnedCycleTraceStdErrors;
    MeanProfiles.TbinnedCycleTraceCount = TbinnedCycleTraceCount;
    
    TbinnedCycleMeanPerNucleusTraces = TbinnedCycleMeanTraces.*TbinnedCycleFractionOn;
    MeanProfiles.TbinnedCycleMeanPerNucleusTraces = TbinnedCycleMeanPerNucleusTraces;
    MeanProfiles.TbinnedCycleNumActiveNuclei = TbinnedCycleActiveNuclei;
    MeanProfiles.TbinnedCycleNumInactiveNuclei =TbinnedCycleInactiveNuclei;
    MeanProfiles.TbinnedCycleNucleiCount = TbinnedCycleNucleiCount;
    MeanProfiles.TbinnedCycleGradedOnNuclei = TbinnedCycleGradedOnNuclei ;
    MeanProfiles.TbinnedCycleGradedOffNuclei =  TbinnedCycleGradedOffNuclei ;
    MeanProfiles.TbinnedCycleGradedTotalNuclei = TbinnedCycleGradedTotalNuclei;
    MeanProfiles.TbinnedCycleGradedFractionOn = TbinnedCycleGradedFractionOn ;
    MeanProfiles.TbinnedCycleGradedActiveNuclei = TbinnedCycleGradedActiveNuclei ;
    MeanProfiles.TbinnedCycleGradedInactiveNuclei = TbinnedCycleGradedInactiveNuclei ;
    MeanProfiles.TbinnedCycleGradedNucleiCount =  TbinnedCycleGradedNucleiCount ;
    
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
    Tbinned3DCycleMeanPerNucleusTraces = Tbinned3DCycleMeanTraces.*Tbinned3DCycleFractionOn;
    MeanProfiles.Tbinned3DCycleMeanPerNucleusTraces = Tbinned3DCycleMeanPerNucleusTraces;
    MeanProfiles.Tbinned3DCycleNumActiveNuclei = Tbinned3DCycleActiveNuclei;
    MeanProfiles.Tbinned3DCycleNumInactiveNuclei =Tbinned3DCycleInactiveNuclei;
    MeanProfiles.Tbinned3DCycleNucleiCount = Tbinned3DCycleNucleiCount;
    MeanProfiles.Tbinned3DCycleGradedOnNuclei = Tbinned3DCycleGradedOnNuclei ;
    MeanProfiles.Tbinned3DCycleGradedOffNuclei =  Tbinned3DCycleGradedOffNuclei ;
    MeanProfiles.Tbinned3DCycleGradedTotalNuclei = Tbinned3DCycleGradedTotalNuclei;
    MeanProfiles.Tbinned3DCycleGradedFractionOn = Tbinned3DCycleGradedFractionOn ;
    MeanProfiles.Tbinned3DCycleGradedActiveNuclei = Tbinned3DCycleGradedActiveNuclei ;
    MeanProfiles.Tbinned3DCycleGradedInactiveNuclei = Tbinned3DCycleGradedInactiveNuclei ;
    MeanProfiles.Tbinned3DCycleGradedNucleiCount =  Tbinned3DCycleGradedNucleiCount ;
    
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













end



