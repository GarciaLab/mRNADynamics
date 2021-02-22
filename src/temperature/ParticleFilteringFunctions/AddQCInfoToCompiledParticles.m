function [CompiledParticles] = AddQCInfoToCompiledParticles(CompiledParticles, schnitzcells, Prefix)

%% Add histone fluorescence information for each schnitz cell z-stack
liveExperiment = LiveExperiment(Prefix);
FrameInfo = getFrameInfo(liveExperiment);
FrameTimes = [FrameInfo(:).Time];
nc_info = [liveExperiment.nc9, liveExperiment.nc10, liveExperiment.nc11,...
    liveExperiment.nc12, liveExperiment.nc13, liveExperiment.nc14, length(FrameInfo)];
PixelSize = liveExperiment.pixelSize_um;
SnippetSize = liveExperiment.snippetSize_px;
nucleusDiameters = zeros(1, 6);
for nc=9:14
    nucleusDiameters(nc-8) = getDefaultParameters(FrameInfo,[['d', num2str(nc)]])/PixelSize; % in pixels
end
xDim = liveExperiment.xDim;
yDim = liveExperiment.yDim;
zDim = liveExperiment.zDim;
schnitzcells = integrateHistoneFluo(Prefix, schnitzcells, FrameInfo);

%  Add spot fluorescence information for each particle z-stack in CompiledParticles


[CompiledParticles] = GetFluoZInfo(liveExperiment, CompiledParticles);

for ChN=1:length(liveExperiment.spotChannels)
    for i =1:length(CompiledParticles{ChN})
        CompiledParticles{ChN}(i).schnitzcell = schnitzcells(CompiledParticles{ChN}(i).schnitz);
        CompiledParticles{ChN}(i).cycle = CompiledParticles{ChN}(i).schnitzcell.cycle;
        CompiledParticles{ChN}(i).FlaggingInfo = {};
    end
    %%
    
    his_zpadding = 1;
    spot_zpadding = 2;
    for i=1:length(CompiledParticles{ChN})%[test_idx]
        p = CompiledParticles{ChN}(i);
        sc = p.schnitzcell;
        if sc.Flag == 6
            p.FlaggingInfo.SickNucleus = true;
        else
            p.FlaggingInfo.SickNucleus = false;
        end
        
        if sc.Approved == 1
            p.FlaggingInfo.ApprovedNucleus = true;
        else
            p.FlaggingInfo.ApprovedNucleus = false;
        end
        % Find True Nuclear start frame by inference or direct observation
        if ~isempty(sc.inferredAnaphaseFrame)
            if sc.inferredAnaphaseFrame == 0
                if ~isempty(sc.anaphaseFrame)
                    if sc.anaphaseFrame >= nc_info(sc.cycle-8)-1
                        p.FlaggingInfo.TrueNucStart = max([sc.anaphaseFrame, 1]);
                        p.FlaggingInfo.TNSinferred = false;
                    else
                        p.FlaggingInfo.TrueNucStart =  max([nc_info(sc.cycle-8), 1]);
                        p.FlaggingInfo.TNSinferred = true;
                    end
                else
                    p.FlaggingInfo.TrueNucStart =  max([1, nc_info(sc.cycle-8)]);
                    p.FlaggingInfo.TNSinferred = true;
                end
            elseif sc.anaphaseFrame >= nc_info(sc.cycle-8)-1
                p.FlaggingInfo.TrueNucStart = max([1, sc.anaphaseFrame]);
                p.FlaggingInfo.TNSinferred = true;
            else
                p.FlaggingInfo.TrueNucStart =  max([1,nc_info(sc.cycle-8)]);
                p.FlaggingInfo.TNSinferred = true;
            end
        elseif ~isempty(sc.anaphaseFrame)
            if sc.anaphaseFrame >= nc_info(sc.cycle-8)-1
                p.FlaggingInfo.TrueNucStart = max([1, sc.anaphaseFrame]);
                p.FlaggingInfo.TNSinferred = false;
            else
                p.FlaggingInfo.TrueNucStart =  max([1, nc_info(sc.cycle-8)]);
                p.FlaggingInfo.TNSinferred = true;
            end
        else
            p.FlaggingInfo.TrueNucStart =  max([1, nc_info(sc.cycle-8)]);
            p.FlaggingInfo.TNSinferred = true;
        end
        % Find True Nuclear end frame by inference or direct observation
        if p.NucEnd >= nc_info(sc.cycle-7)-1
            p.FlaggingInfo.TrueNucEnd = max([1, p.NucEnd]);
            p.FlaggingInfo.TNEinferred =false;
        else
            p.FlaggingInfo.TrueNucEnd = max([1, nc_info(sc.cycle-7)-1]);
            p.FlaggingInfo.TNEinferred =true;
        end
        
        p.FlaggingInfo.TrueFrames =...
            p.FlaggingInfo.TrueNucStart:p.FlaggingInfo.TrueNucEnd;
        
        
        
        % Make Histone and spot Fluorescence traces and z positions
        NFrames = length(p.FlaggingInfo.TrueFrames);
        MaxSpotFluoLevel = NaN(1, NFrames);
        MaxSpotZPos= NaN(1, NFrames);
        SpotZTraceApproved = zeros(1, NFrames);
        DetectedSpot = zeros(1, NFrames);
        MaxHisFluoLevel = NaN(1, NFrames);
        MaxHisZPos= NaN(1, NFrames);
        HisTraceApproved = zeros(1, NFrames);
        DetectedNucleus = zeros(1, NFrames);
        for j = 1:length(sc.frames)
            schnitz_index = find(p.FlaggingInfo.TrueFrames == sc.frames(j), 1);
            DetectedNucleus(schnitz_index) = 1;
            HisVector = sc.HistoneFluo(j,2:zDim+1);
            MaxF = max(HisVector);
            if ~isempty(MaxF)
                MaxHisFluoLevel(schnitz_index) = MaxF;
                MaxZ = find(HisVector(his_zpadding:zDim-his_zpadding) == MaxF, 1);
                if ~isempty(MaxZ)
                    MaxHisZPos(schnitz_index) = MaxZ+his_zpadding;
                    HisTraceApproved(schnitz_index) = 1;
                else
                    MaxZ = find(HisVector == MaxF, 1);
                    if ~isempty(MaxZ)
                        MaxHisZPos(schnitz_index) = MaxZ;
                        HisTraceApproved(schnitz_index) = -1;
                    end
                end
                
            end
            
        end
        p.FlaggingInfo.MaxHisFluoLevel = MaxHisFluoLevel;
        p.FlaggingInfo.MaxHisZPos = MaxHisZPos;
        p.FlaggingInfo.HisTraceApproved = HisTraceApproved;
        
        for j = 1:length(p.Frame)
            particle_index = find(p.FlaggingInfo.TrueFrames == p.Frame(j), 1);
            if isempty(find(sc.frames == p.Frame(j), 1))
                continue
            elseif p.FrameApproved(j) == 0
                continue
            end
            DetectedSpot(particle_index) = 1;
            SpotZSliceVector = p.FluoZInfo(j,:);
            MaxF = max(SpotZSliceVector);
            if ~isempty(MaxF)
                MaxSpotFluoLevel(particle_index) = MaxF;
                MaxZ = find(SpotZSliceVector(spot_zpadding:zDim-spot_zpadding) == MaxF, 1);
                if ~isempty(MaxZ)
                    MaxSpotZPos(particle_index) = MaxZ+spot_zpadding;
                    SpotZTraceApproved(particle_index) = 1;
                else
                    MaxZ = find(SpotZSliceVector == MaxF, 1);
                    if ~isempty(MaxZ)
                        MaxSpotZPos(particle_index) = MaxZ;
                        SpotZTraceApproved(particle_index) = -1;
                    end
                end
                
            end
            
        end
        p.FlaggingInfo.MaxSpotFluoLevel = MaxSpotFluoLevel;
        p.FlaggingInfo.MaxSpotZPos = MaxSpotZPos;
        p.FlaggingInfo.SpotZTraceApproved = SpotZTraceApproved;
        
        p.FlaggingInfo.DetectedNucleus = DetectedNucleus;
        p.FlaggingInfo.DetectedSpot = DetectedSpot;
        % Add spot position deviation velocity flags
        
        
        
        SpotMedianXPos = median(p.xPos);
        SpotMedianYPos = median(p.yPos);
        
        SpotXPosDeviations = NaN(1, NFrames);
        SpotYPosDeviations = NaN(1, NFrames);
        SpotPosDeviations  = NaN(1, NFrames);
        SpotXVelocities = NaN(1, NFrames);
        SpotYVelocities = NaN(1, NFrames);
        SpotVelocities =  NaN(1, NFrames);
        SpotVeloCalcFrameStep = NaN(1,NFrames);
        
        SpotRelativeXPos= NaN(1, NFrames);
        SpotRelativeYPos = NaN(1, NFrames);
        for j = 1:length(p.Frame)
            if isempty(find(sc.frames == p.Frame(j), 1))
                continue
            elseif p.FrameApproved(j) == 0
                continue
            end
            CurrentFrame = p.Frame(j);
            idx = find(p.FlaggingInfo.TrueFrames == CurrentFrame, 1);
            schnitz_idx = find(sc.frames == CurrentFrame, 1);
            if ~isempty(schnitz_idx)
                SpotRelativeXPos(idx) = (p.xPos(j)-sc.cenx(schnitz_idx));
                SpotRelativeYPos(idx) = (p.yPos(j)-sc.ceny(schnitz_idx));
            end
        end
        
        SpotMedianRelativeXPos = nanmedian(SpotRelativeXPos);
        SpotMedianRelativeYPos = nanmedian(SpotRelativeYPos);
        SpotRelativeXPosDeviations = NaN(1, NFrames);
        SpotRelativeYPosDeviations = NaN(1, NFrames);
        SpotRelativePosDeviations  = NaN(1, NFrames);
        SpotRelativeXVelocities = NaN(1, NFrames);
        SpotRelativeYVelocities = NaN(1, NFrames);
        SpotRelativeVelocities =  NaN(1, NFrames);
        for j = 1:length(p.Frame)
            if isempty(find(sc.frames == p.Frame(j), 1))
                continue
            elseif p.FrameApproved(j) == 0
                continue
            end
            CurrentFrame = p.Frame(j);
            idx = find(p.FlaggingInfo.TrueFrames == CurrentFrame, 1);
            if isempty(idx)
                p.FrameApproved(j) = 0;
                continue
            end
            schnitz_idx = find(sc.frames == CurrentFrame, 1);
            SpotXPosDeviations(idx) = abs(p.xPos(j)-SpotMedianXPos);
            SpotXPosDeviations(idx) = abs(p.yPos(j)-SpotMedianYPos);
            SpotPosDeviations(idx) = sqrt((p.xPos(j)-SpotMedianXPos)^2 + (p.yPos(j)-SpotMedianYPos)^2);
            if ~isnan(SpotRelativeXPos(idx)) & ~isnan(SpotRelativeYPos(idx))
                SpotRelativeXPosDeviations(idx) = abs(SpotRelativeXPos(idx)- SpotMedianRelativeXPos);
                SpotRelativeYPosDeviations(idx) = abs(SpotRelativeYPos(idx)- SpotMedianRelativeYPos);
                SpotRelativePosDeviations(idx) = sqrt(SpotRelativeXPosDeviations(idx)^2 + SpotRelativeYPosDeviations(idx)^2);
            end
            if j > 1
                prev_j = find(p.FrameApproved(1:j-1) == 1, 1, 'last');
                if isempty(prev_j)
                    continue
                end
                PreviousFrame = p.Frame(prev_j);
                Tdelta = FrameTimes(CurrentFrame)-FrameTimes(PreviousFrame);
                Xdelta = (p.xPos(j)-p.xPos(prev_j)); % in microns
                Ydelta = (p.yPos(j)-p.yPos(prev_j)); % in microns
                SpotXVelocities(idx) = Xdelta/Tdelta;
                SpotYVelocities(idx) = Ydelta/Tdelta;
                SpotVelocities(idx) = sqrt(SpotXVelocities(idx)^2 + SpotYVelocities(idx)^2);
                SpotVeloCalcFrameStep(idx) = CurrentFrame-PreviousFrame;
                idx2 = find(p.FlaggingInfo.TrueFrames == PreviousFrame, 1);
                if ~isnan(SpotRelativeXPos(idx)) & ~isnan(SpotRelativeYPos(idx)) & ...
                        ~isnan(SpotRelativeXPos(idx2)) & ~isnan(SpotRelativeYPos(idx2))
                    Xdelta = SpotRelativeXPosDeviations(idx)-SpotRelativeXPosDeviations(idx2);
                    Ydelta = SpotRelativeYPosDeviations(idx)-SpotRelativeYPosDeviations(idx2);
                    SpotRelativeXVelocities(idx) =  Xdelta/Tdelta;
                    SpotRelativeYVelocities(idx) = Ydelta/Tdelta;
                    SpotRelativeVelocities(idx) = sqrt(SpotRelativeXVelocities(idx)^2 + SpotRelativeYVelocities(idx)^2);
                end
            end
            
            
            
        end
        
        
        p.FlaggingInfo.SpotMedianXPos = SpotMedianXPos*PixelSize;
        p.FlaggingInfo.SpotMedianYPos = SpotMedianYPos*PixelSize;
        p.FlaggingInfo.SpotXPosDeviations = SpotXPosDeviations*PixelSize;
        p.FlaggingInfo.SpotYPosDeviations = SpotYPosDeviations*PixelSize;
        p.FlaggingInfo.SpotPosDeviations = SpotPosDeviations*PixelSize;
        p.FlaggingInfo.SpotXVelocities = SpotXVelocities*PixelSize;
        p.FlaggingInfo.SpotYVelocities = SpotYVelocities*PixelSize;
        p.FlaggingInfo.SpotVelocities = SpotVelocities*PixelSize;
        p.FlaggingInfo.SpotVeloCalcFrameStep = SpotVeloCalcFrameStep;
        
        p.FlaggingInfo.SpotMedianRelativeXPos = SpotMedianRelativeXPos*PixelSize;
        p.FlaggingInfo.SpotMedianRelativeYPos = SpotMedianRelativeYPos*PixelSize;
        p.FlaggingInfo.SpotRelativeXPosDeviations = SpotRelativeXPosDeviations*PixelSize;
        p.FlaggingInfo.SpotRelativeYPosDeviations = SpotRelativeYPosDeviations*PixelSize;
        p.FlaggingInfo.SpotRelativePosDeviations  = SpotRelativePosDeviations*PixelSize;
        p.FlaggingInfo.SpotRelativeXVelocities = SpotRelativeXVelocities*PixelSize;
        p.FlaggingInfo.SpotRelativeYVelocities = SpotRelativeYVelocities*PixelSize;
        p.FlaggingInfo.SpotRelativeVelocities = SpotRelativeVelocities*PixelSize;
        
        % add nuclear position info
        
        SpotXPos = NaN(1, NFrames);
        SpotYPos = NaN(1, NFrames);
        schnitzXPos = NaN(1, NFrames);
        schnitzYPos = NaN(1, NFrames);
        schnitzXVelo = NaN(1, NFrames);
        schnitzYVelo = NaN(1, NFrames);
        for j=1:NFrames
            schnitz_idx = find(sc.frames == p.FlaggingInfo.TrueFrames(j), 1);
            particle_idx = find(p.Frame == p.FlaggingInfo.TrueFrames(j), 1);
            CurrentFrame = p.FlaggingInfo.TrueFrames(j);
            if ~isempty(schnitz_idx)
                schnitzXPos(j) = sc.cenx(schnitz_idx);
                schnitzYPos(j) = sc.ceny(schnitz_idx);
                if j > 1
                    schnitzXVelo(j) = (schnitzXPos(j) -schnitzXPos(j-1) )/(FrameTimes(CurrentFrame)-FrameTimes(CurrentFrame-1));
                    schnitzYVelo(j) = (schnitzYPos(j) -schnitzYPos(j-1) )/(FrameTimes(CurrentFrame)-FrameTimes(CurrentFrame-1));
                end
            end
            if ~isempty(particle_idx)
                SpotXPos(j) = p.xPos(particle_idx);
                SpotYPos(j) = p.yPos(particle_idx);
            end
        end
        
        p.FlaggingInfo.SpotXPos = SpotXPos;
        p.FlaggingInfo.SpotYPos = SpotYPos;
        p.FlaggingInfo.schnitzXPos = schnitzXPos;
        p.FlaggingInfo.schnitzYPos = schnitzYPos;
        p.FlaggingInfo.schnitzXVelo = schnitzXVelo;
        p.FlaggingInfo.schnitzYVelo = schnitzYVelo;
        
        CompiledParticles{ChN}(i) = p;
        % discontinutity error
    end
    %%
    
    AllRelativeVelos = [];
    AllVelos = [];
    AllFrames = [];
    AllNCs = [];
    AllSpotPosDeviations = [];
    AllSpotRelativeDeviations = [];
    AllSpotFluoDeltas = [];
    AllNormedSpotFluoDeltas = [];
    AllHistoneZs = [];
    AllSpotZs = [];
    for i = 1:length(CompiledParticles{1})
        AllHistoneZs = [AllHistoneZs,  CompiledParticles{ChN}(i).FlaggingInfo.MaxHisZPos];
        AllSpotZs = [AllSpotZs,  CompiledParticles{ChN}(i).FlaggingInfo.MaxSpotZPos];
        AllRelativeVelos = [AllRelativeVelos, CompiledParticles{ChN}(i).FlaggingInfo.SpotRelativeVelocities];
        AllVelos = [AllVelos, CompiledParticles{ChN}(i).FlaggingInfo.SpotVelocities];
        AllFrames = [AllFrames, CompiledParticles{ChN}(i).FlaggingInfo.TrueFrames];
        AllNCs = [AllNCs, repmat(CompiledParticles{ChN}(i).cycle, 1, length(CompiledParticles{ChN}(i).FlaggingInfo.TrueFrames))];
        AllSpotPosDeviations = [AllSpotPosDeviations, CompiledParticles{ChN}(i).FlaggingInfo.SpotPosDeviations];
        AllSpotRelativeDeviations = [AllSpotRelativeDeviations, CompiledParticles{ChN}(i).FlaggingInfo.SpotRelativePosDeviations];
        NFrames = length(CompiledParticles{ChN}(i).FlaggingInfo.TrueFrames);
        FluoDeltas = NaN(1, NFrames);
        NormedFluoDeltas = NaN(1, NFrames);
        for j = 2:NFrames
            if (sum(CompiledParticles{ChN}(i).FlaggingInfo.SpotZTraceApproved(j-1:j)) == 2)
                FluoDeltas(j) =  (CompiledParticles{ChN}(i).FlaggingInfo.MaxSpotFluoLevel(j) - ...
                    CompiledParticles{ChN}(i).FlaggingInfo.MaxSpotFluoLevel(j-1));
                NormedFluoDeltas(j) = FluoDeltas(j)/max(CompiledParticles{ChN}(i).FlaggingInfo.MaxSpotFluoLevel);
            end
        end
        CompiledParticles{ChN}(i).FlaggingInfo.FluoDeltas = FluoDeltas;
        CompiledParticles{ChN}(i).FlaggingInfo.NormedFluoDeltas = NormedFluoDeltas;
        AllSpotFluoDeltas = [AllSpotFluoDeltas, FluoDeltas];
        AllNormedSpotFluoDeltas = [AllNormedSpotFluoDeltas, NormedFluoDeltas];
    end
    
    probs = 0:0.01:1;
    NormedSpotFluoQuantiles = quantile(AllNormedSpotFluoDeltas, probs);
    VeloQuantiles = quantile(AllVelos, probs);
    RelativeVeloQuantiles = quantile(AllRelativeVelos, probs);
    SpotPosDevQuantiles = quantile(AllSpotPosDeviations, probs);
    SpotRelPosDevQuantiles = quantile(AllSpotRelativeDeviations, probs);
    %%
    
    for i=1:length(CompiledParticles{ChN})%
        p = CompiledParticles{ChN}(i);
        sc = p.schnitzcell;
        FI = p.FlaggingInfo;
        NFrames = length(FI.TrueFrames);
        FrameApprovedByFluoValue = zeros(1, NFrames);
        FrameApprovedByZPos = zeros(1, NFrames);
        FrameApprovedByPosition = zeros(1, NFrames);
        FrameApprovedByTiming = zeros(1, NFrames);
        
        if FI.SickNucleus
            p.Approved = false;
            CompiledParticles{ChN}(i) = p;
            continue
        end
        % Checking Fluo Values
        MaxFluo = max(FI.MaxSpotFluoLevel);
        for j=1:NFrames
            if j == 1 | j == NFrames
                FrameApprovedByFluoValue(j) = 1;
                continue
            elseif isnan(FI.MaxSpotFluoLevel(j))
                FrameApprovedByFluoValue(j) = 1;
                continue
            end
            
            if ~isnan(FI.MaxSpotFluoLevel(j-1)) & ~isnan(FI.MaxSpotFluoLevel(j+1))
                if (FI.MaxSpotFluoLevel(j)-FI.MaxSpotFluoLevel(j-1))/MaxFluo < NormedSpotFluoQuantiles(6)
                    if (FI.MaxSpotFluoLevel(j+1) - FI.MaxSpotFluoLevel(j))/MaxFluo > NormedSpotFluoQuantiles(96)
                        FrameApprovedByFluoValue(j) = -1;
                    else
                        FrameApprovedByFluoValue(j) = 1;
                    end
                else
                    FrameApprovedByFluoValue(j) = 1;
                end
            else
                FrameApprovedByFluoValue(j) = 1;
            end
            
        end
        % Checking Z Position Info
        for j=1:NFrames
            if FI.DetectedNucleus(j) ~= 1
                continue
            end
            CurrentFrame = FI.TrueFrames(j);
            if FI.DetectedSpot(j) == 1
                if ~isnan(FI.MaxSpotFluoLevel(j)) & (FI.SpotZTraceApproved(j) == 1)
                    FrameApprovedByZPos(j) = 1;
                else
                    FrameApprovedByZPos(j) = -1;
                end
            elseif FI.HisTraceApproved(j)== 0
                FrameApprovedByZPos(j) = -1;
            else
                if j > 1
                    PreviousDetection = find(FI.DetectedSpot(1:j-1), 1, 'last');
                else
                    PreviousDetection = [];
                end
                if j < NFrames
                    NextDetection = find(FI.DetectedSpot(j+1:end), 1);
                    if ~isempty(NextDetection)
                        NextDetection = NextDetection+j;
                    end
                else
                    NextDetection= [];
                end
                if ~isempty(PreviousDetection)
                    if ~isempty(NextDetection)
                        if FI.SpotZTraceApproved(PreviousDetection) == 1 & ...
                                FI.SpotZTraceApproved(NextDetection)== 1
                            FrameApprovedByZPos(j) = 1;
                        else
                            FrameApprovedByZPos(j) = -1;
                        end
                    elseif FI.SpotZTraceApproved(PreviousDetection) == 1
                        FrameApprovedByZPos(j) = 1;
                    else
                        FrameApprovedByZPos(j) = -1;
                    end
                elseif ~isempty(NextDetection)
                    if FI.SpotZTraceApproved(NextDetection) == 1
                        FrameApprovedByZPos(j) = 1;
                    else
                        FrameApprovedByZPos(j) = -1;
                    end
                end
                
                
                
            end
            
            
            
        end
        
        % Checking Position Info
        for j=1:NFrames
            CurrentFrame = FI.TrueFrames(j);
            if FI.DetectedNucleus(j) ~= 1
                FrameApprovedByPosition(j) = -1;
                continue
            elseif FI.DetectedSpot(j) == 1
                FrameApprovedByPosition(j) = 1;
                continue
            end
            
            xMin_px = FI.schnitzXPos(j) - nucleusDiameters(p.cycle-8)*.75-SnippetSize; % in pixels
            xMax_px = FI.schnitzXPos(j) + nucleusDiameters(p.cycle-8)*.75+SnippetSize; % in pixels
            yMin_px = FI.schnitzYPos(j) - nucleusDiameters(p.cycle-8)*.75-SnippetSize; % in pixels
            yMax_px = FI.schnitzYPos(j) + nucleusDiameters(p.cycle-8)*.75+SnippetSize; % in pixels
            
            if (xMin_px >= 0) & (xMax_px <= xDim) & (yMax_px <= yDim) & (yMin_px > 0)
                FrameApprovedByPosition(j) = 1;
                continue
            end
            if j > 1
                PreviousDetection = find(FI.DetectedSpot(1:j-1), 1, 'last');
            else
                PreviousDetection = [];
            end
            if j < NFrames
                NextDetection = find(FI.DetectedSpot(j+1:end), 1);
                if ~isempty(NextDetection)
                    NextDetection = NextDetection+j;
                end
            else
                NextDetection= [];
            end
            
            if isempty(PreviousDetection)
                FrameApprovedByPosition(j) = -1;
                continue
            end
            
            if isempty(NextDetection)
                PreviousSchnitzXPos = FI.schnitzXPos(PreviousDetection);
                PreviousSchnitzYPos = FI.schnitzYPos(PreviousDetection);
                PreviousSpotXPos = FI.SpotXPos(PreviousDetection);
                PreviousSpotYPos = FI.SpotYPos(PreviousDetection);
                DeltaSchnitzXPos = FI.schnitzXPos(j)-PreviousSchnitzXPos;
                DeltaSchnitzYPos = FI.schnitzYPos(j)-PreviousSchnitzYPos;
                DeltaT = FrameTimes(FI.TrueFrames(j))-FrameTimes(FI.TrueFrames(PreviousDetection));
                RelativeDeltaSpot = (RelativeVeloQuantiles(96)/PixelSize)*DeltaT;
                xMin_px = PreviousSpotXPos + DeltaSchnitzXPos - RelativeDeltaSpot;
                xMax_px = PreviousSpotXPos + DeltaSchnitzXPos + RelativeDeltaSpot;
                yMin_px = PreviousSpotYPos + DeltaSchnitzYPos - RelativeDeltaSpot;
                yMax_px = PreviousSpotYPos + DeltaSchnitzYPos + RelativeDeltaSpot;
                if (xMin_px >= 0) & (xMax_px <= xDim) & (yMax_px <= yDim) & (yMin_px > 0)
                    FrameApprovedByPosition(j) = 1;
                else
                    FrameApprovedByPosition(j) = -1;
                end
            else
                FrameApprovedByPosition(j) = 1;
                %             PreviousSchnitzXPos = FI.schnitzXPos(PreviousDetection);
                %             PreviousSchnitzYPos = FI.schnitzYPos(PreviousDetection);
                %             PreviousSpotXPos = FI.SpotXPos(PreviousDetection);
                %             PreviousSpotYPos = FI.SpotYPos(PreviousDetection);
                %             NextSchnitzXPos = FI.schnitzXPos(NextDetection);
                %             NextSchnitzYPos = FI.schnitzYPos(NextDetection);
                %             NextSpotXPos = FI.SpotXPos(NextDetection);
                %             NextSpotYPos = FI.SpotYPos(NextDetection);
                
            end
        end
        SpotIsOff = 0;
        SchnitzIsInFrame = 5;
        IncludeFrames = true;
        for j = 1:NFrames
            if FI.DetectedNucleus(j) == 1
                SchnitzIsInFrame = SchnitzIsInFrame + 1;
                if FI.DetectedSpot(j) == 1
                    if (SpotIsOff >= 5) & (SchnitzIsInFrame >= 5)
                        FrameApprovedByTiming(j) = 1;
                    else
                        FrameApprovedByTiming(j) = -1;
                    end
                else
                    SpotIsOff = SpotIsOff +1;
                    FrameApprovedByTiming(j) = 1;
                end
            else
                SchnitzIsInFrame = 0;
                SpotIsOff = 0;
                FrameApprovedByTiming(j) = 1;
            end
            
        end
        FI.FrameApprovedByFluoValue = FrameApprovedByFluoValue;
        FI.FrameApprovedByZPos = FrameApprovedByZPos;
        FI.FrameApprovedByPosition = FrameApprovedByPosition;
        FI.FrameApprovedByTiming = FrameApprovedByTiming;
        
        
        
        FrameApprovedFinal = zeros(1, NFrames);
        for j =1:NFrames
            CurrentFrame = FI.TrueFrames(j);
            particle_idx = find(p.Frame == CurrentFrame, 1);
            if ~isempty(particle_idx)
                FrameApprovedCPT = p.FrameApproved(particle_idx);
            else
                FrameApprovedCPT = 1;
            end
            if (FrameApprovedCPT == 1) & (FrameApprovedByFluoValue(j) ~= -1) & ...
                    (FrameApprovedByZPos(j) ~= -1) & (FrameApprovedByPosition(j) ~= -1) & ...
                    (FrameApprovedByTiming(j) ~= -1)
                FrameApprovedFinal(j) = 1;
            end
        end
        
        
        
        FrameStates = zeros(1, NFrames); % 1 means on, -1 means off, and 0 is inconclusive
        for j = 1:NFrames
            if (FrameApprovedFinal(j) == 1) & (FI.DetectedSpot(j) == 0)
                FrameStates(j) = -1;
            elseif (FrameApprovedFinal(j) == 1) & (FI.DetectedSpot(j) == 1)
                FrameStates(j) = 1;
            end
        end
        
        OnFrames = find(FrameStates == 1);
        if ~isempty(OnFrames)
            if length(OnFrames) > 1
                DiffOnFrames = diff(OnFrames);
                IsolatedOffFrames = find(DiffOnFrames == 2);
                for k=IsolatedOffFrames
                    FrameApprovedFinal(OnFrames(k)+1) = 0;
                    FrameStates(OnFrames(k)+1) = 0;
                end
            end
        end
        FrameStates = zeros(1, NFrames); % 1 means on, -1 means off, and 0 is inconclusive
        for j = 1:NFrames
            if (FrameApprovedFinal(j) == 1) & (FI.DetectedSpot(j) == 0)
                FrameStates(j) = -1;
            elseif (FrameApprovedFinal(j) == 1) & (FI.DetectedSpot(j) == 1)
                FrameStates(j) = 1;
            end
        end
        OffFrames = find(FrameStates == -1);
        if ~isempty(OffFrames)
            if length(OffFrames) > 1
                DiffOffFrames = diff(OffFrames);
                IsolatedOnFrames = find(DiffOffFrames == 2);
                for k=IsolatedOnFrames
                    FrameApprovedFinal(OffFrames(k)+1) = 0;
                    FrameStates(OffFrames(k)+1) = 0;
                end
            end
        end
        
        FrameStates = zeros(1, NFrames); % 1 means on, -1 means off, and 0 is inconclusive
        for j = 1:NFrames
            if (FrameApprovedFinal(j) == 1) & (FI.DetectedSpot(j) == 0)
                FrameStates(j) = -1;
            elseif (FrameApprovedFinal(j) == 1) & (FI.DetectedSpot(j) == 1)
                FrameStates(j) = 1;
            end
        end
        
        FI.FrameApprovedFinal = FrameApprovedFinal;
        FI.FrameStates = FrameStates;
        
        if sum(FrameApprovedFinal) == 0
            FI.Approved = 0;
        elseif p.FlaggingInfo.SickNucleus
            FI.Approved = 0;
        else
            FI.Approved = p.Approved;
        end
        
        p.FlaggingInfo = FI;
        CompiledParticles{ChN}(i) = p;
        
        CompiledParticles{ChN}(i).ManualFrameApproved = CompiledParticles{ChN}(i).FrameApproved;
        CompiledParticles{ChN}(i).ManualApproved = CompiledParticles{ChN}(i).Approved;
        CompiledParticles{ChN}(i).Approved = CompiledParticles{ChN}(i).FlaggingInfo.Approved;
        FrameApproved = zeros(1, length(CompiledParticles{ChN}(i).Frame));
        for j=1:length(CompiledParticles{ChN}(i).Frame)
            particle_idx = find(CompiledParticles{ChN}(i).FlaggingInfo.TrueFrames == CompiledParticles{ChN}(i).Frame(j));
            if ~isempty(particle_idx)
                FrameApproved(j) = CompiledParticles{ChN}(i).FlaggingInfo.FrameApprovedFinal(particle_idx);
            end
        end
        CompiledParticles{ChN}(i).FrameApproved = FrameApproved;
        
    end
end
