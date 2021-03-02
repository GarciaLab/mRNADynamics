function CompiledParticles = AddBoundaryPositionInfo(CompiledParticles, ChN, liveExperiment)
%%
xDim = liveExperiment.xDim;
yDim = liveExperiment.xDim;
snippet_size = liveExperiment.snippetSize_px;
PixelSize = liveExperiment.pixelSize_um;

FrameInfo = getFrameInfo(liveExperiment);
FrameTimes = [FrameInfo(:).Time];
nucleusDiameters = zeros(1, 6);
for nc=9:14
    nucleusDiameters(nc-8) = getDefaultParameters(FrameInfo,[['d', num2str(nc)]])/PixelSize; % in pixels
end

for i = 1:length(CompiledParticles{ChN})
    p = CompiledParticles{ChN}(i);
    sc = p.schnitzcell;
    
    SpotMedianXPos = median(p.xPos);
    SpotMedianYPos = median(p.yPos);
    
    NFrames = length(p.FlaggingInfo.TrueFrames);
    
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
            Xdelta = (p.xPos(j)-p.xPos(prev_j)); % in pixels
            Ydelta = (p.yPos(j)-p.yPos(prev_j)); % in pixels
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
    
    % discontinutity error
end

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
AllSpotXs = [];
AllSpotYs = [];

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

VeloCutoff = prctile(AllVelos, 75);
for i = 1:length(CompiledParticles{ChN})
    p = CompiledParticles{ChN}(i); 
    NFrames = length(p.FlaggingInfo.TrueFrames);
    PositionApproved = zeros(1, NFrames);
    SpotX = p.FlaggingInfo.SpotXPos;
    SpotY = p.FlaggingInfo.SpotYPos;
    schnitzX = p.FlaggingInfo.schnitzXPos;
    schnitzY = p.FlaggingInfo.schnitzYPos;
    SchnitzBasedBoundaries = zeros(1, NFrames, 'logical');
    NuclearDiameter = nucleusDiameters(p.cycle-8);
    for j= 1:NFrames
        if (schnitzX(j) > NuclearDiameter) & (schnitzX(j) < xDim-NuclearDiameter) & ...
          (schnitzY(j) > NuclearDiameter) & (schnitzY(j) < yDim-NuclearDiameter)   
            SchnitzBasedBoundaries(j) = true;
        end
    end
    
    FrameHasSpot = ismember(p.FlaggingInfo.TrueFrames, p.Frame);
    
    SpotHasEverBeenOnPreviously = false;
    SpotIsOnInPreviousFrame = false;
    
    FirstOnSpot = find(FrameHasSpot, 1);
    if FirstOnSpot == 1
       PositionApproved(FirstOnSpot) = -1;
    else
    DeltaT = FrameTimes(FirstOnSpot)-FrameTimes(FirstOnSpot-1);
    DeltaS = DeltaT*VeloCutoff;

    if (SpotX(FirstOnSpot) < snippet_size+DeltaS) & ...
            (SpotX(FirstOnSpot) > xDim - snippet_size-DeltaS) & ...
            (SpotY(FirstOnSpot) < snippet_size+DeltaS) & ...
            (SpotY(FirstOnSpot) > yDim - snippet_size-DeltaS)
       PositionApproved(1:FirstOnSpot) = -1;
    else
       NextDeltaS = VeloCutoff*FrameTimes(FirstOnSpot+1)-FrameTimes(FirstOnSpot);
       BadBoundaryPosition = (SpotX(FirstOnSpot) < snippet_size+NextDeltaS) & ...
            (SpotX(FirstOnSpot) > xDim - snippet_size-NextDeltaS) & ...
            (SpotY(FirstOnSpot) < snippet_size+NextDeltaS) & ...
            (SpotY(FirstOnSpot) > yDim - snippet_size-NextDeltaS);
       if FrameHasSpot(FirstOnSpot+1) & ~BadBoundaryPosition
           PositionApproved(FirstOnSpot-1:FirstOnSpot) = 1;
           if FirstOnSpot-1 > 1
               for j = 1:FirstOnSpot-2
                   if SchnitzBasedBoundaries(j) 
                        PositionApproved(j) = 1;
                   else
                       PositionApproved(j) = -1;
                   end
               end
           end
           
       else
           PositionApproved(FirstOnSpot) = -1;
       end
       
       
    end
    end
    
    SpotHasEverBeenOnPreviously = true;
    SpotIsOnInPreviousFrame = true;
    
    if FirstOnSpot + 1 < NFrames
                
        for FrameIndex = (FirstOnSpot+1):NFrames
            DeltaT = FrameTimes(FrameIndex)-FrameTimes(FrameIndex-1);
            DeltaS = DeltaT*VeloCutoff;
            if FrameHasSpot(FrameIndex)
                if SpotIsOnInPreviousFrame
                    PositionApproved(FrameIndex) = PositionApproved(FrameIndex-1);
                elseif (SpotX(FrameIndex) > snippet_size+DeltaS) & ...
                        (SpotX(FrameIndex) < xDim - snippet_size - DeltaS) & ...
                        (SpotY(FrameIndex) > snippet_size + DeltaS) & ...
                        (SpotY(FrameIndex) < yDim - snippet_size - DeltaS)
                    if FrameIndex + 1 >= NFrames 
                        PositionApproved(FrameIndex) = 1;
                    elseif FrameHasSpot(FrameIndex+1)
                        PositionApproved(FrameIndex) = 1;
                    else
                        NextDeltaS = VeloCutoff*FrameTimes(FrameIndex+1)-FrameTimes(FrameIndex);
                        BadBoundaryPosition = (SpotX(FrameIndex) < snippet_size+NextDeltaS) & ...
                            (SpotX(FrameIndex) > xDim - snippet_size-NextDeltaS) & ...
                            (SpotY(FrameIndex) < snippet_size+NextDeltaS) & ...
                            (SpotY(FrameIndex) > yDim - snippet_size-NextDeltaS);
                        
                        if ~BadBoundaryPosition
                            PositionApproved(FrameIndex) = 1;
                        else
                            PositionApproved(FrameIndex) = -1;
                        end
                    end
                end
            elseif SpotIsOnInPreviousFrame
                BadBoundaryPosition = (SpotX(FrameIndex) < snippet_size+DeltaS) & ...
                            (SpotX(FrameIndex) > xDim - snippet_size-DeltaS) & ...
                            (SpotY(FrameIndex) < snippet_size+DeltaS) & ...
                            (SpotY(FrameIndex) > yDim - snippet_size-DeltaS);
                if BadBoundaryPosition 
                    PositionApproved(FrameIndex) = -1;
                else
                    PositionApproved(FrameIndex) = 1;
                end
            else
                PositionApproved(FrameIndex) = PositionApproved(FrameIndex-1);
            end
            if FrameHasSpot(FrameIndex)
                SpotIsOnInPreviousFrame = true;
            else
                SpotIsOnInPreviousFrame = false;
            end
                
        end
    end
    

      
    p.FlaggingInfo.PositionApproved = PositionApproved;
    CompiledParticles{ChN}(i) = p;
end

