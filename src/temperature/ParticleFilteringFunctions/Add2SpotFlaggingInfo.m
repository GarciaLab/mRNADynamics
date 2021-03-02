function CompiledParticles = Add2SpotFlaggingInfo(CompiledParticles, ChN, FluoString, FrameInfo, liveExperiment)


FrameTimes = [FrameInfo(:).Time];
snippet_size = liveExperiment.snippetSize_px;
SchnitzIndexVector = [CompiledParticles{ChN}(:).schnitz];
AllSchnitzesWithSpots = unique(SchnitzIndexVector);
% 
% 
% CloseSchnitzIndices = [493, 498, 609];
% 
% FarSchnitzIndices = [340, 451, 607];

Twin1IDs = [];
Twin2IDs = [];
AllFrames = [];
DistanceBetweenTwins = [];
PrevTwinDistance = [];
NextTwinDistance = [];
Twin1Fluos = [];
Twin2Fluos = [];
DeltaPrevFluoTwin1 = [];
DeltaNextFluoTwin1 = [];
DeltaPrevFluoTwin2 = [];
DeltaNextFluoTwin2 = [];

for SchnitzIndex = 1:length(AllSchnitzesWithSpots)
CurrentTwins = find(SchnitzIndexVector == AllSchnitzesWithSpots(SchnitzIndex));
if length(CurrentTwins) == 2
    Twin1 = CompiledParticles{ChN}(CurrentTwins(1));
    Twin2 = CompiledParticles{ChN}(CurrentTwins(2));
    DistanceToTwin1 = NaN(1, length(Twin1.FlaggingInfo.TrueFrames));
    DistanceToTwin2 = NaN(1, length(Twin2.FlaggingInfo.TrueFrames));
    for Frame = unique([Twin1.FlaggingInfo.TrueFrames Twin2.FlaggingInfo.TrueFrames])
        Twin1IDs = [Twin1IDs, CurrentTwins(1)];
        Twin2IDs = [Twin1IDs, CurrentTwins(2)];
        AllFrames = [AllFrames, Frame];
        Twin1SuperIndex = find(Twin1.FlaggingInfo.TrueFrames == Frame);
        Twin2SuperIndex = find(Twin1.FlaggingInfo.TrueFrames == Frame);
        Twin1Index = find(Twin1.Frame == Frame);
        Twin1PrevIndex = find(Twin1.Frame == Frame-1);
        Twin1NextIndex = find(Twin1.Frame == Frame+1);
        Twin2Index = find(Twin2.Frame == Frame);
        Twin2PrevIndex = find(Twin2.Frame == Frame-1);
        Twin2NextIndex = find(Twin2.Frame == Frame+1);
        PrevDeltaT = FrameTimes(Frame)-FrameTimes(Frame-1);
        if Frame + 1 <= length(FrameTimes)
            NextDeltaT = FrameTimes(Frame+1)-FrameTimes(Frame);
        else
            NextDeltaT = NaN;
        end
        if ~isempty(Twin1Index) & ~isempty(Twin2Index)
            CurrentDistance = sqrt((Twin1.xPos(Twin1Index)-Twin2.xPos(Twin2Index))^2 + (Twin1.yPos(Twin1Index)-Twin2.yPos(Twin2Index))^2);
            DistanceBetweenTwins = [DistanceBetweenTwins, CurrentDistance];
            DistanceToTwin1(Twin1SuperIndex) = CurrentDistance;
            DistanceToTwin2(Twin2SuperIndex) = CurrentDistance;
        else
            DistanceBetweenTwins = [DistanceBetweenTwins, NaN];
        end
        
        if ~isempty(Twin1PrevIndex) & ~isempty(Twin2PrevIndex)
            PreviousDistance = sqrt((Twin1.xPos(Twin1PrevIndex)-Twin2.xPos(Twin2PrevIndex))^2 + (Twin1.yPos(Twin1PrevIndex)-Twin2.yPos(Twin2PrevIndex))^2);
            PrevTwinDistance = [PrevTwinDistance, PreviousDistance];
        else
            PrevTwinDistance = [PrevTwinDistance, NaN];
        end
        
        if ~isempty(Twin1NextIndex) & ~isempty(Twin2NextIndex)
            NextDistance = sqrt((Twin1.xPos(Twin1NextIndex)-Twin2.xPos(Twin2NextIndex))^2 + (Twin1.yPos(Twin1NextIndex)-Twin2.yPos(Twin2NextIndex))^2);
            NextTwinDistance = [NextTwinDistance, NextDistance];
        else
            NextTwinDistance = [NextTwinDistance, NaN];
        end
        
        if ~isempty(Twin1Index)
            Twin1Fluos = [Twin1Fluos, Twin1.(FluoString)(Twin1Index)];
        else
            Twin1Fluos = [Twin1Fluos, NaN];
        end
        
        if ~isempty(Twin2Index)
            Twin2Fluos = [Twin2Fluos, Twin2.(FluoString)(Twin2Index)];
        else
            Twin2Fluos = [Twin2Fluos, NaN];
        end
        
        if ~isempty(Twin1Index) & ~isempty(Twin1PrevIndex)
            DeltaPrevFluoTwin1 = [DeltaPrevFluoTwin1, abs(Twin1.(FluoString)(Twin1Index)-Twin1.(FluoString)(Twin1PrevIndex))/PrevDeltaT];
        else
            DeltaPrevFluoTwin1 = [DeltaPrevFluoTwin1, NaN];
        end
        
        if ~isempty(Twin1Index) & ~isempty(Twin1NextIndex)
            DeltaNextFluoTwin1 = [DeltaNextFluoTwin1, abs(Twin1.(FluoString)(Twin1Index)-Twin1.(FluoString)(Twin1NextIndex))/NextDeltaT];
        else
            DeltaNextFluoTwin1 = [DeltaNextFluoTwin1, NaN];
        end
        
        if ~isempty(Twin2Index) & ~isempty(Twin2PrevIndex)
            DeltaPrevFluoTwin2 = [DeltaPrevFluoTwin2, abs(Twin2.(FluoString)(Twin2Index)-Twin2.(FluoString)(Twin2PrevIndex))/PrevDeltaT];
        else
            DeltaPrevFluoTwin2 = [DeltaPrevFluoTwin2, NaN];
        end
        
        if ~isempty(Twin2Index) & ~isempty(Twin2NextIndex)
            DeltaNextFluoTwin2 = [DeltaNextFluoTwin2, abs(Twin2.(FluoString)(Twin2Index)-Twin2.(FluoString)(Twin2NextIndex))/NextDeltaT];
        else
            DeltaNextFluoTwin2 = [DeltaNextFluoTwin2, NaN];
        end
        
    end
    
    CompiledParticles{ChN}(CurrentTwins(1)).FlaggingInfo.DistanceToTwin = DistanceToTwin1;
    CompiledParticles{ChN}(CurrentTwins(2)).FlaggingInfo.DistanceToTwin = DistanceToTwin2;
else
    CompiledParticles{ChN}(CurrentTwins(1)).FlaggingInfo.DistanceToTwin  = NaN(1, length(CompiledParticles{ChN}(CurrentTwins(1)).FlaggingInfo.TrueFrames)); 

end
end

DeltaFluoCutoff = prctile([DeltaPrevFluoTwin1(DistanceBetweenTwins > nanmedian(DistanceBetweenTwins) & PrevTwinDistance > nanmedian(DistanceBetweenTwins)) DeltaPrevFluoTwin2(DistanceBetweenTwins > nanmedian(DistanceBetweenTwins)& PrevTwinDistance > nanmedian(DistanceBetweenTwins)) ], 99.5);
MaxDeltaFluoCutoff = max([DeltaPrevFluoTwin1(DistanceBetweenTwins > nanmedian(DistanceBetweenTwins) & PrevTwinDistance > nanmedian(DistanceBetweenTwins)) DeltaPrevFluoTwin2(DistanceBetweenTwins > nanmedian(DistanceBetweenTwins)& PrevTwinDistance > nanmedian(DistanceBetweenTwins)) ]);
%%
for cp =1:length(CompiledParticles{ChN})
    Pos2spotFrameApproved = NaN(1, length(CompiledParticles{ChN}(cp).Frame));
    DistanceToTwin = CompiledParticles{ChN}(cp).FlaggingInfo.DistanceToTwin(ismember(CompiledParticles{ChN}(cp).FlaggingInfo.TrueFrames, CompiledParticles{ChN}(cp).Frame));
    if ~isempty(find(DistanceToTwin < snippet_size, 1))
        Pos2spotFrameApproved(DistanceToTwin < snippet_size) = -1;
    end
    if length(CompiledParticles{ChN}(cp).Frame) <= 1
        Pos2spotFrameApproved(1) = 1;
        continue
    end
    for FrameIndex=2:length(CompiledParticles{ChN}(cp).Frame)
        if FrameIndex == 2
            PrevFrameIndex = 1;
        elseif sum(Pos2spotFrameApproved(1:(FrameIndex-1)) == 1) > 0
            PrevFrameIndex = find(Pos2spotFrameApproved(1:(FrameIndex-1)) == 1, 1, 'last');
        elseif sum(Pos2spotFrameApproved(1:(FrameIndex-1)) == 0) > 0
            PrevFrameIndex = find(Pos2spotFrameApproved(1:(FrameIndex-1)) == 0, 1, 'last');
        else
            continue
        end
        DeltaT = FrameTimes(CompiledParticles{ChN}(cp).Frame(FrameIndex))-FrameTimes(CompiledParticles{ChN}(cp).Frame(PrevFrameIndex));
        DeltaFluo = (CompiledParticles{ChN}(cp).(FluoString)(FrameIndex)-CompiledParticles{ChN}(cp).(FluoString)(PrevFrameIndex))/DeltaT;
        if abs(DeltaFluo) <= DeltaFluoCutoff 
            Pos2spotFrameApproved(FrameIndex) = 1;
            if sum(Pos2spotFrameApproved(1:(FrameIndex-1)) == 1) == 0
                Pos2spotFrameApproved(PrevFrameIndex) = 1;
            end
        elseif sum(Pos2spotFrameApproved(1:(FrameIndex-1)) == 1) > 0
            Pos2spotFrameApproved(FrameIndex) = -1;
        elseif FrameIndex == max(CompiledParticles{ChN}(cp).Frame)
            Pos2spotFrameApproved(FrameIndex) = -1;
            Pos2spotFrameApproved(PrevFrameIndex) = -1;
        else
            NextIndex = FrameIndex+1;
            DeltaT2 = FrameTimes(CompiledParticles{ChN}(cp).Frame(NextIndex))-FrameTimes(CompiledParticles{ChN}(cp).Frame(FrameIndex));
            DeltaT3 = FrameTimes(CompiledParticles{ChN}(cp).Frame(NextIndex))-FrameTimes(CompiledParticles{ChN}(cp).Frame(PrevFrameIndex));
            DeltaFluo2 = (CompiledParticles{ChN}(cp).(FluoString)(NextIndex)-CompiledParticles{ChN}(cp).(FluoString)(FrameIndex))/DeltaT2;
            DeltaFluo3 = (CompiledParticles{ChN}(cp).(FluoString)(NextIndex)-CompiledParticles{ChN}(cp).(FluoString)(PrevFrameIndex))/DeltaT3;
            if abs(DeltaFluo2) < abs(DeltaFluo3)
                Pos2spotFrameApproved(FrameIndex) = 1;
                Pos2spotFrameApproved(PrevFrameIndex) = -1;
            else
                Pos2spotFrameApproved(FrameIndex) = -1;
                Pos2spotFrameApproved(PrevFrameIndex) = 1;
            end
        end
    end
    
    
    FullPos2spotFrameApproved = zeros(1, length(CompiledParticles{ChN}(cp).FlaggingInfo.TrueFrames));
    FullPos2spotFrameApproved(find(ismember(CompiledParticles{ChN}(cp).FlaggingInfo.TrueFrames, CompiledParticles{ChN}(cp).Frame))) = ...
        Pos2spotFrameApproved;

    CompiledParticles{ChN}(cp).FlaggingInfo.TwoSpotApproved = FullPos2spotFrameApproved;
end

