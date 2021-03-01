% function CompiledParticles{ChN} =
% Add2SpotFlaggingInfo(CompiledParticles{ChN}, FluoString)
FluoString = 'Fluo';

SchnitzIndexVector = [CompiledParticles{ChN}(:).schnitz];
AllSchnitzesWithSpots = unique(SchnitzIndexVector);


CloseSchnitzIndices = [493, 498, 609];

FarSchnitzIndices = [340, 451, 607];

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
    
    for Frame = unique([Twin1.FlaggingInfo.TrueFrames Twin2.FlaggingInfo.TrueFrames])
        Twin1IDs = [Twin1IDs, CurrentTwins(1)];
        Twin2IDs = [Twin1IDs, CurrentTwins(2)];
        AllFrames = [AllFrames, Frame];
        Twin1Index = find(Twin1.Frame == Frame);
        Twin1PrevIndex = find(Twin1.Frame == Frame-1);
        Twin1NextIndex = find(Twin1.Frame == Frame+1);
        Twin2Index = find(Twin2.Frame == Frame);
        Twin2PrevIndex = find(Twin2.Frame == Frame-1);
        Twin2NextIndex = find(Twin2.Frame == Frame+1);
        
        if ~isempty(Twin1Index) & ~isempty(Twin2Index)
            CurrentDistance = sqrt((Twin1.xPos(Twin1Index)-Twin2.xPos(Twin2Index))^2 + (Twin1.yPos(Twin1Index)-Twin2.yPos(Twin2Index))^2);
            DistanceBetweenTwins = [DistanceBetweenTwins, CurrentDistance];
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
            DeltaPrevFluoTwin1 = [DeltaPrevFluoTwin1, abs(Twin1.(FluoString)(Twin1Index)-Twin1.(FluoString)(Twin1PrevIndex))];
        else
            DeltaPrevFluoTwin1 = [DeltaPrevFluoTwin1, NaN];
        end
        
        if ~isempty(Twin1Index) & ~isempty(Twin1NextIndex)
            DeltaNextFluoTwin1 = [DeltaNextFluoTwin1, abs(Twin1.(FluoString)(Twin1Index)-Twin1.(FluoString)(Twin1NextIndex))];
        else
            DeltaNextFluoTwin1 = [DeltaNextFluoTwin1, NaN];
        end
        
        if ~isempty(Twin2Index) & ~isempty(Twin2PrevIndex)
            DeltaPrevFluoTwin2 = [DeltaPrevFluoTwin2, abs(Twin2.(FluoString)(Twin2Index)-Twin2.(FluoString)(Twin2PrevIndex))];
        else
            DeltaPrevFluoTwin2 = [DeltaPrevFluoTwin2, NaN];
        end
        
        if ~isempty(Twin2Index) & ~isempty(Twin2NextIndex)
            DeltaNextFluoTwin2 = [DeltaNextFluoTwin2, abs(Twin2.(FluoString)(Twin2Index)-Twin2.(FluoString)(Twin2NextIndex))];
        else
            DeltaNextFluoTwin2 = [DeltaNextFluoTwin2, NaN];
        end
        
    end

end
end

DeltaFluoCutoff = prctile([DeltaPrevFluoTwin1(DistanceBetweenTwins > nanmedian(DistanceBetweenTwins) & PrevTwinDistance > nanmedian(DistanceBetweenTwins)) DeltaPrevFluoTwin2(DistanceBetweenTwins > nanmedian(DistanceBetweenTwins)& PrevTwinDistance > nanmedian(DistanceBetweenTwins)) ], 99.5);

