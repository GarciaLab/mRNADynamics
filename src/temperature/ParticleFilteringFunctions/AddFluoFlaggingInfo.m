function CompiledParticles = AddFluoFlaggingInfo(CompiledParticles, ChN, FluoString, FrameInfo, liveExperiment)


FrameTimes = [FrameInfo(:).Time];
snippet_size = liveExperiment.snippetSize_px;
SchnitzIndexVector = [CompiledParticles{ChN}(:).schnitz];
AllSchnitzesWithSpots = unique(SchnitzIndexVector);
% 
% 
% CloseSchnitzIndices = [493, 498, 609];
% 
% FarSchnitzIndices = [340, 451, 607];

DeltaFluo = [];


for CPIndex = 1:length(CompiledParticles{ChN})
    p = CompiledParticles{ChN}(CPIndex);
    if length(p.FlaggingInfo.TrueFrames) == 1
        continue
    end
       
    for FrameIndex = 2:length(p.FlaggingInfo.TrueFrames)
        Frame = p.FlaggingInfo.TrueFrames(FrameIndex);
        pIndex = find(p.Frame == Frame);
        pPrevIndex = find(p.Frame == Frame-1);
        pNextIndex = find(p.Frame == Frame+1);
        PrevDeltaT = FrameTimes(Frame)-FrameTimes(Frame-1);
        if Frame + 1 <= length(FrameTimes)
            NextDeltaT = FrameTimes(Frame+1)-FrameTimes(Frame);
        else
            NextDeltaT = NaN;
        end
        
        
        
     
        if ~isempty(pIndex) & ~isempty(pPrevIndex)
            DeltaFluo = [DeltaFluo, abs(p.(FluoString)(pIndex)-p.(FluoString)(pPrevIndex))/PrevDeltaT];
        else
            DeltaFluo = [DeltaFluo, NaN];
        end
       
    end
    
    CompiledParticles{ChN}(CPIndex).FlaggingInfo.DeltaFluo = DeltaFluo;

end

DeltaFluoCutoff = prctile(DeltaFluo, 99);
%%
for cp =1:length(CompiledParticles{ChN})
    FluoApproved = NaN(1, length(CompiledParticles{ChN}(cp).Frame));
    if length(CompiledParticles{ChN}(cp).Frame) <= 1
        FluoApproved(1) = 1;
        continue
    end
    for FrameIndex=2:length(CompiledParticles{ChN}(cp).Frame)
        if FrameIndex == 2
            PrevFrameIndex = 1;
        elseif sum(FluoApproved(1:(FrameIndex-1)) == 1) > 0
            PrevFrameIndex = find(FluoApproved(1:(FrameIndex-1)) == 1, 1, 'last');
        elseif sum(FluoApproved(1:(FrameIndex-1)) == 0) > 0
            PrevFrameIndex = find(FluoApproved(1:(FrameIndex-1)) == 0, 1, 'last');
        else
            continue
        end
        DeltaT = FrameTimes(CompiledParticles{ChN}(cp).Frame(FrameIndex))-FrameTimes(CompiledParticles{ChN}(cp).Frame(PrevFrameIndex));
        DeltaFluo = (CompiledParticles{ChN}(cp).(FluoString)(FrameIndex)-CompiledParticles{ChN}(cp).(FluoString)(PrevFrameIndex))/DeltaT;
        if abs(DeltaFluo) <= DeltaFluoCutoff 
            FluoApproved(FrameIndex) = 1;
            if sum(FluoApproved(1:(FrameIndex-1)) == 1) == 0
                FluoApproved(PrevFrameIndex) = 1;
            end
        elseif sum(FluoApproved(1:(FrameIndex-1)) == 1) > 0
            FluoApproved(FrameIndex) = -1;
        elseif FrameIndex == max(CompiledParticles{ChN}(cp).Frame)
            FluoApproved(FrameIndex) = -1;
            FluoApproved(PrevFrameIndex) = -1;
        else
            NextIndex = FrameIndex+1;
            DeltaT2 = FrameTimes(CompiledParticles{ChN}(cp).Frame(NextIndex))-FrameTimes(CompiledParticles{ChN}(cp).Frame(FrameIndex));
            DeltaT3 = FrameTimes(CompiledParticles{ChN}(cp).Frame(NextIndex))-FrameTimes(CompiledParticles{ChN}(cp).Frame(PrevFrameIndex));
            DeltaFluo2 = (CompiledParticles{ChN}(cp).(FluoString)(NextIndex)-CompiledParticles{ChN}(cp).(FluoString)(FrameIndex))/DeltaT2;
            DeltaFluo3 = (CompiledParticles{ChN}(cp).(FluoString)(NextIndex)-CompiledParticles{ChN}(cp).(FluoString)(PrevFrameIndex))/DeltaT3;
            if abs(DeltaFluo2) < abs(DeltaFluo3)
                FluoApproved(FrameIndex) = 1;
                FluoApproved(PrevFrameIndex) = -1;
            else
                FluoApproved(FrameIndex) = -1;
                FluoApproved(PrevFrameIndex) = 1;
            end
        end
    end
    
    
    FullFluoApproved = zeros(1, length(CompiledParticles{ChN}(cp).FlaggingInfo.TrueFrames));
    FullFluoApproved(find(ismember(CompiledParticles{ChN}(cp).FlaggingInfo.TrueFrames, CompiledParticles{ChN}(cp).Frame))) = ...
        FluoApproved;

    CompiledParticles{ChN}(cp).FlaggingInfo.FluoApproved = FullFluoApproved;
end

