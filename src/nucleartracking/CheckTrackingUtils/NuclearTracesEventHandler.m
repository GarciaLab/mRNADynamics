function keyInputHandler = NuclearTracesEventHandler(cntState)
    
    function keyInput(cc)
        if cc == 'c'
            [cntState.PreviousNucleus, cntState.schnitzcells] = combineNuclearTraces(cntState.CurrentFrame,...
                cntState.schnitzcells, cntState.Ellipses, cntState.CurrentNucleus);
            
        elseif cc == 'd'
            % Separate traces forward at the current frame.
            [cntState.schnitzcells, cntState.PreviousNucleus] = separateNuclearTraces(cntState.CurrentNucleus, ...
                cntState.CurrentFrame, cntState.schnitzcells);
        elseif cc == 'q'
            % Approve a trace
            if cntState.schnitzcells(cntState.CurrentNucleus).Approved == 1
                cntState.schnitzcells(cntState.CurrentNucleus).Approved = 2;
            elseif cntState.schnitzcells(cntState.CurrentNucleus).Approved == 0
                cntState.schnitzcells(cntState.CurrentNucleus).Approved = 1;
            elseif cntState.schnitzcells(cntState.CurrentNucleus).Approved == 2
                cntState.schnitzcells(cntState.CurrentNucleus).Approved = 0;
            end
            
        elseif cc == 'w'
            % Disapprove a trace
            if cntState.schnitzcells(cntState.CurrentNucleus).Approved == -1
                cntState.schnitzcells(cntState.CurrentNucleus).Approved = 0;
            else
                cntState.schnitzcells(cntState.CurrentNucleus).Approved = -1;
            end

        elseif cc == 'h'
            if cntState.HideApprovedFlag == 0
                cntState.HideApprovedFlag = 1; %Show only non-approved traces
            elseif cntState.HideApprovedFlag == 1
                cntState.HideApprovedFlag = 2; %Show only yellow and red traces
            elseif cntState.HideApprovedFlag == 2
                cntState.HideApprovedFlag = 0;
            end
        end
    end

    keyInputHandler = @keyInput;
end
