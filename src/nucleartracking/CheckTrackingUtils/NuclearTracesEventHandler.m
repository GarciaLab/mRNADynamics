function keyInputHandler = NuclearTracesEventHandler(cntState, FrameInfo, ncFrames)
    % NuclearTracesEventHandler.m
    % author: Gabriella Martini
    % date created: 9/7/20
    % date last modified: 9/13/20
    function keyInput(cc)
        if cc == 'c'
            [cntState.PreviousNucleus, cntState.schnitzcells] = combineNuclearTraces(cntState.CurrentFrame,...
                cntState.schnitzcells, cntState.Ellipses, cntState.CurrentNucleus, FrameInfo, ncFrames);
            
        elseif cc == 'd'
            % Separate traces forward at the current frame.
            if length(cntState.schnitzcells(cntState.CurrentNucleus).frames) > 1
                cntState.schnitzcells = SeparateNuclearTraces(cntState.CurrentNucleus, ...
                    cntState.CurrentFrame, cntState.schnitzcells, FrameInfo, ncFrames);
            end
        elseif cc == 'D'
            
            if length(cntState.schnitzcells(cntState.CurrentNucleus).frames) > 1
                HoldCurrentNucleus = cntState.CurrentNucleus;
                cntState.schnitzcells=SeparateAllFrameNuclearTraces(cntState.CurrentNucleus,...
                    cntState.CurrentFrame, cntState.schnitzcells, FrameInfo, ncFrames);
                [cntState.CurrentNucleus,cntState.CurrentFrame, cntState.ManualZFlag] = ...
                    changeNucleus(HoldCurrentNucleus, cntState.schnitzcells, cntState.numNuclei())
            end

            
            
        elseif cc == 'q'
            % Approve a trace
            
            if (cntState.schnitzcells(cntState.CurrentNucleus).Approved == 1)
                cntState.schnitzcells(cntState.CurrentNucleus).Approved = 2;
            elseif cntState.schnitzcells(cntState.CurrentNucleus).Approved <= 0
                cntState.schnitzcells(cntState.CurrentNucleus).Approved = 1;
            elseif cntState.schnitzcells(cntState.CurrentNucleus).Approved == 2
                cntState.schnitzcells(cntState.CurrentNucleus).Approved = 0;
            end
            cntState.schnitzcells(cntState.CurrentNucleus).Checked = 1;
            
        elseif cc == 'w'
            % Disapprove a trace
            cntState.schnitzcells(cntState.CurrentNucleus).Approved = 0;
            cntState.schnitzcells(cntState.CurrentNucleus).Checked = 1;
        elseif cc == 'f'
            [flag, flag_string]  = chooseFlag;
            disp(flag_string)
            cntState.schnitzcells(cntState.CurrentNucleus).Flag = flag;
            cntState.schnitzcells(cntState.CurrentNucleus).Checked = 1;
        
%         elseif cc == 'h'
%             if cntState.HideApprovedFlag == 0
%                 cntState.HideApprovedFlag = 1; %Show only non-approved traces
%             elseif cntState.HideApprovedFlag == 1
%                 cntState.HideApprovedFlag = 2; %Show only yellow and red traces
%             elseif cntState.HideApprovedFlag == 2
%                 cntState.HideApprovedFlag = 0;
%             end
        end
    end

    keyInputHandler = @keyInput;
end
