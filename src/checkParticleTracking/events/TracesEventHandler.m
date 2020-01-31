function keyInputHandler = TracesEventHandler(cptState)
    
    function keyInput(cc)
        if cc == 'c'
            [cptState.PreviousParticle, cptState.Particles] = combineTraces(cptState.Spots, ...
            cptState.CurrentChannel, cptState.CurrentFrame, cptState.Particles, cptState.CurrentParticle);
            
        elseif cc == 'd'
            % Separate traces forward at the current frame.
            [cptState.Particles, cptState.PreviousParticle] = separateTraces(cptState.Particles, ...
                cptState.CurrentChannel, cptState.CurrentFrame, cptState.CurrentParticle);
        elseif cc == 'q'
            % Approve a trace
            if cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Approved == 1
                cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Approved = 2;
            elseif cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Approved == 0
                cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Approved = 1;
            elseif cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Approved == 2
                cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Approved = 0;
            end
            
        elseif cc == 'w'
            % Disapprove a trace
            if cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Approved == -1
                cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Approved = 0;
            else
                cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Approved = -1;
            end

        elseif cc == 'h'
            if cptState.HideApprovedFlag == 0
                cptState.HideApprovedFlag = 1; %Show only non-approved traces
            elseif cptState.HideApprovedFlag == 1
                cptState.HideApprovedFlag = 2; %Show only yellow and red traces
            elseif cptState.HideApprovedFlag == 2
                cptState.HideApprovedFlag = 0;
            end
        end
    end

    keyInputHandler = @keyInput;
end
