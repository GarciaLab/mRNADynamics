function keyInputHandler = TracesEventHandler(cptState)
    
    function keyInput(cc)
        if cc == 'c'
            warning('NL: this feature has not been tested. Proceed with caution')
            [cptState.PreviousParticle, cptState.Particles] = combineTraces(cptState.Spots, ...
            cptState.CurrentChannelIndex, cptState.CurrentFrame, cptState.Particles, cptState.CurrentParticle);
            
        elseif cc == 'd'
            % Separate traces forward at the current frame.
            cptState = separateTraces(cptState);
              
        elseif cc == 'q'
            % Approve a trace
            aState = cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle).Approved;
            if aState>=0
              % update Particles structure itself
              newState = (aState+1)*(aState~=2);
              cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle).Approved = newState;      
              % update auxiliary particles structures
              cptState = updateAuxiliaryStatus(cptState,newState);
            end
        elseif cc == 'w'
            % Disapprove a trace
            aState = cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle).Approved;
            if aState<=0
              newState = (aState-1)*(aState~=-1);
              cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle).Approved = newState;  
              cptState = updateAuxiliaryStatus(cptState,newState);
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
