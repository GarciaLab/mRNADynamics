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
            oldStatus = cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle).Approved;
            if oldStatus>=0
              % update Particles structure itself
              newStatus = (oldStatus+1)*(oldStatus~=2);
              cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle).Approved = newStatus;      
              % update auxiliary particles structures
              cptState = updateAuxiliaryStatus(cptState,oldStatus,newStatus);
            end
        elseif cc == 'w'
            % Disapprove a trace
            oldStatus = cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle).Approved;
            if oldStatus<=0
              newStatus = (oldStatus-1)*(oldStatus~=-1);
              cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle).Approved = newStatus;  
              cptState = updateAuxiliaryStatus(cptState,oldStatus,newStatus);
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
