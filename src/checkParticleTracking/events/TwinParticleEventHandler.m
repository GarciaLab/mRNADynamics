function keyInputHandler = TwinParticleEventHandler(cptState)
    
    function keyInput(cc)
        
        
        if cc == '7'
            if ~isempty(cptState.TwinParticle)
         
                [cptState.CurrentParticle, ~, cptState.ManualZFlag, cptState.TwinParticle] = ...
                    changeParticle(cptState.TwinParticle, cptState.Particles, cptState.numParticles, cptState.CurrentChannelIndex);
                
                cptState.DisplayRange = [];     
            else
                disp('No twin particle to switch to from current particle.')
            end
        elseif cc == '6'
            if ~isempty(cptState.TwinParticle)
         
                [cptState.CurrentParticle, cptState.CurrentFrame, cptState.ManualZFlag, cptState.TwinParticle] = ...
                    changeParticle(cptState.TwinParticle, cptState.Particles, cptState.numParticles, cptState.CurrentChannelIndex);
                
                cptState.DisplayRange = [];     
            else
                disp('No twin particle to switch to from current particle.')
            end
        elseif cc == '&'
            [cptState.Particles{cptState.CurrentChannelIndex}, cptState.CurrentParticle] = flipTwinTrajectories(cptState);
        elseif cc == '1'
            if cptState.HideSingleSliceTrace
                cptState.HideSingleSliceTrace = false;
            else
                cptState.HideSingleSliceTrace = true;
            end
        end
       
    end

    keyInputHandler = @keyInput;
end
