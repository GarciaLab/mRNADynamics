function keyInputHandler = DeleteSpotEventHandler(cptState)

    function switchedParticlesFlag = doDeleteSpot(shouldQueryUser)
        
        switchedParticlesFlag = false;
        isOnlyFrame = length(cptState.Particles{cptState.CurrentChannelIndex}...
            (cptState.CurrentParticle).Frame) == 1;
        
        disp(['Removing frame ', num2str(cptState.CurrentFrame),...
            ' from particle ' num2str(cptState.CurrentParticle)])
        
        cptState =...
            ...
            removeSpot(cptState, shouldQueryUser);
        

        
        if isOnlyFrame
            
            
            %switch to another particle if this particle stops existing
            cptState.CurrentParticle = cptState.CurrentParticle - 1;
            
            if (cptState.CurrentParticle+1) < cptState.numParticles
                cptState.CurrentParticle = cptState.CurrentParticle+1;
            end
            
            if cptState.numParticles == 1
                cptState.lastParticle = 1;
            end
            
            cptState.CurrentFrame=cptState.Particles{cptState.CurrentChannelIndex}...
                (cptState.CurrentParticle).Frame(1);
            
            switchedParticlesFlag = true;
            
        elseif cptState.CurrentFrame > 1
            cptState.CurrentFrame = getCurrentParticle(cptState).Frame(1);
            cptState.ManualZFlag=0;
        elseif cptState.CurrentFrame < length({cptState.Spots{1}.Fits})
            cptState.CurrentFrame= cptState.CurrentFrame+1;
            cptState.ManualZFlag=0;
        end
        
        if switchedParticlesFlag
            disp('Switched particles.')
        end
    end

    function keyInput(cc)

        if cc == '#' %remove a spot from cptState.Spots and erase its frame in Particles
            
            doDeleteSpot(true);
            
        elseif cc == '^' %remove a whole trace from Spots, Particles, and SpotFilter.
            
            switchedParticlesFlag = false;
            
            while ~switchedParticlesFlag
                
                switchedParticlesFlag = doDeleteSpot(false);
                
            end
            
        end
    end

keyInputHandler = @keyInput;

end