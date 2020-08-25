function keyInputHandler = ParticleChangeEventHandler(cptState)
    
    function keyInput(cc)
        if cc == 'k'
        
            try
                ParticleJump = inputdlg('Particle to jump to:', ...
                    'Move to particle');
                ParticleJump = str2double(ParticleJump{1});
            catch
                ParticleJump = CurrentParticle;
            end
            
            cptState = changeParticle(ParticleJump,cptState);

            cptState.DisplayRange = [];
        
        elseif cc == 'r'
            cptState.Particles = orderParticles(cptState.numParticles(), cptState.CurrentChannelIndex, cptState.Particles);
        
        elseif cc == 'p'
            % Identify a particle. It will also tell you the particle associated with the clicked nucleus.
            identifyParticle(cptState.Spots, cptState.Particles, cptState.CurrentFrame, ...
                cptState.CurrentChannelIndex, cptState.UseHistoneOverlay, cptState.schnitzcells);
        
        elseif cc == '\'
            % Moves to clicked particle.
            [cptState.CurrentParticle, cptState.CurrentFrame, cptState.ManualZFlag] =...
                toNearestParticle(cptState.Spots, ...
                cptState.Particles, cptState.CurrentFrame, cptState.CurrentChannelIndex,...
                cptState.UseHistoneOverlay, cptState.schnitzcells);
        
        elseif cc == 'u'
            [x2, ~, ~] = SpotsXYZ(cptState.Spots{cptState.CurrentChannelIndex}(cptState.CurrentFrame));
            
            if ~isempty(x2)
                ClickedSpot = ginput(1);
                
                UnfilterSpot(cptState.Spots{cptState.CurrentChannelIndex}, cptState.SpotFilter{cptState.CurrentChannelIndex}, ...
                    ClickedSpot, cptState.Particles{cptState.CurrentChannelIndex}, cptState.CurrentFrame)
            end

        elseif (cc == 'm') & (cptState.CurrentParticle < cptState.numParticles())
            cptState = goNextParticle(cptState,0);
        
        elseif (cc == 'n') & (cptState.CurrentParticle > 1)
            cptState = goNextParticle(cptState,1);

        elseif cc == 'i'
            warning(' AR 1/15/18: This is currently deprecated. Talk to HG if you need this function.')
        end
    end

    keyInputHandler = @keyInput;
end
