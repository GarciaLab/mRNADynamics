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
            [cptState.CurrentParticle, cptState.CurrentFrame, cptState.ManualZFlag, cptState.TwinParticle] = ...
                changeParticle(ParticleJump, cptState.Particles, cptState.numParticles, cptState.CurrentChannelIndex, cptState.UseTwinTraces);

            cptState.DisplayRange = [];
        
        elseif cc == 'r'
            cptState.Particles = orderParticles(cptState.numParticles(), cptState.CurrentChannelIndex, cptState.Particles);
        
        elseif cc == 'p'
            % Identify a particle. It will also tell you the particle associated with the clicked nucleus.
            
            identifyParticle(cptState.Spots, cptState.Particles, cptState.CurrentFrame, ...
                cptState.CurrentChannelIndex, cptState.UseHistoneOverlay, cptState.schnitzcells);
        
        elseif cc == '|'
            % Remove nuclear association with current particle
            [cptState.Particles, cptState.CurrentParticle, cptState.CurrentFrame] =...
                removeNucleusAssociatedWithParticle(cptState.Spots, cptState.Particles, cptState.CurrentParticle,...
                cptState.CurrentFrame, cptState.CurrentChannelIndex, cptState.UseHistoneOverlay, cptState.schnitzcells);

        elseif cc == '\'
            % Associates clicked nucleus with current Particle
            [cptState.Particles, cptState.CurrentParticle, cptState.CurrentFrame] =...
                addNucleusAssociatedWithParticle(cptState.Spots, cptState.Particles, cptState.CurrentParticle,...
                cptState.CurrentFrame, cptState.CurrentChannelIndex, cptState.UseHistoneOverlay, cptState.schnitzcells);

        elseif cc == ')'
            % Switches nuclear association to clicked particle and switches to the first frame for that particle.
            % Added by G. Martini on 11/23/20
            [cptState.Particles, cptState.CurrentParticle, cptState.CurrentFrame, cptState.ManualZFlag, cptState.TwinParticle] =...
                changeParticleAssociatedWithNucleus(cptState.Spots, cptState.Particles, cptState.CurrentParticle,...
                cptState.CurrentFrame, cptState.CurrentChannelIndex, cptState.UseHistoneOverlay, cptState.schnitzcells, [], cptState.UseTwinTraces);

        elseif cc == 'u'
            [x2, ~, ~] = SpotsXYZ(cptState.Spots{cptState.CurrentChannelIndex}(cptState.CurrentFrame));
            
            if ~isempty(x2)
                ClickedSpot = ginput(1);
                
                UnfilterSpot(cptState.Spots{cptState.CurrentChannelIndex}, cptState.SpotFilter{cptState.CurrentChannelIndex}, ...
                    ClickedSpot, cptState.Particles{cptState.CurrentChannelIndex}, cptState.CurrentFrame)
            end

        elseif (cc == 'm') & (cptState.CurrentParticle < cptState.numParticles())
  
            [cptState.lineFitted, cptState.CurrentParticle, cptState.CurrentFrame, cptState.ManualZFlag, cptState.DisplayRange, cptState.TwinParticle] = ...
                goNextParticle(cptState.CurrentParticle, cptState.CurrentChannelIndex, cptState.HideApprovedFlag, cptState.Particles, cptState.UseTwinTraces);
        
        elseif (cc == 'n') & (cptState.CurrentParticle > 1)

            [cptState.lineFitted, cptState.CurrentParticle, cptState.CurrentFrame, cptState.ManualZFlag, cptState.DisplayRange, cptState.TwinParticle] = ...
                goPreviousParticle(cptState.CurrentParticle, cptState.CurrentChannelIndex, cptState.HideApprovedFlag, cptState.Particles, cptState.UseTwinTraces);

        elseif (cc == 'v') & (cptState.CurrentParticle < cptState.numParticles())

            [cptState.lineFitted, cptState.CurrentParticle, cptState.CurrentFrame, cptState.ManualZFlag, cptState.DisplayRange, cptState.TwinParticle] = ...
                goNextUncheckedParticle(cptState.CurrentParticle, cptState.CurrentChannelIndex, cptState.HideApprovedFlag, cptState.Particles, cptState.UseTwinTraces);
        
        elseif cc == 'i'
            warning(' AR 1/15/18: This is currently deprecated. Talk to HG if you need this function.')
        end
    end

    keyInputHandler = @keyInput;
end
