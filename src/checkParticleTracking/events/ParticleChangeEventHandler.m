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
            
            [cptState.CurrentParticle, cptState.CurrentFrame, cptState.ManualZFlag] = ...
                changeParticle(ParticleJump, cptState.Particles, cptState.numParticles, cptState.CurrentChannel);

            cptState.DisplayRange = [];
        
        elseif cc == 'r'
            cptState.Particles = orderParticles(cptState.numParticles(), cptState.CurrentChannel, cptState.Particles);
        end
    end

    keyInputHandler = @keyInput;
end
