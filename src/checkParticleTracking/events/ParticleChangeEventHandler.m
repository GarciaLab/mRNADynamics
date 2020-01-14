function [textInputHandler, keyInputHandler] = ParticleChangeEventHandler(cptState, robot, fake_event)
    
    function textInput(particle_num, event)
        figure(ancestor(particle_num, 'figure'));
        
        [cptState.CurrentParticle, cptState.CurrentFrame, cptState.ManualZFlag] = changeParticle(...
            str2double(particle_num.Value), cptState.Particles, cptState.numParticles, cptState.CurrentChannel);

        robot.keyPress(fake_event);
        robot.keyRelease(fake_event);
        disp('Particle Changed.');
    end

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
        end
    end

    textInputHandler = @textInput;
    keyInputHandler = @keyInput;
end
