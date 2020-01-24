function [textInputHandler, keyInputHandler] = DeleteSpotEventHandler(cptState, robot, fake_event)
 
    function doDeleteSpot(frame)
        [cptState.Spots, cptState.SpotFilter, cptState.CurrentFrame, ...
            cptState.CurrentParticle, cptState.Particles, cptState.ManualZFlag, cptState.lastParticle, cptState.PreviousParticle] = ...
            removeSpot(cptState.Frames, frame, ...
            cptState.CurrentChannel, cptState.CurrentParticle, cptState.CurrentParticleIndex, cptState.Particles, cptState.Spots, cptState.SpotFilter);
    end

    function textInput(delete_spot, event)
        Overlay = ancestor(delete_spot, 'figure');
        figure(Overlay);

        cptState.no_clicking = true;

        doDeleteSpot(cptState.CurrentFrame);
        
        robot.keyPress(fake_event);
        robot.keyRelease(fake_event);
        cptState.no_clicking = false;
    end

    function keyInput(cc)
        if cc == '#' %remove a spot from cptState.Spots and erase its frame in Particles
            doDeleteSpot(cptState.CurrentFrame);
            
        elseif cc == '^' %remove a whole trace from cptState.Spots and Particles. AR 7/9/2019 a work in progress
            for frame = 1:length(cptState.Frames)
                doDeleteSpot(frame);
            end
        end
    end

    textInputHandler = @textInput;
    keyInputHandler = @keyInput;
end
