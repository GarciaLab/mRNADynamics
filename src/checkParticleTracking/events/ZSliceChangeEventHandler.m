function [textInputHandler, keyInputHandler] = ZSliceChangeEventHandler(cptState, robot, fake_event)
    
    function textInput(z_num, event)
        figure(ancestor(z_num, 'figure'));
        [cptState.CurrentZ, cptState.ManualZFlag] = changeZSlice(str2double(z_num.Value), cptState.ZSlices);
        robot.keyPress(fake_event);
        robot.keyRelease(fake_event);
        disp('ZSlice Changed.');
    end

    function keyInput(cc)
        if (cc == 'a') & (cptState.CurrentZ < cptState.ZSlices) %Move up in Z
            [cptState.CurrentZ, cptState.ManualZFlag] = changeZSlice(cptState.CurrentZ + 1, cptState.ZSlices);
        elseif (cc == 'z') & (cptState.CurrentZ > 1) %Move down in Z
            [cptState.CurrentZ, cptState.ManualZFlag] = changeZSlice(cptState.CurrentZ - 1, cptState.ZSlices);
        elseif cc == 't'
            
            try
                iJump = inputdlg('z-slice to jump to:', ...
                    'Move to z-slice');
                iJump = str2double(iJump{1});
            catch
                iJump = cptState.CurrentFrame;
            end
            [cptState.CurrentZ, cptState.ManualZFlag] = changeZSlice(iJump, cptState.ZSlices);

            cptState.DisplayRange = [];
       end
    end

    textInputHandler = @textInput;
    keyInputHandler = @keyInput;
end
