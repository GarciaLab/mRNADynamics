function [textInputHandler, keyInputHandler] = FrameChangeEventHandler(cptState, Spots, robot, fake_event)
    numValidFrames = length({Spots{1}.Fits}); %check handle of spots is updated when is needed down the road
 
    function textInput(frame_num, event)
        figure(ancestor(frame_num, 'figure'));
        [cptState.CurrentFrame, cptState.ManualZFlag] = changeFrame(str2double(frame_num.Value), numValidFrames);
        robot.keyPress(fake_event);
        robot.keyRelease(fake_event);
        disp('Frame Changed.');
    end

    function keyInput(cc)
        if cc == '.' %Move forward one frame
            cptState.PreviousFrame = cptState.CurrentFrame;
            [cptState.CurrentFrame, cptState.ManualZFlag] = changeFrame(cptState.CurrentFrame + 1, numValidFrames);

        elseif (cc == ',') %Move backward one frame
            cptState.PreviousFrame = cptState.CurrentFrame;
            [cptState.CurrentFrame, cptState.ManualZFlag] = changeFrame(cptState.CurrentFrame - 1, numValidFrames);

        elseif (cc == '>') %Move forward five frames
            cptState.PreviousFrame = cptState.CurrentFrame;
            [cptState.CurrentFrame, cptState.ManualZFlag] = changeFrame(cptState.CurrentFrame + 5, numValidFrames);

        elseif (cc == '<') %#ok<*AND2>%Move backward five frames
            cptState.PreviousFrame = cptState.CurrentFrame;
            [cptState.CurrentFrame, cptState.ManualZFlag] = changeFrame(cptState.CurrentFrame - 5, numValidFrames);
 
        elseif cc == 'j'
            cptState.PreviousFrame = cptState.CurrentFrame;
            try
                iJump = inputdlg('Frame to jump to:', ...
                    'Move to frame');
                iJump = str2double(iJump{1});
            catch
                iJump = cptState.CurrentFrame;
            end

            [cptState.CurrentFrame, cptState.ManualZFlag] = changeFrame(iJump, numValidFrames);
        end
    end

    textInputHandler = @textInput;
    keyInputHandler = @keyInput;
end
