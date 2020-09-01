function keyInputHandler = FrameChangeEventHandler(cptState)
 
    function keyInput(cc)
        CC = cptState.CurrentChannelIndex;
        CP = cptState.CurrentParticle;
        
        numValidFrames = length({cptState.Spots{CC}.Fits}); %check handle of spots is updated when is needed down the road
        
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

            cptState.DisplayRange = [];

        elseif (cc == '''') & (cptState.CurrentFrame < cptState.numValidFrames())
            % Move to the next skipped frame within the particle
            cptState.PreviousFrame = cptState.CurrentFrame;
            cptState.CurrentFrame = nextSkippedFrame(cptState.Particles, CC, ...
                CP, cptState.CurrentFrame);
        
        elseif (cc == ';') & (cptState.CurrentFrame > 1)
            % Move to the previous skipped frame within the particle
            cptState.PreviousFrame = cptState.CurrentFrame;
            cptState.CurrentFrame = previousSkippedFrame(cptState.Particles, CC, ...
                CP, cptState.CurrentFrame);
        
        elseif cc == 'e'
            % Approve/Disapprove a frame within a trace
            FrameFilter = cptState.Particles{CC}(CP).Frame == cptState.CurrentFrame;
            cptState.Particles{CC}(CP).FrameApproved(FrameFilter) = ...                
                ~cptState.Particles{CC}(CP).FrameApproved(...
                FrameFilter);
            if cptState.Particles{CC}(CP).FrameApproved(FrameFilter)  
              for v = 1:length(cptState.qcFlagFields)
                if cptState.Particles{CC}(CP).(cptState.qcFlagFields{v})(FrameFilter) >= 1
                  cptState.Particles{CC}(CP).(cptState.qcFlagFields{v})(FrameFilter) = 0;
                end              
              end
            end           
        end
    end

    keyInputHandler = @keyInput;
end
