function keyInputHandler = NuclearFrameChangeEventHandler(cntState)
 
    function keyInput(cc)
        numValidFrames = length(cntState.Ellipses); %check handle of spots is updated when is needed down the road

        if cc == '.' %Move forward one frame
            cntState.PreviousFrame = cntState.CurrentFrame;
            [cntState.CurrentFrame, ~] = changeFrame(cntState.CurrentFrame + 1, numValidFrames);

        elseif (cc == ',') %Move backward one frame
            cntState.PreviousFrame = cntState.CurrentFrame;
            [cntState.CurrentFrame, cntState.ManualZFlag] = changeFrame(cntState.CurrentFrame - 1, numValidFrames);

         
        elseif (cc == '>') %Move forward five frames
            cntState.PreviousFrame = cntState.CurrentFrame;
            [cntState.CurrentFrame, cntState.ManualZFlag] = changeFrame(cntState.CurrentFrame + 5, numValidFrames);

        elseif (cc == '<') %#ok<*AND2>%Move backward five frames
            cntState.PreviousFrame = cntState.CurrentFrame;
            [cntState.CurrentFrame, cntState.ManualZFlag] = changeFrame(cntState.CurrentFrame - 5, numValidFrames);
 
        elseif cc == 'j'
            cntState.PreviousFrame = cntState.CurrentFrame;
            try
                iJump = inputdlg('Frame to jump to:', ...
                    'Move to frame');
                iJump = str2double(iJump{1});
            catch
                iJump = cntState.CurrentFrame;
            end

            [cntState.CurrentFrame, cntState.ManualZFlag] = changeFrame(iJump, numValidFrames);

            cntState.DisplayRange = [];

        elseif (cc == ':') & (cntState.CurrentFrame < cntState.numValidFrames())
            % Move to the next skipped frame within the particle
            cntState.PreviousFrame = cntState.CurrentFrame;
            cntState.CurrentFrame = nextNuclearUnapprovedFrame(cntState.schnitzcells, ...
                cntState.CurrentNucleus, cntState.CurrentFrame);
        
        elseif (cc == ';') & (cntState.CurrentFrame > 1)
            % Move to the previous skipped frame within the particle
            cntState.PreviousFrame = cntState.CurrentFrame;
            cntState.CurrentFrame = previousNuclearUnapprovedFrame(cntState.schnitzcells, ...
                cntState.CurrentNucleus, cntState.CurrentFrame);
        
        elseif cc == 'e'
            % Approve/Disapprove a frame within a trace
            cntState.schnitzcells...
                (cntState.CurrentNucleus).FrameApproved(cntState.schnitzcells...
                (cntState.CurrentNucleus).frames == cntState.CurrentFrame) = ...
                ...
                ~cntState.schnitzcells(cntState.CurrentNucleus).FrameApproved(...
                cntState.schnitzcells(cntState.CurrentNucleus).frames == cntState.CurrentFrame);
        end
    end

    keyInputHandler = @keyInput;
end
