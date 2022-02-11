function keyInputHandler = NuclearTracesEventHandler(cntState, FrameInfo, ncFrames)
    % NuclearTracesEventHandler.m
    % author: Gabriella Martini
    % date created: 9/7/20
    % date last modified: 9/13/20
    function keyInput(cc)
        if cc == 'c'
            [cntState.CurrentNucleus, cntState.schnitzcells] = combineNuclearTraces(cntState.CurrentFrame,...
                cntState.schnitzcells, cntState.Ellipses, cntState.CurrentNucleus, FrameInfo, ncFrames);
        elseif cc == 'v'
            [cntState.CurrentNucleus, cntState.schnitzcells] = backwardsCombineNuclearTraces(cntState.CurrentFrame,...
                cntState.schnitzcells, cntState.Ellipses, cntState.CurrentNucleus, FrameInfo, ncFrames);
            
        elseif cc == 'd'
            % Separate traces forward at the current frame.
            frame_idx = find(cntState.schnitzcells(cntState.CurrentNucleus).frames == cntState.CurrentFrame);
            if ~isempty(frame_idx)
                if frame_idx(1) > 1 
                    cntState.schnitzcells = SeparateNuclearTraces(cntState.CurrentNucleus, ...
                        cntState.CurrentFrame, cntState.schnitzcells, FrameInfo, ncFrames);
                end
            end
%         elseif cc == 'D'
%             
%             if length(cntState.schnitzcells(cntState.CurrentNucleus).frames) > 1
%                 HoldCurrentNucleus = cntState.CurrentNucleus;
%                 cntState.schnitzcells=SeparateAllFrameNuclearTraces(cntState.CurrentNucleus,...
%                     cntState.CurrentFrame, cntState.schnitzcells, FrameInfo, ncFrames);
%                 [cntState.CurrentNucleus,cntState.CurrentFrame, cntState.ManualZFlag] = ...
%                     changeNucleus(HoldCurrentNucleus, cntState.schnitzcells, cntState.numNuclei())
%             end
% 
%             
            
        elseif cc == 'q'
            % Approve a trace
            
            if (cntState.schnitzcells(cntState.CurrentNucleus).Approved == 1)
                cntState.schnitzcells(cntState.CurrentNucleus).Approved = 2;
            elseif cntState.schnitzcells(cntState.CurrentNucleus).Approved <= 0
                cntState.schnitzcells(cntState.CurrentNucleus).Approved = 1;
            elseif cntState.schnitzcells(cntState.CurrentNucleus).Approved == 2
                cntState.schnitzcells(cntState.CurrentNucleus).Approved = 0;
            end
            cntState.schnitzcells(cntState.CurrentNucleus).Checked = 1;
            
        elseif cc == 'w'
            % Disapprove a trace
            cntState.schnitzcells(cntState.CurrentNucleus).Approved = 0;
            cntState.schnitzcells(cntState.CurrentNucleus).Checked = 1;
        elseif cc == 'f'
            try
                [flag, flag_string]  = chooseFlag;
                disp(flag_string)
                cntState.schnitzcells(cntState.CurrentNucleus).Flag = flag;
                cntState.schnitzcells(cntState.CurrentNucleus).Checked = 1;
            catch
                disp('No Flag Selected')
            end
        elseif cc == '1'
            cntState.schnitzcells(cntState.CurrentNucleus).Flag = 1;
            cntState.schnitzcells(cntState.CurrentNucleus).Checked = 1;
        elseif cc == '2'
            cntState.schnitzcells(cntState.CurrentNucleus).Flag = 2;
            cntState.schnitzcells(cntState.CurrentNucleus).Checked = 1;
        elseif cc == '3'
            cntState.schnitzcells(cntState.CurrentNucleus).Flag = 3;
            cntState.schnitzcells(cntState.CurrentNucleus).Checked = 1;
        elseif cc == '4'
            cntState.schnitzcells(cntState.CurrentNucleus).Flag = 4;
            cntState.schnitzcells(cntState.CurrentNucleus).Checked = 1;
        elseif cc == '5'
            cntState.schnitzcells(cntState.CurrentNucleus).Flag = 5;
            cntState.schnitzcells(cntState.CurrentNucleus).Checked = 1;
        elseif cc == '6'
            cntState.schnitzcells(cntState.CurrentNucleus).Flag = 6;
            cntState.schnitzcells(cntState.CurrentNucleus).Checked = 1;
        elseif cc == '7'
            cntState.schnitzcells(cntState.CurrentNucleus).Flag = 7;
            cntState.schnitzcells(cntState.CurrentNucleus).Checked = 1;
        elseif cc == '0'
            cntState.schnitzcells(cntState.CurrentNucleus).Flag = 0;
            cntState.schnitzcells(cntState.CurrentNucleus).Checked = 1;
        elseif cc == 'a'
            % Indicate anaphase frame for nucleus 
            CurrentNucleusInFrame = sum(find(cntState.schnitzcells(cntState.CurrentNucleus).frames == cntState.CurrentFrame, 1));
            if CurrentNucleusInFrame
                if ~isfield(cntState.schnitzcells, 'anaphaseFrame')
                   for i=1:cntState.numNuclei()
                       cntState.schnitzcells(i).anaphaseFrame = [];
                   end
                end
                if ~isfield(cntState.schnitzcells, 'inferredAnaphaseFrame')
                   for i=1:cntState.numNuclei()
                       cntState.schnitzcells(i).anaphaseFrame = false;
                   end
                end
                CurrentNucleusInAnyPreviousFrame  = sum(find(cntState.schnitzcells(cntState.CurrentNucleus).frames < cntState.CurrentFrame, 1));
                if CurrentNucleusInAnyPreviousFrame
                    cntState.schnitzcells = SeparateNuclearTraces(cntState.CurrentNucleus, ...
                        cntState.CurrentFrame, cntState.schnitzcells, FrameInfo, ncFrames);
                    
                    cntState.schnitzcells(cntState.CurrentNucleus).anaphaseFrame = [];
                    cntState.schnitzcells(cntState.CurrentNucleus).inferredAnaphaseFrame = false;
                    
                    cntState.CurrentNucleus = cntState.CurrentNucleus+1;
                end
                
                cntState.schnitzcells(cntState.CurrentNucleus).anaphaseFrame = cntState.CurrentFrame;
                
                cntState.schnitzcells(cntState.CurrentNucleus).inferredAnaphaseFrame = false;
                
                
                if isfield(cntState.schnitzcells,'timeSinceAnaphase')
                    disp('Updating anaphase frame information')
                    ncFrames(ncFrames==0) = 1;
                    ind = find(isnan(ncFrames));
                    ncFrames(ind) = ncFrames(ind-1);
                    time = [FrameInfo.Time]/60; %frame times in minutes 
                    cntState.schnitzcells(cntState.CurrentNucleus).timeSinceAnaphase = time(cntState.schnitzcells(cntState.CurrentNucleus).frames) - time(cntState.schnitzcells(cntState.CurrentNucleus).anaphaseFrame);
                end
            end

        end
    end

    keyInputHandler = @keyInput;
end
