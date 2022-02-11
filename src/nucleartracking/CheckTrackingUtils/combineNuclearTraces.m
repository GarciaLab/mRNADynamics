function [CurrentNucleus, schnitzcells] =...
    combineNuclearTraces(CurrentFrame, schnitzcells, Ellipses, CurrentNucleus, FrameInfo, ncFrames)
%COMBINETRACES Summary of this function goes here
%   Detailed explanation goes here

PreviousNucleus=0;
debug_mode = false;
exitConnectFlag = 0;

nucleiExistInFrame = size(Ellipses{CurrentFrame}, 1);
currentNucleusExistsInPreviousFrame =...
    sum(schnitzcells(CurrentNucleus).frames<CurrentFrame);

if currentNucleusExistsInPreviousFrame && nucleiExistInFrame
    
    [ConnectPositionx,ConnectPositiony]=ginput(1);
    ConnectPosition = [ConnectPositionx,ConnectPositiony];
    
    if ~isempty(ConnectPosition)
        % find index of the particle we want to add (a.k.a output particle) to current
        % particle (current particle)
        [NucleusOutput,~]=FindClickedNucleus(ConnectPosition,CurrentFrame,...
            schnitzcells);
        if CurrentNucleus ~= NucleusOutput
            if debug_mode
                disp(['Current Nucleus: ', num2str(CurrentNucleus), ', Clicked Nucleus: ',num2str(NucleusOutput)])
                curr_frames = schnitzcells(CurrentNucleus).frames;
                fmt = ['Current Nucleus Frames: [', repmat('%g, ', 1, numel(curr_frames)-1), '%g]\n'];
                fprintf(fmt, curr_frames)
                clicked_frames = schnitzcells(NucleusOutput).frames;
                fmt = ['Clicked Nucleus Frames: [', repmat('%g, ', 1, numel(clicked_frames)-1), '%g]\n'];
                fprintf(fmt, clicked_frames)
            end
            %Check that the clicked particle doesn't exist in a previous
            %frame, that there is no overlap of frames. If it does
            %exist in a previous frame we will have to disconnect it.
            
            
            clickedNucleusExistsInAnyPreviousFrame = sum(schnitzcells(NucleusOutput).frames<CurrentFrame);
            
            
            if clickedNucleusExistsInAnyPreviousFrame
                schnitzcells = SeparateNuclearTraces(NucleusOutput, ...
                    CurrentFrame, schnitzcells, FrameInfo, ncFrames);
                schnitzcells(NucleusOutput).Approved = 1;
                schnitzcells(NucleusOutput).Checked = 0;
                schnitzcells(NucleusOutput).FirstPass = 1;
                if NucleusOutput <= CurrentNucleus
                    CurrentNucleus = CurrentNucleus+1;
                end
                NucleusOutput = NucleusOutput + 1;
                if debug_mode
                    disp('Clicked Nucleus Exists in a previous frame. Separating clicked nucleus previous frames from current frame');
                    disp(['New Current Nucleus Index: ', num2str(CurrentNucleus), ', New Clicked Nucleus Index: ',num2str(NucleusOutput)])
                    clicked_frames = schnitzcells(NucleusOutput).frames;
                    fmt = ['Clicked Nucleus Frames: [', repmat('%g, ', 1, numel(clicked_frames)-1), '%g]\n'];
                    fprintf(fmt, clicked_frames)
                end
            end
            currentNucleusExistsInAnyFollowingFrame = sum(schnitzcells(CurrentNucleus).frames>=CurrentFrame);
            if currentNucleusExistsInAnyFollowingFrame
                approved_status = schnitzcells(CurrentNucleus).Approved;
                schnitzcells = SeparateNuclearTraces(CurrentNucleus, ...
                    CurrentFrame, schnitzcells, FrameInfo, ncFrames);
                schnitzcells(CurrentNucleus+1).Approved = approved_status;
                schnitzcells(CurrentNucleus+1).FirstPass = 1;
                schnitzcells(NucleusOutput).Checked = 0;
                if NucleusOutput >= CurrentNucleus
                    NucleusOutput = NucleusOutput+1;
                end
                if debug_mode
                    disp('Current Nucleus Exists in current or future frame. Separating current nucleus current frames from previous frames');
                    disp(['New Current Nucleus Index: ', num2str(CurrentNucleus), ', New Clicked Nucleus Index: ',num2str(NucleusOutput)])
                    curr_frames = schnitzcells(CurrentNucleus).frames;
                    fmt = ['Current Nucleus Frames: [', repmat('%g, ', 1, numel(curr_frames)-1), '%g]\n'];
                    fprintf(fmt, curr_frames)
                    
                end
                
            end
            
            
            schnitzcells=JoinNuclearTraces(CurrentNucleus,NucleusOutput,schnitzcells, FrameInfo, ncFrames);
            %Deals with the indexing changing because of the removal of
            %the old particle.
            
            if NucleusOutput<CurrentNucleus
                CurrentNucleus=CurrentNucleus-1;
            end
            schnitzcells(CurrentNucleus).FirstPass = 1;
            if debug_mode
                disp(['Current Nucleus: ', num2str(CurrentNucleus), ', Clicked Nucleus: ',num2str(NucleusOutput)])
                curr_frames = schnitzcells(CurrentNucleus).frames;
                fmt = ['Current Nucleus Frames: [', repmat('%g, ', 1, numel(curr_frames)-1), '%g]\n'];
                fprintf(fmt, curr_frames)
            end
        end
    end
    
else
    
    disp('Cannnot connect to two schnitz cells!')
    
    
end
end

