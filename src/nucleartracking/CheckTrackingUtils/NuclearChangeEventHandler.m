function keyInputHandler = NuclearChangeEventHandler(cntState)
    
    function keyInput(cc)
        if cc == 'k'
        
            try
                NucleusJump = inputdlg('Schnitz cell to jump to:', ...
                    'Move to schnitzcell');
                NucleusJump = str2double(NucleusJump{1});
            catch
                NucleusJump = cntState.CurrentNucleus;
            end
            if NucleusJump ~= cntState.CurrentNucleus
                cntState.schnitzcells(cntState.CurrentNucleus).Checked = 1;
                
            end
            [cntState.CurrentNucleus, cntState.CurrentFrame, cntState.ManualZFlag] = ...
                changeNucleus(NucleusJump, cntState.schnitzcells, cntState.numNuclei);


            cntState.DisplayRange = [];
        
        elseif cc == 'r'
            cntState.schnitzcells = orderNuclei(cntState.numNuclei, cntState.schnitzcells);
        
        elseif cc == 'p'
            % Identify a Nucleus.
            identifyNucleus(cntState.schnitzcells, cntState.CurrentFrame, ...
                cntState.UseHistoneOverlay, cntState);
        
        elseif cc == '\'
            % Moves to clicked nucleus.

            cntState.schnitzcells(cntState.CurrentNucleus).Checked = 1;

            [cntState.CurrentNucleus, cntState.CurrentFrame, cntState.ManualZFlag] =...
                toNearestNucleus(cntState.schnitzcells, ...
                cntState.CurrentFrame, cntState.UseHistoneOverlay, cntState);
        

        elseif (cc == 'm') & (cntState.CurrentNucleus < cntState.numNuclei())
            cntState.schnitzcells(cntState.CurrentNucleus).Checked = 1;
            [cntState.CurrentNucleus, cntState.CurrentFrame, cntState.ManualZFlag, cntState.DisplayRange] = ...
                goNextNucleus(cntState.CurrentNucleus, cntState.HideApprovedFlag, cntState.schnitzcells);
        
        elseif (cc == 'n') & (cntState.CurrentNucleus > 1)
            cntState.schnitzcells(cntState.CurrentNucleus).Checked = 1;
            [cntState.CurrentNucleus, cntState.CurrentFrame, cntState.ManualZFlag, cntState.DisplayRange] = ...
                goPreviousNucleus(cntState.CurrentNucleus, cntState.HideApprovedFlag, cntState.schnitzcells);

        end
    end

    keyInputHandler = @keyInput;
end
