function keyInputHandler = NuclearChangeEventHandler(cntState)
    
    function keyInput(cc)
        if cc == 'k'
        
            try
                NucleusJump = inputdlg('Particle to jump to:', ...
                    'Move to particle');
                NucleusJump = str2double(NucleusJump{1});
            catch
                NucleusJump = CurrentNucleus;
            end
            
            [cntState.CurrentNucleus, cntState.CurrentFrame, cntState.ManualZFlag] = ...
                changeNucleus(NucleusJump, cntState.schnitzcells, cntState.numNuclei);


            cntState.DisplayRange = [];
        
        elseif cc == 'r'
            cntState.schnitzcells = orderNuclei(cntState.numNuclei, cntState.schnitzcells);
        
        elseif cc == 'p'
            % Identify a particle. It will also tell you the particle associated with the clicked nucleus.
            identifyNucleus(cntState.schnitzcells, cntState.CurrentFrame, ...
                cntState.UseHistoneOverlay);
        
        elseif cc == '\'
            % Moves to clicked particle.
            [cntState.CurrentNucleus, cntState.CurrentFrame, cntState.ManualZFlag] =...
                toNearestNucleus(cntState.schnitzcells, ...
                cntState.CurrentFrame, cntState.UseHistoneOverlay);
        

        elseif (cc == 'm') & (cntState.CurrentNucleus < cntState.numNuclei())
            [cntState.CurrentNucleus, cntState.CurrentFrame, cntState.ManualZFlag, cntState.DisplayRange] = ...
                goNextNucleus(cntState.CurrentNucleus, cntState.HideApprovedFlag, cntState.schnitzcells);
        
        elseif (cc == 'n') & (cntState.CurrentNucleus > 1)
            [cntState.CurrentNucleus, cntState.CurrentFrame, cntState.ManualZFlag, cntState.DisplayRange] = ...
                goPreviousNucleus(cntState.CurrentNucleus, cntState.HideApprovedFlag, cntState.schnitzcells);

        elseif cc == 'i'
            warning(' AR 1/15/18: This is currently deprecated. Talk to HG if you need this function.')
        end
    end

    keyInputHandler = @keyInput;
end
