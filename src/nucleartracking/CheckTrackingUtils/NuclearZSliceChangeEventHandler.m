function keyInputHandler = NuclearZSliceChangeEventHandler(cntState)
    
    function keyInput(cc)
        if (cc == 'a') & (cntState.CurrentZ < cntState.ZSlices) %Move up in Z
            [cntState.CurrentZ, cntState.ManualZFlag] = changeZSlice(cntState.CurrentZ + 1, cntState.ZSlices);
        elseif (cc == 'z') & (cntState.CurrentZ > 1) %Move down in Z
            [cntState.CurrentZ, cntState.ManualZFlag] = changeZSlice(cntState.CurrentZ - 1, cntState.ZSlices);
        elseif cc == 't'
            
            try
                iJump = inputdlg('z-slice to jump to:', ...
                    'Move to z-slice');
                iJump = str2double(iJump{1});
            catch
                iJump = cntState.CurrentFrame;
            end
            [cntState.CurrentZ, cntState.ManualZFlag] = changeZSlice(iJump, cntState.ZSlices);

            cntState.DisplayRange = [];
       end
    end

    keyInputHandler = @keyInput;
end
