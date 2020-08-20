function keyInputHandler = HistoneContrastChangeEventHandler(cptState)
 
    function keyInput(cc)
        if cc == 'g' & cptState.UseHistoneOverlay %Increase histone channel contrast
        
            if isempty(cptState.DisplayRange)'
                cptState.DisplayRange = [min(min(cptState.ImageHis)), max(max(cptState.ImageHis)) / 1.5];
            else
                cptState.DisplayRange = [cptState.DisplayRange(1), cptState.DisplayRange(2) / 1.5];
            end
            
            disp('increased nuclear contrast');

        elseif cc == 'b' & cptState.UseHistoneOverlay %Decrease histone channel contrast
            cptState.DisplayRange = [min(min(cptState.ImageHis)), max(max(cptState.ImageHis)) * 1.5];
            disp('decreased nuclear contrast');
        end
    end

    keyInputHandler = @keyInput;
end
