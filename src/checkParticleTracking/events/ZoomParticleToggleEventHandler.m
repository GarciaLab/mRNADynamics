 function keyInputHandler = ZoomParticleToggleEventHandler(cptState)
    
    function keyInput(cc)
     	if cc == 'o'
	        if ~cptState.GlobalZoomMode
	            cptState.ZoomMode = ~cptState.ZoomMode;
	        elseif cptState.GlobalZoomMode & ~cptState.ZoomMode
	            cptState.GlobalZoomMode = false;
	        end
	    end
    end

    keyInputHandler = @keyInput;
end


