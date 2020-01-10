 function keyInputHandler = ZoomAnywhereEventHandler(cptState)
    
    function keyInput(cc)
 		if cc == '+'
	        if ~cptState.ZoomMode & ~cptState.GlobalZoomMode
	                [ConnectPositionx, ConnectPositiony] = ginput(1);
	                cptState.xForZoom = round(ConnectPositionx);
	                cptState.yForZoom = round(ConnectPositiony);
	            cptState.GlobalZoomMode = true;
	        elseif cptState.ZoomMode & ~cptState.GlobalZoomMode
	            cptState.ZoomMode = false;
	        elseif ~cptState.ZoomMode & cptState.GlobalZoomMode
	            cptState.GlobalZoomMode = false;
	        end
	    end
    end

    keyInputHandler = @keyInput;
end


