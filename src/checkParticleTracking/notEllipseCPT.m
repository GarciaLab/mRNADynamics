function ellipseHandle = notEllipseCPT(cptState, color, pointsToDraw, overlayAxes,...
        xTrace, yTrace, ZoomRange) % last 3 arguments added by G. Martini 2/13/21 to facilitate faster use of zoom mode
	ra = cptState.Ellipses{cptState.CurrentFrame}(:,3);
	rb = cptState.Ellipses{cptState.CurrentFrame}(:,4);
	ang = pi - cptState.Ellipses{cptState.CurrentFrame}(:,5);
	x0 = cptState.Ellipses{cptState.CurrentFrame}(:,1) + 1;
	y0 = cptState.Ellipses{cptState.CurrentFrame}(:,2) + 1;
    
    % Speed up plotting by excluding ellipses outside of ZoomRange wherever
    % possible
    if cptState.ZoomMode & exist('xTrace', 'var') & ~isempty(xTrace)
        EllipsesInZoomFrame = (abs(x0-xTrace)-ra <= ZoomRange) & ...
            (abs(y0-yTrace)-ra <= ZoomRange);
        ra = ra(EllipsesInZoomFrame);
        rb = rb(EllipsesInZoomFrame);
        ang = ang(EllipsesInZoomFrame);
        x0 = x0(EllipsesInZoomFrame);
        y0 = y0(EllipsesInZoomFrame);
        
    end
    

	hold(overlayAxes, 'on')
	ellipseHandle = notEllipse(ra, rb, ang, x0, y0, color, pointsToDraw, overlayAxes);
	hold(overlayAxes, 'off')
end
