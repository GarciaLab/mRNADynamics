function ellipseHandle = notEllipseCNT(cntState, color, pointsToDraw, overlayAxes)
	ra = cntState.Ellipses{cntState.CurrentFrame}(:,3);
	rb = cntState.Ellipses{cntState.CurrentFrame}(:,4);
	ang = pi - cntState.Ellipses{cntState.CurrentFrame}(:,5);
	x0 = cntState.Ellipses{cntState.CurrentFrame}(:,1) + 1;
	y0 = cntState.Ellipses{cntState.CurrentFrame}(:,2) + 1;

	hold(overlayAxes, 'on')
	ellipseHandle = notEllipse(ra, rb, ang, x0, y0, color, pointsToDraw, overlayAxes);
	hold(overlayAxes, 'off')
end
