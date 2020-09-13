function ellipseHandle = ellipseCellCNT(cntState, aCell, color, pointsToDraw, overlayAxes)
	ra = cntState.Ellipses{cntState.CurrentFrame}(aCell, 3);
	rb = cntState.Ellipses{cntState.CurrentFrame}(aCell, 4);
	ang = pi - cntState.Ellipses{cntState.CurrentFrame}(aCell, 5);
	x0 = cntState.Ellipses{cntState.CurrentFrame}(aCell, 1) + 1;
	y0 = cntState.Ellipses{cntState.CurrentFrame}(aCell, 2) + 1;

	hold(overlayAxes, 'on')
	ellipseHandle = ellipse(ra, rb, ang, x0, y0, color, pointsToDraw, overlayAxes);
	hold(overlayAxes, 'off')
end
