function ellipseHandle = ellipseCellCPT(cptState, aCell, color, pointsToDraw, overlayAxes)
	ra = cptState.Ellipses{cptState.CurrentFrame}(aCell, 3);
	rb = cptState.Ellipses{cptState.CurrentFrame}(aCell, 4);
	ang = pi - cptState.Ellipses{cptState.CurrentFrame}(aCell, 5);
	x0 = cptState.Ellipses{cptState.CurrentFrame}(aCell, 1) + 1;
	y0 = cptState.Ellipses{cptState.CurrentFrame}(aCell, 2) + 1;

	hold(overlayAxes, 'on')
	ellipseHandle = ellipse(ra, rb, ang, x0, y0, color, pointsToDraw, overlayAxes);
	hold(overlayAxes, 'off')
end
