function ellipseHandle = notEllipseCellCPT(cptState, schnitzCellNo, color, pointsToDraw, overlayAxes)
	ra = cptState.Ellipses{cptState.CurrentFrame}(schnitzCellNo, 3);
	rb = cptState.Ellipses{cptState.CurrentFrame}(schnitzCellNo, 4);
	ang = cptState.Ellipses{cptState.CurrentFrame}(schnitzCellNo, 5);
	x0 = cptState.Ellipses{cptState.CurrentFrame}(schnitzCellNo, 1) + 1;
	y0 = cptState.Ellipses{cptState.CurrentFrame}(schnitzCellNo, 2) + 1;

	hold(overlayAxes, 'on')
	ellipseHandle = notEllipse(ra, rb, ang, x0, y0, color, pointsToDraw, overlayAxes);
	hold(overlayAxes, 'off')
end
