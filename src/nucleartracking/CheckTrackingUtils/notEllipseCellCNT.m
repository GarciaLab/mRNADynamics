function ellipseHandle = notEllipseCellCNT(cntState, schnitzCellNo, color, pointsToDraw, overlayAxes)
	ra = cntState.Ellipses{cntState.CurrentFrame}(schnitzCellNo, 3);
	rb = cntState.Ellipses{cntState.CurrentFrame}(schnitzCellNo, 4);
	ang = cntState.Ellipses{cntState.CurrentFrame}(schnitzCellNo, 5);
	x0 = cntState.Ellipses{cntState.CurrentFrame}(schnitzCellNo, 1) + 1;
	y0 = cntState.Ellipses{cntState.CurrentFrame}(schnitzCellNo, 2) + 1;

	hold(overlayAxes, 'on')
	ellipseHandle = notEllipse(ra, rb, ang, x0, y0, color, pointsToDraw, overlayAxes);
	hold(overlayAxes, 'off')
end
