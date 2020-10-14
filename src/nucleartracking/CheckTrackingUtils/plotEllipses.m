function ellipseHandle = plotEllipses(Ellipses, CurrentFrame,  schnitzCellNo, color, pointsToDraw, overlayAxes)
	ra = Ellipses{CurrentFrame}(schnitzCellNo, 3);
	rb = Ellipses{CurrentFrame}(schnitzCellNo, 4);
	ang = Ellipses{CurrentFrame}(schnitzCellNo, 5);
	x0 = Ellipses{CurrentFrame}(schnitzCellNo, 1) + 1;
	y0 = Ellipses{CurrentFrame}(schnitzCellNo, 2) + 1;

	hold(overlayAxes, 'on')
	ellipseHandle = notEllipseThick(ra, rb, ang, x0, y0, color, pointsToDraw, overlayAxes);
	hold(overlayAxes, 'off')
end
