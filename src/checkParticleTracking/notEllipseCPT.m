function ellipseHandle = notEllipseCPT(cptState, color, pointsToDraw, overlayAxes)
	ra = cptState.Ellipses{cptState.CurrentFrame}(:,3);
	rb = cptState.Ellipses{cptState.CurrentFrame}(:,4);
	ang = cptState.Ellipses{cptState.CurrentFrame}(:,5);
	x0 = cptState.Ellipses{cptState.CurrentFrame}(:,1) + 1;
	y0 = cptState.Ellipses{cptState.CurrentFrame}(:,2) + 1;

	hold(overlayAxes, 'on')
	ellipseHandle = notEllipse(ra, rb, ang, x0, y0, color, pointsToDraw, overlayAxes);
	hold(overlayAxes, 'off')
end
