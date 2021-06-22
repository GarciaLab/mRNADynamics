function ellipseHandle = notEllipseCPT(cptState, color, pointsToDraw, overlayAxes)

    ra = cptState.Ellipses{cptState.CurrentFrame}(:,3);
    rb = cptState.Ellipses{cptState.CurrentFrame}(:,4);
    ang = pi - cptState.Ellipses{cptState.CurrentFrame}(:,5);
    x0 = cptState.Ellipses{cptState.CurrentFrame}(:,1) + 1;
    y0 = cptState.Ellipses{cptState.CurrentFrame}(:,2) + 1;

    if cptState.displayOnlyCurrentZEllipses
        disp("Displaying Ellipses only for the current Z slice.")
        zCenterSlices = cptState.Ellipses{cptState.CurrentFrame}(:,10);
        
        currentZ = cptState.CurrentZ;
        T = table(ra, rb, ang, x0, y0, zCenterSlices);
        
        CurrentEllipsesTable = T(T.zCenterSlices == currentZ,:);
        
        ra = CurrentEllipsesTable.ra;
        rb = CurrentEllipsesTable.rb;
        ang = CurrentEllipsesTable.ang;
        x0 = CurrentEllipsesTable.x0;
        y0 = CurrentEllipsesTable.y0;
    end
    
    hold(overlayAxes, 'on')
	ellipseHandle = notEllipse(ra, rb, ang, x0, y0, color, pointsToDraw, overlayAxes);
	hold(overlayAxes, 'off')
end
