function ellipseFrame = findNuclearXYthenAdjustContours(...
    FrameInfo, hisImage, pixelSize_um, nuclearCycle)

diameters_um = [nan(1,8), 7.26, 7.26, 6.16, 5.72, 4.84, 3.96]; %nuclear diameters from ncs 1 to 14.
diameters_px = diameters_um / pixelSize_um;

ellipseFrame = findNuclei(FrameInfo, hisImage,...
    diameters_um(nuclearCycle) );

%we need some rough circles here before refinement
ellipseFrame(:, 3) = diameters_px(nuclearCycle);
ellipseFrame(:, 4) = diameters_px(nuclearCycle);
ellipseFrame(:, 5) = 0; %ellipse orientation angle

ellipseFrame = adjustNuclearContours(ellipseFrame, hisImage, pixelSize_um);


end