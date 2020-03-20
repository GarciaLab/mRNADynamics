function ellipsesFrame = registerEllipses(ellipsesFrame, movingImage, fixedImage)


[movingReg.DisplacementField,~] = imregdemons(movingImage, fixedImage,100,'AccumulatedFieldSmoothing',1.0,'PyramidLevels',3);
Dx = movingReg.DisplacementField(:, :, 1); %displacement left to right
Dy = movingReg.DisplacementField(:, :, 2); %displacement top to bottom


for i = 1:length(ellipsesFrame)
    xOld = ellipsesFrame(i, 1);
    xOldSub = round(xOld);
    yOld = ellipsesFrame(i, 2);
    yOldSub = round(yOld);
    ellipsesFrame(i, 1) = xOld + Dx(yOldSub, xOldSub);
    ellipsesFrame(i, 2) = yOld + Dy(yOldSub, xOldSub);
end

end