function ellipsesFrame = registerEllipses(ellipsesFrame, movingImage, fixedImage)

[movingReg.DisplacementField,~] =...
    imregdemons(gpuArray(movingImage), gpuArray(fixedImage),100,'AccumulatedFieldSmoothing',1.0,'PyramidLevels',3);
Dx = movingReg.DisplacementField(:, :, 1); %displacement left to right
Dy = movingReg.DisplacementField(:, :, 2); %displacement top to bottom


for i = 1:length(ellipsesFrame)
    xOld = ellipsesFrame(i, 1);
    xOldSub = round(xOld);
    yOld = ellipsesFrame(i, 2);
    yOldSub = round(yOld);
    ellipsesFrame(i, 1) = xOld + gather(Dx(yOldSub, xOldSub));
    ellipsesFrame(i, 2) = yOld + gather(Dy(yOldSub, xOldSub));
end

end