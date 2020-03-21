function ellipsesFrame = registerEllipses(ellipsesFrame, movingImage, fixedImage)

xDim = size(movingImage, 2);
yDim = size(movingImage, 1);

[movingReg.DisplacementField,~] =...
    imregdemons(gpuArray(movingImage), gpuArray(fixedImage),100,'AccumulatedFieldSmoothing',1.0,'PyramidLevels',3);
Dx = movingReg.DisplacementField(:, :, 1); %displacement left to right
Dy = movingReg.DisplacementField(:, :, 2); %displacement top to bottom


for i = 1:length(ellipsesFrame)
    xOld = ellipsesFrame(i, 1);
    xOldSub = min(max(round(xOld), 1), xDim);
    yOld = ellipsesFrame(i, 2);
    yOldSub = min(max(round(yOld), 1), yDim);
    ellipsesFrame(i, 1) = xOld + gather(Dx(yOldSub, xOldSub));
    ellipsesFrame(i, 2) = yOld + gather(Dy(yOldSub, xOldSub));
end

end