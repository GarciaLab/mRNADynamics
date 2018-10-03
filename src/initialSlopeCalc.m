function [coefficientsOfFittedLine,indexOfMax] =...
    initialSlopeCalc(frames,intensities)

[~, indexOfMax]= max(intensities);
regionOfRiseX = frames(1:indexOfMax);
regionOfRiseY = intensities(1:indexOfMax);

if length(regionOfRiseX) == 1
    coefficientsOfFittedLine = NaN;
else
    coefficientsOfFittedLine = polyfit(regionOfRiseX,regionOfRiseY,1);
end

end